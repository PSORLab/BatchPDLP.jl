
function main_loop_rHalpern_kernel(
    solutions,                        # [n_LPs, n_vars−1] = Holds solutions to LPs
    objectives,                       # [n_LPs] = Holds dual objectives
    original_variable_lower_bounds,   # [n_LPs, n_vars] = Used for getting unscaled convergence info
    original_variable_upper_bounds,   # [n_LPs, n_vars] = ""
    original_right_hand_side,         # [n_LPs × total_LP_length] = RHS h before scaling
    original_objective_vector,        # [n_LPs, n_vars] = Cost vector c before scaling of each LP objective
    original_objective_constant,      # [n_LPs] = the offset constant in the objective function of each LP
    # =============================================================================================================================#
    # The original_* versions are only used to compute unscaled convergence statistics
    # the algorithm runs entirely on the scaled_* versions.
    # =============================================================================================================================#
    scaled_variable_lower_bounds,     # [n_LPs, n_vars] = Lower bounds after initial rescaling using Ruiz/Pock-Chambolle rescaling
    scaled_variable_upper_bounds,     # [n_LPs, n_vars] = Upper bounds after initial rescaling using Ruiz/Pock-Chambolle rescaling
    scaled_constraint_matrix,         # [n_LPS × total_LP_length, n_vars] = Constraint matrix G~ after initial rescaling
    scaled_right_hand_side,           # [n_LPs × total_LP_length] = Right-hand h~ side after initial rescaling
    scaled_objective_vector,          # [n_LPs, n_vars] = Objective vector c~ after initial rescaling
    nz_count,                         # Int = Total number of nonzeros in the constraint matrix
    nz_rows,                          # [nz_count] = Row indices containing nonzeros in the constraint matrix
    nz_cols,                          # [nz_count] = Column indices containing nonzeros in the constraint matrix
    active_constraint,                # [n_LPs × total_LP_length] = Flag to indicate whether a constraint is in use.
    variable_rescaling,               # [n_LPs, n_vars] = Total value used for variable rescaling - Per-variable scale factor Dx so that unscaled_x = scaled_x / D_x
    constraint_rescaling,             # [n_LPs × total_LP_length] = Per-constraint scale factor Dc - ​Total value used for constraint rescaling
    # input_current_primal_solution,    # Usually 0s, but used for hot-starts (if desired in the future)
    # input_current_dual_solution,      # Usually 0s, but used for hot-starts (if desired in the future)
    # input_current_primal_product,     # Usually 0s, but used for hot-starts (if desired in the future)
    # input_current_dual_product,       # Usually 0s, but used for hot-starts (if desired in the future)
    # Last Restart Solutions
    # Current Iterate State - These are updated every PDLP step and represent where the algorithm currently is
    current_primal_solution,          # [n_LPs, n_vars] = Current x~ (scaled) - Currently active solution state
    current_dual_solution,            # [n_LPs × total_LP_length] = Current y~ (scaled)
    current_dual_product,             # [n_LPs, n_vars] = Cached G~^T * y~ to avoid recomputing
    current_primal_product,           # [n_LPs × total_LP_length] = Cached G~ * x~ to avoid recomputing (scaled)
    current_primal_gradient,           # [n_LPs, n_vars] = 
    initial_primal_solution,
    initial_dual_solution,
    pdhg_primal_solution,
    pdhg_dual_solution,
    reflected_primal_solution, 
    reflected_dual_solution,
    dual_slack,
    residual_primal_product,
    residual_dual_product,
    primal_residual,                  # [n_LPs × total_LP_length]]
    primal_slack,            # [n_LPs × total_LP_length]]
    dual_residual,                    # [n_LPs , n_vars]]
    # Unscaled Solution for convergence checking
    original_primal_solution,         # [n_LPs, n_vars] = Unscaled solution state information, used to calculate actual infeasibility
    original_primal_gradient,         # [n_LPs, n_vars] = Unscaled primal gradient
    original_dual_solution,           # [n_LPs × total_LP_length] = Unscaled dual solution y
    original_primal_product,          # [n_LPs × total_LP_length] = G~ * x~ unscaled (== G * x)
    # Step Computation Buffers (scratch space per PDLP step)
    buffer_kkt_primal_solution,       # [n_LPs, n_vars] = normalized x~ with primal_ray_norm for infeasibility ray checks Intermediate storage for calculating infeasibility
    buffer_kkt_primal_product,        # [n_LPs × total_LP_length] = normalized original_primal_product using primal_ray_norm (G~ * buffer_kkt_primal_solution)
    buffer_kkt_lower_variable_violation, # [n_LPs, n_vars] = max(l - x, 0) : Infeasible variable/constraint violation storage
    buffer_kkt_upper_variable_violation, # [n_LPs, n_vars] = max(x - u, 0) : Infeasible variable/constraint violation storage
    buffer_kkt_reduced_costs,         # [n_LPs, n_vars] = 
    delta_primal,                     # [n_LPs, n_vars] = Δx~ = x'~ - x~ : the value used in calculating PDLP steps
    delta_primal_product,             # [n_LPs × total_LP_length] = G~ * Δx~
    delta_dual,                       # [n_LPs × total_LP_length] = Δy~ = y'~ - y~ : the value used in calculating PDLP steps
    # Parameters for Control & Termination
    input_primal_weight,              # [n_LPs] = The original primal weight (could be moved internally if needed)
    input_step_size,                  # [n_LPs] = The original step size η (could be moved internally if needed)
    termination_reason,               # [n_LPs] = TerminationReason enum per LP : Field to let the user know why the LP terminated
    current_LP_length,                # Int = The current length of each individual LP (i.e., number of constraints)
    total_LP_length,                  # Int = The total possible size of each individual LP (i.e., max allowed number of constraints)
    n_LPs,                            # Int = The current number of LPs being solved (fewer than the max may be used during B&B)
    n_vars,                           # Int = The number of primal variables in the problem
    iteration_limit,                  # Int = The maximum number of PDLP steps allowed before termination
    kkt_matrix_pass_limit,            # Float = The maximum number of KKT matrix passes (generally unused)
    termination_evaluation_frequency, # Number of PDLP steps to take before checking termination criteria (default: 200)
    necessary_reduction_for_restart,  # Float = β_necessary = Necessary KKT error decrease factor needed for a restart if no progress is being made (default: 0.5)
    sufficient_reduction_for_restart, # Float = β_sufficient = Sufficient KKT error decrease factor that would trigger a restart immediately (default: 0.2)
    artificial_ratio_for_restart,     # Float = β_artificial = artificial restart to avoid long inner loop (default = )
    extrapolation_coefficient,        # Float = Value used in PDLP steps (default: 1.0)
    reflection_coefficient,           # Float = Reflection coefficient γ ∈ [0,1] at reflection step (default: 1.0)
    pid_KP,                           # Float = Primal Weight Update PID controller proportional coefficient (default: 0.99)
    pid_KI,                           # Float = Primal Weight Update PID controller integral coefficient (default: 0.01)
    pid_KD,                           # Float = Primal Weight Update PID controller derivative coefficient (default: 0)
    i_smooth,                         # Float = smoothing factor for restart error sum
    abs_tol,                          # Float = Absolute tolerance for termination
    rel_tol,                          # Float = Relative tolerance for termination
    eps_primal_infeasible,            # Float = Primal infeasibility tolerance
    eps_dual_infeasible,              # Float = Dual infeasibility tolerance
    return_code,                      # Int = Indicator for returning primal obj (1), dual obj (2), or both (3)
    global_upper_bound,               # Float = Information about the B&B upper bound (PDLP terminates if a dual feasible solution is above this value)
    skip_hard_problems,               # Bool = Flag to skip problems with too many iterations
    global_counter,                   # [1] (atomic Int32) = Count of LPs being solved - Running count of successfully completed LPs (shared across blocks via atomics)
    iteration_counter,                # [1] (atomic Int32) = Count of total iterations for solved LPs - Total iterations over all completed LPs (used to compute average for skip_hard_problems)
    skip_flag,                        # [n_LPs] = Flag to completely skip an individual LP - Per-LP flag to skip entirely (e.g., already pruned in B&B)
    iterations,                       # [n_LPs] = The final number of iterations needed for each LP
    temp,
    temp_dual,
    residual_delta_dual,
    )
    # In this kernel, assume that the number of blocks is equal to the number
    # of LPs, so that all threads are working on one LP
    LP = blockIdx().x
    idx = threadIdx().x
    
    block_stride = blockDim().x # total number of threads in the block
    grid_stride = gridDim().x # total number of blocks in the grid

    # Calculate strides for parallel reductions as the largest power of 2
    # less than n_vars (or current_LP_length)
    var_stride = Int32(1) << floor(Int32, log2(n_vars))
    len_stride = Int32(1) << floor(Int32, log2(current_LP_length))

    # Set up dynamic shared space
    shared_space = @cuDynamicSharedMem(Float64, max(n_vars, current_LP_length))


    while LP <= n_LPs
          
        # Check if we're supposed to skip this LP
        if skip_flag[LP]
            LP += grid_stride
            continue
        end

        # Initialize basic information for this LP
        ## information bookkeeping
        ## Only thread 1 of this block will store these information 
        ## the information below does not need to be stored per thread
        if idx==1
            cumulative_kkt_passes = 0.0

            # Initialize temporary Float64 values needed for calculation
            restart_primal_distance = 0.0
            restart_dual_distance = 0.0

            CI_primal_objective = 0.0 # Needed for optimality termination check
            CI_dual_objective = 0.0
            CI_l2_primal_residual = 0.0 # Needed for optimality termination check
            CI_l2_dual_residual = 0.0

            
            last_fixed_point_error = Inf
            restart_error = 0.0
            sum_restart_error = 0.0
            last_restart_error = 0.0

            best_primal_dual_residual_gap = Inf
            primal_dual_residual_gap = Inf

            restart_choice = RESTART_CHOICE_NO_RESTART
            cross_term = 0
            squared_delta_dual = 0
            squared_delta_primal = 0
            anchor_cross_term = 0
            anchor_squared_delta_primal = 0
            anchor_squared_delta_dual = 0
            initial_fixed_error = 0.0
            candidate_fixed_error = 0.0

            relative_dual_residual = 0.0
            relative_primal_residual = 0.0
            objective_gap = 0.0
            relative_objective_gap = 0.0

            ratio_infeas = Inf
        end

        # Other information that every thread needs (stored in the Registery memory)
        ## The values are always identical across threads — it is replicated rather than shared for performance purposes.
        
        total_iterations = Int32(0) # total iteration counter (T)
        inner_iterations = Int32(0)

        
        # Information that is much easier to save as static shared memory (L1 Block memory)
        do_restart = @cuStaticSharedMem(Bool, 1)
        primal_weight = @cuStaticSharedMem(Float64, 1)
        best_primal_weight = @cuStaticSharedMem(Float64, 1)
        numerical_error = @cuStaticSharedMem(Bool, 1)
        step_size = @cuStaticSharedMem(Float64, 1)

        if idx==1
            do_restart[1] = false
            primal_weight[1] = input_primal_weight[LP]
            best_primal_weight[1] = primal_weight[1]
            step_size[1] = input_step_size[LP]
            numerical_error[1] = false
        end

        
        # Set up the starting row for this LP (minus 1, so that the first
        # row to consider is `active_row + 1`)
        active_row = (LP-Int32(1)) * total_LP_length 


        # Set up current values to match inputs
        while idx <= n_vars
            current_primal_gradient[LP, idx] = scaled_objective_vector[LP, idx]
            idx += block_stride
        end
        # resetting idx to its original thread id after tweaking it to work on different variables above
        idx = threadIdx().x
        

        # Set up objective vector and RHS norms
        while idx <= n_vars
            shared_space[idx] = original_objective_vector[LP, idx]^2
            idx += block_stride
        end
        idx = threadIdx().x
        parallel_sum(shared_space, block_stride, var_stride, n_vars)

        if idx==1
            cache_l2_norm_primal_linear_objective = sqrt(shared_space[1])
        end

        while idx <= current_LP_length
            shared_space[idx] = original_right_hand_side[active_row + idx]^2
            idx += block_stride
        end
        idx = threadIdx().x
        parallel_sum(shared_space, block_stride, len_stride, current_LP_length)

        if idx==1
            cache_l2_norm_primal_right_hand_side = sqrt(shared_space[1])
        end

        # Begin the main loop
        while total_iterations <= iteration_limit

                ##########################################################################################
                #                       Taking Halpern Reflected PDHG steps                              #
                ##########################################################################################

                epoch_iterations = Int32(0) # counter for inner loop for computation purposes

                while epoch_iterations < termination_evaluation_frequency

                    local_primal_weight = unsafe_load(CUDA.pointer(primal_weight, 1))
                    local_step_size = unsafe_load(CUDA.pointer(step_size, 1))

                    # Step 1) Add one to the total iterations tracker and cumulative kkt passes
                    
                    if idx==1
                        cumulative_kkt_passes += Int32(1)
                    end
                        
                    # Step 2) Loop over nonzeros to update current_dual_product: 
                    while idx <= current_LP_length
                        shared_space[idx] = 0.0
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    sync_threads()

                    while idx <= nz_count
                        if active_constraint[active_row + nz_rows[idx]]
                            CUDA.atomic_add!(CUDA.pointer(shared_space, nz_cols[idx]), scaled_constraint_matrix[active_row + nz_rows[idx], nz_cols[idx]] * current_dual_solution[active_row + nz_rows[idx]]) # adds one to kkt passes
                        end
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    sync_threads()
                    while idx <= n_vars
                        current_dual_product[LP, idx] = shared_space[idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # Step 3) Calculate the next primal solution information                    

                    while idx <= n_vars
                        temp[LP, idx] = current_primal_solution[LP, idx] - (local_step_size/local_primal_weight) * (scaled_objective_vector[LP, idx] - current_dual_product[LP, idx])
                        pdhg_primal_solution[LP, idx] = max(scaled_variable_lower_bounds[LP, idx], min(scaled_variable_upper_bounds[LP, idx], temp[LP, idx]))
                        dual_slack[LP, idx] = (pdhg_primal_solution[LP, idx] - temp[LP, idx]) / (local_step_size/local_primal_weight)
                        reflected_primal_solution[LP, idx] = 2.0 * pdhg_primal_solution[LP, idx] - current_primal_solution[LP, idx]
                        current_primal_solution[LP, idx] = ((inner_iterations + 1)/ (inner_iterations + 2)) * (reflection_coefficient * reflected_primal_solution[LP, idx] + 
                                                            (1 - reflection_coefficient) * current_primal_solution[LP, idx]) + (1 - ((inner_iterations + 1)/ (inner_iterations + 2))) * initial_primal_solution[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # Step 4) Loop over nonzeros to update current_primal_product, based on the reflected primal solution (A(2x(k+1)-xk)) [reflected]
                    while idx <= current_LP_length
                        shared_space[idx] = 0.0
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    sync_threads()

                    while idx <= nz_count
                        if active_constraint[active_row + nz_rows[idx]]
                            CUDA.atomic_add!(CUDA.pointer(shared_space, nz_rows[idx]), scaled_constraint_matrix[active_row + nz_rows[idx], nz_cols[idx]] * reflected_primal_solution[LP, nz_cols[idx]]) # adds one to kkt passes
                        end
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    sync_threads()
                    while idx <= current_LP_length
                        current_primal_product[active_row + idx] = shared_space[idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x


                    # Step 5) Compute the next dual solution

                    while idx <= current_LP_length

                        temp_dual[active_row + idx] = (current_dual_solution[active_row + idx] / (local_primal_weight*local_step_size) - current_primal_product[active_row + idx])
                        pdhg_dual_solution[active_row + idx] = (temp_dual[active_row + idx] - min(-scaled_right_hand_side[active_row + idx], temp_dual[active_row + idx])) * (local_primal_weight*local_step_size)
                        reflected_dual_solution[active_row + idx] = 2.0 * pdhg_dual_solution[active_row + idx] - current_dual_solution[active_row + idx]
                        current_dual_solution[active_row + idx] =  ((inner_iterations + 1)/ (inner_iterations + 2)) * (reflection_coefficient * reflected_dual_solution[active_row + idx] + (1 - reflection_coefficient) * current_dual_solution[active_row + idx]) +
                                                                    (1.0 -  ((inner_iterations + 1)/ (inner_iterations + 2))) * initial_dual_solution[active_row + idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    sync_threads() 


                    if inner_iterations == Int32(0)

                        # Compute Δx = reflected_primal - pdhg_primal, and ||Δx||²
                        while idx <= n_vars
                            delta_primal[LP, idx] = reflected_primal_solution[LP, idx] - pdhg_primal_solution[LP, idx]
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        while idx <= n_vars 
                            
                            shared_space[idx] = delta_primal[LP, idx] ^ 2
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        parallel_sum(shared_space, block_stride, var_stride, n_vars)
                        if idx == 1
                            anchor_squared_delta_primal = shared_space[1]
                        end

                        # Compute Δy = reflected_dual - pdhg_dual, and ||Δy||²
                        while idx <= current_LP_length
                            delta_dual[active_row + idx] = reflected_dual_solution[active_row + idx] - pdhg_dual_solution[active_row + idx]
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        while idx <= current_LP_length
                            shared_space[idx] = delta_dual[active_row + idx] ^ 2
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                        if idx == 1
                            anchor_squared_delta_dual = shared_space[1]
                        end

                        # Compute A^T Δy (store in shared_space, which has n_vars width)
                        while idx <= n_vars
                            shared_space[idx] = 0.0
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        sync_threads()

                        while idx <= nz_count
                            if active_constraint[active_row + nz_rows[idx]]
                                CUDA.atomic_add!(CUDA.pointer(shared_space, nz_cols[idx]), scaled_constraint_matrix[active_row + nz_rows[idx], nz_cols[idx]] * delta_dual[active_row + nz_rows[idx]])
                            end
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        sync_threads()

                        # Compute cross_term = dot(A^T Δy, Δx)
                        while idx <= n_vars
                            shared_space[idx] = shared_space[idx] * delta_primal[LP, idx]
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        parallel_sum(shared_space, block_stride, var_stride, n_vars)
                        if idx == 1
                            anchor_cross_term = shared_space[1]
                        end

                        if unsafe_load(CUDA.pointer(do_restart, 1)) == true 
                            if idx == 1
                                initial_fixed_error = sqrt( (local_primal_weight) * anchor_squared_delta_primal + (1 / (local_primal_weight)) * anchor_squared_delta_dual + 2 * local_step_size * anchor_cross_term )
                                do_restart[1] = false
                            end
                        end

                    end

                    inner_iterations += Int32(1)
                    epoch_iterations += Int32(1)     
                    sync_threads()   
                end     

                if idx == 1
                    cumulative_kkt_passes += Int32(2)
                end

                while idx <= n_vars
                    current_primal_gradient[LP, idx] = scaled_objective_vector[LP, idx] - current_dual_product[LP, idx]
                    idx += block_stride
                end
                idx = threadIdx().x

                
                ##########################################################################################
                #                            Compute Fixed-point error                                   #
                ##########################################################################################
                sync_threads()
                
                while idx <= n_vars
                    shared_space[idx] = (reflected_primal_solution[LP, idx] - pdhg_primal_solution[LP, idx]) ^ 2
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx == 1
                    squared_delta_primal = shared_space[1]
                end
                while idx <= current_LP_length
                    shared_space[idx] = 0.0
                    idx += block_stride
                end
                idx = threadIdx().x
                sync_threads()

                while idx <= current_LP_length
                    shared_space[idx] = (reflected_dual_solution[active_row + idx] - pdhg_dual_solution[active_row + idx]) ^ 2
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                if idx == 1
                    squared_delta_dual = shared_space[1]
                end

                
                while idx <= current_LP_length
                    residual_delta_dual[active_row + idx] = reflected_dual_solution[active_row + idx] - pdhg_dual_solution[active_row + idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                sync_threads()

                while idx <= n_vars
                    shared_space[idx] = 0.0
                    idx += block_stride
                end
                idx = threadIdx().x
                sync_threads()

                while idx <= nz_count
                    if active_constraint[active_row + nz_rows[idx]] 
                        CUDA.atomic_add!(CUDA.pointer(shared_space, nz_cols[idx]), scaled_constraint_matrix[active_row + nz_rows[idx], nz_cols[idx]] * residual_delta_dual[active_row + nz_rows[idx]])
                    end
                    idx += block_stride
                end
                idx = threadIdx().x
                sync_threads()

                while idx <= n_vars
                    shared_space[idx] *= (reflected_primal_solution[LP, idx] - pdhg_primal_solution[LP, idx])
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx == 1
                    cross_term = shared_space[1]
                end

                if idx==1
                    candidate_fixed_error = sqrt(primal_weight[1] * squared_delta_primal + squared_delta_dual / (primal_weight[1]) + 2.0 * step_size[1] * cross_term)
                end

                
                ##########################################################################################
                #                            Compute Residuals                                           #
                ##########################################################################################
                sync_threads()
                
                # ++++++++++++++++++++++++++++++++++++++++++++++++ Updating Primal Product and Dual Product  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
                # Recompute primal product as G * pdhg_primal
                while idx <= current_LP_length
                    shared_space[idx] = 0.0
                    idx += block_stride
                end
                idx = threadIdx().x
                sync_threads()

                while idx <= nz_count
                    if active_constraint[active_row + nz_rows[idx]]
                        CUDA.atomic_add!(CUDA.pointer(shared_space, nz_rows[idx]), 
                            scaled_constraint_matrix[active_row + nz_rows[idx], nz_cols[idx]] * pdhg_primal_solution[LP, nz_cols[idx]])
                    end
                    idx += block_stride
                end
                idx = threadIdx().x
                sync_threads()
                while idx <= current_LP_length
                    residual_primal_product[active_row + idx] = shared_space[idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                # Dual Product as G^T * pdhg_dual
                while idx <= n_vars
                    shared_space[idx] = 0.0
                    idx += block_stride
                end
                idx = threadIdx().x
                sync_threads()
                
                # calculate residual dual product
                while idx <= nz_count
                    if active_constraint[active_row + nz_rows[idx]]
                        CUDA.atomic_add!(CUDA.pointer(shared_space, nz_cols[idx]), scaled_constraint_matrix[active_row + nz_rows[idx], nz_cols[idx]] * pdhg_dual_solution[active_row + nz_rows[idx]])
                    end
                    idx += block_stride
                end
            
                idx = threadIdx().x
                sync_threads()
                while idx <= n_vars
                    residual_dual_product[LP, idx] = shared_space[idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                
                sync_threads()

                while idx <= current_LP_length
                    clamped_val = max(residual_primal_product[active_row + idx], scaled_right_hand_side[active_row + idx])
                    primal_residual[active_row + idx] = (residual_primal_product[active_row + idx] - clamped_val) * constraint_rescaling[active_row + idx]
                    primal_slack[active_row + idx] = max(pdhg_dual_solution[active_row + idx], 0.0) * scaled_right_hand_side[active_row + idx]
                    idx += block_stride
                end
                idx = threadIdx().x

                # calculting dual_residual
                while idx <= n_vars
                    dual_residual[LP, idx] = (scaled_objective_vector[LP, idx] - residual_dual_product[LP, idx] - dual_slack[LP, idx]) * variable_rescaling[LP, idx]
                    idx += block_stride
                end
                idx = threadIdx().x

                # cuPDLPx has L ∞ as default, I ignored it and used L2 norm only 

                while idx <= current_LP_length
                    shared_space[idx] = primal_residual[active_row + idx] ^ 2
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, len_stride, current_LP_length) # this sum all squared terms
                # taking its sqrt and assiging as primal residual norm 
                if idx == 1
                    CI_l2_primal_residual = sqrt(shared_space[1])
                end

                #TODO: if bound objective rescaling is set to true we need to uncomment the line below to rescale back primal residual norm 
                # CI_l2_primal_residual /= constraint_bound_norm

                while idx <= n_vars
                    shared_space[idx] = dual_residual[LP, idx] ^ 2
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars) # this sum all squared terms
                # taking its sqrt and assiging as primal residual norm 
                if idx == 1
                    CI_l2_dual_residual = sqrt(shared_space[1])
                end

                #TODO: if bound objective rescaling is set to true we need to uncomment the line below to rescale back primal residual norm 
                # CI_l2_dual_residual /= objective_vector_norm

                # Compute the primal objective 
                while idx <= n_vars
                    shared_space[idx] = scaled_objective_vector[LP, idx] * pdhg_primal_solution[LP, idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    CI_primal_objective = original_objective_constant[LP] + shared_space[1]
                end

                # Compute the dual objective
                while idx <= n_vars
                    shared_space[idx] = dual_slack[LP, idx] * pdhg_primal_solution[LP, idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx == 1
                    CI_dual_objective = shared_space[1]
                end

                while idx <= current_LP_length
                    shared_space[idx] = primal_slack[active_row + idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                if idx == 1
                    CI_dual_objective += shared_space[1] + original_objective_constant[LP]
                end

                if idx == 1
                    # #CUDA.@cuprint("constraint bound norm = $constraint_bound_norm")
                    constraint_bound_norm = cache_l2_norm_primal_right_hand_side
                    relative_primal_residual = CI_l2_primal_residual / Float64(1.0 + constraint_bound_norm)
                    objective_vector_norm = cache_l2_norm_primal_linear_objective
                    relative_dual_residual = CI_l2_dual_residual / Float64(1.0  + objective_vector_norm)
                    objective_gap = abs(CI_primal_objective - CI_dual_objective)
                    relative_objective_gap = objective_gap / (1.0 + abs(CI_primal_objective) + abs(CI_dual_objective))
                end
                
                # +++++++++++++++++++++++++++++++++ Infeasibility Detection +++++++++++++++++++++++++++++++++++++++++ #
                #TODO: finish feasiblity detection, cuPDLPx does not do it, but we need it in subproblems
                
                total_iterations += Int32(termination_evaluation_frequency)

                ##########################################################################################
                #                            Check Termination Criteria                                  #
                ##########################################################################################
                sync_threads()

                # Check optimality criteria
                if idx==1

                    # Check if the dual objective value is above the B&B global upper bound and we
                    # satisfy the tolerance for dual feasibility
                    # (if we're past the first 10 iterations)
                    if (CI_dual_objective > global_upper_bound + abs_tol) &&
                        (CI_l2_dual_residual < abs_tol + rel_tol*cache_l2_norm_primal_linear_objective) 
                        termination_reason[LP] = TERMINATION_REASON_GLOBAL_UPPER_BOUND_HIT
                    end

                    # If we want to skip harder-than-average problems, and we've already solved at least
                    # 100 problems, check if the current number of iterations is over two times the running
                    # average
                    if skip_hard_problems && (unsafe_load(CUDA.pointer(global_counter, 1)) > 100)
                        if total_iterations > 2*(unsafe_load(CUDA.pointer(iteration_counter, 1))/unsafe_load(CUDA.pointer(global_counter, 1)))
                            termination_reason[LP] = TERMINATION_REASON_IMPATIENCE
                        end
                    end

                    # Check iteration limit
                    if total_iterations >= iteration_limit
                        
                        termination_reason[LP] = TERMINATION_REASON_ITERATION_LIMIT
                    end

                    # Check KKT matrix pass limit
                    if cumulative_kkt_passes >= kkt_matrix_pass_limit
                        termination_reason[LP] = TERMINATION_REASON_KKT_MATRIX_PASS_LIMIT
                    end

                    # Check for numerical errors
                    if numerical_error[1]
                        termination_reason[LP] = TERMINATION_REASON_NUMERICAL_ERROR
                    end

                    # Check if we're within the tolerances for primal and dual infeasibility, and that there's
                    # a sufficiently small gap between the primal and dual objective values.

                    if (relative_dual_residual < rel_tol) &&
                        (relative_primal_residual < rel_tol) &&
                        (relative_objective_gap < rel_tol)

                        termination_reason[LP] = TERMINATION_REASON_OPTIMAL
                    end
                    
                    # TODO: add checking on infeasibility after convergence failure
                end
                

                # If the LP is finished, update the solutions and objective(s), and
                # then break out of the iteration loop and move on to the next LP
                sync_threads()

                reason = unsafe_load(CUDA.pointer(termination_reason, LP))
                if reason == TERMINATION_REASON_OPTIMAL
                    while idx <= n_vars 
                        solutions[LP, idx] = current_primal_solution[LP, idx] / variable_rescaling[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    if idx==1
                        if return_code==Int32(1)
                            objectives[LP] = CI_primal_objective
                        elseif return_code==Int32(2)
                            objectives[LP] = CI_dual_objective
                        else
                            objectives[LP, Int32(1)] = CI_primal_objective
                            objectives[LP, Int32(2)] = CI_dual_objective
                        end
                        iterations[LP] = total_iterations
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(total_iterations))
                    end
                    
                    LP += grid_stride
                    break
                elseif (reason == TERMINATION_REASON_PRIMAL_INFEASIBLE) || 
                    (reason == TERMINATION_REASON_DUAL_INFEASIBLE)
                    while idx <= n_vars 
                        solutions[LP, idx] = current_primal_solution[LP, idx] / variable_rescaling[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    if idx==1
                        if return_code != Int32(3)
                            objectives[LP] = Inf
                        else
                            objectives[LP, Int32(1)] = Inf
                            objectives[LP, Int32(2)] = Inf
                        end
                        iterations[LP] = total_iterations
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(total_iterations))
                    end
                    LP += grid_stride
                    break
                elseif (reason == TERMINATION_REASON_ITERATION_LIMIT) || 
                    (reason == TERMINATION_REASON_KKT_MATRIX_PASS_LIMIT) ||
                    (reason == TERMINATION_REASON_NUMERICAL_ERROR)
                    while idx <= n_vars 
                        solutions[LP, idx] = current_primal_solution[LP, idx] / variable_rescaling[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    if idx==1
                        if return_code != Int32(3)
                            objectives[LP] = -Inf
                        else
                            objectives[LP, Int32(1)] = -Inf
                            objectives[LP, Int32(2)] = -Inf
                        end
                        iterations[LP] = total_iterations
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(total_iterations))
                    end

                    LP += grid_stride
                    break
                elseif (reason == TERMINATION_REASON_GLOBAL_UPPER_BOUND_HIT)
                    while idx <= n_vars
                        solutions[LP, idx] = current_primal_solution[LP, idx] / variable_rescaling[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    if idx==1
                        # Special case: This checks for a feasible dual objective value
                        # above the global upper bound in a B&B algorithm. The primal
                        # objective value may not be meaningful/feasible.
                        if return_code==Int32(1)
                            objectives[LP] = CI_primal_objective
                        elseif return_code==Int32(2)
                            objectives[LP] = CI_dual_objective
                        else
                            objectives[LP, Int32(1)] = CI_primal_objective
                            objectives[LP, Int32(2)] = CI_dual_objective
                        end
                        iterations[LP] = total_iterations
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(total_iterations))
                    end
                    LP += grid_stride
                    break
                elseif (reason == TERMINATION_REASON_IMPATIENCE)
                    while idx < n_vars
                        solutions[LP, idx] = current_primal_solution[LP, idx+1] / variable_rescaling[LP, idx+1]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    if idx==1
                        if return_code != Int32(3)
                            objectives[LP] = -Inf
                        else
                            objectives[LP, Int32(1)] = -Inf
                            objectives[LP, Int32(2)] = -Inf
                        end
                        iterations[LP] = total_iterations
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(total_iterations))
                    end
                    LP += grid_stride
                    break
                end

                ##########################################################################################
                #                            Check Restart Conditions                                    #
                ##########################################################################################
                sync_threads()
                # Check if we need to do a restart
                if idx == 1
                    if total_iterations == termination_evaluation_frequency

                        do_restart[1] = true 

                    elseif total_iterations > termination_evaluation_frequency

                        if candidate_fixed_error < sufficient_reduction_for_restart * initial_fixed_error

                            do_restart[1] = true
                        end

                        if candidate_fixed_error < necessary_reduction_for_restart * initial_fixed_error
                            if candidate_fixed_error > last_fixed_point_error

                                do_restart[1] = true
                            end
                        end
                        

                        if inner_iterations >= artificial_ratio_for_restart * (total_iterations)
                            do_restart[1] = true
                        end
                    end

                    last_fixed_point_error = candidate_fixed_error
                    
                end

                ##########################################################################################
                #                            Perform Restart if do_restart = true                        #
                ##########################################################################################
                sync_threads()
                
                if unsafe_load(CUDA.pointer(do_restart, 1)) == true 
                    
                    # Primal distance
                    while idx <= n_vars 
                        shared_space[idx] = (pdhg_primal_solution[LP, idx] - initial_primal_solution[LP, idx]) ^ 2
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    parallel_sum(shared_space, block_stride, var_stride, n_vars)

                    if idx == 1
                        restart_primal_distance = sqrt(shared_space[1])
                    end
                    
                    # Dual distance
                    while idx <= current_LP_length
                        shared_space[idx] = (pdhg_dual_solution[active_row + idx] - initial_dual_solution[active_row + idx]) ^ 2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, len_stride, current_LP_length)

                    if idx == 1
                        restart_dual_distance = sqrt(shared_space[1])
                        ratio_infeas = relative_dual_residual / relative_primal_residual
                    end
                    
                    # initialize the primal anchor 
                    while idx <= n_vars
                        initial_primal_solution[LP, idx] = pdhg_primal_solution[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # initialize the dual anchor
                    while idx <= current_LP_length
                        initial_dual_solution[active_row + idx] = pdhg_dual_solution[active_row + idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # updating the current_primal_solution and current_dual_solution 
                    while idx <= n_vars 
                        current_primal_solution[LP, idx] = pdhg_primal_solution[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # initialize the dual anchor
                    while idx <= current_LP_length
                        current_dual_solution[active_row + idx] = pdhg_dual_solution[active_row + idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # resetting step itertation back to 0 
                    inner_iterations = Int32(0)
                    
                    # we do not restart to average in Halpern Reflection Scheme
                    if idx==1
                        last_fixed_point_error = Inf
                        restart_choice = RESTART_CHOICE_LAST_ITERATE_RESET
                    end
                else
                    if idx==1
                        restart_choice = RESTART_CHOICE_NO_RESTART
                    end
                end
                


                ###################################################################
                ##### Compute a New Primal Weight (if a restart was used) 
                ###################################################################
                sync_threads()
                if idx==1
                    if restart_choice == RESTART_CHOICE_LAST_ITERATE_RESET

                        if restart_primal_distance > 1e-16 && restart_dual_distance > 1e-16 && restart_primal_distance < 1e12 && restart_dual_distance < 1e12 && ratio_infeas > 1e-8 && ratio_infeas < 1e8
                            restart_error = log(restart_dual_distance) - log(restart_primal_distance) - log(primal_weight[1])
                            sum_restart_error *= i_smooth
                            sum_restart_error += restart_error

                            
                            primal_weight[1] *= exp(pid_KP * restart_error + pid_KI * sum_restart_error + pid_KD * (restart_error - last_restart_error))
                                                        
                            # update last restart error 
                            last_restart_error = restart_error
                        else 
                            primal_weight[1] = best_primal_weight[1]
                                                        
                            sum_restart_error = 0.0
                            last_restart_error = 0.0
                        end

                        # Update best primal weight if residual gap improved
                        primal_dual_residual_gap = abs(log10(relative_dual_residual / relative_primal_residual))

                        if primal_dual_residual_gap < best_primal_dual_residual_gap
                            best_primal_dual_residual_gap = primal_dual_residual_gap
                            best_primal_weight[1] = primal_weight[1]
                        end
                    end
                end

                sync_threads()
        end
        
    end

    return nothing
end

# A quick function to calculate a parallel sum over the first max_len elements in the shared space
function parallel_sum(shared, block_stride, reduction_stride, maxlen)
    sync_threads()
    idx = threadIdx().x
    while reduction_stride > 0
        while idx <= reduction_stride && idx + reduction_stride <= maxlen
            shared[idx] += shared[idx + reduction_stride]
            idx += block_stride
        end
        idx = threadIdx().x
        sync_threads()
        reduction_stride >>>= 1
    end
    return nothing
end

# A quick function to calculate a parallel max over the first max_len elements in the shared space
function parallel_max(shared, block_stride, reduction_stride, maxlen)
    sync_threads()
    idx = threadIdx().x
    while reduction_stride > 0
        while idx <= reduction_stride && idx + reduction_stride <= maxlen
            shared[idx] = max(shared[idx], shared[idx + reduction_stride])
            idx += block_stride
        end
        idx = threadIdx().x
        sync_threads()
        reduction_stride >>>= 1
    end
    return nothing
end