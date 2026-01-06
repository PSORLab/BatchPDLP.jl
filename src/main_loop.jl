
function main_loop_kernel(
    solutions,                        # Holds solutions to LPs
    objectives,                       # Holds dual objectives
    original_variable_lower_bounds,   # Used for getting unscaled convergence info
    original_variable_upper_bounds,   # ""
    original_right_hand_side,         # ""
    original_objective_vector,        # ""
    original_objective_constant,      # ""
    scaled_variable_lower_bounds,     # Lower bounds after initial rescaling
    scaled_variable_upper_bounds,     # Upper bounds after initial rescaling
    scaled_constraint_matrix,         # Constraint matrix after initial rescaling
    scaled_right_hand_side,           # Right-hand side after initial rescaling
    scaled_objective_vector,          # Objective vector after initial rescaling
    scaled_objective_constant,        # Objective constant after initial rescaling
    nz_count,                         # Total number of nonzeros in the constraint matrix
    nz_rows,                          # Rows containing nonzeros in the constraint matrix
    nz_cols,                          # Columns containing nonzeros in the constraint matrix
    active_constraint,                # Flag to indicate whether a constraint is in use.
    variable_rescaling,               # Total value used for variable rescaling
    constraint_rescaling,             # Total value used for constraint rescaling
    # input_current_primal_solution,    # Usually 0s, but used for hot-starts (if desired in the future)
    # input_current_dual_solution,      # Usually 0s, but used for hot-starts (if desired in the future)
    # input_current_primal_product,     # Usually 0s, but used for hot-starts (if desired in the future)
    # input_current_dual_product,       # Usually 0s, but used for hot-starts (if desired in the future)
    last_restart_primal_solution,     # Solution state the last time a restart occurred
    last_restart_primal_gradient,     # ""
    last_restart_dual_solution,       # ""
    last_restart_primal_product,      # ""
    current_primal_solution,          # Currently active solution state
    current_dual_solution,            # ""
    current_dual_product,             # ""
    current_primal_product,           # ""
    buffer_primal_gradient,           # Generally the avg_primal_gradient, but used for some other purposes
    avg_primal_solution,              # Running average solution state information, kept since last restart
    avg_primal_gradient,              # ""
    avg_dual_solution,                # ""
    avg_primal_product,               # ""
    sum_primal_solutions,             # Running sum solution state information, used to calculate averages, since last restart
    sum_dual_solutions,               # ""
    sum_primal_product,               # ""
    sum_dual_product,                 # ""
    original_primal_solution,         # Unscaled solution state information, used to calculate actual infeasibility
    original_primal_gradient,         # ""
    original_dual_solution,           # ""
    original_primal_product,          # ""
    buffer_kkt_primal_solution,       # Intermediate storage for calculating infeasibility
    buffer_kkt_primal_product,        # ""
    buffer_kkt_lower_variable_violation, # Infeasible variable/constraint violation storage
    buffer_kkt_upper_variable_violation, # ""
    buffer_kkt_reduced_costs,         # Storage used for calculating, e.g., dual objective values
    delta_primal,                     # Value used in calculating PDLP steps
    delta_primal_product,             # ""
    delta_dual,                       # ""
    input_primal_weight,              # The original primal weight (could be moved internally if needed)
    input_step_size,                  # The original step size (could be moved internally if needed)
    termination_reason,               # Field to let the user know why the LP terminated
    current_LP_length,                # The current length of each individual LP (i.e., number of constraints)
    total_LP_length,                  # The total possible size of each individual LP (i.e., max allowed number of constraints)
    n_LPs,                            # The current number of LPs being solved (fewer than the max may be used during B&B)
    n_vars,                           # The number of variables in the problem
    iteration_limit,                  # The maximum number of PDLP steps allowed before termination
    kkt_matrix_pass_limit,            # The maximum number of KKT matrix passes (generally unused)
    termination_evaluation_frequency, # Number of PDLP steps to take before checking termination criteria (default: 64)
    necessary_reduction_for_restart,  # KKT error decrease factor needed for a restart if no progress is being made (default: 0.8)
    sufficient_reduction_for_restart, # KKT error decrease factor that would trigger a restart immediately (default: 0.2)
    extrapolation_coefficient,        # Value used in PDLP steps (default: 1.0)
    reduction_exponent,               # Step size adjustment based on step iteration count (default: 0.3)
    growth_exponent,                  # Step size adjustment based on step iteration count (default: 0.6)
    abs_tol,                          # Absolute tolerance for termination
    rel_tol,                          # Relative tolerance for termination
    eps_primal_infeasible,            # Primal infeasibility tolerance
    eps_dual_infeasible,              # Dual infeasibility tolerance
    return_code,                      # Indicator for returning primal obj (1), dual obj (2), or both (3)
    global_upper_bound,               # Information about the B&B upper bound (PDLP terminates if a dual feasible solution is above this value)
    skip_hard_problems,               # Flag to skip problems with too many iterations
    global_counter,                   # Count of LPs being solved
    iteration_counter,                # Count of total iterations for solved LPs
    skip_flag,                        # Flag to completely skip an individual LP
    iterations,                       # The final number of iterations needed for each LP
    )
    # In this kernel, assume that the number of blocks is equal to the number
    # of LPs, so that all threads are working on one LP
    idx = threadIdx().x
    LP = blockIdx().x
    block_stride = blockDim().x
    grid_stride = gridDim().x

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
        if idx==1
            cumulative_kkt_passes = 0.5

            # Initialize temporary Float64 values needed for calculation
            last_restart_primal_distance_moved = 0.0
            last_restart_dual_distance_moved = 0.0
            buffer_kkt_dual_objective = 0.0
            buffer_kkt_dual_res_inf = 0.0
            CI_primal_objective = 0.0 # Needed for optimality termination check
            CI_dual_objective = 0.0
            CI_l2_primal_residual = 0.0 # Needed for optimality termination check
            CI_l2_dual_residual = 0.0
            last_restart_length = 1.0
            last_reduction_ratio = 1.0
            avg_kkt_residual = 0.0
            current_kkt_residual = 0.0
            candidate_kkt_residual = 0.0
            last_kkt_residual = 0.0
            kkt_reduction_ratio = 0.0
            restart_choice = RESTART_CHOICE_NO_RESTART
        end

        # Other information that every thread needs
        step_iterations = Int32(0)
        iteration = Int32(1)
        sum_primal_solution_weights = 0.0
        sum_dual_solution_weights = 0.0
        primal_ray_norm = 0.0
        sum_primal_solutions_count = Int32(0)
        sum_dual_solutions_count = Int32(0)

        # Information that is much easier to save as static shared memory
        do_restart = @cuStaticSharedMem(Bool, 1)
        reset_to_average = @cuStaticSharedMem(Bool, 1)
        primal_weight = @cuStaticSharedMem(Float64, 1)
        step_size = @cuStaticSharedMem(Float64, 1)
        step_size_limit = @cuStaticSharedMem(Float64, 1)
        numerical_error = @cuStaticSharedMem(Bool, 1)
        if idx==1
            do_restart[1] = false
            reset_to_average[1] = false
            primal_weight[1] = input_primal_weight[LP]
            step_size[1] = input_step_size[LP]
            step_size_limit[1] = 0.0
            numerical_error[1] = false
        end
        
        # Set up the starting row for this LP (minus 1, so that the first
        # row to consider is `active_row + 1`)
        active_row = (LP-Int32(1))*total_LP_length 

        # Set up current values to match inputs
        while idx <= n_vars
            buffer_primal_gradient[LP, idx] = scaled_objective_vector[LP, idx]
            idx += block_stride
        end
        idx = threadIdx().x
        
        # Set up last restart info
        while idx <= n_vars
            last_restart_primal_solution[LP, idx] = current_primal_solution[LP, idx]
            last_restart_primal_gradient[LP, idx] = buffer_primal_gradient[LP, idx]
            idx += block_stride
        end
        idx = threadIdx().x
        while idx <= current_LP_length
            last_restart_dual_solution[active_row + idx] = current_dual_solution[active_row + idx]
            last_restart_primal_product[active_row + idx] = current_primal_product[active_row + idx]
            idx += block_stride
        end
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
        while true
            # Check termination occasionally
            if iteration <= 10 || mod(iteration - 1, termination_evaluation_frequency) == 0 || (iteration >= iteration_limit)
                # Add to the cumulative kkt pass count
                if idx==1
                    cumulative_kkt_passes += 2.0
                end

                # Update the average primal/dual solutions, average primal gradient,
                # and average primal product
                sync_threads()
                if iszero(sum_primal_solutions_count) || iszero(sum_dual_solutions_count) || unsafe_load(CUDA.pointer(numerical_error, 1))
                    while idx <= n_vars
                        avg_primal_solution[LP, idx] = current_primal_solution[LP, idx]
                        avg_primal_gradient[LP, idx] = buffer_primal_gradient[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    while idx <= current_LP_length
                        avg_dual_solution[active_row + idx] = current_dual_solution[active_row + idx]
                        avg_primal_product[active_row + idx] = current_primal_product[active_row + idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                else
                    while idx <= n_vars
                        avg_primal_solution[LP, idx] = sum_primal_solutions[LP, idx] / sum_primal_solution_weights
                        avg_primal_gradient[LP, idx] = (-sum_dual_product[LP, idx] / sum_dual_solution_weights) + scaled_objective_vector[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    while idx <= current_LP_length
                        avg_dual_solution[active_row + idx] = sum_dual_solutions[active_row + idx] / sum_dual_solution_weights
                        avg_primal_product[active_row + idx] = sum_primal_product[active_row + idx] / sum_primal_solution_weights
                        idx += block_stride
                    end
                    idx = threadIdx().x
                end


                ###################################################################
                ##### Begin "evaluate_unscaled_iteraiton_stats"
                ###################################################################
                while idx <= n_vars
                    original_primal_solution[LP, idx] = avg_primal_solution[LP, idx] / variable_rescaling[LP, idx]
                    original_primal_gradient[LP, idx] = avg_primal_gradient[LP, idx] * variable_rescaling[LP, idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                while idx <= current_LP_length
                    original_dual_solution[active_row + idx] = avg_dual_solution[active_row + idx] / constraint_rescaling[active_row + idx]
                    original_primal_product[active_row + idx] = avg_primal_product[active_row + idx] * constraint_rescaling[active_row + idx]
                    idx += block_stride
                end
                idx = threadIdx().x


                ###################################################################
                ##### Compute Convergence Information
                ###################################################################
                # Compute primal variable/constraint violations
                while idx <= n_vars
                    buffer_kkt_lower_variable_violation[LP, idx] = max(original_variable_lower_bounds[LP, idx] - original_primal_solution[LP, idx], 0.0)
                    buffer_kkt_upper_variable_violation[LP, idx] = max(original_primal_solution[LP, idx] - original_variable_upper_bounds[LP, idx], 0.0)
                    idx += block_stride
                end
                idx = threadIdx().x

                # Compute the primal objective and residual
                while idx <= n_vars
                    shared_space[idx] = original_objective_vector[LP, idx] * original_primal_solution[LP, idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    CI_primal_objective = original_objective_constant[LP] + shared_space[1]
                end

                # CI_l_inf_primal_residual
                while idx <= n_vars
                    shared_space[idx] = abs(buffer_kkt_lower_variable_violation[LP, idx])^2 +
                                        abs(buffer_kkt_upper_variable_violation[LP, idx])^2
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    CI_l2_primal_residual = shared_space[1]
                end
                while idx <= current_LP_length
                    shared_space[idx] = abs(max(original_right_hand_side[active_row + idx] - original_primal_product[active_row + idx], 0.0))^2
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    CI_l2_primal_residual = sqrt(CI_l2_primal_residual + shared_space[1])
                end
                
                # Compute dual stats
                while idx <= n_vars
                    buffer_kkt_reduced_costs[LP, idx] = max(original_primal_gradient[LP, idx], 0.0) * isfinite(original_variable_lower_bounds[LP, idx]) + 
                                                        min(original_primal_gradient[LP, idx], 0.0) * isfinite(original_variable_upper_bounds[LP, idx])
                    idx += block_stride
                end
                idx = threadIdx().x

                # Calculate the dual objective
                while idx <= n_vars
                    if buffer_kkt_reduced_costs[LP, idx] > 0
                        shared_space[idx] = original_variable_lower_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                    elseif buffer_kkt_reduced_costs[LP, idx] < 0.0
                        shared_space[idx] = original_variable_upper_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                    else
                        shared_space[idx] = 0.0
                    end
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    CI_dual_objective = original_objective_constant[LP] + shared_space[1]
                end
                while idx <= current_LP_length
                    shared_space[idx] = original_right_hand_side[active_row + idx] * original_dual_solution[active_row + idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    CI_dual_objective += shared_space[1]
                end

                # Calculate the l2 dual residual
                while idx <= n_vars
                    shared_space[idx] = abs(original_primal_gradient[LP, idx] - buffer_kkt_reduced_costs[LP, idx])^2
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    CI_l2_dual_residual = shared_space[1]
                end
                while idx <= current_LP_length
                    shared_space[idx] = abs(max(-original_dual_solution[active_row + idx], 0.0))^2
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    CI_l2_dual_residual = sqrt(CI_l2_dual_residual + shared_space[1])
                end


                ###################################################################
                ##### Compute Infeasibility Information
                ###################################################################
                while idx <= n_vars
                    shared_space[idx] = abs(original_primal_solution[LP, idx])
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, var_stride, n_vars)
                # Everyone needs to know the primal ray norm
                primal_ray_norm = shared_space[1]

                # Calculate infeasibility primal solution and primal product
                if !iszero(primal_ray_norm)
                    while idx <= n_vars
                        buffer_kkt_primal_solution[LP, idx] = original_primal_solution[LP, idx] / primal_ray_norm
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    while idx <= current_LP_length
                        buffer_kkt_primal_product[active_row + idx] = original_primal_product[active_row + idx] / primal_ray_norm
                        idx += block_stride
                    end
                    idx = threadIdx().x
                else
                    while idx <= n_vars
                        buffer_kkt_primal_solution[LP, idx] = original_primal_solution[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    while idx <= current_LP_length
                        buffer_kkt_primal_product[active_row + idx] = original_primal_product[active_row + idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                end

                # Compute infeasible variable/constraint violations
                while idx <= n_vars
                    if isfinite(original_variable_lower_bounds[LP, idx])
                        buffer_kkt_lower_variable_violation[LP, idx] = max(-buffer_kkt_primal_solution[LP, idx], 0.0)
                    else
                        buffer_kkt_lower_variable_violation[LP, idx] = 0.0
                    end
                    if isfinite(original_variable_upper_bounds[LP, idx])
                        buffer_kkt_upper_variable_violation[LP, idx] = max(buffer_kkt_primal_solution[LP, idx], 0.0)
                    else
                        buffer_kkt_upper_variable_violation[LP, idx] = 0.0
                    end
                    idx += block_stride
                end
                idx = threadIdx().x

                # Calculate the max primal ray infeasibility
                while idx <= n_vars
                    shared_space[idx] = max(abs(buffer_kkt_lower_variable_violation[LP, idx]), 
                                            abs(buffer_kkt_upper_variable_violation[LP, idx]))
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    II_max_primal_ray_infeasibility = shared_space[1]
                end
                while idx <= current_LP_length
                    shared_space[idx] = abs(max(-buffer_kkt_primal_product[active_row + idx], 0.0))
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    II_max_primal_ray_infeasibility = max(shared_space[1], II_max_primal_ray_infeasibility)
                end

                # Calculate the primal ray linear objective
                while idx <= n_vars
                    shared_space[idx] = original_objective_vector[LP, idx] * buffer_kkt_primal_solution[LP, idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    II_primal_ray_linear_objective = shared_space[1]
                end

                # Compute reduced costs and reduced costs violation
                while idx <= n_vars
                    buffer_kkt_reduced_costs[LP, idx] = max(original_primal_gradient[LP, idx] - original_objective_vector[LP, idx], 0.0) * isfinite(original_variable_lower_bounds[LP, idx]) + 
                                                        min(original_primal_gradient[LP, idx] - original_objective_vector[LP, idx], 0.0) * isfinite(original_variable_upper_bounds[LP, idx])
                    idx += block_stride
                end
                idx = threadIdx().x

                # Compute the dual residual
                while idx <= n_vars
                    shared_space[idx] = abs(original_primal_gradient[LP, idx] - original_objective_vector[LP, idx] - buffer_kkt_reduced_costs[LP, idx])
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    buffer_kkt_dual_res_inf = shared_space[1]
                end
                while idx <= current_LP_length
                    shared_space[idx] = abs(max(-original_dual_solution[active_row + idx], 0.0))
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    buffer_kkt_dual_res_inf = max(shared_space[1], buffer_kkt_dual_res_inf)
                end

                # Compute the dual objective
                while idx <= n_vars
                    if buffer_kkt_reduced_costs[LP, idx] > 0.0
                        shared_space[idx] = original_variable_lower_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                    elseif buffer_kkt_reduced_costs[LP, idx] < 0.0
                        shared_space[idx] = original_variable_upper_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                    else
                        shared_space[idx] = 0.0
                    end
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    buffer_kkt_dual_objective = shared_space[1] + original_objective_constant[LP]
                end
                while idx <= current_LP_length
                    shared_space[idx] = original_right_hand_side[active_row + idx] * original_dual_solution[active_row + idx]
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    buffer_kkt_dual_objective += shared_space[1]
                end

                # Compute infeasibility information using a scaling factor
                while idx <= n_vars
                    shared_space[idx] = abs(buffer_kkt_reduced_costs[LP, idx])
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, var_stride, n_vars)
                if idx==1
                    scaling_factor = shared_space[1]
                end
                while idx <= current_LP_length
                    shared_space[idx] = abs(original_dual_solution[active_row + idx])
                    idx += block_stride
                end
                idx = threadIdx().x
                parallel_max(shared_space, block_stride, len_stride, current_LP_length)
                if idx==1
                    scaling_factor = max(shared_space[1], scaling_factor)
                    if scaling_factor==0.0
                        II_max_dual_ray_infeasibility = 0.0
                        II_dual_ray_objective = 0.0
                    else
                        II_max_dual_ray_infeasibility = buffer_kkt_dual_res_inf / scaling_factor
                        II_dual_ray_objective = buffer_kkt_dual_objective / scaling_factor
                    end
                end

                
                ###################################################################
                ##### Set up Primal and Dual Norm Params
                ###################################################################
                if idx==1
                    primal_norm_params = 1 / step_size[1] * primal_weight[1]
                    dual_norm_params = 1 / step_size[1] / primal_weight[1]
                end


                ###################################################################
                ##### Check Termination Criteria
                ###################################################################
                # Check optimality criteria
                if idx==1
                    # Check if the dual objective value is above the B&B global upper bound and we
                    # satisfy the tolerance for dual feasibility
                    # (if we're past the first 10 iterations)
                    if (CI_dual_objective > global_upper_bound + abs_tol) &&
                        (CI_l2_dual_residual < abs_tol + rel_tol*cache_l2_norm_primal_linear_objective) &&
                        (iteration >= 10)
                        termination_reason[LP] = TERMINATION_REASON_GLOBAL_UPPER_BOUND_HIT
                    end

                    # If we want to skip harder-than-average problems, and we've already solved at least
                    # 100 problems, check if the current number of iterations is over two times the running
                    # average
                    if skip_hard_problems && (unsafe_load(CUDA.pointer(global_counter, 1)) > 100)
                        if iteration > 2*(unsafe_load(CUDA.pointer(iteration_counter, 1))/unsafe_load(CUDA.pointer(global_counter, 1)))
                            termination_reason[LP] = TERMINATION_REASON_IMPATIENCE
                        end
                    end

                    # Check if we're within the tolerances for primal and dual infeasibility, and that there's
                    # a sufficiently small gap between the primal and dual objective values.
                    if (CI_l2_dual_residual < abs_tol + rel_tol*cache_l2_norm_primal_linear_objective) &&
                        (CI_l2_primal_residual < abs_tol + rel_tol*cache_l2_norm_primal_right_hand_side) &&
                        (abs(CI_primal_objective - CI_dual_objective) < abs_tol + rel_tol*(abs(CI_primal_objective)+abs(CI_dual_objective)))
                        termination_reason[LP] = TERMINATION_REASON_OPTIMAL
                    end
                    
                    # Check primal infeasibility (if we're past the first 10 iterations)
                    if (II_dual_ray_objective > 0.0) &&
                        ((II_max_dual_ray_infeasibility / II_dual_ray_objective) <= eps_primal_infeasible) &&
                        (iteration >= 10)
                        termination_reason[LP] = TERMINATION_REASON_PRIMAL_INFEASIBLE
                    end

                    # Check dual infeasibility (if we're past the first 10 iterations)
                    if (II_primal_ray_linear_objective < 0.0) && 
                        ((II_max_primal_ray_infeasibility / (-II_primal_ray_linear_objective)) <= eps_dual_infeasible) &&
                        (iteration >= 10)
                        termination_reason[LP] = TERMINATION_REASON_DUAL_INFEASIBLE
                    end

                    # Check iteration limit
                    if iteration >= iteration_limit
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
                end


                ###################################################################
                ##### Update Solutions
                ###################################################################

                # If the LP is finished, update the solutions and objective(s), and
                # then break out of the iteration loop and move on to the next LP
                sync_threads()
                reason = unsafe_load(CUDA.pointer(termination_reason, LP))
                if reason == TERMINATION_REASON_OPTIMAL
                    while idx < n_vars # Intentionally skipping the epigraph variable
                        solutions[LP, idx] = avg_primal_solution[LP, idx+1] / variable_rescaling[LP, idx+1]
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
                        iterations[LP] = iteration
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(iteration))
                    end
                    LP += grid_stride
                    break
                elseif (reason == TERMINATION_REASON_PRIMAL_INFEASIBLE) || 
                    (reason == TERMINATION_REASON_DUAL_INFEASIBLE)
                    while idx < n_vars # Intentionally skipping the epigraph variable
                        solutions[LP, idx] = avg_primal_solution[LP, idx+1] / variable_rescaling[LP, idx+1]
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
                        iterations[LP] = iteration
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(iteration))
                    end
                    LP += grid_stride
                    break
                elseif (reason == TERMINATION_REASON_ITERATION_LIMIT) || 
                    (reason == TERMINATION_REASON_KKT_MATRIX_PASS_LIMIT) ||
                    (reason == TERMINATION_REASON_NUMERICAL_ERROR)
                    while idx < n_vars # Intentionally skipping the epigraph variable
                        solutions[LP, idx] = avg_primal_solution[LP, idx+1] / variable_rescaling[LP, idx+1]
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
                        iterations[LP] = iteration
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(iteration))
                    end
                    LP += grid_stride
                    break
                elseif (reason == TERMINATION_REASON_GLOBAL_UPPER_BOUND_HIT)
                    while idx < n_vars
                        solutions[LP, idx] = avg_primal_solution[LP, idx+1] / variable_rescaling[LP, idx+1]
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
                        iterations[LP] = iteration
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(iteration))
                    end
                    LP += grid_stride
                    break
                elseif (reason == TERMINATION_REASON_IMPATIENCE)
                    while idx < n_vars
                        solutions[LP, idx] = avg_primal_solution[LP, idx+1] / variable_rescaling[LP, idx+1]
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
                        iterations[LP] = iteration
                        CUDA.atomic_add!(CUDA.pointer(global_counter, 1), Int32(1))
                        CUDA.atomic_add!(CUDA.pointer(iteration_counter, 1), Int32(iteration))
                    end
                    LP += grid_stride
                    break
                end
                

                ###################################################################
                ##### Update Buffer Primal Gradient
                ###################################################################
                while idx <= n_vars
                    buffer_primal_gradient[LP, idx] = scaled_objective_vector[LP, idx] - current_dual_product[LP, idx]
                    idx += block_stride
                end
                idx = threadIdx().x


                ###################################################################
                ##### Run the Restart Scheme
                ###################################################################
                if (sum_primal_solutions_count <= 0) || (sum_dual_solutions_count <= 0)
                    if idx==1
                        restart_choice = RESTART_CHOICE_NO_RESTART
                    end
                else
                    if idx==1
                        # Initialize do_restart
                        do_restart[1] = false

                        if sum_primal_solutions_count >= 0.36 * (iteration - 1)
                            do_restart[1] = true
                        end
                    end

                    # Check CURRENT residual

                    # Compute the current primal residual
                    while idx <= n_vars
                        buffer_kkt_lower_variable_violation[LP, idx] = max(scaled_variable_lower_bounds[LP, idx] - current_primal_solution[LP, idx], 0.0)
                        buffer_kkt_upper_variable_violation[LP, idx] = max(current_primal_solution[LP, idx] - scaled_variable_upper_bounds[LP, idx], 0.0)
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # Compute the current primal objective and l2 primal residual
                    while idx <= n_vars
                        shared_space[idx] = scaled_objective_vector[LP, idx] * current_primal_solution[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, var_stride, n_vars)
                    if idx==1
                        primal_objective_storage = shared_space[1] + scaled_objective_constant[LP]
                    end

                    while idx <= n_vars
                        shared_space[idx] = abs(buffer_kkt_lower_variable_violation[LP, idx])^2 +
                                            abs(buffer_kkt_upper_variable_violation[LP, idx])^2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, var_stride, n_vars)
                    if idx==1
                        l2_primal_residual = shared_space[1]
                    end
                    while idx <= current_LP_length
                        shared_space[idx] = abs(max(scaled_right_hand_side[active_row + idx] - current_primal_product[active_row + idx], 0.0))^2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                    if idx==1
                        l2_primal_residual = sqrt(shared_space[1] + l2_primal_residual)
                    end

                    # Compute the current dual residual
                    while idx <= n_vars
                        buffer_kkt_reduced_costs[LP, idx] = max(buffer_primal_gradient[LP, idx], 0.0) * isfinite(scaled_variable_lower_bounds[LP, idx]) + 
                                                            min(buffer_primal_gradient[LP, idx], 0.0) * isfinite(scaled_variable_upper_bounds[LP, idx])
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # Calculate the current dual objective
                    while idx <= n_vars
                        if buffer_kkt_reduced_costs[LP, idx] > 0.0
                            shared_space[idx] = scaled_variable_lower_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                        elseif buffer_kkt_reduced_costs[LP, idx] < 0.0
                            shared_space[idx] = scaled_variable_upper_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                        else
                            shared_space[idx] = 0.0
                        end
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, var_stride, n_vars)
                    if idx==1
                        buffer_kkt_dual_objective = scaled_objective_constant[LP] + shared_space[1]
                    end
                    while idx <= current_LP_length
                        shared_space[idx] = scaled_right_hand_side[active_row + idx] * current_dual_solution[active_row + idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                    if idx==1
                        buffer_kkt_dual_objective += shared_space[1]
                    end
                    
                    # Calculate the l2 dual residual
                    while idx <= n_vars
                        shared_space[idx] = abs(buffer_primal_gradient[LP, idx] - buffer_kkt_reduced_costs[LP, idx])^2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, var_stride, n_vars)
                    if idx==1
                        l2_dual_residual = shared_space[1]
                    end
                    while idx <= current_LP_length
                        shared_space[idx] = abs(max(-current_dual_solution[active_row + idx], 0.0))^2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                    if idx==1
                        l2_dual_residual += shared_space[1]
                    end

                    # Calculate the "current" kkt residual
                    if idx==1
                        current_kkt_residual = sqrt(primal_weight[1] * l2_primal_residual^2 + 
                                                    (1/primal_weight[1]) * l2_dual_residual^2 + 
                                                    abs(primal_objective_storage - buffer_kkt_dual_objective)^2)
                    end

                    # Check AVERAGE residual

                    # Compute the current primal residual
                    while idx <= n_vars
                        buffer_kkt_lower_variable_violation[LP, idx] = max(scaled_variable_lower_bounds[LP, idx] - avg_primal_solution[LP, idx], 0.0)
                        buffer_kkt_upper_variable_violation[LP, idx] = max(avg_primal_solution[LP, idx] - scaled_variable_upper_bounds[LP, idx], 0.0)
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # Compute the current primal objective and l2 primal residual
                    while idx <= n_vars
                        shared_space[idx] = scaled_objective_vector[LP, idx] * avg_primal_solution[LP, idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, var_stride, n_vars)
                    if idx==1
                        primal_objective_storage = shared_space[1] + scaled_objective_constant[LP]
                    end

                    while idx <= n_vars
                        shared_space[idx] = abs(buffer_kkt_lower_variable_violation[LP, idx])^2 + 
                                            abs(buffer_kkt_upper_variable_violation[LP, idx])^2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, var_stride, n_vars)
                    if idx==1
                        l2_primal_residual = shared_space[1]
                    end
                    while idx <= current_LP_length
                        shared_space[idx] = abs(max(scaled_right_hand_side[active_row + idx] - avg_primal_product[active_row + idx], 0.0))^2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                    if idx==1
                        l2_primal_residual = sqrt(shared_space[1] + l2_primal_residual)
                    end

                    # Compute the current dual residual
                    while idx <= n_vars
                        buffer_kkt_reduced_costs[LP, idx] = max(avg_primal_gradient[LP, idx], 0.0) * isfinite(scaled_variable_lower_bounds[LP, idx]) + 
                                                            min(avg_primal_gradient[LP, idx], 0.0) * isfinite(scaled_variable_upper_bounds[LP, idx])
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # Calculate the current dual objective
                    while idx <= n_vars
                        if buffer_kkt_reduced_costs[LP, idx] > 0.0
                            shared_space[idx] = scaled_variable_lower_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                        elseif buffer_kkt_reduced_costs[LP, idx] < 0.0
                            shared_space[idx] = scaled_variable_upper_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                        else
                            shared_space[idx] = 0.0
                        end
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, var_stride, n_vars)
                    if idx==1
                        buffer_kkt_dual_objective = scaled_objective_constant[LP] + shared_space[1]
                    end
                    while idx <= current_LP_length
                        shared_space[idx] = scaled_right_hand_side[active_row + idx] * avg_dual_solution[active_row + idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                    if idx==1
                        buffer_kkt_dual_objective += shared_space[1]
                    end

                    # Calculate the l2 dual residual
                    while idx <= n_vars
                        shared_space[idx] = abs(avg_primal_gradient[LP, idx] - buffer_kkt_reduced_costs[LP, idx])^2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, var_stride, n_vars)
                    if idx==1
                        l2_dual_residual = shared_space[1]
                    end
                    while idx <= current_LP_length
                        shared_space[idx] = abs(max(-avg_dual_solution[active_row + idx], 0.0))^2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                    if idx==1
                        l2_dual_residual += shared_space[1]
                    end

                    # Calculate the "current" kkt residual
                    if idx==1
                        avg_kkt_residual = sqrt(primal_weight[1] * l2_primal_residual^2 + 
                                                    (1/primal_weight[1]) * l2_dual_residual^2 + 
                                                    abs(primal_objective_storage - buffer_kkt_dual_objective)^2)
                    end

                    # Determine if we should reset to average
                    if idx==1
                        if avg_kkt_residual < current_kkt_residual
                            reset_to_average[1] = true
                            candidate_kkt_residual = avg_kkt_residual
                        else
                            reset_to_average[1] = false
                            candidate_kkt_residual = current_kkt_residual
                        end
                    end

                    sync_threads()
                    if unsafe_load(CUDA.pointer(do_restart, 1))==false
                        # Check LAST RESTART residual

                        # Compute the current primal residual
                        while idx <= n_vars
                            buffer_kkt_lower_variable_violation[LP, idx] = max(scaled_variable_lower_bounds[LP, idx] - last_restart_primal_solution[LP, idx], 0.0)
                            buffer_kkt_upper_variable_violation[LP, idx] = max(last_restart_primal_solution[LP, idx] - scaled_variable_upper_bounds[LP, idx], 0.0)
                            idx += block_stride
                        end
                        idx = threadIdx().x

                        # Compute the current primal objective and l2 primal residual
                        while idx <= n_vars
                            shared_space[idx] = scaled_objective_vector[LP, idx] * last_restart_primal_solution[LP, idx]
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        parallel_sum(shared_space, block_stride, var_stride, n_vars)
                        if idx==1
                            primal_objective_storage = scaled_objective_constant[LP] + shared_space[1]
                        end

                        while idx <= n_vars
                            shared_space[idx] = abs(buffer_kkt_lower_variable_violation[LP, idx])^2 +
                                                abs(buffer_kkt_upper_variable_violation[LP, idx])^2
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        parallel_sum(shared_space, block_stride, var_stride, n_vars)
                        if idx==1
                            l2_primal_residual = shared_space[1]
                        end
                        while idx <= current_LP_length
                            shared_space[idx] = abs(max(scaled_right_hand_side[active_row + idx] - last_restart_primal_product[active_row + idx], 0.0))^2
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                        if idx==1
                            l2_primal_residual = sqrt(shared_space[1] + l2_primal_residual)
                        end

                        # Compute the current dual residual
                        while idx <= n_vars
                            buffer_kkt_reduced_costs[LP, idx] = max(last_restart_primal_gradient[LP, idx], 0.0) * isfinite(scaled_variable_lower_bounds[LP, idx]) + 
                                                                min(last_restart_primal_gradient[LP, idx], 0.0) * isfinite(scaled_variable_upper_bounds[LP, idx])
                            idx += block_stride
                        end
                        idx = threadIdx().x
                    
                        # Calculate the current dual objective
                        while idx <= n_vars
                            if buffer_kkt_reduced_costs[LP, idx] > 0.0
                                shared_space[idx] = scaled_variable_lower_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                            elseif buffer_kkt_reduced_costs[LP, idx] < 0.0
                                shared_space[idx] = scaled_variable_upper_bounds[LP, idx] * buffer_kkt_reduced_costs[LP, idx]
                            else
                                shared_space[idx] = 0.0
                            end
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        parallel_sum(shared_space, block_stride, var_stride, n_vars)
                        if idx==1
                            buffer_kkt_dual_objective = scaled_objective_constant[LP] + shared_space[1]
                        end
                        while idx <= current_LP_length
                            shared_space[idx] = scaled_right_hand_side[active_row + idx] * last_restart_dual_solution[active_row + idx]
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                        if idx==1
                            buffer_kkt_dual_objective += shared_space[1]
                        end

                        # Calculate the l2 dual residual
                        while idx <= n_vars
                            shared_space[idx] = abs(last_restart_primal_gradient[LP, idx] - buffer_kkt_reduced_costs[LP, idx])^2
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        parallel_sum(shared_space, block_stride, var_stride, n_vars)
                        if idx==1
                            l2_dual_residual = shared_space[1]
                        end
                        while idx <= current_LP_length
                            shared_space[idx] = abs(max(-last_restart_dual_solution[active_row + idx], 0.0))^2
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                        if idx==1
                            l2_dual_residual += shared_space[1]
                        end

                        # Calculate the "current" kkt residual
                        if idx==1
                            last_kkt_residual = sqrt(primal_weight[1] * l2_primal_residual^2 + 
                                                        (1/primal_weight[1]) * l2_dual_residual^2 + 
                                                        abs(primal_objective_storage - buffer_kkt_dual_objective)^2)
                            kkt_reduction_ratio = candidate_kkt_residual / last_kkt_residual

                            # Check if we need to do a restart
                            if kkt_reduction_ratio < necessary_reduction_for_restart
                                if kkt_reduction_ratio < sufficient_reduction_for_restart
                                    do_restart[1] = true
                                elseif kkt_reduction_ratio > last_reduction_ratio
                                    do_restart[1] = true
                                end
                            end
                            last_reduction_ratio = kkt_reduction_ratio
                        end
                    end


                #>>>###################################################################
                #>>>##### Apply the Reset (if do_restart==true)
                #>>>###################################################################
                    sync_threads()
                    if unsafe_load(CUDA.pointer(do_restart, 1)) == true
                        if unsafe_load(CUDA.pointer(reset_to_average, 1)) == true
                            while idx <= n_vars
                                current_primal_solution[LP, idx] = avg_primal_solution[LP, idx]
                                current_dual_product[LP, idx] = scaled_objective_vector[LP, idx] - avg_primal_gradient[LP, idx]
                                buffer_primal_gradient[LP, idx] = avg_primal_gradient[LP, idx]
                                idx += block_stride
                            end
                            idx = threadIdx().x
                            while idx <= current_LP_length
                                current_dual_solution[active_row + idx] = avg_dual_solution[active_row + idx]
                                current_primal_product[active_row + idx] = avg_primal_product[active_row + idx]
                                idx += block_stride
                            end
                            idx = threadIdx().x
                        end

                        # Reset sum variables
                        sum_primal_solutions_count = Int32(0)
                        sum_dual_solutions_count = Int32(0)
                        sum_primal_solution_weights = 0.0
                        sum_dual_solution_weights = 0.0
                        while idx <= n_vars
                            sum_primal_solutions[LP, idx] = 0.0
                            sum_dual_product[LP, idx] = 0.0
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        while idx <= current_LP_length
                            sum_dual_solutions[active_row + idx] = 0.0
                            sum_primal_product[active_row + idx] = 0.0
                            idx += block_stride
                        end
                        idx = threadIdx().x


                #>>>>>>>###################################################################
                #>>>>>>>##### Update Information About the Last Restart
                #>>>>>>>###################################################################
                        # Primal distance
                        while idx <= n_vars
                            shared_space[idx] = (avg_primal_solution[LP, idx] - last_restart_primal_solution[LP, idx])^2
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        parallel_sum(shared_space, block_stride, var_stride, n_vars)
                        if idx==1
                            last_restart_primal_distance_moved = sqrt(primal_norm_params)*sqrt(shared_space[1])/sqrt(primal_weight[1])
                        end
                        
                        # Dual distance
                        while idx <= current_LP_length
                            shared_space[idx] = (avg_dual_solution[active_row + idx] - last_restart_dual_solution[active_row + idx])^2
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                        if idx==1
                            last_restart_dual_distance_moved = sqrt(dual_norm_params)*sqrt(shared_space[1])*sqrt(primal_weight[1])
                        end

                        # Primal solution
                        while idx <= n_vars
                            last_restart_primal_solution[LP, idx] = current_primal_solution[LP, idx]
                            idx += block_stride
                        end
                        idx = threadIdx().x

                        # Dual solution
                        while idx <= current_LP_length
                            last_restart_dual_solution[active_row + idx] = current_dual_solution[active_row + idx]
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        
                        # Primal product
                        while idx <= current_LP_length
                            last_restart_primal_product[active_row + idx] = current_primal_product[active_row + idx]
                            idx += block_stride
                        end
                        idx = threadIdx().x

                        # Primal gradient
                        while idx <= n_vars
                            last_restart_primal_gradient[LP, idx] = buffer_primal_gradient[LP, idx]
                            idx += block_stride
                        end
                        idx = threadIdx().x

                        if idx==1
                            if reset_to_average[1]==true
                                restart_choice = RESTART_CHOICE_RESTART_TO_AVERAGE
                            else
                                restart_choice = RESTART_CHOICE_WEIGHTED_AVERAGE_RESET
                            end
                        end
                    else
                        if idx==1
                            restart_choice = RESTART_CHOICE_NO_RESTART
                        end
                    end
                end


                ###################################################################
                ##### Compute a New Primal Weight (if a restart was used)
                ###################################################################
                if idx==1
                    if restart_choice != RESTART_CHOICE_NO_RESTART
                        if (last_restart_primal_distance_moved > eps()) &&
                        (last_restart_dual_distance_moved > eps())
                            primal_weight[1] = exp(0.5 * log(last_restart_dual_distance_moved / last_restart_primal_distance_moved) +
                                                (1 - 0.5) * log(primal_weight[1]))
                        end
                    end
                end
                sync_threads()
            end


            ###################################################################
            ##### Take n Steps
            ###################################################################
            if iteration < 10
                n_steps = Int32(1)
            else
                n_steps = Int32(min(termination_evaluation_frequency - mod(iteration, termination_evaluation_frequency) + 1,
                                iteration_limit - iteration))
            end

            # Store the primal weight locally in each thread
            sync_threads()
            local_primal_weight = unsafe_load(CUDA.pointer(primal_weight, 1))
            local_step_size = unsafe_load(CUDA.pointer(step_size, 1))

            step = Int32(1)
            while step <= n_steps
                # Prepare a `done` flag for this particular step
                done = false
                iter = 0
                
                while !done
                    iter += 1
                    if iter > 100000 # Prevent an endless loop if no progress is made finding a step size
                        if idx==1 
                            numerical_error[1] = true
                        end
                        break
                    end
                    # Store the current step size locally in each thread
                    sync_threads()
                    local_step_size = unsafe_load(CUDA.pointer(step_size, 1))

                    # Step 1) Add one to the total iterations tracker and cumulative kkt passes
                    step_iterations += Int32(1)
                    if idx==1
                        cumulative_kkt_passes += Int32(1)
                    end
                    
                    # Step 2) Operate on matrices to:
                    # A) Update delta primal
                    # B) Calculate the norm of delta primal
                    while idx <= n_vars
                        delta_primal[LP, idx] = min(scaled_variable_upper_bounds[LP, idx], max(scaled_variable_lower_bounds[LP, idx], current_primal_solution[LP, idx] - 
                                                    (local_step_size/local_primal_weight) * (scaled_objective_vector[LP, idx] - current_dual_product[LP, idx]))) - 
                                                    current_primal_solution[LP, idx]
                        shared_space[idx] = delta_primal[LP, idx]^2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, var_stride, n_vars)
                    if idx==1
                        norm_delta_primal_sq = shared_space[1]
                    end

                    # Step 3) Perform the following steps for each row in the current LP:
                    # A) Update delta primal product using LP-segmented matrix-vector multiplication (size is (n_constraints,))
                    #   (delta_primal_product[row] = sum(constraint_matrix[row, :] .* delta_primal[LP, :])
                    # B) Update delta dual
                    # C) Calculate the norm of delta_dual
                    # D) Compute the interaction and movement terms

                    # Reset delta_primal_product and the shared space
                    while idx <= current_LP_length
                        delta_primal_product[active_row + idx] = 0.0
                        shared_space[idx] = 0.0
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    sync_threads()

                    # Loop over nonzeros to update delta_primal_product
                    while idx <= nz_count
                        if active_constraint[active_row + nz_rows[idx]]
                            CUDA.atomic_add!(CUDA.pointer(shared_space, nz_rows[idx]), scaled_constraint_matrix[active_row + nz_rows[idx], nz_cols[idx]] * delta_primal[LP, nz_cols[idx]])
                        end
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    sync_threads()
                    while idx <= current_LP_length
                        delta_primal_product[active_row + idx] += shared_space[idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x

                    # Compute other terms normally
                    while idx <= current_LP_length
                        delta_dual[active_row + idx] = max(0.0, (current_dual_solution[active_row + idx] + 
                                            (local_primal_weight * local_step_size) * (scaled_right_hand_side[active_row + idx] -
                                            ((Int32(1) + extrapolation_coefficient) * 
                                            delta_primal_product[active_row + idx]) - 
                                            (extrapolation_coefficient * current_primal_product[active_row + idx])))) - 
                                            current_dual_solution[active_row + idx]
                        shared_space[idx] = delta_primal_product[active_row + idx] * delta_dual[active_row + idx]
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                    if idx==1
                        interaction = abs(shared_space[1])
                    end
                    sync_threads()
                    while idx <= current_LP_length
                        shared_space[idx] = delta_dual[active_row + idx]^2
                        idx += block_stride
                    end
                    idx = threadIdx().x
                    parallel_sum(shared_space, block_stride, len_stride, current_LP_length)
                    if idx==1
                        norm_delta_dual_sq = shared_space[1]
                        movement = 0.5 * local_primal_weight * norm_delta_primal_sq + (0.5 / local_primal_weight) * norm_delta_dual_sq
                    end

                    # Step 4) Update the step size limit, and check for numerical errors
                    if idx==1
                        if interaction > 0.0
                            step_size_limit[1] = movement/interaction
                            if iszero(movement)
                                numerical_error[1] = true
                            end
                        else
                            step_size_limit[1] = Inf
                        end
                    end
                    
                    # Step 5) Update the solution (if we're done)
                    sync_threads()
                    if local_step_size <= unsafe_load(CUDA.pointer(step_size_limit, 1))
                        done = true # This will be the final iteration of this step
                        while idx <= n_vars
                            current_primal_solution[LP, idx] += delta_primal[LP, idx]
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        while idx <= current_LP_length
                            current_primal_product[active_row + idx] += delta_primal_product[active_row + idx]
                            current_dual_solution[active_row + idx] += delta_dual[active_row + idx]
                            idx += block_stride
                        end
                        idx = threadIdx().x

                        # Reset current_dual_product and the shared space
                        while idx <= n_vars
                            current_dual_product[LP, idx] = 0.0
                            shared_space[idx] = 0.0
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        sync_threads()

                        # Loop over nonzeros to update current_dual_product
                        while idx <= nz_count
                            if active_constraint[active_row + nz_rows[idx]]
                                CUDA.atomic_add!(CUDA.pointer(shared_space, nz_cols[idx]), scaled_constraint_matrix[active_row + nz_rows[idx], nz_cols[idx]] * current_dual_solution[active_row + nz_rows[idx]])
                            end
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        sync_threads()
                        while idx <= n_vars
                            current_dual_product[LP, idx] += shared_space[idx]
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        sync_threads()
                        
                        # Also update solution weighted average
                        while idx <= n_vars
                            sum_primal_solutions[LP, idx] += current_primal_solution[LP, idx] * local_step_size
                            sum_dual_product[LP, idx] += current_dual_product[LP, idx] * local_step_size
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        sum_primal_solutions_count += Int32(1)
                        sum_primal_solution_weights += local_step_size
                        while idx <= current_LP_length
                            sum_primal_product[active_row + idx] += current_primal_product[active_row + idx] * local_step_size
                            sum_dual_solutions[active_row + idx] += current_dual_solution[active_row + idx] * local_step_size
                            idx += block_stride
                        end
                        idx = threadIdx().x
                        sum_dual_solutions_count += Int32(1)
                        sum_dual_solution_weights += local_step_size
                    end


                    if idx==1
                        step_size[1] = min(((1.0 - 1.0/(step_iterations + 1.0)^reduction_exponent)) * step_size_limit[1],
                                            ((1.0 + 1.0/(step_iterations + 1.0)^growth_exponent)) * step_size[1])
                    end
                    sync_threads()
                end

                step += Int32(1)
            end
            iteration += n_steps
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