
function PDLP(
    PDLP_data::PDLPData; 
    solutions::CuArray{Float64}=CuArray{Float64}(undef, PDLP_data.dims.n_LPs, PDLP_data.dims.n_vars), 
    objectives::CuArray{Float64}=CuArray{Float64}(undef, PDLP_data.dims.n_LPs),
    return_dual_obj::Bool=false,
    return_both_obj::Bool=false,
    global_upper_bound::Float64=Inf
    )

    # Check return conditions and verify that `objectives` storage is correctly sized
    if return_dual_obj + return_both_obj == 2
        error("Only one return condition allowed")
    end

    if !return_dual_obj && !return_both_obj
        return_code = Int32(1)
        if size(objectives, 2) != 1
            error("Objective storage sized incorrectly")
        end
    elseif return_dual_obj
        return_code = Int32(2)
        if size(objectives, 2) != 1
            error("Objective storage sized incorrectly")
        end
    elseif return_both_obj
        return_code = Int32(3)
        if size(objectives, 2) != 2
            error("Objective storage sized incorrectly")
        end
    end

    # Validate the LP data we've been given to make sure the numbers are all valid and the dimensions
    # of participating matrices are correct
    validate(PDLP_data)

    # Reset all fields relevant to problem status
    reset_all_fields!(PDLP_data)

    # Perform rescaling. Note that if hot-starting is to be added in the future, you should save all
    # primal/dual information (rather than only primal solutions and dual objectives), and then here,
    # instead of rescaling the problem immediately, first un-scale the primal/dual statuses of the
    # previous run, then do re-scaling, then scale the primal/dual statuses again. This will also
    # require adding storage for these values, and un-commenting some lines in `main_loop.jl`
    # to allow hot-starting to impact the main PDLP algorithm. 
    rescale_problem(
        PDLP_data.original_problem, 
        PDLP_data.scaled_problem, 
        PDLP_data.variable_rescaling, 
        PDLP_data.constraint_rescaling, 
        PDLP_data.dims,
        PDLP_data.parameters
        )

    # Scale the primal weight if desired (otherwise it should be 1.0)
    # (Could also put this inside the main kernel)
    if PDLP_data.parameters.scale_initial_primal_weight
        select_initial_primal_weight(PDLP_data.primal_weight, PDLP_data.scaled_problem, PDLP_data.dims)
    else
        PDLP_data.primal_weight .= 1.0
    end

    # Come up with a starting step size (Could also put this inside the kernel)
    update_step_size(PDLP_data.scaled_problem, PDLP_data.step_size, PDLP_data.dims)

    # Run the main loop kernel
    max_size = max(PDLP_data.dims.n_vars, PDLP_data.dims.current_LP_length)
    max_req = Int32(min(256, max(32, ceil(Int, max_size/32)*32)))

    # Reset total solve and iteration number counters
    PDLP_data.global_counter .= Int32(0)
    PDLP_data.iteration_counter .= Int32(0)

    # Call the main PDLP kernel
    CUDA.@sync @cuda blocks=PDLP_data.dims.n_LPs threads=max_req shmem=max_size*sizeof(Float64) main_loop_kernel(
            solutions,
            objectives,
            PDLP_data.original_problem.variable_lower_bounds,
            PDLP_data.original_problem.variable_upper_bounds,
            PDLP_data.original_problem.right_hand_side,
            PDLP_data.original_problem.objective_vector,
            PDLP_data.original_problem.objective_constant,
            PDLP_data.scaled_problem.variable_lower_bounds,
            PDLP_data.scaled_problem.variable_upper_bounds,
            PDLP_data.scaled_problem.constraint_matrix,
            PDLP_data.scaled_problem.right_hand_side,
            PDLP_data.scaled_problem.objective_vector,
            PDLP_data.scaled_problem.objective_constant,
            PDLP_data.sparsity.nz_count[PDLP_data.dims.current_LP_length],
            PDLP_data.sparsity.nz_rows,
            PDLP_data.sparsity.nz_cols,
            PDLP_data.active_constraint,
            PDLP_data.variable_rescaling, 
            PDLP_data.constraint_rescaling, 
            PDLP_data.kernel_storage.last_restart_primal_solution,
            PDLP_data.kernel_storage.last_restart_primal_gradient,
            PDLP_data.kernel_storage.last_restart_dual_solution,
            PDLP_data.kernel_storage.last_restart_primal_product,
            PDLP_data.kernel_storage.current_primal_solution,
            PDLP_data.kernel_storage.current_dual_solution,
            PDLP_data.kernel_storage.current_dual_product,
            PDLP_data.kernel_storage.current_primal_product,
            PDLP_data.kernel_storage.buffer_primal_gradient,
            PDLP_data.kernel_storage.avg_primal_solution,
            PDLP_data.kernel_storage.avg_primal_gradient,
            PDLP_data.kernel_storage.avg_dual_solution,
            PDLP_data.kernel_storage.avg_primal_product,
            PDLP_data.kernel_storage.sum_primal_solutions,
            PDLP_data.kernel_storage.sum_dual_solutions,
            PDLP_data.kernel_storage.sum_primal_product,
            PDLP_data.kernel_storage.sum_dual_product,
            PDLP_data.kernel_storage.original_primal_solution,
            PDLP_data.kernel_storage.original_primal_gradient,
            PDLP_data.kernel_storage.original_dual_solution,
            PDLP_data.kernel_storage.original_primal_product,
            PDLP_data.kernel_storage.buffer_kkt_primal_solution,
            PDLP_data.kernel_storage.buffer_kkt_primal_product,
            PDLP_data.kernel_storage.buffer_kkt_lower_variable_violation,
            PDLP_data.kernel_storage.buffer_kkt_upper_variable_violation,
            PDLP_data.kernel_storage.buffer_kkt_reduced_costs,
            PDLP_data.kernel_storage.delta_primal,
            PDLP_data.kernel_storage.delta_primal_product,
            PDLP_data.kernel_storage.delta_dual,
            PDLP_data.primal_weight,
            PDLP_data.step_size,
            PDLP_data.termination_reason,
            PDLP_data.dims.current_LP_length,
            PDLP_data.dims.total_LP_length,
            PDLP_data.dims.n_LPs,
            PDLP_data.dims.n_vars,
            PDLP_data.parameters.iteration_limit,
            PDLP_data.parameters.kkt_matrix_pass_limit,
            PDLP_data.parameters.termination_evaluation_frequency,
            PDLP_data.parameters.necessary_reduction_for_restart,
            PDLP_data.parameters.sufficient_reduction_for_restart,
            PDLP_data.parameters.extrapolation_coefficient,
            PDLP_data.parameters.reduction_exponent,
            PDLP_data.parameters.growth_exponent,
            PDLP_data.parameters.termination_criteria.eps_optimal_absolute,
            PDLP_data.parameters.termination_criteria.eps_optimal_relative,
            PDLP_data.parameters.termination_criteria.eps_primal_infeasible,
            PDLP_data.parameters.termination_criteria.eps_dual_infeasible,
            return_code,
            global_upper_bound,
            PDLP_data.parameters.skip_hard_problems,
            PDLP_data.global_counter,
            PDLP_data.iteration_counter,
            PDLP_data.skip_flag,
            PDLP_data.iterations,
            )
    return nothing
end

# Check that all the LPs we've created are valid and that the sizes of objects are correct
function validate(PDLP_data::PDLPData)
    # Make sure all entries in the constraint matrix and right-hand side are valid
    if any(!isfinite, PDLP_data.original_problem.constraint_matrix)
        error("Something in the constraint matrix is NaN or Inf")
    end
    if any(!isfinite, PDLP_data.original_problem.right_hand_side)
        error("Something in the right-hand side is NaN or Inf")
    end

    # Make sure the dimensions all line up
    n_LPs = PDLP_data.dims.n_LPs
    n_vars = PDLP_data.dims.n_vars
    tot_len = PDLP_data.dims.total_LP_length

    if size(PDLP_data.original_problem.variable_lower_bounds, 1) < n_LPs || 
       size(PDLP_data.original_problem.variable_lower_bounds, 2) != n_vars
        @show size(PDLP_data.original_problem.variable_lower_bounds)
        @show (n_LPs, n_vars)
        error("Lower bound matrix is the wrong size, or the listed number of LPs/variables is incorrect.")
    end
    if size(PDLP_data.original_problem.variable_upper_bounds, 1) < n_LPs || 
       size(PDLP_data.original_problem.variable_upper_bounds, 2) != n_vars
        @show size(PDLP_data.original_problem.variable_upper_bounds)
        @show (n_LPs, n_vars)
        error("Lower bound matrix is the wrong size, or the listed number of LPs/variables is incorrect.")
    end
    if size(PDLP_data.original_problem.constraint_matrix, 1) < n_LPs*tot_len || 
       size(PDLP_data.original_problem.constraint_matrix, 2) != n_vars
        @show size(PDLP_data.original_problem.constraint_matrix)
        @show (n_LPs*tot_len, n_vars)
        error("Constraint matrix is the wrong size, or the listed number of LPs/variables is incorrect.")
    end
    if size(PDLP_data.original_problem.right_hand_side, 1) < n_LPs*tot_len
        @show size(PDLP_data.original_problem.right_hand_side)
        @show (n_LPs*tot_len)
        error("Right-hand side is the wrong size, or the listed number of LPs/variables is incorrect.")
    end
    if size(PDLP_data.original_problem.objective_vector, 1) < n_LPs || 
       size(PDLP_data.original_problem.objective_vector, 2) != n_vars
       @show size(PDLP_data.original_problem.objective_vector)
       @show (n_LPs, n_vars)
        error("Objective vector is the wrong size, or the listed number of LPs/variables is incorrect.")
    end
    if size(PDLP_data.original_problem.objective_constant, 1) < n_LPs
        @show size(PDLP_data.original_problem.objective_constant)
        @show n_LPs
        error("Objective constant is the wrong size, or the listed number of LPs/variables is incorrect.")
    end
    return nothing
end

function rescale_problem(
    original_problem::LinearProgramSet, 
    scaled_problem::LinearProgramSet, 
    variable_rescaling::CuArray{Float64}, 
    constraint_rescaling::CuArray{Float64},  
    dims::PDLPDims,
    params::PDLPParams
    )

    # Copy problem info to the scaled problem so as to not overwrite the original
    copyto!(scaled_problem, original_problem)

    # Reset rescaling values to be 1.0.
    constraint_rescaling .= 1.0
    variable_rescaling .= 1.0

    # Perform L_inf Ruiz rescaling for `ruiz_iterations` iterations
    if params.ruiz_iterations > 0 # Default is 10
        ruiz_rescaling(
            scaled_problem, 
            params.ruiz_iterations, 
            variable_rescaling, 
            constraint_rescaling, 
            dims,
            )
    end

    # cuPDLP also has a section for l2 norm rescaling, but it's off by default, so
    # I haven't included it here

    # Perform Pock-Chambolle rescaling if alpha isn't nothing
    if !isnothing(params.pock_chambolle_alpha) # Default is 1.0
        pock_chambolle_rescaling(
            scaled_problem, 
            params.pock_chambolle_alpha, 
            variable_rescaling, 
            constraint_rescaling,  
            dims,
            )
    end
    
    return nothing
end

function reset_all_fields!(PDLP_data::PDLPData)
    # Reset all kernel storage
    for field in fieldnames(KernelStorage)
        CUDA.fill!(getfield(PDLP_data.kernel_storage, field), 0.0)
    end

    # Reset the termination indicator
    CUDA.fill!(PDLP_data.termination_reason, TERMINATION_REASON_UNSPECIFIED)
end