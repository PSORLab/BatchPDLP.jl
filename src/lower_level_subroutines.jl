
# Rescaling and update functions called from `primal_subroutines.jl`. These
# are generally based on cuPDLP.jl.
function ruiz_rescaling(
    problem::LinearProgramSet, 
    n_iterations::Int, 
    variable_rescaling::CuArray{Float64}, 
    constraint_rescaling::CuArray{Float64},  
    dims::PDLPDims,
    )

    # Set up temporary variables for intermediate scaling results
    temp_variable_rescaling = CuArray{Float64}(undef, size(variable_rescaling))
    temp_constraint_rescaling = CuArray{Float64}(undef, size(constraint_rescaling))

    for _ in 1:n_iterations
        # Variable rescaling. sqrt of the maximum value of each column in each LP
        CUDA.@sync @cuda blocks=GPU_blocks threads=512 ruiz_variable_kernel(
            temp_variable_rescaling, 
            problem.constraint_matrix, 
            dims.current_LP_length, 
            dims.total_LP_length
            )

        # Constraint resscaling. sqrt of the maximum value of each row of the constraint matrix
        CUDA.@sync @cuda blocks=GPU_blocks threads=512 ruiz_constraint_kernel(
            temp_constraint_rescaling, 
            problem.constraint_matrix, 
            dims.current_LP_length, 
            dims.total_LP_length
            )

        # Use the variable and constraint rescaling values to scale the problem
        scale_problem(
            problem, 
            temp_variable_rescaling, 
            temp_constraint_rescaling,  
            dims,
            )

        # Update the overall rescaling variables
        variable_rescaling .*= temp_variable_rescaling
        constraint_rescaling .*= temp_constraint_rescaling
    end

    # Free up temporary variables
    CUDA.unsafe_free!(temp_variable_rescaling)
    CUDA.unsafe_free!(temp_constraint_rescaling)
    return nothing
end

function pock_chambolle_rescaling(
    problem::LinearProgramSet, 
    alpha::Float64, 
    variable_rescaling::CuArray{Float64}, 
    constraint_rescaling::CuArray{Float64},  
    dims::PDLPDims,
    )

    # Preallocate space for intermediate rescaling terms
    temp_variable_rescaling = CuArray{Float64}(undef, size(variable_rescaling))
    temp_constraint_rescaling = CuArray{Float64}(undef, size(constraint_rescaling))

    # Call the kernels to determine variable and constraint rescaling
    CUDA.@sync @cuda blocks=GPU_blocks threads=512 pock_chambolle_variable_kernel(
        temp_variable_rescaling, 
        problem.constraint_matrix, 
        dims.current_LP_length, 
        dims.total_LP_length, 
        alpha
        )
    CUDA.@sync @cuda blocks=GPU_blocks threads=512 pock_chambolle_constraint_kernel(
        temp_constraint_rescaling, 
        problem.constraint_matrix, 
        dims.current_LP_length, 
        dims.total_LP_length,
        alpha, 
        )

    # Apply the rescaling values to the problem
    scale_problem(
        problem, 
        temp_variable_rescaling, 
        temp_constraint_rescaling,  
        dims
        )

    # Update the overall rescaling variables
    variable_rescaling .*= temp_variable_rescaling
    constraint_rescaling .*= temp_constraint_rescaling
    
    # Free up temporary variables
    CUDA.unsafe_free!(temp_variable_rescaling)
    CUDA.unsafe_free!(temp_constraint_rescaling)
    return nothing
end

function scale_problem(
    problem::LinearProgramSet, 
    variable_rescaling::CuArray{Float64}, 
    constraint_rescaling::CuArray{Float64},  
    dims::PDLPDims,
    )

    # Perform the following steps:
    # 1) problem.objective_vector = problem.objective_vector ./ variable_rescaling
    # 2) problem.variable_lower_bound = problem.variable_lower_bound .* variable_rescaling
    # 3) problem.variable_upper_bound = problem.variable_upper_bound .* variable_rescaling
    # 4) problem.right_hand_side = problem.right_hand_side ./ constraint_rescaling
    # 5) problem.constraint_matrix...
    #       Multiply each ROW by 1/constraint_rescaling (because there's that many constraints)
    #       Multiply each COL by 1/variable_rescaling (because there's that many variables)

    # Handle the first 3 steps in one kernel
    CUDA.@sync @cuda blocks=GPU_blocks threads=512 scaling_part1_kernel(
        problem.objective_vector, 
        problem.variable_lower_bounds,                            
        problem.variable_upper_bounds, 
        variable_rescaling
        )
    
    # Handle the final 2 steps in another kernel. It may be faster, though a bit messier, to
    # combine this with the previous kernel
    CUDA.@sync @cuda blocks=GPU_blocks threads=512 scaling_part2_kernel(
        problem.right_hand_side,
        problem.constraint_matrix, 
        dims.current_LP_length, 
        dims.total_LP_length,  
        variable_rescaling, 
        constraint_rescaling
        )
    return nothing
end

function select_initial_primal_weight(
    primal_weight::CuArray{Float64},
    problem::LinearProgramSet,
    dims::PDLPDims
    )
    # Theoretically the primal importance can change, but the default in the MOI_wrapper
    # is to set it to 1.0. The other parameters are un-settable in cuPDLP but theoretically
    # could be changed as well
    CUDA.@sync @cuda blocks=GPU_blocks threads=512 primal_weight_kernel(
        primal_weight, 
        problem.objective_vector, 
        problem.right_hand_side, 
        dims.n_LPs, 
        dims.total_LP_length, 
        dims.current_LP_length
        )
    return nothing
end

function update_step_size(problem::LinearProgramSet, step_size::CuArray{Float64}, dims::PDLPDims)
    CUDA.@sync @cuda blocks=GPU_blocks threads=512 group_max_kernel(
        step_size, 
        problem.constraint_matrix, 
        dims.n_LPs, 
        dims.total_LP_length, 
        dims.current_LP_length
        )
    return nothing
end

# Various functions for adding constraints to the LP set's constraint matrix and RHS storage
# based on relaxations of the objective (or other constraints) calculated via SourceCodeMcCormick.jl.

# This adds an LP constraint based on a relaxation of the objective function. I.e., it adds a
# constraint that the epigraph variable is GEQ the subtangent hyperplane of a convex relaxation
# of the objective function evaluated at `eval_points`. The relaxations themselves are input
# in `obj_result`, and there should be one relaxation per LP (as one constraint gets added
# per LP).
function add_LP_objective_constraint(
    LPs::LinearProgramSet, 
    obj_result, 
    eval_points,
    active_constraint::CuArray{Bool},
    dims::PDLPDims
    )

    # Add the constraint
    CUDA.@sync @cuda blocks=GPU_blocks threads=256 add_LP_constraint_kernel(
        LPs.constraint_matrix, 
        LPs.right_hand_side, 
        obj_result, 
        eval_points, 
        active_constraint,
        dims.n_vars, 
        dims.n_LPs,
        dims.total_LP_length,
        dims.current_LP_length + Int32(1), # Place this constraint one after the current number of constraints
        false, # GEQ flag. For relaxations of the objective function, leave this as `false`
        true # Objective function relaxation indicator. Always `true` for this function
        )
    
    # Add 1 to the number of rows for each LP
    dims.current_LP_length += Int32(1)
    return nothing
end

# This is a more general version of `add_LP_objective_constraint`, in that it handles adding
# LP constraints based on either LEQ or GEQ constraints in the original global optimization
# problem. Depending on the status of the `geq` flag, this uses either subtangent hyperplanes
# derived from the convex or concave relaxation stored in `cons_result`.
function add_LP_constraint(
    LPs::LinearProgramSet, 
    cons_result, 
    eval_points,
    active_constraint,
    dims::PDLPDims;
    geq::Bool=true
    )

    # Add the constraint
    CUDA.@sync @cuda blocks=GPU_blocks threads=256 add_LP_constraint_kernel(
        LPs.constraint_matrix, 
        LPs.right_hand_side, 
        cons_result, 
        eval_points, 
        active_constraint, 
        dims.n_vars, 
        dims.n_LPs,
        dims.total_LP_length,
        dims.current_LP_length + Int32(1), # Place this constraint one after the current number of constraints
        geq, 
        false # Objective function relaxation indicator. Always `false` for this function
        )

    # Add 1 to the number of rows for each LP
    dims.current_LP_length += Int32(1)
    return nothing
end

# A relatively simple version of adding an LP constraint. This adds a lower bound to the
# epigraph variable as a constraint, and uses the lower bound of the inclusion monotonic
# interval extension from SourceCodeMcCormick.jl, stored in `obj_result`.
function add_LP_lower_bound(
    LPs::LinearProgramSet, 
    obj_result,
    active_constraint::CuArray{Bool},
    dims::PDLPDims
    )

    # Add the constraint
    CUDA.@sync @cuda blocks=GPU_blocks threads=256 add_LP_lower_bound_kernel(
        LPs.constraint_matrix, 
        LPs.right_hand_side, 
        obj_result, 
        active_constraint,
        dims.n_LPs, 
        dims.total_LP_length,
        dims.current_LP_length + Int32(1), # Place this constraint one after the current number of constraints
        )
    
    # Add 1 to the number of rows for each LP
    dims.current_LP_length += Int32(1)
    return nothing
end

# The `add_LP_objective_constraint` function assumed there was one point per LP, but
# there may be cases where multiple points are calculated per LP, and we want to add
# all of them as constraints at once. This handles that case, with the `n_points`
# value to indicate how many points to expected per LP.
# DEPRECATION WARNING: This function is outdated and does not change the `active_constraint`
#                      field. This field alerts the main PDLP kernel that the constraint is
#                      in use. If `active_constraint` is not adjusted, the constraint will be
#                      ignored.
# function add_multiple_LP_objective_constraint(
#     LPs::LinearProgramSet, 
#     obj_result, 
#     eval_points,
#     dims::PDLPDims,
#     n_points::Int32,
#     )

#     # Add the constraint
#     CUDA.@sync @cuda blocks=GPU_blocks threads=256 add_multiple_LP_constraint_kernel(
#         LPs.constraint_matrix, 
#         LPs.right_hand_side, 
#         obj_result, 
#         eval_points, 
#         dims.n_vars, 
#         dims.n_LPs,
#         dims.total_LP_length,
#         dims.current_LP_length + Int32(1), # Place this constraint one after the current number of constraints
#         false, # GEQ flag. For relaxations of the objective function, leave this as `false`
#         true, # Objective function relaxation indicator. Always `true` for this function
#         n_points,
#         )
    
#     # Add n_points to the number of rows for each LP
#     dims.current_LP_length += n_points
#     return nothing
# end

# Much like `add_multiple_LP_objective_constraint`, this function allows multiple relaxations
# per LP to be added as constraints, with `n_points` indicating how many points per LP to
# expect.
# DEPRECATION WARNING: This function is outdated and does not change the `active_constraint`
#                      field. This field alerts the main PDLP kernel that the constraint is
#                      in use. If `active_constraint` is not adjusted, the constraint will be
#                      ignored.
# function add_multiple_LP_constraint(
#     LPs::LinearProgramSet, 
#     cons_result, 
#     eval_points,
#     dims::PDLPDims,
#     n_points::Int32;
#     geq::Bool=true,
#     )

#     # Add the constraint
#     CUDA.@sync @cuda blocks=GPU_blocks threads=256 add_multiple_LP_constraint_kernel(
#         LPs.constraint_matrix, 
#         LPs.right_hand_side, 
#         cons_result, 
#         eval_points, 
#         dims.n_vars, 
#         dims.n_LPs,
#         dims.total_LP_length,
#         dims.current_LP_length + Int32(1), # Place this constraint one after the current number of constraints
#         geq, 
#         false, # Objective function relaxation indicator. Always `false` for this function
#         n_points,
#         )

#     # Add n_points to the number of rows for each LP
#     dims.current_LP_length += n_points
#     return nothing
# end

# This function assumes there are `n_points` relaxations per LP, but that we only want
# to add a select few as constraints. This takes relaxations of the objective function
# as input in `obj_result`, the points where these relaxations were evaluated in
# `eval_points`, a vector of the weighted sum of convex relaxation subgradient values squared
# in `comparison_vector`, and a value related to the signs of convex relaxation subgradient
# elements in `subgradient_checksum`. The "best" point in each pass is selected as the one
# that has the flattest overall slope (i.e., the lowest value in `comparison_vector`). 
# Once a point is selected, all points that have the same signs on all subgradient coefficients
# (i.e., the same value in `subgradient_checksum`) are removed from consideration in
# future passes. This essentially tries to enforce the idea that the multiple points
# being added provide varied information.
function add_best_obj_LP_constraints(
    LPs::LinearProgramSet,
    obj_result,
    eval_points,
    comparison_vector,
    subgradient_checksum,
    dims::PDLPDims,
    active_constraint::CuArray{Bool},
    n_points::Int32,
    num_linearizations::Int,
    )

    # Add the best constraint, based on the comparison vector, num_linearizations times
    CUDA.@sync @cuda blocks=GPU_blocks threads=256 add_best_obj_LP_constraints_kernel(
        LPs.constraint_matrix,
        LPs.right_hand_side,
        obj_result,
        eval_points,
        comparison_vector,
        subgradient_checksum,
        dims.n_vars,
        dims.n_LPs,
        dims.total_LP_length,
        dims.current_LP_length, # Placement begins after the current number of constraints
        active_constraint,
        n_points,
        Int32(num_linearizations),
        )
    
    # Add num_linearizations to the number of rows for each LP
    dims.current_LP_length += num_linearizations
    return nothing
end

# This function is similar to `add_best_obj_LP_constraints`, but instead of relaxations
# of the objective function, it adds relaxations of LEQ or GEQ constraints. Additionally,
# instead of `comparison_vector` containing the weighted sum of squared subgradient
# elements, here it should contain either the convex or concave relaxation values evaluated
# at `eval_points`. I.e., the "best" constraints are those whose evaluations occur closest
# to 0, as these points are likely to give the most accurate information about the zero
# level-set of the convex relaxation itself.
function add_best_cons_LP_constraints(
    LPs::LinearProgramSet,
    cons_result,
    eval_points,
    comparison_vector,
    subgradient_checksum,
    dims::PDLPDims,
    active_constraint::CuArray{Bool},
    n_points::Int32,
    num_linearizations::Int;
    geq::Bool=true,
    )

    # Add the best constraint, based on the comparison vector, num_linearizations times
    CUDA.@sync @cuda blocks=GPU_blocks threads=256 add_best_cons_LP_constraints_kernel(
        LPs.constraint_matrix,
        LPs.right_hand_side,
        cons_result,
        eval_points,
        comparison_vector,
        subgradient_checksum,
        dims.n_vars,
        dims.n_LPs,
        dims.total_LP_length,
        dims.current_LP_length, # Placement begins after the current number of constraints
        active_constraint,
        geq, # GEQ indicator
        n_points,
        Int32(num_linearizations),
        )
    
    # Add num_linearizations to the number of rows for each LP
    dims.current_LP_length += num_linearizations
    return nothing
end

# Even if multiple constraints are being added in general per LP, and multiple relaxations
# are calculated per LP, we still only want to add one lower bound for the epigraph variable,
# since the inclusion monotonic interval extensions calculated by SourceCodeMcCormick.jl
# do not depend on where in the domain the evaluation takes place. I.e., all the interval
# extensions will be the same for any specific domain, so we only need to add one lower bound
# constraint. This is different from `add_LP_lower_bound` in that it recognizes that `n_points`
# relaxations exist for every LP, as opposed to mapping relaxations in a 1:1 manner.
function add_multiple_LP_lower_bound(
    LPs::LinearProgramSet, 
    obj_result,
    dims::PDLPDims,
    active_constraint::CuArray{Bool},
    n_points::Int32,
    )

    # Add the constraint
    CUDA.@sync @cuda blocks=GPU_blocks threads=256 add_multiple_LP_lower_bound_kernel(
        LPs.constraint_matrix, 
        LPs.right_hand_side, 
        obj_result, 
        dims.n_LPs, 
        dims.total_LP_length,
        dims.current_LP_length + Int32(1), # Place this constraint one after the current number of constraints
        active_constraint,
        n_points,
        )
    
    # Add 1 to the number of rows for each LP (only one lower bound is actually used, from n_points)
    dims.current_LP_length += Int32(1)
    return nothing
end


