
# Structs that are needed for BatchPDLP to function
"""
TerminationReason explains why the solver stopped. See termination.jl for the
precise criteria used to check termination.

# Values

- `TERMINATION_REASON_UNSPECIFIED`: Default value.
- `TERMINATION_REASON_OPTIMAL`
- `TERMINATION_REASON_PRIMAL_INFEASIBLE`: Note in this situation the dual could
be either unbounded or infeasible.
- `TERMINATION_REASON_DUAL_INFEASIBLE`: Note in this situation the primal could
either unbounded or infeasible.
- `TERMINATION_REASON_TIME_LIMIT`
- `TERMINATION_REASON_KKT_MATRIX_PASS_LIMIT`
- `TERMINATION_REASON_NUMERICAL_ERROR`
- `TERMINATION_REASON_INVALID_PROBLEM`: Indicates that the solver detected
invalid problem data, e.g., inconsistent bounds.
- `TERMINATION_REASON_OTHER`
"""
@enum TerminationReason begin
    TERMINATION_REASON_UNSPECIFIED
    TERMINATION_REASON_OPTIMAL
    TERMINATION_REASON_PRIMAL_INFEASIBLE
    TERMINATION_REASON_DUAL_INFEASIBLE
    TERMINATION_REASON_TIME_LIMIT
    TERMINATION_REASON_ITERATION_LIMIT
    TERMINATION_REASON_KKT_MATRIX_PASS_LIMIT
    TERMINATION_REASON_NUMERICAL_ERROR
    TERMINATION_REASON_INVALID_PROBLEM
    TERMINATION_REASON_OTHER
    TERMINATION_REASON_GLOBAL_UPPER_BOUND_HIT
    TERMINATION_REASON_IMPATIENCE
    DO_NOT_TERMINATE
end

"""
RestartChoice specifies whether a restart was performed on a given iteration.

# Values

- `RESTART_CHOICE_UNSPECIFIED`: Default value.
- `RESTART_CHOICE_NO_RESTART`: No restart on this iteration.
- `RESTART_CHOICE_WEIGHTED_AVERAGE_RESET`: The weighted average of iterates is
  cleared and reset to the current point. Note that from a mathematical
  perspective this can be equivalently viewed as restarting the algorithm but
  picking the restart point to be the current iterate.
- `RESTART_CHOICE_RESTART_TO_AVERAGE`: The algorithm is restarted at the average
  of iterates since the last restart.
"""
@enum RestartChoice begin
    RESTART_CHOICE_UNSPECIFIED
    RESTART_CHOICE_NO_RESTART
    RESTART_CHOICE_WEIGHTED_AVERAGE_RESET
    RESTART_CHOICE_RESTART_TO_AVERAGE
end

mutable struct LinearProgramSet
    variable_lower_bounds::CuArray{Float64}
    variable_upper_bounds::CuArray{Float64}
    constraint_matrix::CuArray{Float64}
    right_hand_side::CuArray{Float64}
    objective_vector::CuArray{Float64}
    objective_constant::CuArray{Float64}
end
import Base.copyto!
function copyto!(to::LinearProgramSet, from::LinearProgramSet)
    copy!(to.variable_lower_bounds, from.variable_lower_bounds)
    copy!(to.variable_upper_bounds, from.variable_upper_bounds)
    copy!(to.constraint_matrix, from.constraint_matrix)
    copy!(to.right_hand_side, from.right_hand_side)
    copy!(to.objective_vector, from.objective_vector)
    copy!(to.objective_constant, from.objective_constant)
    return nothing
end


mutable struct TerminationCriteria
    # Let p correspond to the norm we are using as specified by optimality_norm.
    # If the algorithm terminates with termination_reason =
    # TERMINATION_REASON_OPTIMAL then the following hold:
    # | primal_objective - dual_objective | <= eps_optimal_absolute +
    #  eps_optimal_relative * ( | primal_objective | + | dual_objective | )
    # norm(primal_residual, p) <= eps_optimal_absolute + eps_optimal_relative *
    #  norm(right_hand_side, p)
    # norm(dual_residual, p) <= eps_optimal_absolute + eps_optimal_relative *
    #   norm(objective_vector, p)
    # It is possible to prove that a solution satisfying the above conditions
    # also satisfies SCS's optimality conditions (see link above) with ϵ_pri =
    # ϵ_dual = ϵ_gap = eps_optimal_absolute = eps_optimal_relative. (ϵ_pri,
    # ϵ_dual, and ϵ_gap are SCS's parameters).

    """
    Absolute tolerance on the duality gap, primal feasibility, and dual
    feasibility.
    """
    eps_optimal_absolute::Float64

    """
    Relative tolerance on the duality gap, primal feasibility, and dual
    feasibility.
    """
    eps_optimal_relative::Float64

    """
    If the following two conditions hold we say that we have obtained an
    approximate dual ray, which is an approximate certificate of primal
    infeasibility.
    (1) dual_ray_objective > 0.0,
    (2) max_dual_ray_infeasibility / dual_ray_objective <=
        eps_primal_infeasible.
    """
    eps_primal_infeasible::Float64

    """
    If the following three conditions hold we say we have obtained an
    approximate primal ray, which is an approximate certificate of dual
    infeasibility.
    (1) primal_ray_linear_objective < 0.0,
    (2) max_primal_ray_infeasibility / (-primal_ray_linear_objective) <=
        eps_dual_infeasible,
    (3) primal_ray_quadratic_norm / (-primal_ray_linear_objective) <=
        eps_dual_infeasible.
    """
    eps_dual_infeasible::Float64

    """
    If termination_reason = TERMINATION_REASON_TIME_LIMIT then the solver has
    taken at least time_sec_limit time. (NOTE: currently not used)
    """
    time_sec_limit::Float64 # NOTE: currently not used

    """
    If termination_reason = TERMINATION_REASON_ITERATION_LIMIT then the solver has taken at least iterations_limit iterations.
    """
    iteration_limit::Int32

    """
    If termination_reason = TERMINATION_REASON_KKT_MATRIX_PASS_LIMIT then
    cumulative_kkt_matrix_passes is at least kkt_pass_limit.
    """
    kkt_matrix_pass_limit::Float64
end

mutable struct KernelStorage
    last_restart_primal_solution::CuArray{Float64}
    last_restart_primal_gradient::CuArray{Float64}
    last_restart_dual_solution::CuArray{Float64}
    last_restart_primal_product::CuArray{Float64}
    current_primal_solution::CuArray{Float64}
    current_dual_solution::CuArray{Float64}
    current_dual_product::CuArray{Float64}
    current_primal_product::CuArray{Float64}
    buffer_primal_gradient::CuArray{Float64}
    avg_primal_solution::CuArray{Float64}
    avg_primal_gradient::CuArray{Float64}
    avg_dual_solution::CuArray{Float64}
    avg_primal_product::CuArray{Float64}
    sum_primal_solutions::CuArray{Float64}
    sum_dual_solutions::CuArray{Float64}
    sum_primal_product::CuArray{Float64}
    sum_dual_product::CuArray{Float64}
    original_primal_solution::CuArray{Float64}
    original_primal_gradient::CuArray{Float64}
    original_dual_solution::CuArray{Float64}
    original_primal_product::CuArray{Float64}
    buffer_kkt_primal_solution::CuArray{Float64}
    buffer_kkt_primal_product::CuArray{Float64}
    buffer_kkt_lower_variable_violation::CuArray{Float64}
    buffer_kkt_upper_variable_violation::CuArray{Float64}
    buffer_kkt_reduced_costs::CuArray{Float64}
    delta_primal::CuArray{Float64}
    delta_primal_product::CuArray{Float64}
    delta_dual::CuArray{Float64}
end


mutable struct PDLPParams
    ruiz_iterations::Int
    pock_chambolle_alpha::Union{Nothing,Float64}
    scale_initial_primal_weight::Bool
    termination_evaluation_frequency::Int32
    extrapolation_coefficient::Float64
    reduction_exponent::Float64
    growth_exponent::Float64
    kkt_matrix_pass_limit::Float64
    necessary_reduction_for_restart::Float64
    sufficient_reduction_for_restart::Float64
    iteration_limit::Int32
    skip_hard_problems::Bool
    termination_criteria::TerminationCriteria
end

mutable struct PDLPDims
    current_LP_length::Int32
    total_LP_length::Int32
    n_LPs::Int32
    n_vars::Int32
end

mutable struct PDLPSparsity
    nz_count::Vector{Int}
    nz_rows::CuArray{Int32}
    nz_cols::CuArray{Int32}
end

function sparse_constructor(mat::Matrix{Bool})
    nz_count = zeros(Int, size(mat, 1))
    nz_rows = Int32[]
    nz_cols = Int32[]
    for row = 1:size(mat,1)
        if row != 1
            nz_count[row] = nz_count[row-1]
        end
        for col = 1:size(mat,2)
            if mat[row, col]
                nz_count[row] += 1
                push!(nz_rows, Int32(row))
                push!(nz_cols, Int32(col))
            end
        end
    end
    return nz_count, CuArray(nz_rows), CuArray(nz_cols)
end

mutable struct PDLPData
    original_problem::LinearProgramSet
    scaled_problem::LinearProgramSet
    sparsity::PDLPSparsity
    variable_rescaling::CuArray{Float64}
    constraint_rescaling::CuArray{Float64}
    kernel_storage::KernelStorage
    # initial_primal_solution::CuArray{Float64}
    # initial_dual_solution::CuArray{Float64}
    # initial_primal_product::CuArray{Float64}
    # initial_dual_product::CuArray{Float64}
    parameters::PDLPParams
    dims::PDLPDims
    primal_weight::CuArray{Float64}
    step_size::CuArray{Float64}
    termination_reason::CuArray{TerminationReason}
    skip_flag::CuArray{Bool}
    iterations::CuArray{Int32}
    active_constraint::CuArray{Bool}
    global_counter::CuArray{Int32}
    iteration_counter::CuArray{Int32}
end

function PDLPData(
    n_LPs::Int, 
    n_vars::Int, 
    total_LP_length::Int;
    sparsity::Matrix{Bool} = fill(true, total_LP_length, n_vars),
    iteration_limit::Int = Int(typemax(Int32)),
    termination_evaluation_frequency::Int = 64,
    extrapolation_coefficient::Float64 = 1.0,
    reduction_exponent::Float64 = 0.3,
    growth_exponent::Float64 = 0.6,
    kkt_matrix_pass_limit::Float64 = Inf,
    necessary_reduction_for_restart::Float64 = 0.8,
    sufficient_reduction_for_restart::Float64 = 0.2,
    abs_tol::Float64 = 1.0E-8,
    rel_tol::Float64 = 1.0E-8,
    skip_hard_problems::Bool = false,
    )
    # Call the sparse constructor to get sparsity information
    nz, nz_rows, nz_cols = sparse_constructor(sparsity)
    return PDLPData(
        LinearProgramSet( # Original LPs
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, total_LP_length * n_LPs, n_vars),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, n_LPs)),
        LinearProgramSet( # Scaled LPs
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, total_LP_length * n_LPs, n_vars),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, n_LPs)),
        PDLPSparsity(
            nz,
            nz_rows,
            nz_cols),
        CuArray{Float64}(undef, n_LPs, n_vars), # Variable rescaling
        CuArray{Float64}(undef, n_LPs * total_LP_length), # Constraint rescaling
        # CuArray{Float64}(undef, n_LPs, n_vars), # Initial primal solution
        # CuArray{Float64}(undef, n_LPs * total_LP_length), # Initial dual solution
        # CuArray{Float64}(undef, n_LPs * total_LP_length), # Initial primal product
        # CuArray{Float64}(undef, n_LPs, n_vars), # Initial dual product
        KernelStorage(
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, n_LPs, n_vars),
            CUDA.zeros(Float64, total_LP_length * n_LPs),
            CUDA.zeros(Float64, total_LP_length * n_LPs)
        ),
        PDLPParams( # Parameters
            10,                               # Iterations for Ruiz rescaling
            1.0,                              # Alpha for Pock Chambolle rescaling
            true,                             # Scale initial primal weight flag
            termination_evaluation_frequency, # Termination evaluation frequency
            extrapolation_coefficient,        # Extrapolation coefficient used for taking steps
            reduction_exponent,               # Reduction exponent (for step size updates)
            growth_exponent,                  # Growth coefficient (for step size updates)
            kkt_matrix_pass_limit,            # Limit for KKT matrix passes (default Inf)
            necessary_reduction_for_restart,  # Necessary reduction for restart (default 0.8)
            sufficient_reduction_for_restart, # Sufficient reduction for restart (default 0.2)
            iteration_limit,                  # Iteration limit (default typemax(Int32))
            skip_hard_problems,               # Flag to skip problems with too many iterations
            TerminationCriteria(
                abs_tol,          # Absolute tolerance epsilon (default 1E-8)
                rel_tol,          # Relative tolerance epsilon (default 1E-8)
                1.0e-8,           # Primal infeasible epsilon
                1.0e-8,           # Dual infeasible epsilon
                Inf,              # Time limit (s)
                iteration_limit,  # Iteration limit
                kkt_matrix_pass_limit # KKT matrix pass limit,
            )
        ),
        PDLPDims(   
            Int32(0),               # Current LP length
            Int32(total_LP_length), # Total LP length
            Int32(n_LPs),           # Number of LPs
            Int32(n_vars),          # Number of variables
        ),
        CuArray{Float64}(undef, n_LPs), # Primal weights
        CuArray{Float64}(undef, n_LPs), # Step sizes
        CUDA.fill(TERMINATION_REASON_UNSPECIFIED, n_LPs), # Termination reason
        CUDA.fill(false, n_LPs), #Skip flag
        CuArray{Int32}(undef, n_LPs), # Iteration counter
        CuArray{Bool}(undef, total_LP_length * n_LPs), # Active LP constraints
        CuArray{Int32}(undef, 1), # Global counter
        CuArray{Int32}(undef, 1), # Iteration counter
    )
end
