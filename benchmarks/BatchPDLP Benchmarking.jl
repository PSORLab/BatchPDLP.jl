
# Use the necessary packages
using BatchPDLP, SourceCodeMcCormick, CUDA, JuMP, GLPK, Gurobi, HiGHS, Dates, Random
const GRB_ENV = Gurobi.Env()

# Include all files from the master list
include("master_list.jl")

# Using the number of problems and the dimensionality, create a partition of the 
# space of the original problem
function case_generator(lbd, ubd, n, m)
    all_lvbs = zeros(m, n)
    all_uvbs = zeros(m, n)
    all_lvbs[1,:] .= lbd
    all_uvbs[1,:] .= ubd

    # Deal with +/-Inf by setting them to +/- 1e6
    for i in eachindex(lbd)
        if all_lvbs[1,i]==Inf
            all_lvbs[1,i] = 1E6
        elseif all_lvbs[1,i]==-Inf
            all_lvbs[1,i] = -1E6
        end
        if all_uvbs[1,i]==Inf
            all_uvbs[1,i] = 1E6
        elseif all_uvbs[1,i]==-Inf
            all_uvbs[1,i] = -1E6
        end
    end

    var = 1
    curr_len = 1
    while curr_len < m
        # Only make up to m nodes in total, m_new at a time
        m_new = min(m - curr_len, curr_len)
        active_len = curr_len

        for i = 1:m_new
            # Create a copy of the i'th set of bounds
            all_lvbs[curr_len+i,:] .= all_lvbs[i,:]
            all_uvbs[curr_len+i,:] .= all_uvbs[i,:]

            # Get the midpoint of this subdomain in the `var` dimension
            midpoint = (all_uvbs[i,var] + all_lvbs[i,var])/2

            # Adjust the original and copy of the i'th node to contain
            # the two halves of variable `var`, split at the midpoint
            all_uvbs[i,var] = midpoint
            all_lvbs[curr_len+i,var] = midpoint
            
            # Increment the current length tracker
            active_len += 1
        end
        curr_len = active_len

        # Increment the variable index being split, up to varmax
        var < n ? (var += 1) : (var = 1)
    end
    return all_lvbs, all_uvbs
end

# Now, given a problem, I want to create the right objects to calculate
# relaxations, create the optimization problems, and pass them to both
# PDLP and, e.g., GLPK. 
function run_example(example::LoadedProblem, n_LPs::Int, n_cuts::Int)
    ##############################################################################
    ######### Step 1) Get Bounds        
    ##############################################################################

    # Based on the problem, partition the domain into n_LPs portions.
    lvbs, uvbs = case_generator(example.lvbs, example.uvbs, example.nvars, n_LPs)

    # Set up the bounds as CuArrays
    lvbs_d = CuArray(lvbs)
    uvbs_d = CuArray(uvbs)

    leq_len = length(example.leq_cons)
    geq_len = length(example.geq_cons)
    eq_len = length(example.eq_cons)

    ##############################################################################
    ######### Step 2) Set up the PDLP struct
    ##############################################################################

    # Calculate the sparsity for PDLP
    cut_height = (1 + leq_len + geq_len + 2*eq_len)
    sparsity = zeros(Bool, 1 + n_cuts*cut_height, example.nvars+1)

    # Lower bound, epigraph variable is always 1, other variables are 0
    sparsity[1, 1] = true

    # Objective function
    for i = 1:n_cuts
        sparsity[2 + (i-1)*cut_height, :] .= [true, example.obj_sp...]
    end

    # LEQ constraints
    for i = 1:n_cuts
        for j = 1:leq_len
            sparsity[2 + j + (i-1)*cut_height, :] .= [false, example.leq_sp[j]...]
        end
    end

    # GEQ constraints
    for i = 1:n_cuts
        for j = 1:geq_len
            sparsity[2 + leq_len + j + (i-1)*cut_height, :] .= [false, example.geq_sp[j]...]
        end
    end

    # EQ constraints
    for i = 1:n_cuts
        for j = 1:eq_len
            sparsity[2 + leq_len + geq_len + (2j-1) + (i-1)*cut_height, :] .= [false, example.eq_sp[j]...]
            sparsity[2 + leq_len + geq_len + (2j-0) + (i-1)*cut_height, :] .= [false, example.eq_sp[j]...]
        end
    end

    # Feed in data to the PDLPData struct
    PDLP_data = PDLPData(n_LPs, example.nvars+1, 1+n_cuts*cut_height, sparsity=sparsity, iteration_limit=1000000)

    
    ##############################################################################
    ######### Step 3) Evaluate terms using SCMC
    ##############################################################################

    # Fill in eval_points very simply, just doing a straight line through all
    # the variables based on the number of cuts
    eval_points = CUDA.zeros(Float64, n_LPs*n_cuts, example.nvars)
    
    points = collect(0:(1/(n_cuts+1)):1)
    for i=1:n_cuts
        @. eval_points[(i-1)*n_LPs + 1 : i*n_LPs, :] = (lvbs_d + uvbs_d).*points[i+1]
    end

    # The objective function and constraints take variables as inputs, where each
    # variable needs its eval point, the lower bound, and the upper bound. Create those.
    input_storage = [CuArray{Float64}(undef, n_LPs*n_cuts, 3) for i=1:example.nvars]
    for i = 1:example.nvars
        input_storage[i] = hcat(eval_points[:,i], repeat(lvbs_d[:,i], n_cuts), repeat(uvbs_d[:,i], n_cuts))
    end

    # Calculate relaxations for the objective function and constraints
    result_width = 4 + 2*example.nvars
    obj_result_storage = CUDA.zeros(Float64, n_LPs*n_cuts, result_width)
    leq_result_storage = [CUDA.zeros(Float64, n_LPs*n_cuts, result_width) for i = 1:leq_len]
    geq_result_storage = [CUDA.zeros(Float64, n_LPs*n_cuts, result_width) for i = 1:geq_len]
    eq_result_storage = [CUDA.zeros(Float64, n_LPs*n_cuts, result_width) for i = 1:eq_len]

    # Objective function
    CUDA.@sync example.obj_fun(obj_result_storage, [input_storage[j] for j=1:example.nvars]...)

    # LEQ constraints
    for i in 1:leq_len
        CUDA.@sync example.leq_cons[i](leq_result_storage[i], [input_storage[j] for j=1:example.nvars]...)
    end
    
    # GEQ constraints
    for i in 1:geq_len
        CUDA.@sync example.geq_cons[i](geq_result_storage[i], [input_storage[j] for j=1:example.nvars]...)
    end
    
    # EQ constraints
    for i in 1:eq_len
        CUDA.@sync example.eq_cons[i](eq_result_storage[i], [input_storage[j] for j=1:example.nvars]...)
    end

    # Since the PDLP_GPU struct already has an allocated field for PDLP data, we
    # only need to update that field with the new LP to solve. The main PDLP
    # algorithm already resets all fields except for the original LP, so to prepare
    # for this next PDLP run, all we need to do is update the original_problem
    # field.
    LPs = PDLP_data.original_problem

    LPs.variable_lower_bounds .= hcat(CUDA.fill(-Inf, n_LPs), lvbs_d) # The left-most column is for the epigraph variable
    LPs.variable_upper_bounds .= hcat(CUDA.fill( Inf, n_LPs), uvbs_d)
    CUDA.fill!(LPs.constraint_matrix, 0.0)
    CUDA.fill!(LPs.right_hand_side, 0.0)
    LPs.objective_vector .= hcat(CUDA.ones(Float64, n_LPs), CUDA.zeros(Float64, n_LPs, example.nvars)) # Minimize the epigraph variable
    CUDA.fill!(LPs.objective_constant, 0.0)
    CUDA.fill!(PDLP_data.active_constraint, false)

    PDLP_data.dims.current_LP_length = Int32(0) # Zero, until we add constraints
    # PDLP_data.dims.total_LP_length #unchanged
    PDLP_data.dims.n_LPs = Int32(n_LPs)
    # PDLP_data.dims.n_vars = Int32(varmax) #unchanged


    # Now that the original_problem field is fully reset, we can begin filling in constraints.

    # First, add the lower bound from the inclusion monotonic interval extension as a constraint. 
    # This has a specialized (simpler) function call
    add_LP_lower_bound(LPs, obj_result_storage[1:n_LPs,:], PDLP_data.active_constraint, PDLP_data.dims)

    # For each of the cuts...
    for i = 1:n_cuts
        # Next, add in the objective function's relaxation at the midpoint as a constraint
        add_LP_objective_constraint(LPs, obj_result_storage[(i-1)*n_LPs+1:i*n_LPs,:], eval_points[(i-1)*n_LPs+1:i*n_LPs,:], PDLP_data.active_constraint, PDLP_data.dims)
        
        # Add in each of the regular constraints
        # LEQ constraints
        for j in 1:leq_len
            add_LP_constraint(LPs, leq_result_storage[j][(i-1)*n_LPs+1:i*n_LPs,:], eval_points[(i-1)*n_LPs+1:i*n_LPs,:], PDLP_data.active_constraint, PDLP_data.dims, geq=false)
        end

        # GEQ constraints
        for j in 1:geq_len
            add_LP_constraint(LPs, geq_result_storage[j][(i-1)*n_LPs+1:i*n_LPs,:], eval_points[(i-1)*n_LPs+1:i*n_LPs,:], PDLP_data.active_constraint, PDLP_data.dims, geq=true)
        end
        
        # EQ constraints
        for j in 1:eq_len
            add_LP_constraint(LPs, eq_result_storage[j][(i-1)*n_LPs+1:i*n_LPs,:], eval_points[(i-1)*n_LPs+1:i*n_LPs,:], PDLP_data.active_constraint, PDLP_data.dims, geq=false)
            add_LP_constraint(LPs, eq_result_storage[j][(i-1)*n_LPs+1:i*n_LPs,:], eval_points[(i-1)*n_LPs+1:i*n_LPs,:], PDLP_data.active_constraint, PDLP_data.dims, geq=true)
        end

    end

    PDLP_setup_time = time() - start_time

    # Solve the problem set with PDLP
    PDLP_solutions = CUDA.zeros(Float64, n_LPs, example.nvars+1)
    PDLP_objectives = CUDA.zeros(Float64, n_LPs, 1)
    PDLP_lowres_solutions = CUDA.zeros(Float64, n_LPs, example.nvars+1)
    PDLP_lowres_objectives = CUDA.zeros(Float64, n_LPs, 1)
    

    # display(sparsity)

    println("Example: $(example.name)")
    println("Base problem: $(example.nvars) variables, $(example.ncons) constraints")
    println("There are $leq_len LEQ constraints, $geq_len GEQ constraints, and $eq_len EQ constraints.")
    println("Doing $n_cuts cuts, so each LP has $(PDLP_data.dims.total_LP_length) constraints")

    # Solve it once for compilation
    try
        PDLP(PDLP_data, solutions=PDLP_solutions, objectives=PDLP_objectives)
    catch
        solving_time = NaN
    end
    device!(0)

    # Solve the problem once at higher tolerances (default is 1E-8 for abs,rel,primal/dual infeas)
    println("Starting PDLP")
    GC.gc()
    GC.enable(false)

    # Timers
    PDLP_solving_times = zeros(3)
    PDLP_lowres_solving_times = zeros(3)
    PDLP_95pct_solving_times = zeros(3)

    
    # Solve for the first time, but only do 95% of problems 
    # (then reset and run all the problems as normal)
    iterations = Array(PDLP_data.iterations)
    stop_here = 0
    for i in sort(unique(iterations))
        if count(iterations .<= i) >= 9500
            stop_here = i
            break
        end
    end
    PDLP_data.parameters.iteration_limit = Int32(stop_here)
    for attempt = 1:3
        try
            PDLP_95pct_solving_times[attempt] = @elapsed PDLP(PDLP_data, solutions=PDLP_solutions, objectives=PDLP_objectives)
        catch
            PDLP_95pct_solving_times[attempt] = NaN
        end
    end
    PDLP_data.parameters.iteration_limit = Int32(1000000)

    PDLP_correct = 0
    PDLP_lowres_correct = 0
    PDLP_iteration_limit_reached = 0
    for attempt = 1:3
        try
            device!(0)
            PDLP_solving_times[attempt] = @elapsed PDLP(PDLP_data, solutions=PDLP_solutions, objectives=PDLP_objectives)
        catch
            PDLP_solving_times[attempt] = NaN
        end
    end
    device!(0)

    # Check the termination reasons. If the iteration limit was reached, add one to that counter
    PDLP_iteration_limit_reached += sum(Array(PDLP_data.termination_reason).==TERMINATION_REASON_ITERATION_LIMIT)

    # Save solutions as regular arrays for comparisons
    PDLP_termination_array = Array(PDLP_data.termination_reason)
    PDLP_objectives_array = Array(PDLP_objectives)

    # Now solve it again using much lower tolerances
    PDLP_data.parameters.termination_criteria.eps_optimal_absolute = 1E-4
    PDLP_data.parameters.termination_criteria.eps_optimal_relative = 1E-4
    PDLP_data.parameters.termination_criteria.eps_primal_infeasible= 1E-8
    PDLP_data.parameters.termination_criteria.eps_dual_infeasible  = 1E-8

    # Solve again with lower tolerances
    for attempt = 1:3
        try
            PDLP_lowres_solving_times[attempt] = @elapsed PDLP(PDLP_data, solutions=PDLP_lowres_solutions, objectives=PDLP_lowres_objectives)
        catch
            PDLP_lowres_solving_times[attempt] = NaN
        end
    end
    GC.enable(true)

    PDLP_solving_time = minimum(PDLP_solving_times)
    PDLP_lowres_solving_time = minimum(PDLP_lowres_solving_times)
    PDLP_95pct_solving_time = minimum(PDLP_95pct_solving_times)

    println("And how do times compare:")
    println("Single GPU:  | $(round.(PDLP_solving_times, digits=6))")
    println("==============================")
    println("What if we only solve 95% of problems:")
    println("Single GPU:  | $(round.(PDLP_95pct_solving_times, digits=6))")

    # Save constraints and the RHS for CPU solvers
    constraint_mat = Array(LPs.constraint_matrix)
    rhs = Array(LPs.right_hand_side)

    # Save lowres solutions as regular arrays for comparisons
    PDLP_lowres_termination_array = Array(PDLP_data.termination_reason)
    PDLP_lowres_objectives_array = Array(PDLP_lowres_objectives)

    #############################################################
    ##################### Other Solvers #########################
    #############################################################
    n_solvers = 11
    
    # Now, we wish to create optimization problems and time how long it takes
    # to solve them using the commercial optimizers.
    setup_times = zeros(Float64, n_solvers)
    hires_solving_times = zeros(Float64, n_solvers, 3)
    hires_termination_statuses = fill(OPTIMIZE_NOT_CALLED, n_solvers)
    hires_objectives = zeros(Float64, n_solvers)
    hires_correct = zeros(Int, n_solvers)
    lowres_solving_times = zeros(Float64, n_solvers, 3)
    lowres_termination_statuses = fill(OPTIMIZE_NOT_CALLED, n_solvers)
    lowres_objectives = zeros(Float64, n_solvers)
    lowres_correct = zeros(Int, n_solvers)

    println("Starting Other Solvers")
    
    # Create continue flags, in case solvers run into errors
    hires_continue = ones(Bool, n_solvers)
    lowres_continue = ones(Bool, n_solvers)

    # Initialize all solvers
    setup_times[1] += @elapsed m_GLPK_pSimplex = GLPK.Optimizer(method=GLPK.SIMPLEX)
    setup_times[1] += @elapsed MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tol_bnd"), 1E-8) # Primal feasibility tolerance
    setup_times[1] += @elapsed MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tol_dj"), 1E-8) # Dual feasibility tolerance
    setup_times[1] += @elapsed MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("meth"), GLPK.GLP_PRIMAL)
    setup_times[1] += @elapsed MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time limit
    MOI.set(m_GLPK_pSimplex, MOI.Silent(), true)

    setup_times[2] += @elapsed m_GLPK_dSimplex = GLPK.Optimizer(method=GLPK.SIMPLEX)
    setup_times[2] += @elapsed MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tol_bnd"), 1E-8) # Primal feasibility tolerance
    setup_times[2] += @elapsed MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tol_dj"), 1E-8) # Dual feasibility tolerance
    setup_times[2] += @elapsed MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("meth"), GLPK.GLP_DUAL)
    setup_times[2] += @elapsed MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time 
    MOI.set(m_GLPK_dSimplex, MOI.Silent(), true)

    setup_times[3] += @elapsed m_GLPK_Interior = GLPK.Optimizer(method=GLPK.INTERIOR)
    setup_times[3] += @elapsed MOI.set(m_GLPK_Interior, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time limit
    MOI.set(m_GLPK_Interior, MOI.Silent(), true)
 
    setup_times[4] += @elapsed m_Gurobi_pSimplex = Gurobi.Optimizer(GRB_ENV)
    setup_times[4] += @elapsed MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    setup_times[4] += @elapsed MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    setup_times[4] += @elapsed MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    setup_times[4] += @elapsed MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("Method"), 0) # Primal simplex
    setup_times[4] += @elapsed MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    MOI.set(m_Gurobi_pSimplex, MOI.Silent(), true)

    setup_times[5] += @elapsed m_Gurobi_dSimplex = Gurobi.Optimizer(GRB_ENV)
    setup_times[5] += @elapsed MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    setup_times[5] += @elapsed MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    setup_times[5] += @elapsed MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    setup_times[5] += @elapsed MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("Method"), 1) # Dual simplex
    setup_times[5] += @elapsed MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    MOI.set(m_Gurobi_dSimplex, MOI.Silent(), true)

    setup_times[6] += @elapsed m_Gurobi_Interior = Gurobi.Optimizer(GRB_ENV)
    setup_times[6] += @elapsed MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    setup_times[6] += @elapsed MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("BarConvTol"), 1E-8) # Relative convergence tolerance for barrier method
    setup_times[6] += @elapsed MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    setup_times[6] += @elapsed MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    setup_times[6] += @elapsed MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("Method"), 2) # Barrier
    setup_times[6] += @elapsed MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    MOI.set(m_Gurobi_Interior, MOI.Silent(), true)

    setup_times[7] += @elapsed m_Gurobi_PDLP = Gurobi.Optimizer(GRB_ENV)
    setup_times[7] += @elapsed MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    setup_times[7] += @elapsed MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    setup_times[7] += @elapsed MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    setup_times[7] += @elapsed MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("Method"), 6) # PDLP (not in documentation)
    setup_times[7] += @elapsed MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    MOI.set(m_Gurobi_PDLP, MOI.Silent(), true)

    setup_times[8] += @elapsed m_HiGHS_pSimplex = HiGHS.Optimizer()
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("solver"), "simplex") # Simplex in general. See next option.
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("simplex_strategy"), 4) # Primal simplex
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-8) 
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-8) 
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    MOI.set(m_HiGHS_pSimplex, MOI.Silent(), true)

    setup_times[9] += @elapsed m_HiGHS_dSimplex = HiGHS.Optimizer()
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("solver"), "simplex") # Simplex in general. See next option.
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("simplex_strategy"), 1) # Dual simplex
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-8) 
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-8) 
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    MOI.set(m_HiGHS_dSimplex, MOI.Silent(), true)

    setup_times[10] += @elapsed m_HiGHS_Interior = HiGHS.Optimizer()
    setup_times[10] += @elapsed MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("solver"), "ipm") # Interior point method
    setup_times[10] += @elapsed MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-8) 
    setup_times[10] += @elapsed MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-8) 
    setup_times[10] += @elapsed MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    setup_times[10] += @elapsed MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    setup_times[10] += @elapsed MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    MOI.set(m_HiGHS_Interior, MOI.Silent(), true)

    setup_times[11] += @elapsed m_HiGHS_PDLP = HiGHS.Optimizer()
    setup_times[11] += @elapsed MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("solver"), "pdlp") # PDLP
    setup_times[11] += @elapsed MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-8) 
    setup_times[11] += @elapsed MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-8) 
    setup_times[11] += @elapsed MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    setup_times[11] += @elapsed MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    setup_times[11] += @elapsed MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    MOI.set(m_HiGHS_PDLP, MOI.Silent(), true)
                         
                                      

    # List of all solvers being used
    solver_list = [m_GLPK_pSimplex, m_GLPK_dSimplex, m_GLPK_Interior,
                   m_Gurobi_pSimplex, m_Gurobi_dSimplex, m_Gurobi_Interior, m_Gurobi_PDLP,
                   m_HiGHS_pSimplex, m_HiGHS_dSimplex, m_HiGHS_Interior, m_HiGHS_PDLP]

    # Compile solvers
    for i in eachindex(solver_list)
        try
            m = solver_list[i]
            epi = MOI.add_variable(m)
            vi = Vector{MOI.VariableIndex}(undef, example.nvars)
            for j = 1:example.nvars
                vi[j], (_,_) = MOI.add_constrained_variable(m, (MOI.GreaterThan(lvbs[1,j]), MOI.LessThan(uvbs[1,j])))
            end

            for j = 1:PDLP_data.dims.current_LP_length
                start = (1-1)*PDLP_data.dims.total_LP_length
                @views MOI.add_constraint(m, constraint_mat[start+j,1]*epi + 
                                            sum(constraint_mat[start+j,2:end].*vi[1:example.nvars]),
                                            MOI.GreaterThan(rhs[start+j]))
            end

            MOI.set(m, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), 0.0+epi)

            MOI.optimize!(m)
        catch
            hires_continue[i] = false
            lowres_continue[i] = false
        end
    end
    # println("Solvers compiled $(time() - big_time)")

    # Now solve all the LPs
    for attempt = 1:3
        for i = 1:n_LPs
            # Reset statuses
            hires_termination_statuses = fill(OPTIMIZE_NOT_CALLED, n_solvers)
            hires_objectives = zeros(Float64, n_solvers)
            lowres_termination_statuses = fill(OPTIMIZE_NOT_CALLED, n_solvers)
            lowres_objectives = zeros(Float64, n_solvers)
            
            for solver = 1:length(solver_list)
                current_solver = solver_list[solver]
                if hires_continue[solver]
                    try
                        # Empty the solver of previous run information
                        MOI.empty!(current_solver)

                        # Add the epigraph variable
                        setup_times[solver] += @elapsed epi = MOI.add_variable(current_solver)

                        # Add other problem variables
                        setup_times[solver] += @elapsed begin
                            vi = Vector{MOI.VariableIndex}(undef, example.nvars)
                            for j = 1:example.nvars
                                vi[j], (_,_) = MOI.add_constrained_variable(current_solver, (MOI.GreaterThan(lvbs[i,j]), MOI.LessThan(uvbs[i,j])))
                            end
                        end

                        # Add constraints
                        for j = 1:PDLP_data.dims.current_LP_length
                            start = (i-1)*PDLP_data.dims.total_LP_length
                            constraint_section = constraint_mat[start+j,1:end]
                            rhs_section = rhs[start+j]
                            setup_times[solver] += @elapsed MOI.add_constraint(current_solver, constraint_section[1]*epi + 
                                                        sum(constraint_section[2:end].*vi[1:example.nvars]),
                                                        MOI.GreaterThan(rhs_section))
                        end

                        # Add the objective function (always already in epigraph form)
                        MOI.set(current_solver, MOI.ObjectiveSense(), MOI.MIN_SENSE)
                        MOI.set(current_solver, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), 0.0+epi)

                    
                        # Optimize, and add the time to the solving time
                        hires_solving_times[solver, attempt] += @elapsed(MOI.optimize!(current_solver))

                        # If it solved, save the termination status and objective value
                        hires_termination_statuses[solver] = MOI.get(current_solver, MOI.TerminationStatus())
                        if hires_termination_statuses[solver] == OPTIMAL
                            hires_objectives[solver] = MOI.get(current_solver, MOI.ObjectiveValue())
                        elseif hires_termination_statuses[solver] == OTHER_ERROR ||
                            hires_termination_statuses[solver] == NUMERICAL_ERROR ||
                            hires_termination_statuses[solver] == ITERATION_LIMIT ||
                            hires_termination_statuses[solver] == TIME_LIMIT
                            hires_objectives[solver] = NaN
                        else
                            hires_objectives[solver] = Inf
                        end
                    catch
                        hires_continue[solver] = false
                        hires_solving_times[solver, attempt] = NaN
                        hires_termination_statuses[solver] = OTHER_ERROR
                        hires_objectives[solver] = NaN
                    end
                else
                    hires_solving_times[solver, attempt] = NaN
                    hires_termination_statuses[solver] = OTHER_ERROR
                    hires_objectives[solver] = NaN
                end

                # If the total time is above an hour, just stop doing this solver
                if hires_solving_times[solver, attempt] > 3600.0
                    hires_solving_times[solver, attempt] = 3600.0
                    hires_continue[solver] = false
                end
            end

            # Determine what the "correct" answer is and count up corrects
            if attempt==1
                PDLP_correct += answer_verification(
                    hires_correct, 
                    hires_termination_statuses, 
                    hires_objectives, 
                    PDLP_termination_array[i],
                    PDLP_objectives_array[i],
                    atol=1e-6, 
                    rtol=1e-6,
                    )
            end
        end
    end

    println("Starting Other Solvers (lowres)")
    # No way to set GLPK abs/rel convergence tolerance!
    m_GLPK_pSimplex = GLPK.Optimizer(method=GLPK.SIMPLEX)
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tol_bnd"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tol_dj"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("meth"), GLPK.GLP_PRIMAL)
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time limit
    m_GLPK_dSimplex = GLPK.Optimizer(method=GLPK.SIMPLEX)
    MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tol_bnd"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tol_dj"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("meth"), GLPK.GLP_DUAL)
    MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time limit
    m_GLPK_Interior = GLPK.Optimizer(method=GLPK.INTERIOR)
    MOI.set(m_GLPK_Interior, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time limit
 
    # Gurobi also doesn't really have ways to set convergence tolerances
    # NOTE: It should be possible to set the absolute/relative/convergence tolerances for PDLP, but
    #       it is currently undocumented. Lu2025 (cuPDLPx) suggests that the parameters "GURO_PAR_PDHGABSTOL",
    #       "GURO_PAR_PDHGCONVTOL", and "GURO_PAR_PDHGRELTOL" exist for this purpose, but I can't figure
    #       out how to set them.
    m_Gurobi_pSimplex = Gurobi.Optimizer(GRB_ENV)
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("Method"), 0) # Primal simplex
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    m_Gurobi_dSimplex = Gurobi.Optimizer(GRB_ENV)
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("Method"), 1) # Dual simplex
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    m_Gurobi_Interior = Gurobi.Optimizer(GRB_ENV)
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("BarConvTol"), 1E-4) # Relative convergence tolerance for barrier method
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("Method"), 2) # Barrier
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    m_Gurobi_PDLP = Gurobi.Optimizer(GRB_ENV)
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("Method"), 6) # PDLP (not in documentation)
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)

    # HiGHS does have a way to set primal/dual residual tolerances (Not sure if absolute or relative)
    m_HiGHS_pSimplex = HiGHS.Optimizer()
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("solver"), "simplex") # Simplex in general. See next option.
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("simplex_strategy"), 4) # Primal simplex
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    m_HiGHS_dSimplex = HiGHS.Optimizer()
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("solver"), "simplex") # Simplex in general. See next option.
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("simplex_strategy"), 1) # Dual simplex
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    m_HiGHS_Interior = HiGHS.Optimizer()
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("solver"), "ipm") # Interior point method
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    m_HiGHS_PDLP = HiGHS.Optimizer()
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("solver"), "pdlp") # PDLP
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    


    # Same thing but for lower tolerance
    for attempt = 1:3
        for i = 1:n_LPs
            # Reset statuses
            # setup_times = zeros(Float64, n_solvers)
            # hires_solving_times = zeros(Float64, n_solvers)
            hires_termination_statuses = fill(OPTIMIZE_NOT_CALLED, n_solvers)
            hires_objectives = zeros(Float64, n_solvers)
            # hires_correct = zeros(Int, n_solvers)
            # lowres_solving_times = zeros(Float64, n_solvers)
            lowres_termination_statuses = fill(OPTIMIZE_NOT_CALLED, n_solvers)
            lowres_objectives = zeros(Float64, n_solvers)
            # lowres_correct = zeros(Int, n_solvers)

            for solver = 1:length(solver_list)
                current_solver = solver_list[solver]
                if lowres_continue[solver]
                    try
                        # Empty the solver of previous run information
                        MOI.empty!(current_solver)

                        # Add the epigraph variable
                        setup_times[solver] += @elapsed epi = MOI.add_variable(current_solver)

                        # Add other problem variables
                        setup_times[solver] += @elapsed begin
                            vi = Vector{MOI.VariableIndex}(undef, example.nvars)
                            for j = 1:example.nvars
                                vi[j], (_,_) = MOI.add_constrained_variable(current_solver, (MOI.GreaterThan(lvbs[i,j]), MOI.LessThan(uvbs[i,j])))
                            end
                        end

                        # Add constraints
                        for j = 1:PDLP_data.dims.current_LP_length
                            start = (i-1)*PDLP_data.dims.total_LP_length
                            constraint_section = constraint_mat[start+j,1:end]
                            rhs_section = rhs[start+j]
                            setup_times[solver] += @elapsed MOI.add_constraint(current_solver, constraint_section[1]*epi + 
                                                        sum(constraint_section[2:end].*vi[1:example.nvars]),
                                                        MOI.GreaterThan(rhs_section))
                        end

                        # Add the objective function (always already in epigraph form)
                        MOI.set(current_solver, MOI.ObjectiveSense(), MOI.MIN_SENSE)
                        MOI.set(current_solver, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), 0.0+epi)
                        
                        # Optimize, and add the time to the solving time
                        lowres_solving_times[solver, attempt] += @elapsed(MOI.optimize!(current_solver))

                        # If it solved, save the termination status and objective value
                        lowres_termination_statuses[solver] = MOI.get(current_solver, MOI.TerminationStatus())
                        if lowres_termination_statuses[solver] == OPTIMAL
                            lowres_objectives[solver] = MOI.get(current_solver, MOI.ObjectiveValue())
                        elseif lowres_termination_statuses[solver] == OTHER_ERROR ||
                            lowres_termination_statuses[solver] == NUMERICAL_ERROR ||
                            lowres_termination_statuses[solver] == ITERATION_LIMIT ||
                            lowres_termination_statuses[solver] == TIME_LIMIT
                            lowres_objectives[solver] = NaN
                        else
                            lowres_objectives[solver] = Inf
                        end
                    catch
                        lowres_continue[solver] = false
                        lowres_solving_times[solver, attempt] = NaN
                        lowres_termination_statuses[solver] = OTHER_ERROR
                        lowres_objectives[solver] = NaN
                    end
                else
                    lowres_solving_times[solver, attempt] = NaN
                    lowres_termination_statuses[solver] = OTHER_ERROR
                    lowres_objectives[solver] = NaN
                end

                # If the total time is above an hour, just stop doing this solver
                if lowres_solving_times[solver, attempt] > 3600.0
                    lowres_solving_times[solver, attempt] = 3600.0
                    lowres_continue[solver] = false
                end
                # println("Solver $solver done (lowres) ($(time() - big_time))")
            end

            # Determine what the "correct" answer is and count up corrects
            if attempt==1
                PDLP_lowres_correct += answer_verification(
                    lowres_correct, 
                    lowres_termination_statuses, 
                    lowres_objectives, 
                    PDLP_lowres_termination_array[i],
                    PDLP_lowres_objectives_array[i],
                    atol=1e-2, 
                    rtol=1e-2,
                    )
            end
            # println("Answers verified (lowres) ($(time() - big_time))")
        end
    end

    println("")
    return hcat(PDLP_setup_time, PDLP_correct, PDLP_iteration_limit_reached, PDLP_solving_time, PDLP_95pct_solving_time, PDLP_lowres_correct, PDLP_lowres_solving_time,  
                setup_times[1],  hires_correct[1], minimum(hires_solving_times[1,:]), lowres_correct[1], minimum(lowres_solving_times[1,:]),
                setup_times[2],  hires_correct[2], minimum(hires_solving_times[2,:]), lowres_correct[2], minimum(lowres_solving_times[2,:]),
                setup_times[3],  hires_correct[3], minimum(hires_solving_times[3,:]), lowres_correct[3], minimum(lowres_solving_times[3,:]),
                setup_times[4],  hires_correct[4], minimum(hires_solving_times[4,:]), lowres_correct[4], minimum(lowres_solving_times[4,:]),
                setup_times[5],  hires_correct[5], minimum(hires_solving_times[5,:]), lowres_correct[5], minimum(lowres_solving_times[5,:]),
                setup_times[6],  hires_correct[6], minimum(hires_solving_times[6,:]), lowres_correct[6], minimum(lowres_solving_times[6,:]),
                setup_times[7],  hires_correct[7], minimum(hires_solving_times[7,:]), lowres_correct[7], minimum(lowres_solving_times[7,:]),
                setup_times[8],  hires_correct[8], minimum(hires_solving_times[8,:]), lowres_correct[8], minimum(lowres_solving_times[8,:]),
                setup_times[9],  hires_correct[9], minimum(hires_solving_times[9,:]), lowres_correct[9], minimum(lowres_solving_times[9,:]),
                setup_times[10], hires_correct[10], minimum(hires_solving_times[10,:]), lowres_correct[10], minimum(lowres_solving_times[10,:]),
                setup_times[11], hires_correct[11], minimum(hires_solving_times[11,:]), lowres_correct[11], minimum(lowres_solving_times[11,:]),
                )
end

# A version of run_example that uses 2 GPUs
function run_example_2GPU(example::LoadedProblem, n_LPs::Int, n_cuts::Int)
    ##############################################################################
    ######### Step 1) Get Bounds        
    ##############################################################################

    # Based on the problem, partition the domain into n_LPs portions.
    lvbs, uvbs = case_generator(example.lvbs, example.uvbs, example.nvars, n_LPs)

    # Set up the bounds as CuArrays
    lvbs_d = CuArray(lvbs)
    uvbs_d = CuArray(uvbs)

    leq_len = length(example.leq_cons)
    geq_len = length(example.geq_cons)
    eq_len = length(example.eq_cons)

    ##############################################################################
    ######### Step 2) Set up the PDLP struct
    ##############################################################################

    # Calculate the sparsity for PDLP
    cut_height = (1 + leq_len + geq_len + 2*eq_len)
    sparsity = zeros(Bool, 1 + n_cuts*cut_height, example.nvars+1)

    # Lower bound, epigraph variable is always 1, other variables are 0
    sparsity[1, 1] = true

    # Objective function
    for i = 1:n_cuts
        sparsity[2 + (i-1)*cut_height, :] .= [true, example.obj_sp...]
    end

    # LEQ constraints
    for i = 1:n_cuts
        for j = 1:leq_len
            sparsity[2 + j + (i-1)*cut_height, :] .= [false, example.leq_sp[j]...]
        end
    end

    # GEQ constraints
    for i = 1:n_cuts
        for j = 1:geq_len
            sparsity[2 + leq_len + j + (i-1)*cut_height, :] .= [false, example.geq_sp[j]...]
        end
    end

    # EQ constraints
    for i = 1:n_cuts
        for j = 1:eq_len
            sparsity[2 + leq_len + geq_len + (2j-1) + (i-1)*cut_height, :] .= [false, example.eq_sp[j]...]
            sparsity[2 + leq_len + geq_len + (2j-0) + (i-1)*cut_height, :] .= [false, example.eq_sp[j]...]
        end
    end

    # Feed in data to the PDLPData struct
    device!(0)
    PDLP_data = PDLPData(n_LPs, example.nvars+1, 1+n_cuts*cut_height, sparsity=sparsity, iteration_limit=1000000)

    # Make two more PDLPData structs for the 2-GPU setup
    PDLP_data_A = PDLPData(Int(n_LPs/2), example.nvars+1, 1+n_cuts*cut_height, sparsity=sparsity, iteration_limit=1000000)
    device!(1)
    PDLP_data_B = PDLPData(Int(n_LPs/2), example.nvars+1, 1+n_cuts*cut_height, sparsity=sparsity, iteration_limit=1000000)

    # Switch back to the first GPU
    device!(0)
    
    ##############################################################################
    ######### Step 3) Evaluate terms using SCMC
    ##############################################################################

    # Fill in eval_points very simply, just doing a straight line through all
    # the variables based on the number of cuts
    eval_points = CUDA.zeros(Float64, n_LPs*n_cuts, example.nvars)
    
    points = collect(0:(1/(n_cuts+1)):1)
    for i=1:n_cuts
        @. eval_points[(i-1)*n_LPs + 1 : i*n_LPs, :] = (lvbs_d + uvbs_d).*points[i+1]
    end

    # The objective function and constraints take variables as inputs, where each
    # variable needs its eval point, the lower bound, and the upper bound. Create those.
    input_storage = [CuArray{Float64}(undef, n_LPs*n_cuts, 3) for i=1:example.nvars]
    for i = 1:example.nvars
        input_storage[i] = hcat(eval_points[:,i], repeat(lvbs_d[:,i], n_cuts), repeat(uvbs_d[:,i], n_cuts))
    end

    # Calculate relaxations for the objective function and constraints
    result_width = 4 + 2*example.nvars
    obj_result_storage = CUDA.zeros(Float64, n_LPs*n_cuts, result_width)
    leq_result_storage = [CUDA.zeros(Float64, n_LPs*n_cuts, result_width) for i = 1:leq_len]
    geq_result_storage = [CUDA.zeros(Float64, n_LPs*n_cuts, result_width) for i = 1:geq_len]
    eq_result_storage = [CUDA.zeros(Float64, n_LPs*n_cuts, result_width) for i = 1:eq_len]

    # Objective function
    CUDA.@sync example.obj_fun(obj_result_storage, [input_storage[j] for j=1:example.nvars]...)

    # LEQ constraints
    for i in 1:leq_len
        CUDA.@sync example.leq_cons[i](leq_result_storage[i], [input_storage[j] for j=1:example.nvars]...)
    end
    
    # GEQ constraints
    for i in 1:geq_len
        CUDA.@sync example.geq_cons[i](geq_result_storage[i], [input_storage[j] for j=1:example.nvars]...)
    end
    
    # EQ constraints
    for i in 1:eq_len
        CUDA.@sync example.eq_cons[i](eq_result_storage[i], [input_storage[j] for j=1:example.nvars]...)
    end

    # Since the PDLP_GPU struct already has an allocated field for PDLP data, we
    # only need to update that field with the new LP to solve. The main PDLP
    # algorithm already resets all fields except for the original LP, so to prepare
    # for this next PDLP run, all we need to do is update the original_problem
    # field.
    LPs = PDLP_data.original_problem

    LPs.variable_lower_bounds .= hcat(CUDA.fill(-Inf, n_LPs), lvbs_d) # The left-most column is for the epigraph variable
    LPs.variable_upper_bounds .= hcat(CUDA.fill( Inf, n_LPs), uvbs_d)
    CUDA.fill!(LPs.constraint_matrix, 0.0)
    CUDA.fill!(LPs.right_hand_side, 0.0)
    LPs.objective_vector .= hcat(CUDA.ones(Float64, n_LPs), CUDA.zeros(Float64, n_LPs, example.nvars)) # Minimize the epigraph variable
    CUDA.fill!(LPs.objective_constant, 0.0)
    CUDA.fill!(PDLP_data.active_constraint, false)

    PDLP_data.dims.current_LP_length = Int32(0) # Zero, until we add constraints
    # PDLP_data.dims.total_LP_length #unchanged
    PDLP_data.dims.n_LPs = Int32(n_LPs)
    # PDLP_data.dims.n_vars = Int32(varmax) #unchanged

    # Fill in the 2-GPU PDLP structs similarly
    PDLP_data_A.original_problem.variable_lower_bounds .= hcat(CUDA.fill(-Inf, Int(n_LPs/2)), lvbs_d[1:Int(n_LPs/2),:])
    PDLP_data_A.original_problem.variable_upper_bounds .= hcat(CUDA.fill( Inf, Int(n_LPs/2)), uvbs_d[1:Int(n_LPs/2),:])
    CUDA.fill!(PDLP_data_A.original_problem.constraint_matrix, 0.0)
    CUDA.fill!(PDLP_data_A.original_problem.right_hand_side, 0.0)
    PDLP_data_A.original_problem.objective_vector .= hcat(CUDA.ones(Float64, Int(n_LPs/2)), CUDA.zeros(Float64, Int(n_LPs/2), example.nvars))
    CUDA.fill!(PDLP_data_A.original_problem.objective_constant, 0.0)
    CUDA.fill!(PDLP_data_A.active_constraint, 0.0)
    device!(1)
    PDLP_data_B.original_problem.variable_lower_bounds .= hcat(CUDA.fill(-Inf, Int(n_LPs/2)), lvbs_d[Int(n_LPs/2)+1:end,:])
    PDLP_data_B.original_problem.variable_upper_bounds .= hcat(CUDA.fill( Inf, Int(n_LPs/2)), uvbs_d[Int(n_LPs/2)+1:end,:])
    CUDA.fill!(PDLP_data_B.original_problem.constraint_matrix, 0.0)
    CUDA.fill!(PDLP_data_B.original_problem.right_hand_side, 0.0)
    PDLP_data_B.original_problem.objective_vector .= hcat(CUDA.ones(Float64, Int(n_LPs/2)), CUDA.zeros(Float64, Int(n_LPs/2), example.nvars))
    CUDA.fill!(PDLP_data_B.original_problem.objective_constant, 0.0)
    CUDA.fill!(PDLP_data_B.active_constraint, 0.0)
    device!(0)

    # Now that the original_problem field is fully reset, we can begin filling in constraints.

    # First, add the lower bound from the inclusion monotonic interval extension as a constraint. 
    # This has a specialized (simpler) function call
    add_LP_lower_bound(LPs, obj_result_storage[1:n_LPs,:], PDLP_data.active_constraint, PDLP_data.dims)

    # For each of the cuts...
    for i = 1:n_cuts
        # Next, add in the objective function's relaxation at the midpoint as a constraint
        add_LP_objective_constraint(LPs, obj_result_storage[(i-1)*n_LPs+1:i*n_LPs,:], eval_points[(i-1)*n_LPs+1:i*n_LPs,:], PDLP_data.active_constraint, PDLP_data.dims)
        
        # Add in each of the regular constraints
        # LEQ constraints
        for j in 1:leq_len
            add_LP_constraint(LPs, leq_result_storage[j][(i-1)*n_LPs+1:i*n_LPs,:], eval_points[(i-1)*n_LPs+1:i*n_LPs,:], PDLP_data.active_constraint, PDLP_data.dims, geq=false)
        end

        # GEQ constraints
        for j in 1:geq_len
            add_LP_constraint(LPs, geq_result_storage[j][(i-1)*n_LPs+1:i*n_LPs,:], eval_points[(i-1)*n_LPs+1:i*n_LPs,:], PDLP_data.active_constraint, PDLP_data.dims, geq=true)
        end
        
        # EQ constraints
        for j in 1:eq_len
            add_LP_constraint(LPs, eq_result_storage[j][(i-1)*n_LPs+1:i*n_LPs,:], eval_points[(i-1)*n_LPs+1:i*n_LPs,:], PDLP_data.active_constraint, PDLP_data.dims, geq=false)
            add_LP_constraint(LPs, eq_result_storage[j][(i-1)*n_LPs+1:i*n_LPs,:], eval_points[(i-1)*n_LPs+1:i*n_LPs,:], PDLP_data.active_constraint, PDLP_data.dims, geq=true)
        end

    end

    # Move information into 2-GPU structs
    L = size(PDLP_data_A.original_problem.constraint_matrix, 1)
    device!(0)
    PDLP_data_A.original_problem.constraint_matrix .= copy(PDLP_data.original_problem.constraint_matrix)[1:L,:]
    PDLP_data_A.original_problem.right_hand_side .= copy(PDLP_data.original_problem.right_hand_side)[1:L]
    PDLP_data_A.dims.current_LP_length = PDLP_data.dims.current_LP_length
    PDLP_data_A.dims.n_LPs = Int32(n_LPs/2)
    PDLP_data_A.active_constraint .= copy(PDLP_data.active_constraint[1:L])
    device!(1)
    PDLP_data_B.original_problem.constraint_matrix .= copy(PDLP_data.original_problem.constraint_matrix)[L+1:end,:]
    PDLP_data_B.original_problem.right_hand_side .= copy(PDLP_data.original_problem.right_hand_side)[L+1:end]
    PDLP_data_B.dims.current_LP_length = PDLP_data.dims.current_LP_length
    PDLP_data_B.dims.n_LPs = Int32(n_LPs/2)
    PDLP_data_B.active_constraint .= copy(PDLP_data.active_constraint[L+1:end])
    device!(0)

    PDLP_setup_time = time() - start_time

    # Solve the problem set with PDLP
    PDLP_solutions = CUDA.zeros(Float64, n_LPs, example.nvars+1)
    PDLP_objectives = CUDA.zeros(Float64, n_LPs, 1)
    PDLP_lowres_solutions = CUDA.zeros(Float64, n_LPs, example.nvars+1)
    PDLP_lowres_objectives = CUDA.zeros(Float64, n_LPs, 1)
    
    device!(0)
    PDLP_A_solutions = CUDA.zeros(Float64, Int(n_LPs/2), example.nvars+1)
    PDLP_A_objectives = CUDA.zeros(Float64, Int(n_LPs/2), 1)
    PDLP_A_lowres_solutions = CUDA.zeros(Float64, Int(n_LPs/2), example.nvars+1)
    PDLP_A_lowres_objectives = CUDA.zeros(Float64, Int(n_LPs/2), 1)
    device!(1)
    PDLP_B_solutions = CUDA.zeros(Float64, Int(n_LPs/2), example.nvars+1)
    PDLP_B_objectives = CUDA.zeros(Float64, Int(n_LPs/2), 1)
    PDLP_B_lowres_solutions = CUDA.zeros(Float64, Int(n_LPs/2), example.nvars+1)
    PDLP_B_lowres_objectives = CUDA.zeros(Float64, Int(n_LPs/2), 1)
    device!(0)

    # display(sparsity)

    println("Example: $(example.name)")
    println("Base problem: $(example.nvars) variables, $(example.ncons) constraints")
    println("There are $leq_len LEQ constraints, $geq_len GEQ constraints, and $eq_len EQ constraints.")
    println("Doing $n_cuts cuts, so each LP has $(PDLP_data.dims.total_LP_length) constraints")

    # Solve it once for compilation
    try
        device!(0)
        PDLP(PDLP_data, solutions=PDLP_solutions, objectives=PDLP_objectives)
        @sync begin
            @async begin
                device!(0)
                PDLP(PDLP_data_A, solutions=PDLP_A_solutions, objectives=PDLP_A_objectives)
            end
            @async begin
                device!(1)
                PDLP(PDLP_data_B, solutions=PDLP_B_solutions, objectives=PDLP_B_objectives)
            end
        end
    catch
        solving_time = NaN
    end
    device!(0)

    # Solve the problem once at higher tolerances (default is 1E-8 for abs,rel,primal/dual infeas)
    println("Starting PDLP")
    GC.gc()
    GC.enable(false)

    # Timers
    PDLP_solving_times = zeros(3)
    PDLP_multi_solving_times = zeros(3)
    PDLP_lowres_solving_times = zeros(3)
    PDLP_lowres_multi_solving_times = zeros(3)

    PDLP_95pct_solving_times = zeros(3)
    PDLP_95pct_multi_solving_times = zeros(3)

    
    # Solve for the first time, but only do 95% of problems 
    # (then reset and run all the problems as normal)
    iterations = Array(PDLP_data.iterations)
    stop_here = 0
    for i in sort(unique(iterations))
        if count(iterations .<= i) >= 9500
            stop_here = i
            break
        end
    end
    PDLP_data.parameters.iteration_limit = Int32(stop_here)
    PDLP_data_A.parameters.iteration_limit = Int32(stop_here)
    PDLP_data_B.parameters.iteration_limit = Int32(stop_here)
    for attempt = 1:3
        try
            device!(0)
            PDLP_95pct_solving_times[attempt] = @elapsed PDLP(PDLP_data, solutions=PDLP_solutions, objectives=PDLP_objectives)
        catch
            PDLP_95pct_solving_times[attempt] = NaN
        end
        try 
            device!(1)
            PDLP_95pct_multi_solving_times[attempt] = @elapsed begin
                @sync begin
                    @async begin
                        device!(0)
                        PDLP(PDLP_data_A, solutions=PDLP_A_solutions, objectives=PDLP_A_objectives)
                    end
                    @async begin
                        device!(1)
                        PDLP(PDLP_data_B, solutions=PDLP_B_solutions, objectives=PDLP_B_objectives)
                    end
                end
            end
        catch
            println("Error in one of the other PDLP modes")
        end
    end
    PDLP_data.parameters.iteration_limit = Int32(1000000)
    PDLP_data_A.parameters.iteration_limit = Int32(1000000)
    PDLP_data_B.parameters.iteration_limit = Int32(1000000)

    PDLP_correct = 0
    PDLP_lowres_correct = 0
    PDLP_iteration_limit_reached = 0
    for attempt = 1:3
        try
            device!(0)
            PDLP_solving_times[attempt] = @elapsed PDLP(PDLP_data, solutions=PDLP_solutions, objectives=PDLP_objectives)
        catch
            PDLP_solving_times[attempt] = NaN
        end
        try 
            device!(1)
            PDLP_multi_solving_times[attempt] = @elapsed begin
                @sync begin
                    @async begin
                        device!(0)
                        PDLP(PDLP_data_A, solutions=PDLP_A_solutions, objectives=PDLP_A_objectives)
                    end
                    @async begin
                        device!(1)
                        PDLP(PDLP_data_B, solutions=PDLP_B_solutions, objectives=PDLP_B_objectives)
                    end
                end
            end
        catch
            println("Error in one of the other PDLP modes")
        end
    end
    device!(0)

    # Check the termination reasons. If the iteration limit was reached, add one to that counter
    PDLP_iteration_limit_reached += sum(Array(PDLP_data.termination_reason).==TERMINATION_REASON_ITERATION_LIMIT)

    # Save solutions as regular arrays for comparisons
    PDLP_termination_array = Array(PDLP_data.termination_reason)
    PDLP_objectives_array = Array(PDLP_objectives)

    # Now solve it again using much lower tolerances
    PDLP_data.parameters.termination_criteria.eps_optimal_absolute = 1E-4
    PDLP_data.parameters.termination_criteria.eps_optimal_relative = 1E-4
    PDLP_data.parameters.termination_criteria.eps_primal_infeasible= 1E-8
    PDLP_data.parameters.termination_criteria.eps_dual_infeasible  = 1E-8
    device!(0)
    PDLP_data_A.parameters.termination_criteria.eps_optimal_absolute = 1E-4
    PDLP_data_A.parameters.termination_criteria.eps_optimal_relative = 1E-4
    PDLP_data_A.parameters.termination_criteria.eps_primal_infeasible= 1E-8
    PDLP_data_A.parameters.termination_criteria.eps_dual_infeasible  = 1E-8
    device!(1)
    PDLP_data_B.parameters.termination_criteria.eps_optimal_absolute = 1E-4
    PDLP_data_B.parameters.termination_criteria.eps_optimal_relative = 1E-4
    PDLP_data_B.parameters.termination_criteria.eps_primal_infeasible= 1E-8
    PDLP_data_B.parameters.termination_criteria.eps_dual_infeasible  = 1E-8
    device!(0)

    # Solve again with lower tolerances
    for attempt = 1:3
        try
            device!(0)
            PDLP_lowres_solving_times[attempt] = @elapsed PDLP(PDLP_data, solutions=PDLP_lowres_solutions, objectives=PDLP_lowres_objectives)
        catch
            PDLP_lowres_solving_times[attempt] = NaN
        end
        try 
            device!(1)
            PDLP_lowres_multi_solving_times[attempt] = @elapsed begin
                @sync begin
                    @async begin
                        device!(0)
                        PDLP(PDLP_data_A, solutions=PDLP_A_lowres_solutions, objectives=PDLP_A_lowres_objectives)
                    end
                    @async begin
                        device!(1)
                        PDLP(PDLP_data_B, solutions=PDLP_B_lowres_solutions, objectives=PDLP_B_lowres_objectives)
                    end
                end
            end
        catch
            println("Error in one of the other PDLP modes (lowres)")
        end
    end
    GC.enable(true)
    device!(0)

    PDLP_solving_time = minimum(PDLP_solving_times)
    PDLP_multi_solving_time = minimum(PDLP_multi_solving_times)
    PDLP_lowres_solving_time = minimum(PDLP_lowres_solving_times)
    PDLP_lowres_multi_solving_time = minimum(PDLP_lowres_multi_solving_times)
    PDLP_95pct_solving_time = minimum(PDLP_95pct_solving_times)
    PDLP_95pct_multi_solving_time = minimum(PDLP_95pct_multi_solving_times)

    println("And how do times compare:")
    println("Single GPU:  | $(round.(PDLP_solving_times, digits=6))")
    println("Multi GPU:   | $(round.(PDLP_multi_solving_times, digits=6))")
    println("==============================")
    println("What if we only solve 95% of problems:")
    println("Single GPU:  | $(round.(PDLP_95pct_solving_times, digits=6))")
    println("Multi GPU:   | $(round.(PDLP_95pct_multi_solving_times, digits=6))")

    # Save constraints and the RHS for CPU solvers
    constraint_mat = Array(LPs.constraint_matrix)
    rhs = Array(LPs.right_hand_side)

    # Save lowres solutions as regular arrays for comparisons
    PDLP_lowres_termination_array = Array(PDLP_data.termination_reason)
    PDLP_lowres_objectives_array = Array(PDLP_lowres_objectives)

    #############################################################
    ##################### Other Solvers #########################
    #############################################################
    n_solvers = 11
    
    # Now, we wish to create optimization problems and time how long it takes
    # to solve them using the commercial optimizers.
    setup_times = zeros(Float64, n_solvers)
    hires_solving_times = zeros(Float64, n_solvers, 3)
    hires_termination_statuses = fill(OPTIMIZE_NOT_CALLED, n_solvers)
    hires_objectives = zeros(Float64, n_solvers)
    hires_correct = zeros(Int, n_solvers)
    lowres_solving_times = zeros(Float64, n_solvers, 3)
    lowres_termination_statuses = fill(OPTIMIZE_NOT_CALLED, n_solvers)
    lowres_objectives = zeros(Float64, n_solvers)
    lowres_correct = zeros(Int, n_solvers)

    println("Starting Other Solvers")
    
    # Create continue flags, in case solvers run into errors
    hires_continue = ones(Bool, n_solvers)
    lowres_continue = ones(Bool, n_solvers)

    # Initialize all solvers
    setup_times[1] += @elapsed m_GLPK_pSimplex = GLPK.Optimizer(method=GLPK.SIMPLEX)
    setup_times[1] += @elapsed MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tol_bnd"), 1E-8) # Primal feasibility tolerance
    setup_times[1] += @elapsed MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tol_dj"), 1E-8) # Dual feasibility tolerance
    setup_times[1] += @elapsed MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("meth"), GLPK.GLP_PRIMAL)
    setup_times[1] += @elapsed MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time limit
    MOI.set(m_GLPK_pSimplex, MOI.Silent(), true)

    setup_times[2] += @elapsed m_GLPK_dSimplex = GLPK.Optimizer(method=GLPK.SIMPLEX)
    setup_times[2] += @elapsed MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tol_bnd"), 1E-8) # Primal feasibility tolerance
    setup_times[2] += @elapsed MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tol_dj"), 1E-8) # Dual feasibility tolerance
    setup_times[2] += @elapsed MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("meth"), GLPK.GLP_DUAL)
    setup_times[2] += @elapsed MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time 
    MOI.set(m_GLPK_dSimplex, MOI.Silent(), true)

    setup_times[3] += @elapsed m_GLPK_Interior = GLPK.Optimizer(method=GLPK.INTERIOR)
    setup_times[3] += @elapsed MOI.set(m_GLPK_Interior, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time limit
    MOI.set(m_GLPK_Interior, MOI.Silent(), true)
 
    setup_times[4] += @elapsed m_Gurobi_pSimplex = Gurobi.Optimizer(GRB_ENV)
    setup_times[4] += @elapsed MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    setup_times[4] += @elapsed MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    setup_times[4] += @elapsed MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    setup_times[4] += @elapsed MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("Method"), 0) # Primal simplex
    setup_times[4] += @elapsed MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    MOI.set(m_Gurobi_pSimplex, MOI.Silent(), true)

    setup_times[5] += @elapsed m_Gurobi_dSimplex = Gurobi.Optimizer(GRB_ENV)
    setup_times[5] += @elapsed MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    setup_times[5] += @elapsed MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    setup_times[5] += @elapsed MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    setup_times[5] += @elapsed MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("Method"), 1) # Dual simplex
    setup_times[5] += @elapsed MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    MOI.set(m_Gurobi_dSimplex, MOI.Silent(), true)

    setup_times[6] += @elapsed m_Gurobi_Interior = Gurobi.Optimizer(GRB_ENV)
    setup_times[6] += @elapsed MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    setup_times[6] += @elapsed MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("BarConvTol"), 1E-8) # Relative convergence tolerance for barrier method
    setup_times[6] += @elapsed MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    setup_times[6] += @elapsed MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    setup_times[6] += @elapsed MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("Method"), 2) # Barrier
    setup_times[6] += @elapsed MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    MOI.set(m_Gurobi_Interior, MOI.Silent(), true)

    setup_times[7] += @elapsed m_Gurobi_PDLP = Gurobi.Optimizer(GRB_ENV)
    setup_times[7] += @elapsed MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    setup_times[7] += @elapsed MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    setup_times[7] += @elapsed MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    setup_times[7] += @elapsed MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("Method"), 6) # PDLP (not in documentation)
    setup_times[7] += @elapsed MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    MOI.set(m_Gurobi_PDLP, MOI.Silent(), true)

    setup_times[8] += @elapsed m_HiGHS_pSimplex = HiGHS.Optimizer()
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("solver"), "simplex") # Simplex in general. See next option.
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("simplex_strategy"), 4) # Primal simplex
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-8) 
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-8) 
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    setup_times[8] += @elapsed MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    MOI.set(m_HiGHS_pSimplex, MOI.Silent(), true)

    setup_times[9] += @elapsed m_HiGHS_dSimplex = HiGHS.Optimizer()
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("solver"), "simplex") # Simplex in general. See next option.
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("simplex_strategy"), 1) # Dual simplex
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-8) 
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-8) 
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    setup_times[9] += @elapsed MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    MOI.set(m_HiGHS_dSimplex, MOI.Silent(), true)

    setup_times[10] += @elapsed m_HiGHS_Interior = HiGHS.Optimizer()
    setup_times[10] += @elapsed MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("solver"), "ipm") # Interior point method
    setup_times[10] += @elapsed MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-8) 
    setup_times[10] += @elapsed MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-8) 
    setup_times[10] += @elapsed MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    setup_times[10] += @elapsed MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    setup_times[10] += @elapsed MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    MOI.set(m_HiGHS_Interior, MOI.Silent(), true)

    setup_times[11] += @elapsed m_HiGHS_PDLP = HiGHS.Optimizer()
    setup_times[11] += @elapsed MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("solver"), "pdlp") # PDLP
    setup_times[11] += @elapsed MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-8) 
    setup_times[11] += @elapsed MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-8) 
    setup_times[11] += @elapsed MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    setup_times[11] += @elapsed MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    setup_times[11] += @elapsed MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    MOI.set(m_HiGHS_PDLP, MOI.Silent(), true)
                         
                                      

    # List of all solvers being used
    solver_list = [m_GLPK_pSimplex, m_GLPK_dSimplex, m_GLPK_Interior,
                   m_Gurobi_pSimplex, m_Gurobi_dSimplex, m_Gurobi_Interior, m_Gurobi_PDLP,
                   m_HiGHS_pSimplex, m_HiGHS_dSimplex, m_HiGHS_Interior, m_HiGHS_PDLP]

    # Compile solvers
    for i in eachindex(solver_list)
        try
            m = solver_list[i]
            epi = MOI.add_variable(m)
            vi = Vector{MOI.VariableIndex}(undef, example.nvars)
            for j = 1:example.nvars
                vi[j], (_,_) = MOI.add_constrained_variable(m, (MOI.GreaterThan(lvbs[1,j]), MOI.LessThan(uvbs[1,j])))
            end

            for j = 1:PDLP_data.dims.current_LP_length
                start = (1-1)*PDLP_data.dims.total_LP_length
                @views MOI.add_constraint(m, constraint_mat[start+j,1]*epi + 
                                            sum(constraint_mat[start+j,2:end].*vi[1:example.nvars]),
                                            MOI.GreaterThan(rhs[start+j]))
            end

            MOI.set(m, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), 0.0+epi)

            MOI.optimize!(m)
        catch
            hires_continue[i] = false
            lowres_continue[i] = false
        end
    end
    # println("Solvers compiled $(time() - big_time)")

    # Now solve all the LPs
    for attempt = 1:3
        for i = 1:n_LPs
            # Reset statuses
            hires_termination_statuses = fill(OPTIMIZE_NOT_CALLED, n_solvers)
            hires_objectives = zeros(Float64, n_solvers)
            lowres_termination_statuses = fill(OPTIMIZE_NOT_CALLED, n_solvers)
            lowres_objectives = zeros(Float64, n_solvers)
            
            for solver = 1:length(solver_list)
                current_solver = solver_list[solver]
                if hires_continue[solver]
                    try
                        # Empty the solver of previous run information
                        MOI.empty!(current_solver)

                        # Add the epigraph variable
                        setup_times[solver] += @elapsed epi = MOI.add_variable(current_solver)

                        # Add other problem variables
                        setup_times[solver] += @elapsed begin
                            vi = Vector{MOI.VariableIndex}(undef, example.nvars)
                            for j = 1:example.nvars
                                vi[j], (_,_) = MOI.add_constrained_variable(current_solver, (MOI.GreaterThan(lvbs[i,j]), MOI.LessThan(uvbs[i,j])))
                            end
                        end

                        # Add constraints
                        for j = 1:PDLP_data.dims.current_LP_length
                            start = (i-1)*PDLP_data.dims.total_LP_length
                            constraint_section = constraint_mat[start+j,1:end]
                            rhs_section = rhs[start+j]
                            setup_times[solver] += @elapsed MOI.add_constraint(current_solver, constraint_section[1]*epi + 
                                                        sum(constraint_section[2:end].*vi[1:example.nvars]),
                                                        MOI.GreaterThan(rhs_section))
                        end

                        # Add the objective function (always already in epigraph form)
                        MOI.set(current_solver, MOI.ObjectiveSense(), MOI.MIN_SENSE)
                        MOI.set(current_solver, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), 0.0+epi)

                    
                        # Optimize, and add the time to the solving time
                        hires_solving_times[solver, attempt] += @elapsed(MOI.optimize!(current_solver))

                        # If it solved, save the termination status and objective value
                        hires_termination_statuses[solver] = MOI.get(current_solver, MOI.TerminationStatus())
                        if hires_termination_statuses[solver] == OPTIMAL
                            hires_objectives[solver] = MOI.get(current_solver, MOI.ObjectiveValue())
                        elseif hires_termination_statuses[solver] == OTHER_ERROR ||
                            hires_termination_statuses[solver] == NUMERICAL_ERROR ||
                            hires_termination_statuses[solver] == ITERATION_LIMIT ||
                            hires_termination_statuses[solver] == TIME_LIMIT
                            hires_objectives[solver] = NaN
                        else
                            hires_objectives[solver] = Inf
                        end
                    catch
                        hires_continue[solver] = false
                        hires_solving_times[solver, attempt] = NaN
                        hires_termination_statuses[solver] = OTHER_ERROR
                        hires_objectives[solver] = NaN
                    end
                else
                    hires_solving_times[solver, attempt] = NaN
                    hires_termination_statuses[solver] = OTHER_ERROR
                    hires_objectives[solver] = NaN
                end

                # If the total time is above an hour, just stop doing this solver
                if hires_solving_times[solver, attempt] > 3600.0
                    hires_solving_times[solver, attempt] = 3600.0
                    hires_continue[solver] = false
                end
            end

            # Determine what the "correct" answer is and count up corrects
            if attempt==1
                PDLP_correct += answer_verification(
                    hires_correct, 
                    hires_termination_statuses, 
                    hires_objectives, 
                    PDLP_termination_array[i],
                    PDLP_objectives_array[i],
                    atol=1e-6, 
                    rtol=1e-6,
                    )
            end
        end
    end

    println("Starting Other Solvers (lowres)")
    # No way to set GLPK abs/rel convergence tolerance!
    m_GLPK_pSimplex = GLPK.Optimizer(method=GLPK.SIMPLEX)
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tol_bnd"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tol_dj"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("meth"), GLPK.GLP_PRIMAL)
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time limit
    m_GLPK_dSimplex = GLPK.Optimizer(method=GLPK.SIMPLEX)
    MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tol_bnd"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tol_dj"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("meth"), GLPK.GLP_DUAL)
    MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time limit
    m_GLPK_Interior = GLPK.Optimizer(method=GLPK.INTERIOR)
    MOI.set(m_GLPK_Interior, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time limit
 
    # Gurobi also doesn't really have ways to set convergence tolerances
    # NOTE: It should be possible to set the absolute/relative/convergence tolerances for PDLP, but
    #       it is currently undocumented. Lu2025 (cuPDLPx) suggests that the parameters "GURO_PAR_PDHGABSTOL",
    #       "GURO_PAR_PDHGCONVTOL", and "GURO_PAR_PDHGRELTOL" exist for this purpose, but I can't figure
    #       out how to set them.
    m_Gurobi_pSimplex = Gurobi.Optimizer(GRB_ENV)
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("Method"), 0) # Primal simplex
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    m_Gurobi_dSimplex = Gurobi.Optimizer(GRB_ENV)
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("Method"), 1) # Dual simplex
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    m_Gurobi_Interior = Gurobi.Optimizer(GRB_ENV)
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("BarConvTol"), 1E-4) # Relative convergence tolerance for barrier method
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("Method"), 2) # Barrier
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    m_Gurobi_PDLP = Gurobi.Optimizer(GRB_ENV)
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("Method"), 6) # PDLP (not in documentation)
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)

    # HiGHS does have a way to set primal/dual residual tolerances (Not sure if absolute or relative)
    m_HiGHS_pSimplex = HiGHS.Optimizer()
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("solver"), "simplex") # Simplex in general. See next option.
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("simplex_strategy"), 4) # Primal simplex
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    m_HiGHS_dSimplex = HiGHS.Optimizer()
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("solver"), "simplex") # Simplex in general. See next option.
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("simplex_strategy"), 1) # Dual simplex
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    m_HiGHS_Interior = HiGHS.Optimizer()
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("solver"), "ipm") # Interior point method
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    m_HiGHS_PDLP = HiGHS.Optimizer()
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("solver"), "pdlp") # PDLP
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-4) 
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    


    # Same thing but for lower tolerance
    for attempt = 1:3
        for i = 1:n_LPs
            # Reset statuses
            # setup_times = zeros(Float64, n_solvers)
            # hires_solving_times = zeros(Float64, n_solvers)
            hires_termination_statuses = fill(OPTIMIZE_NOT_CALLED, n_solvers)
            hires_objectives = zeros(Float64, n_solvers)
            # hires_correct = zeros(Int, n_solvers)
            # lowres_solving_times = zeros(Float64, n_solvers)
            lowres_termination_statuses = fill(OPTIMIZE_NOT_CALLED, n_solvers)
            lowres_objectives = zeros(Float64, n_solvers)
            # lowres_correct = zeros(Int, n_solvers)

            for solver = 1:length(solver_list)
                current_solver = solver_list[solver]
                if lowres_continue[solver]
                    try
                        # Empty the solver of previous run information
                        MOI.empty!(current_solver)

                        # Add the epigraph variable
                        setup_times[solver] += @elapsed epi = MOI.add_variable(current_solver)

                        # Add other problem variables
                        setup_times[solver] += @elapsed begin
                            vi = Vector{MOI.VariableIndex}(undef, example.nvars)
                            for j = 1:example.nvars
                                vi[j], (_,_) = MOI.add_constrained_variable(current_solver, (MOI.GreaterThan(lvbs[i,j]), MOI.LessThan(uvbs[i,j])))
                            end
                        end

                        # Add constraints
                        for j = 1:PDLP_data.dims.current_LP_length
                            start = (i-1)*PDLP_data.dims.total_LP_length
                            constraint_section = constraint_mat[start+j,1:end]
                            rhs_section = rhs[start+j]
                            setup_times[solver] += @elapsed MOI.add_constraint(current_solver, constraint_section[1]*epi + 
                                                        sum(constraint_section[2:end].*vi[1:example.nvars]),
                                                        MOI.GreaterThan(rhs_section))
                        end

                        # Add the objective function (always already in epigraph form)
                        MOI.set(current_solver, MOI.ObjectiveSense(), MOI.MIN_SENSE)
                        MOI.set(current_solver, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), 0.0+epi)
                        
                        # Optimize, and add the time to the solving time
                        lowres_solving_times[solver, attempt] += @elapsed(MOI.optimize!(current_solver))

                        # If it solved, save the termination status and objective value
                        lowres_termination_statuses[solver] = MOI.get(current_solver, MOI.TerminationStatus())
                        if lowres_termination_statuses[solver] == OPTIMAL
                            lowres_objectives[solver] = MOI.get(current_solver, MOI.ObjectiveValue())
                        elseif lowres_termination_statuses[solver] == OTHER_ERROR ||
                            lowres_termination_statuses[solver] == NUMERICAL_ERROR ||
                            lowres_termination_statuses[solver] == ITERATION_LIMIT ||
                            lowres_termination_statuses[solver] == TIME_LIMIT
                            lowres_objectives[solver] = NaN
                        else
                            lowres_objectives[solver] = Inf
                        end
                    catch
                        lowres_continue[solver] = false
                        lowres_solving_times[solver, attempt] = NaN
                        lowres_termination_statuses[solver] = OTHER_ERROR
                        lowres_objectives[solver] = NaN
                    end
                else
                    lowres_solving_times[solver, attempt] = NaN
                    lowres_termination_statuses[solver] = OTHER_ERROR
                    lowres_objectives[solver] = NaN
                end

                # If the total time is above an hour, just stop doing this solver
                if lowres_solving_times[solver, attempt] > 3600.0
                    lowres_solving_times[solver, attempt] = 3600.0
                    lowres_continue[solver] = false
                end
                # println("Solver $solver done (lowres) ($(time() - big_time))")
            end

            # Determine what the "correct" answer is and count up corrects
            if attempt==1
                PDLP_lowres_correct += answer_verification(
                    lowres_correct, 
                    lowres_termination_statuses, 
                    lowres_objectives, 
                    PDLP_lowres_termination_array[i],
                    PDLP_lowres_objectives_array[i],
                    atol=1e-2, 
                    rtol=1e-2,
                    )
            end
            # println("Answers verified (lowres) ($(time() - big_time))")
        end
    end

    println("")
    return hcat(PDLP_setup_time, PDLP_correct, PDLP_iteration_limit_reached, PDLP_solving_time,  PDLP_multi_solving_time, PDLP_95pct_solving_time, PDLP_95pct_multi_solving_time, PDLP_lowres_correct, PDLP_lowres_solving_time, PDLP_lowres_multi_solving_time, 
                setup_times[1],  hires_correct[1], minimum(hires_solving_times[1,:]), lowres_correct[1], minimum(lowres_solving_times[1,:]),
                setup_times[2],  hires_correct[2], minimum(hires_solving_times[2,:]), lowres_correct[2], minimum(lowres_solving_times[2,:]),
                setup_times[3],  hires_correct[3], minimum(hires_solving_times[3,:]), lowres_correct[3], minimum(lowres_solving_times[3,:]),
                setup_times[4],  hires_correct[4], minimum(hires_solving_times[4,:]), lowres_correct[4], minimum(lowres_solving_times[4,:]),
                setup_times[5],  hires_correct[5], minimum(hires_solving_times[5,:]), lowres_correct[5], minimum(lowres_solving_times[5,:]),
                setup_times[6],  hires_correct[6], minimum(hires_solving_times[6,:]), lowres_correct[6], minimum(lowres_solving_times[6,:]),
                setup_times[7],  hires_correct[7], minimum(hires_solving_times[7,:]), lowres_correct[7], minimum(lowres_solving_times[7,:]),
                setup_times[8],  hires_correct[8], minimum(hires_solving_times[8,:]), lowres_correct[8], minimum(lowres_solving_times[8,:]),
                setup_times[9],  hires_correct[9], minimum(hires_solving_times[9,:]), lowres_correct[9], minimum(lowres_solving_times[9,:]),
                setup_times[10], hires_correct[10], minimum(hires_solving_times[10,:]), lowres_correct[10], minimum(lowres_solving_times[10,:]),
                setup_times[11], hires_correct[11], minimum(hires_solving_times[11,:]), lowres_correct[11], minimum(lowres_solving_times[11,:]),
                )
end

# Verify answers by comparing solutions against one another
function answer_verification(correct, status, objective, PDLP_status, PDLP_objective; atol=1e-5, rtol=1e-5)
    len = length(correct)

    # Given `len` solvers' answers, how do we know which is right? First,
    # check to see if multiple (non-PDLP) solvers think the answer is OPTIMAL 
    # or some kind of INFEASIBLE.
    optimal_count = count(status .== OPTIMAL)
    infeas_count = count(status .== INFEASIBLE .|| status .== INFEASIBLE_OR_UNBOUNDED)
    dualinf_count = count(status .== DUAL_INFEASIBLE)

    # If multiple solvers think the answer is OPTIMAL, do they have roughly the
    # same solution values?
    consensus = -Inf
    for i = 1:len
        if consensus == -Inf
            if status[i] == OPTIMAL
                for j = i+1:len
                    if status[j] == OPTIMAL
                        if isapprox(objective[i], objective[j], atol=atol, rtol=rtol)
                            consensus = objective[i]
                        end
                    end
                end
            end
        end
    end

    # If there's no consensus, and the solvers all disagree or have errors, see if
    # PDLP agrees with any.
    if (consensus == -Inf) && (infeas_count < 2) && (dualinf_count < 2)
        if PDLP_status == TERMINATION_REASON_OPTIMAL
            for i = 1:len
                if status[i] == OPTIMAL
                    if isapprox(PDLP_objective, objective[i], atol=atol, rtol=rtol)
                        consensus = objective[i]
                    end
                end
            end
        elseif (PDLP_status == TERMINATION_REASON_PRIMAL_INFEASIBLE) && (infeas_count == 1)
            infeas_count += 1
        elseif (PDLP_status == TERMINATION_REASON_DUAL_INFEASIBLE) && (dualinf_count == 1)
            dualinf_count += 1
        end
    end

    # If we've reached a consensus, then the answer is OPTIMAL with value `consensus`
    if consensus != -Inf
        for i = 1:len
            correct[i] += (status[i] == OPTIMAL) && (isapprox(objective[i], consensus, atol=atol, rtol=rtol))
        end
        PDLP_correct = (PDLP_status == TERMINATION_REASON_OPTIMAL) && (isapprox(PDLP_objective, consensus, atol=atol, rtol=rtol))
        return PDLP_correct
    else
        # No consensus was reached on optimality. Is the answer INFEASIBLE?
        if infeas_count > 1
            # 2 or more solvers think this is infeasible.
            for i = 1:len
                correct[i] += (status[i] == INFEASIBLE) || (status[i] == INFEASIBLE_OR_UNBOUNDED)
            end
            PDLP_correct = (PDLP_status == TERMINATION_REASON_PRIMAL_INFEASIBLE)
            return PDLP_correct
        elseif dualinf_count > 1
            # 2 or more solvers think this is DUAL infeasible
            for i = 1:len
                correct[i] += (status[i] == DUAL_INFEASIBLE)
            end
            PDLP_correct = (PDLP_status == TERMINATION_REASON_DUAL_INFEASIBLE)
            return PDLP_correct
        else
            # There's no consensus at all. Let's just mark them as correct if they
            # didn't get an error.
            for i = 1:len
                correct[i] += (status[i] == OPTIMAL) ||
                             (status[i] == INFEASIBLE) ||
                             (status[i] == INFEASIBLE_OR_UNBOUNDED) ||
                             (status[i] == DUAL_INFEASIBLE)
            end
            PDLP_correct = (PDLP_status == TERMINATION_REASON_OPTIMAL) ||
                           (PDLP_status == TERMINATION_REASON_PRIMAL_INFEASIBLE) ||
                           (PDLP_status == TERMINATION_REASON_DUAL_INFEASIBLE)
            return PDLP_correct
        end
    end
    return nothing
end

# Run examples in a loop and save the results in a .txt file
overall_start_time = time()
file_path = "PDLP Test Results.txt"
open(file_path, "a") do file
    # If using 2 GPUs:
    # write(file, "ID,name,n_LPs,n_cuts,PDLP_setup_time,PDLP_hires_correct,PDLP_hires_iteration_limit_reached,PDLP_hires_solving_time,TwoGPU_PDLP_hires_solving_time,PDLP_hires_95pct_solving_time,TwoGPU_PDLP_hires_95pct_solving_time,PDLP_lowres_correct,PDLP_lowres_solving_time,TwoGPU_PDLP_lowres_solving_time,GLPK_pSimplex_setup_time,GLPK_pSimplex_hires_correct,GLPK_pSimplex_hires_solving_time,GLPK_pSimplex_lowres_correct,GLPK_pSimplex_lowres_solving_time,GLPK_dSimplex_setup_time,GLPK_dSimplex_hires_correct,GLPK_dSimplex_hires_solving_time,GLPK_dSimplex_lowres_correct,GLPK_dSimplex_lowres_solving_time,GLPK_IPM_setup_time,GLPK_IPM_hires_correct,GLPK_IPM_hires_solving_time,GLPK_IPM_lowres_correct,GLPK_IPM_lowres_solving_time,Gurobi_pSimplex_setup_time,Gurobi_pSimplex_hires_correct,Gurobi_pSimplex_hires_solving_time,Gurobi_pSimplex_lowres_correct,Gurobi_pSimplex_lowres_solving_time,Gurobi_dSimplex_setup_time,Gurobi_dSimplex_hires_correct,Gurobi_dSimplex_hires_solving_time,Gurobi_dSimplex_lowres_correct,Gurobi_dSimplex_lowres_solving_time,Gurobi_IPM_setup_time,Gurobi_IPM_hires_correct,Gurobi_IPM_hires_solving_time,Gurobi_IPM_lowres_correct,Gurobi_IPM_lowres_solving_time,Gurobi_PDLP_setup_time,Gurobi_PDLP_hires_correct,Gurobi_PDLP_hires_solving_time,Gurobi_PDLP_lowres_correct,Gurobi_PDLP_lowres_solving_time,HiGHS_pSimplex_setup_time,HiGHS_pSimplex_hires_correct,HiGHS_pSimplex_hires_solving_time,HiGHS_pSimplex_lowres_correct,HiGHS_pSimplex_lowres_solving_time,HiGHS_dSimplex_setup_time,HiGHS_dSimplex_hires_correct,HiGHS_dSimplex_hires_solving_time,HiGHS_dSimplex_lowres_correct,HiGHS_dSimplex_lowres_solving_time,HiGHS_IPM_setup_time,HiGHS_IPM_hires_correct,HiGHS_IPM_hires_solving_time,HiGHS_IPM_lowres_correct,HiGHS_IPM_lowres_solving_time,HiGHS_PDLP_setup_time,HiGHS_PDLP_hires_correct,HiGHS_PDLP_hires_solving_time,HiGHS_PDLP_lowres_correct,HiGHS_PDLP_lowres_solving_time\n")

    # If using 1 GPU:
    write(file, "ID,name,n_LPs,n_cuts,PDLP_setup_time,PDLP_hires_correct,PDLP_hires_iteration_limit_reached,PDLP_hires_solving_time,PDLP_hires_95pct_solving_time,PDLP_lowres_correct,PDLP_lowres_solving_time,GLPK_pSimplex_setup_time,GLPK_pSimplex_hires_correct,GLPK_pSimplex_hires_solving_time,GLPK_pSimplex_lowres_correct,GLPK_pSimplex_lowres_solving_time,GLPK_dSimplex_setup_time,GLPK_dSimplex_hires_correct,GLPK_dSimplex_hires_solving_time,GLPK_dSimplex_lowres_correct,GLPK_dSimplex_lowres_solving_time,GLPK_IPM_setup_time,GLPK_IPM_hires_correct,GLPK_IPM_hires_solving_time,GLPK_IPM_lowres_correct,GLPK_IPM_lowres_solving_time,Gurobi_pSimplex_setup_time,Gurobi_pSimplex_hires_correct,Gurobi_pSimplex_hires_solving_time,Gurobi_pSimplex_lowres_correct,Gurobi_pSimplex_lowres_solving_time,Gurobi_dSimplex_setup_time,Gurobi_dSimplex_hires_correct,Gurobi_dSimplex_hires_solving_time,Gurobi_dSimplex_lowres_correct,Gurobi_dSimplex_lowres_solving_time,Gurobi_IPM_setup_time,Gurobi_IPM_hires_correct,Gurobi_IPM_hires_solving_time,Gurobi_IPM_lowres_correct,Gurobi_IPM_lowres_solving_time,Gurobi_PDLP_setup_time,Gurobi_PDLP_hires_correct,Gurobi_PDLP_hires_solving_time,Gurobi_PDLP_lowres_correct,Gurobi_PDLP_lowres_solving_time,HiGHS_pSimplex_setup_time,HiGHS_pSimplex_hires_correct,HiGHS_pSimplex_hires_solving_time,HiGHS_pSimplex_lowres_correct,HiGHS_pSimplex_lowres_solving_time,HiGHS_dSimplex_setup_time,HiGHS_dSimplex_hires_correct,HiGHS_dSimplex_hires_solving_time,HiGHS_dSimplex_lowres_correct,HiGHS_dSimplex_lowres_solving_time,HiGHS_IPM_setup_time,HiGHS_IPM_hires_correct,HiGHS_IPM_hires_solving_time,HiGHS_IPM_lowres_correct,HiGHS_IPM_lowres_solving_time,HiGHS_PDLP_setup_time,HiGHS_PDLP_hires_correct,HiGHS_PDLP_hires_solving_time,HiGHS_PDLP_lowres_correct,HiGHS_PDLP_lowres_solving_time\n")

end
compilation_example = LoadedProblem(ex8_1_4)
run_example(compilation_example, 16, 3)

for i in 1:size(included, 1)
    elapsed_time = time() - overall_start_time
    days, rem1 = divrem(elapsed_time, 86400)
    hours, rem2 = divrem(rem1, 3600)
    minutes, seconds = divrem(rem2, 60)
    D = isone(days) ? "day" : "days"
    H = isone(hours) ? "hour" : "hours"
    M = isone(minutes) ? "minute" : "minutes"
    S = isone(seconds) ? "second" : "seconds"
    if days > 0
        time_string = "$(Int(days)) $D, $(Int(hours)) $H, $(Int(minutes)) $M, $(round(seconds, digits=2)) $S"
    elseif hours > 0
        time_string = "$(Int(hours)) $H, $(Int(minutes)) $M, $(round(seconds, digits=2)) $S"
    elseif minutes > 0
        time_string = "$(Int(minutes)) $M, $(round(seconds, digits=2)) $S"
    else
        time_string = "$(round(seconds, digits=2)) $S"
    end
    curr = Dates.Libc.TmStruct(Libc.TimeVal().sec)
    ampm = (curr.hour < 12 || curr.hour==24) ? "AM" : "PM"
    min_zeros = "0"^(2-length(string(curr.min)))
    sec_zeros = "0"^(2-length(string(curr.sec)))

    println("$i of $(size(included,1)) (Elapsed time: $time_string) [$(curr.year+1900)-$(curr.month+1)-$(curr.mday) $(mod(curr.hour-1,12)+1):$min_zeros$(curr.min):$sec_zeros$(curr.sec) $ampm]")
    ex = LoadedProblem(eval(Symbol(included[i,1])))
    open(file_path, "a") do file
        @show ex
        write(file, "$i,$(included[i,1]),10000,3,$(join(run_example(ex, 10000, 3),","))\n")
    end
end

