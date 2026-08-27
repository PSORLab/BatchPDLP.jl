
# Use the necessary packages
using Revise
using CSV, Tables, DataFrames
using BatchPDLP, SourceCodeMcCormick, CUDA, Dates, Random
# const GRB_ENV = Gurobi.Env()

# Include all files from the master list
include("./examples/benchmarks/Original Version Benchmarks/only_optimal_master_list.jl")

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

function run_example(example::LoadedProblem, n_LPs::Int, n_cuts::Int; run_BatchPDLP::Bool=true)
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

    # Feed in data to the PDLPData struct
    PDLP_data = PDLPData(n_LPs, example.nvars+1, 1+n_cuts*cut_height, sparsity=sparsity) 

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


    ##############################################################################
    ######### Step 4) Print Problem Data
    ##############################################################################

    # Print out problem information
    CUDA.memory_status()
    println("Example: $(example.name)")
    println("Base problem: $(example.nvars) variables, $(example.ncons) constraints")
    println("There are $leq_len LEQ constraints, $geq_len GEQ constraints, and $eq_len EQ constraints.")
    println("Doing $n_cuts cuts, so each LP has $(PDLP_data.dims.total_LP_length) constraints")

    
    ##############################################################################
    ######### Step 5) Solve Using BatchPDLP
    ##############################################################################

    # Allocate storage for solutions and objectives
    PDLP_solutions = CUDA.zeros(Float64, n_LPs, example.nvars+1)
    PDLP_objectives = CUDA.zeros(Float64, n_LPs, 2)
    PDLP_lowres_solutions = CUDA.zeros(Float64, n_LPs, example.nvars+1)
    PDLP_lowres_objectives = CUDA.zeros(Float64, n_LPs, 2)
    
    # Solve the problems once for compilation
    if run_BatchPDLP
        try
            PDLP(PDLP_data, solutions=PDLP_solutions, objectives=PDLP_objectives, return_both_obj=true)
        catch
            nothing
        end
    end
    
    # Solve the problem once at higher tolerances (default is 1E-8 for abs,rel,primal/dual infeas)
    GC.gc()
    GC.enable(false) # Disable garbage collection temporarily to not impact results

    # Set up storage values
    PDLP_solving_time = zeros(Float64, n_LPs)
    PDLP_lowres_solving_time = zeros(Float64, n_LPs)
    PDLP_95pct_solving_time = zeros(Float64, n_LPs)

    # Solve for the first time, but only do (at least) 95% of problems 
    # (then reset and run all the problems as normal)
    if run_BatchPDLP
        println("Starting BatchPDLPx")
        iterations = Array(PDLP_data.iterations)

        stop_here = 0
        for i in sort(unique(iterations))
            
            if count(iterations .<= i) >= 0.95*n_LPs
                stop_here = i
                break
            end
        end

        
        PDLP_data.parameters.iteration_limit = Int32(stop_here)
        try
            PDLP_95pct_solving_time[1] = @elapsed PDLP(PDLP_data, solutions=PDLP_solutions, objectives=PDLP_objectives, version=:rHalpern, return_both_obj=true)
        catch
            PDLP_95pct_solving_time[1] = NaN
        end
    end

    
    # Save objectives and termination status
    PDLP_95pct_objectives_array = Array(PDLP_objectives)
    PDLP_95pct_termination_array = Array(PDLP_data.termination_reason)
    PDLP_95pct_iterations_array = Array(PDLP_data.iterations)
    
    # Reset the iteration limit for future runs
    PDLP_data.parameters.iteration_limit = Int32(1000000)

    # Run PDLP with an iteration limit of 1E6 and tolerances of 1E-8
    println("-------------------------------------------------------------------------------------------------------------------------------")
    CUDA.memory_status()
    println("Solving with tolerance of 1E-8")
    if run_BatchPDLP
        try
            PDLP_solving_time[1] = @elapsed PDLP(PDLP_data, solutions=PDLP_solutions, objectives=PDLP_objectives, version=:rHalpern, return_both_obj=true)
        catch
            PDLP_solving_time[1] = NaN
        end
    end

    # # Save objectives and termination statuses as CPU arrays
    PDLP_objectives_array = Array(PDLP_objectives)
    PDLP_termination_array = Array(PDLP_data.termination_reason)
    PDLP_iterations_array = Array(PDLP_data.iterations)

    println("-------------------------------------------------------------------------------------------------------------------------------")
    CUDA.memory_status()
    println("Solving with tolerance of 1E-4")
    # Re-run PDLP using lower tolerances (1E-4 instead of 1E-8)
    # (Keep infeasibility tolerances at 1E-8)
    PDLP_data.parameters.termination_criteria.eps_optimal_absolute = 1E-4
    PDLP_data.parameters.termination_criteria.eps_optimal_relative = 1E-4
    PDLP_data.parameters.termination_criteria.eps_primal_infeasible= 1E-8
    PDLP_data.parameters.termination_criteria.eps_dual_infeasible  = 1E-8
    # Solve again with lower tolerances
    if run_BatchPDLP
        try
            PDLP_lowres_solving_time[1] = @elapsed PDLP(PDLP_data, solutions=PDLP_lowres_solutions, objectives=PDLP_lowres_objectives, version=:rHalpern, return_both_obj=true)
        catch
            PDLP_lowres_solving_time[1] = NaN
        end
    end


    GC.enable(true) # PDLP runs are finished; re-enable garbage collection
    GC.gc()

    # Print out some PDLP solve times
    if run_BatchPDLP
        println("BatchPDLPx Solve Times:")
        println("============================================================")
        println("All problems (high-res):     | $(round.(PDLP_solving_time[1], digits=6))")
        println("All problems (low-res):      | $(round.(PDLP_lowres_solving_time[1], digits=6))")
        println(">95% of problems (high-res): | $(round.(PDLP_95pct_solving_time[1], digits=6))")
        println("===========================================================")
    end


    # Save lowres solutions as regular arrays
    PDLP_lowres_objectives_array = Array(PDLP_lowres_objectives)
    PDLP_lowres_termination_array = Array(PDLP_data.termination_reason)
    PDLP_lowres_iterations_array = Array(PDLP_data.iterations)

    res = hcat(PDLP_termination_array, 
                PDLP_iterations_array,
                PDLP_objectives_array,
                PDLP_solving_time,
                PDLP_lowres_termination_array,
                PDLP_lowres_iterations_array,
                PDLP_lowres_objectives_array,
                PDLP_lowres_solving_time,
                PDLP_95pct_termination_array,
                PDLP_95pct_iterations_array,
                PDLP_95pct_objectives_array,
                PDLP_95pct_solving_time
                )
    return res
    
end

# Load and run the problem
ex = LoadedProblem(ex4_1_7, overwrite=true)

data = run_example(ex, 10000, 3)
open("./examples/benchmarks/benchmarking_results/rundata_$(ex.name).csv", "w") do io
    CSV.write(io, Tables.table(data), header=[
            "BatchPDLPx (1xGPU) Hi-Res - Term Status",
            "BatchPDLPx (1xGPU) Hi-Res - Iterations",
            "BatchPDLPx (1xGPU) Hi-Res - Primal Obj",
            "BatchPDLPx (1xGPU) Hi-Res - Dual Obj",
            "BatchPDLPx (1xGPU) Hi-Res - Time (s)",
            "BatchPDLPx (1xGPU) Lo-Res - Term Status",
            "BatchPDLPx (1xGPU) Lo-Res - Iterations",
            "BatchPDLPx (1xGPU) Lo-Res - Primal Obj",
            "BatchPDLPx (1xGPU) Lo-Res - Dual Obj",
            "BatchPDLPx (1xGPU) Lo-Res - Time (s)",
            "BatchPDLPx (1xGPU) Hi-Res 95pct - Term Status",
            "BatchPDLPx (1xGPU) Hi-Res 95pct - Iterations",
            "BatchPDLPx (1xGPU) Hi-Res 95pct - Primal Obj",
            "BatchPDLPx (1xGPU) Hi-Res 95pct - Dual Obj",
            "BatchPDLPx (1xGPU) Hi-Res 95pct - Time (s)",
            ],)
end