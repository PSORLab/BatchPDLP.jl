
# Use the necessary packages
using CSV, Tables
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
# relaxations, create the optimization problems, and pass them to PDLP
# and other comparison solvers (GLPK, Gurobi, HiGHS)
function run_example(example::LoadedProblem, n_LPs::Int, n_cuts::Int; run_GLPK::Bool=true, run_Gurobi::Bool=true, run_HiGHS::Bool=true, run_2xGPU::Bool=false)
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
    try
        PDLP(PDLP_data, solutions=PDLP_solutions, objectives=PDLP_objectives, return_both_obj=true)
    catch
        nothing
    end

    # Solve the problem once at higher tolerances (default is 1E-8 for abs,rel,primal/dual infeas)
    println("Starting BatchPDLP")
    GC.gc()
    GC.enable(false) # Disable garbage collection temporarily to not impact results

    # Set up storage values
    PDLP_solving_time = zeros(Float64, n_LPs)
    PDLP_lowres_solving_time = zeros(Float64, n_LPs)
    PDLP_95pct_solving_time = zeros(Float64, n_LPs)

    # Solve for the first time, but only do (at least) 95% of problems 
    # (then reset and run all the problems as normal)
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
        PDLP_95pct_solving_time[1] = @elapsed PDLP(PDLP_data, solutions=PDLP_solutions, objectives=PDLP_objectives, return_both_obj=true)
    catch
        PDLP_95pct_solving_time[1] = NaN
    end

    # Save objectives and termination status
    PDLP_95pct_objectives_array = Array(PDLP_objectives)
    PDLP_95pct_termination_array = Array(PDLP_data.termination_reason)
    PDLP_95pct_iterations_array = Array(PDLP_data.iterations)
    PDLP_95pct_2xGPU_objectives_array = zeros(Float64, n_LPs, 2)
    PDLP_95pct_2xGPU_termination_array = fill("", n_LPs)
    
    # Reset the iteration limit for future runs
    PDLP_data.parameters.iteration_limit = Int32(1000000)

    # Run PDLP with an iteration limit of 1E6 and tolerances of 1E-8
    try
        PDLP_solving_time[1] = @elapsed PDLP(PDLP_data, solutions=PDLP_solutions, objectives=PDLP_objectives, return_both_obj=true)
    catch
        PDLP_solving_time[1] = NaN
    end

    # Save objectives and termination statuses as CPU arrays
    PDLP_objectives_array = Array(PDLP_objectives)
    PDLP_termination_array = Array(PDLP_data.termination_reason)
    PDLP_iterations_array = Array(PDLP_data.iterations)

    # Re-run PDLP using lower tolerances (1E-4 instead of 1E-8)
    # (Keep infeasibility tolerances at 1E-8)
    PDLP_data.parameters.termination_criteria.eps_optimal_absolute = 1E-4
    PDLP_data.parameters.termination_criteria.eps_optimal_relative = 1E-4
    PDLP_data.parameters.termination_criteria.eps_primal_infeasible= 1E-8
    PDLP_data.parameters.termination_criteria.eps_dual_infeasible  = 1E-8

    # Solve again with lower tolerances
    try
        PDLP_lowres_solving_time[1] = @elapsed PDLP(PDLP_data, solutions=PDLP_lowres_solutions, objectives=PDLP_lowres_objectives, return_both_obj=true)
    catch
        PDLP_lowres_solving_time[1] = NaN
    end

    # Initialize 2xGPU fields
    PDLP_2xGPU_solving_time = zeros(Float64, n_LPs)
    PDLP_2xGPU_objectives_array = zeros(Float64, n_LPs, 2)
    PDLP_2xGPU_termination_array = fill(BatchPDLP.TERMINATION_REASON_UNSPECIFIED, n_LPs)
    PDLP_2xGPU_iterations_array = zeros(Int64, n_LPs)

    PDLP_lowres_2xGPU_solving_time = zeros(Float64, n_LPs)
    PDLP_lowres_2xGPU_objectives_array = zeros(Float64, n_LPs, 2)
    PDLP_lowres_2xGPU_termination_array = fill(BatchPDLP.TERMINATION_REASON_UNSPECIFIED, n_LPs)
    PDLP_lowres_2xGPU_iterations_array = zeros(Int64, n_LPs)

    PDLP_95pct_2xGPU_solving_time = zeros(Float64, n_LPs)
    PDLP_95pct_2xGPU_objectives_array = zeros(Float64, n_LPs, 2)
    PDLP_95pct_2xGPU_termination_array = fill(BatchPDLP.TERMINATION_REASON_UNSPECIFIED, n_LPs)
    PDLP_95pct_2xGPU_iterations_array = zeros(Int64, n_LPs)

    # Run BatchPDLP with 2xGPU if selected
    
    # Feed in data to the PDLPData struct
    PDLP_data = PDLPData(n_LPs, example.nvars+1, 1+n_cuts*cut_height, sparsity=sparsity, iteration_limit=1000000)
    if run_2xGPU
        # Make two more PDLPData structs for the 2-GPU setup
        device!(0)
        PDLP_data_A = PDLPData(Int(n_LPs/2), example.nvars+1, 1+n_cuts*cut_height, sparsity=sparsity, iteration_limit=1000000)
        device!(1)
        PDLP_data_B = PDLPData(Int(n_LPs/2), example.nvars+1, 1+n_cuts*cut_height, sparsity=sparsity, iteration_limit=1000000)

        # Switch back to the first GPU
        device!(0)

        # Fill in the 2-GPU PDLP structs similarly to the original PDLP_data struct
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

        # Fill in matrices appropriately based on the original struct
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

        # Setup storage fields
        device!(0)
        PDLP_A_solutions = CUDA.zeros(Float64, Int(n_LPs/2), example.nvars+1)
        PDLP_A_objectives = CUDA.zeros(Float64, Int(n_LPs/2), 2)
        PDLP_A_lowres_solutions = CUDA.zeros(Float64, Int(n_LPs/2), example.nvars+1)
        PDLP_A_lowres_objectives = CUDA.zeros(Float64, Int(n_LPs/2), 2)
        device!(1)
        PDLP_B_solutions = CUDA.zeros(Float64, Int(n_LPs/2), example.nvars+1)
        PDLP_B_objectives = CUDA.zeros(Float64, Int(n_LPs/2), 2)
        PDLP_B_lowres_solutions = CUDA.zeros(Float64, Int(n_LPs/2), example.nvars+1)
        PDLP_B_lowres_objectives = CUDA.zeros(Float64, Int(n_LPs/2), 2)
        device!(0)

        # Compilation run
        try
            @sync begin
                @async begin
                    device!(0)
                    PDLP(PDLP_data_A, solutions=PDLP_A_solutions, objectives=PDLP_A_objectives, return_both_obj=true)
                end
                @async begin
                    device!(1)
                    PDLP(PDLP_data_B, solutions=PDLP_B_solutions, objectives=PDLP_B_objectives, return_both_obj=true)
                end
            end
        catch
            nothing
        end

        # 95% run
        PDLP_data_A.parameters.iteration_limit = Int32(stop_here)
        PDLP_data_B.parameters.iteration_limit = Int32(stop_here)
        PDLP_95pct_2xGPU_solving_time[1] = @elapsed begin
            @sync begin
                @async begin
                    device!(0)
                    PDLP(PDLP_data_A, solutions=PDLP_A_solutions, objectives=PDLP_A_objectives, return_both_obj=true)
                end
                @async begin
                    device!(1)
                    PDLP(PDLP_data_B, solutions=PDLP_B_solutions, objectives=PDLP_B_objectives, return_both_obj=true)
                end
            end
        end
        device!(0)
        PDLP_95pct_2xGPU_objectives_array .= vcat(Array(PDLP_A_objectives), Array(PDLP_B_objectives))
        PDLP_95pct_2xGPU_termination_array .= vcat(Array(PDLP_data_A.termination_reason), Array(PDLP_data_B.termination_reason))
        PDLP_95pct_2xGPU_iterations_array .= vcat(Array(PDLP_data_A.iterations), Array(PDLP_data_B.iterations))

        # Full run
        PDLP_data_A.parameters.iteration_limit = Int32(1000000)
        PDLP_data_B.parameters.iteration_limit = Int32(1000000)
        try 
            PDLP_2xGPU_solving_time[1] = @elapsed begin
                @sync begin
                    @async begin
                        device!(0)
                        PDLP(PDLP_data_A, solutions=PDLP_A_solutions, objectives=PDLP_A_objectives, return_both_obj=true)
                    end
                    @async begin
                        device!(1)
                        PDLP(PDLP_data_B, solutions=PDLP_B_solutions, objectives=PDLP_B_objectives, return_both_obj=true)
                    end
                end
            end
        catch
            PDLP_2xGPU_solving_time[1] = NaN
        end
        device!(0)
        PDLP_2xGPU_objectives_array .= vcat(Array(PDLP_A_objectives), Array(PDLP_B_objectives))
        PDLP_2xGPU_termination_array .= vcat(Array(PDLP_data_A.termination_reason), Array(PDLP_data_B.termination_reason))
        PDLP_2xGPU_iterations_array .= vcat(Array(PDLP_data_A.iterations), Array(PDLP_data_B.iterations))

        # Full run, lower tolerance
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
        try 
            PDLP_lowres_2xGPU_solving_time[1] = @elapsed begin
                @sync begin
                    @async begin
                        device!(0)
                        PDLP(PDLP_data_A, solutions=PDLP_A_lowres_solutions, objectives=PDLP_A_lowres_objectives, return_both_obj=true)
                    end
                    @async begin
                        device!(1)
                        PDLP(PDLP_data_B, solutions=PDLP_B_lowres_solutions, objectives=PDLP_B_lowres_objectives, return_both_obj=true)
                    end
                end
            end
        catch
            PDLP_lowres_2xGPU_solving_time[1] = NaN
        end
        PDLP_lowres_2xGPU_objectives_array .= vcat(Array(PDLP_A_lowres_objectives), Array(PDLP_B_lowres_objectives))
        PDLP_lowres_2xGPU_termination_array .= vcat(Array(PDLP_data_A.termination_reason), Array(PDLP_data_B.termination_reason))
        PDLP_lowres_2xGPU_iterations_array .= vcat(Array(PDLP_data_A.iterations), Array(PDLP_data_B.iterations))
        device!(0)
    end

    GC.enable(true) # PDLP runs are finished; re-enable garbage collection
    GC.gc()

    # Print out some PDLP solve times
    println("BatchPDLP Solve Times:")
    println("============================================================")
    println("All problems (high-res):     | $(round.(PDLP_solving_time[1], digits=6))")
    println("All problems (low-res):      | $(round.(PDLP_lowres_solving_time[1], digits=6))")
    println(">95% of problems (high-res): | $(round.(PDLP_95pct_solving_time[1], digits=6))")
    println("===========================================================")
    if run_2xGPU
        println("All problems (high-res) (2xGPU):     | $(round.(PDLP_2xGPU_solving_time[1], digits=6))")
        println("All problems (low-res) (2xGPU):      | $(round.(PDLP_lowres_2xGPU_solving_time[1], digits=6))")
        println(">95% of problems (high-res) (2xGPU): | $(round.(PDLP_95pct_2xGPU_solving_time[1], digits=6))")
        println("===========================================================")
    end

    # Save constraints and the RHS for CPU solvers
    constraint_mat = Array(LPs.constraint_matrix)
    rhs = Array(LPs.right_hand_side)

    # Save lowres solutions as regular arrays
    PDLP_lowres_objectives_array = Array(PDLP_lowres_objectives)
    PDLP_lowres_termination_array = Array(PDLP_data.termination_reason)
    PDLP_lowres_iterations_array = Array(PDLP_data.iterations)

    #############################################################
    ##################### Other Solvers #########################
    #############################################################

    # Number of alternate solvers
    n_solvers = 11
    
    # Now, we wish to create optimization problems and time how long it takes
    # to solve them using the commercial optimizers.
    solving_times = zeros(Float64, n_LPs, n_solvers)
    termination_statuses = fill(OPTIMIZE_NOT_CALLED, n_LPs, n_solvers)
    primal_objectives = zeros(Float64, n_LPs, n_solvers)
    dual_objectives = zeros(Float64, n_LPs, n_solvers)
    
    # Initialize all solvers
    m_GLPK_pSimplex = GLPK.Optimizer(method=GLPK.SIMPLEX)
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tol_bnd"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tol_dj"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("meth"), GLPK.GLP_PRIMAL)
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time limit
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("msg_lev"), GLPK.GLP_MSG_OFF) # Terminal output
    MOI.set(m_GLPK_pSimplex, MOI.Silent(), true)

    m_GLPK_dSimplex = GLPK.Optimizer(method=GLPK.SIMPLEX)
    MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tol_bnd"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tol_dj"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("meth"), GLPK.GLP_DUAL)
    MOI.set(m_GLPK_dSimplex, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time 
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("msg_lev"), GLPK.GLP_MSG_OFF) # Terminal output
    MOI.set(m_GLPK_dSimplex, MOI.Silent(), true)

    m_GLPK_Interior = GLPK.Optimizer(method=GLPK.INTERIOR)
    MOI.set(m_GLPK_Interior, MOI.RawOptimizerAttribute("tm_lim"), 1.0) # Time limit
    MOI.set(m_GLPK_pSimplex, MOI.RawOptimizerAttribute("msg_lev"), GLPK.GLP_MSG_OFF) # Terminal output
    MOI.set(m_GLPK_Interior, MOI.Silent(), true)
 
    m_Gurobi_pSimplex = Gurobi.Optimizer(GRB_ENV)
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("Method"), 0) # Primal simplex
    MOI.set(m_Gurobi_pSimplex, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    MOI.set(m_Gurobi_pSimplex, MOI.Silent(), true)

    m_Gurobi_dSimplex = Gurobi.Optimizer(GRB_ENV)
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("Method"), 1) # Dual simplex
    MOI.set(m_Gurobi_dSimplex, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    MOI.set(m_Gurobi_dSimplex, MOI.Silent(), true)

    m_Gurobi_Interior = Gurobi.Optimizer(GRB_ENV)
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("BarConvTol"), 1E-8) # Relative convergence tolerance for barrier method
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("Method"), 2) # Barrier
    MOI.set(m_Gurobi_Interior, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    MOI.set(m_Gurobi_Interior, MOI.Silent(), true)

    m_Gurobi_PDLP = Gurobi.Optimizer(GRB_ENV)
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("OutputFlag"), 0)
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("FeasibilityTol"), 1E-8) # Primal feasibility tolerance
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("OptimalityTol"), 1E-8) # Dual feasibility tolerance
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("Method"), 6) # PDLP (not in documentation)
    MOI.set(m_Gurobi_PDLP, MOI.RawOptimizerAttribute("TimeLimit"), 1.0)
    MOI.set(m_Gurobi_PDLP, MOI.Silent(), true)

    m_HiGHS_pSimplex = HiGHS.Optimizer()
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("solver"), "simplex") # Simplex in general. See next option.
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("simplex_strategy"), 4) # Primal simplex
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_pSimplex, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    MOI.set(m_HiGHS_pSimplex, MOI.Silent(), true)

    m_HiGHS_dSimplex = HiGHS.Optimizer()
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("solver"), "simplex") # Simplex in general. See next option.
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("simplex_strategy"), 1) # Dual simplex
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_dSimplex, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    MOI.set(m_HiGHS_dSimplex, MOI.Silent(), true)

    m_HiGHS_Interior = HiGHS.Optimizer()
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("solver"), "ipm") # Interior point method
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_Interior, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    MOI.set(m_HiGHS_Interior, MOI.Silent(), true)

    m_HiGHS_PDLP = HiGHS.Optimizer()
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("solver"), "pdlp") # PDLP
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("primal_residual_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("dual_residual_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("primal_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("dual_feasibility_tolerance"), 1E-8) 
    MOI.set(m_HiGHS_PDLP, MOI.RawOptimizerAttribute("time_limit"), 1.0)
    MOI.set(m_HiGHS_PDLP, MOI.Silent(), true)
                         
    # List of all solvers being used
    solver_list = [m_GLPK_pSimplex, m_GLPK_dSimplex, m_GLPK_Interior,
                   m_Gurobi_pSimplex, m_Gurobi_dSimplex, m_Gurobi_Interior, m_Gurobi_PDLP,
                   m_HiGHS_pSimplex, m_HiGHS_dSimplex, m_HiGHS_Interior, m_HiGHS_PDLP]

    # Compile solvers
    for i in eachindex(solver_list)
        if !run_GLPK && i <= 3
            solving_times[:,i] .= NaN
            primal_objectives[i,:] .= NaN
            dual_objectives[i,:] .= NaN
            continue
        elseif ~run_Gurobi && (i >= 4 && i <= 8)
            solving_times[:,i] .= NaN
            primal_objectives[i,:] .= NaN
            dual_objectives[i,:] .= NaN
            continue
        elseif ~run_HiGHS && (i >= 9 && i <= 12)
            solving_times[:,i] .= NaN
            primal_objectives[i,:] .= NaN
            dual_objectives[i,:] .= NaN
            continue
        end
        
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
            nothing
        end
    end

    # Now solve all the LPs
    for solver = eachindex(solver_list)
        if !run_GLPK && solver <= 3
            continue
        elseif ~run_Gurobi && (solver >= 4 && solver <= 7)
            continue
        elseif ~run_HiGHS && (solver >= 8 && solver <= 11)
            continue
        end
        current_solver = solver_list[solver]
        if solver==1
            println("Starting GLPK")
        elseif solver==4
            println("Starting Gurobi")
        elseif solver==8
            println("Starting HiGHS")
        end
        for i = 1:n_LPs
            try
                # Empty the solver of previous run information
                MOI.empty!(current_solver)

                # Add the epigraph variable
                epi = MOI.add_variable(current_solver)

                # Add other problem variables
                vi = Vector{MOI.VariableIndex}(undef, example.nvars)
                for j = 1:example.nvars
                    vi[j], (_,_) = MOI.add_constrained_variable(current_solver, (MOI.GreaterThan(lvbs[i,j]), MOI.LessThan(uvbs[i,j])))
                end

                # Add constraints
                for j = 1:PDLP_data.dims.current_LP_length
                    start = (i-1)*PDLP_data.dims.total_LP_length
                    constraint_section = constraint_mat[start+j,1:end]
                    rhs_section = rhs[start+j]
                    MOI.add_constraint(current_solver, constraint_section[1]*epi + 
                                                sum(constraint_section[2:end].*vi[1:example.nvars]),
                                                MOI.GreaterThan(rhs_section))
                end

                # Add the objective function (always already in epigraph form)
                MOI.set(current_solver, MOI.ObjectiveSense(), MOI.MIN_SENSE)
                MOI.set(current_solver, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), 0.0+epi)

                # Optimize, and add the time to the solving time
                solving_times[i, solver] = @elapsed(MOI.optimize!(current_solver))

                # If it solved, save the termination status and objective value
                termination_statuses[i, solver] = MOI.get(current_solver, MOI.TerminationStatus())
                if termination_statuses[i, solver] == OPTIMAL
                    primal_objectives[i, solver] = MOI.get(current_solver, MOI.ObjectiveValue())
                    dual_objectives[i, solver] = MOI.get(current_solver, MOI.DualObjectiveValue())
                elseif termination_statuses[i, solver] == OTHER_ERROR ||
                    termination_statuses[i, solver] == NUMERICAL_ERROR ||
                    termination_statuses[i, solver] == ITERATION_LIMIT ||
                    termination_statuses[i, solver] == TIME_LIMIT
                    primal_objectives[i, solver] = NaN
                    dual_objectives[i, solver] = NaN
                else
                    primal_objectives[i, solver] = Inf
                    dual_objectives[i, solver] = Inf
                end

                # Empty the solver
                MOI.empty!(current_solver)
            catch
                solving_times[i, solver] = NaN
                termination_statuses[i, solver] = OTHER_ERROR
                primal_objectives[i, solver] = NaN
                dual_objectives[i, solver] = NaN
            end
        end
        if solver==3
            println("GLPK Solve Times:")
            println("============================================================")
            println("Primal Simplex: | $(round(sum(solving_times[:,1]), digits=6))")
            println("Dual Simplex:   | $(round(sum(solving_times[:,2]), digits=6))")
            println("Interior Point: | $(round(sum(solving_times[:,3]), digits=6))")
            println("===========================================================")
        elseif solver==7
            println("Gurobi Solve Times:")
            println("============================================================")
            println("Primal Simplex: | $(round(sum(solving_times[:,4]), digits=6))")
            println("Dual Simplex:   | $(round(sum(solving_times[:,5]), digits=6))")
            println("Interior Point: | $(round(sum(solving_times[:,6]), digits=6))")
            println("PDLP:           | $(round(sum(solving_times[:,7]), digits=6))")
            println("===========================================================")
        elseif solver==11
            println("HiGHS Solve Times:")
            println("============================================================")
            println("Primal Simplex: | $(round(sum(solving_times[:,8]), digits=6))")
            println("Dual Simplex:   | $(round(sum(solving_times[:,9]), digits=6))")
            println("Interior Point: | $(round(sum(solving_times[:,10]), digits=6))")
            println("PDLP:           | $(round(sum(solving_times[:,11]), digits=6))")
            println("===========================================================")
            println("")
        end

        # Remove reference to the current solver and force garbage collection
        current_solver = nothing
        GC.gc()
    end

    return hcat(PDLP_termination_array, 
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
                PDLP_95pct_solving_time,
                PDLP_2xGPU_termination_array, 
                PDLP_2xGPU_iterations_array, 
                PDLP_2xGPU_objectives_array, 
                PDLP_2xGPU_solving_time,
                PDLP_lowres_2xGPU_termination_array,
                PDLP_lowres_2xGPU_iterations_array, 
                PDLP_lowres_2xGPU_objectives_array,
                PDLP_lowres_2xGPU_solving_time,
                PDLP_95pct_2xGPU_termination_array,
                PDLP_95pct_2xGPU_iterations_array, 
                PDLP_95pct_2xGPU_objectives_array,
                PDLP_95pct_2xGPU_solving_time,
                termination_statuses[:,1],
                primal_objectives[:,1],
                dual_objectives[:,1],
                solving_times[:,1],
                termination_statuses[:,2],
                primal_objectives[:,2],
                dual_objectives[:,2],
                solving_times[:,2],
                termination_statuses[:,3],
                primal_objectives[:,3],
                dual_objectives[:,3],
                solving_times[:,3],
                termination_statuses[:,4],
                primal_objectives[:,4],
                dual_objectives[:,4],
                solving_times[:,4],
                termination_statuses[:,5],
                primal_objectives[:,5],
                dual_objectives[:,5],
                solving_times[:,5],
                termination_statuses[:,6],
                primal_objectives[:,6],
                dual_objectives[:,6],
                solving_times[:,6],
                termination_statuses[:,7],
                primal_objectives[:,7],
                dual_objectives[:,7],
                solving_times[:,7],
                termination_statuses[:,8],
                primal_objectives[:,8],
                dual_objectives[:,8],
                solving_times[:,8],
                termination_statuses[:,9],
                primal_objectives[:,9],
                dual_objectives[:,9],
                solving_times[:,9],
                termination_statuses[:,10],
                primal_objectives[:,10],
                dual_objectives[:,10],
                solving_times[:,10],
                termination_statuses[:,11],
                primal_objectives[:,11],
                dual_objectives[:,11],
                solving_times[:,11],
                )
end

# Run examples in a loop and save the results in separate CSV files
overall_start_time = time()
compilation_example = LoadedProblem(ex8_1_4)
res = run_example(compilation_example, 16, 3, run_2xGPU=true);
for i = 1:size(included, 1)
    # Print timing and problem information
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
    
    # Load and run the problem
    ex = LoadedProblem(eval(Symbol(included[i,1])))
    data = run_example(ex, 10000, 3, run_2xGPU=true)

    # Write data to a CSV
    open("./benchmarking_results/rundata_$(included[i,1]).csv", "w") do io
        CSV.write(io, Tables.table(data), header=[
                "BatchPDLP (1xGPU) Hi-Res - Term Status",
                "BatchPDLP (1xGPU) Hi-Res - Iterations",
                "BatchPDLP (1xGPU) Hi-Res - Primal Obj",
                "BatchPDLP (1xGPU) Hi-Res - Dual Obj",
                "BatchPDLP (1xGPU) Hi-Res - Time (s)",
                "BatchPDLP (1xGPU) Lo-Res - Term Status",
                "BatchPDLP (1xGPU) Lo-Res - Iterations",
                "BatchPDLP (1xGPU) Lo-Res - Primal Obj",
                "BatchPDLP (1xGPU) Lo-Res - Dual Obj",
                "BatchPDLP (1xGPU) Lo-Res - Time (s)",
                "BatchPDLP (1xGPU) Hi-Res 95pct - Term Status",
                "BatchPDLP (1xGPU) Hi-Res 95pct - Iterations",
                "BatchPDLP (1xGPU) Hi-Res 95pct - Primal Obj",
                "BatchPDLP (1xGPU) Hi-Res 95pct - Dual Obj",
                "BatchPDLP (1xGPU) Hi-Res 95pct - Time (s)",
                "BatchPDLP (2xGPU) Hi-Res - Term Status",
                "BatchPDLP (2xGPU) Hi-Res - Iterations",
                "BatchPDLP (2xGPU) Hi-Res - Primal Obj",
                "BatchPDLP (2xGPU) Hi-Res - Dual Obj",
                "BatchPDLP (2xGPU) Hi-Res - Time (s)",
                "BatchPDLP (2xGPU) Lo-Res - Term Status",
                "BatchPDLP (2xGPU) Lo-Res - Iterations",
                "BatchPDLP (2xGPU) Lo-Res - Primal Obj",
                "BatchPDLP (2xGPU) Lo-Res - Dual Obj",
                "BatchPDLP (2xGPU) Lo-Res - Time (s)",
                "BatchPDLP (2xGPU) Hi-Res 95pct - Term Status",
                "BatchPDLP (2xGPU) Hi-Res 95pct - Iterations",
                "BatchPDLP (2xGPU) Hi-Res 95pct - Primal Obj",
                "BatchPDLP (2xGPU) Hi-Res 95pct - Dual Obj",
                "BatchPDLP (2xGPU) Hi-Res 95pct - Time (s)",
                "GLPK Primal Simplex - Term Status",
                "GLPK Primal Simplex - Primal Obj",
                "GLPK Primal Simplex - Dual Obj",
                "GLPK Primal Simplex - Time (s)",
                "GLPK Dual Simplex - Term Status",
                "GLPK Dual Simplex - Primal Obj",
                "GLPK Dual Simplex - Dual Obj",
                "GLPK Dual Simplex - Time (s)",
                "GLPK IPM - Term Status",
                "GLPK IPM - Primal Obj",
                "GLPK IPM - Dual Obj",
                "GLPK IPM - Time (s)",
                "Gurobi Primal Simplex - Term Status",
                "Gurobi Primal Simplex - Primal Obj",
                "Gurobi Primal Simplex - Dual Obj",
                "Gurobi Primal Simplex - Time (s)",
                "Gurobi Dual Simplex - Term Status",
                "Gurobi Dual Simplex - Primal Obj",
                "Gurobi Dual Simplex - Dual Obj",
                "Gurobi Dual Simplex - Time (s)",
                "Gurobi IPM - Term Status",
                "Gurobi IPM - Primal Obj",
                "Gurobi IPM - Dual Obj",
                "Gurobi IPM - Time (s)",
                "Gurobi PDLP - Term Status",
                "Gurobi PDLP - Primal Obj",
                "Gurobi PDLP - Dual Obj",
                "Gurobi PDLP - Time (s)",
                "HiGHS Primal Simplex - Term Status",
                "HiGHS Primal Simplex - Primal Obj",
                "HiGHS Primal Simplex - Dual Obj",
                "HiGHS Primal Simplex - Time (s)",
                "HiGHS Dual Simplex - Term Status",
                "HiGHS Dual Simplex - Primal Obj",
                "HiGHS Dual Simplex - Dual Obj",
                "HiGHS Dual Simplex - Time (s)",
                "HiGHS IPM - Term Status",
                "HiGHS IPM - Primal Obj",
                "HiGHS IPM - Dual Obj",
                "HiGHS IPM - Time (s)",
                "HiGHS PDLP - Term Status",
                "HiGHS PDLP - Primal Obj",
                "HiGHS PDLP - Dual Obj",
                "HiGHS PDLP - Time (s)",
                ],)
    end
end