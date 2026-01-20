# BatchPDLP.jl

This package applies the PDLP algorithm [[1](#references)] to solve batches of small, structurally similar linear programs (LPs) simultaneously using NVIDIA GPU resources through CUDA.jl [[2](#references)]. It assumes that LPs are provided (or can be generated) in GPU memory, and only stores solution information in GPU memory. I.e., this package does not manage the transfer of LP information to or from the CPU. `BatchPDLP.jl` is meant to be used within a global optimizer such as `EAGO.jl` [[3](#references)] to solve large numbers of LPs of sizes typically seen within global optimization. `BatchPDLP.jl` was tested on batches of LPs with up to roughly 100 variables and 1000 constraints each.


## Usage Note

`BatchPDLP.jl` is primarily designed to work with `SourceCodeMcCormick.jl` [[4](#references)], which calculates McCormick relaxations and their subgradients for factorable expressions. These relaxations and their subgradients can then be converted into LP constraints through `BatchPDLP.jl`, after which the PDLP algorithm can be run. Internally, `BatchPDLP.jl` works by launching a kernel with the number of blocks set equal to the number of LPs being solved, meaning the number of LPs per batch should generally be larger than the number of streaming multiprocessors (SMs) in the GPU. 

More specifically, `BatchPDLP.jl` is not meant to handle single LP instances. It is explicitly not designed to handle the typical "large" LPs seen in, e.g., LP benchmark test sets, which makes it distinct from many other GPU-accelerated LP solvers such as `cuPDLP.jl` [[5](#references)]. `BatchPDLP.jl` is only meant to be performant when solving hundreds, thousands, or tens of thousands of structurally similar LPs.

## Basic Functionality

The main way to interact with `BatchPDLP.jl` is through its `PDLPData` struct. The struct is instantiated with the maximum number of LPs being solved, the number of variables in each LP (i.e., the width of the constraint matrix, i.e., the epigraph variable should be included), and the maximum number of constraints allotted for each LP. For example, to create an empty `PDLPData` struct to solve 10000 simultaneous LPs at a time, where the LPs come from a 2-variable NLP with an objective function and one GEQ constraint, one could do the following:
```julia
using BatchPDLP

# Total LPs per batch (needed to preallocate space)
num_LPs = 10000

# Number of variables is the number of variables in 
# the NLP plus one for the epigraph variable
actual_vars = 2 + 1

# The maximum number of constraints depends on how many
# you intend to add. E.g., if we're adding a lower bound
# constraint for the epigraph variable, one constraint
# from a relaxation of the objective function, and one
# constraint from a relaxation of the one GEQ constraint,
# we could set max_constraints to 3:
max_constraints = 1 + 1 + 1

PDLP_data = PDLPData(num_LPs, actual_vars, max_constraints)
```

One of the strengths of `BatchPDLP.jl` is that in its main intended use case, all of the LPs will be derived from an original NLP, which means the sparsity pattern of each LP will be (mostly) identical. `BatchPDLP.jl` saves memory space by only storing one copy of the sparsity pattern, which it applies to every LP during the PDLP algorithm to skip known zeros in the constraint matrix. Including this sparsity information speeds up `BatchPDLP.jl` and should be added if it is known. If no sparsity information is added, `BatchPDLP.jl` will assume that any element of the constraint matrix could be non-zero (i.e., it relaxes the requirement that all LPs must be structurally similar, but each LP has a time penalty due to many more unnecessary calculations). This sparsity pattern should be calculated separately from `BatchPDLP.jl`, but is provided in the form of a `Matrix{Bool}` to the `PDLPData` struct for simplicity of use. Internally, this `Matrix{Bool}` is immediately converted into `PDLPSparsity` information that stores the same data on the GPU. The sparsity matrix is passed as an argument to `PDLPData` as follows:
```julia
using BatchPDLP

# Sparsity information for a 2-dimensional NLP like the following:
# min x^2 + y^2
# s.t. sqrt(x) > 1.0
#      sqrt(y) > 1.5
#      x + y > 3.5

sparsity = Bool[1 0 0; # Lower bound of the epigraph variable (which is always first)
                1 1 1;  # Objective function, which becomes an LP constraint that depends on {epi, x, y}
                0 1 0;  # First constraint, only depends on x
                0 0 1;  # Second constraint, only depends on y
                0 1 1]  # Final constraint, depends on x and y

PDLP_data = PDLPData(10000, 3, 3, sparsity=sparsity)
```

Various parameters can also be adjusted when creating the `PDLPData` struct, and/or can be directly modified later from within either the `PDLPData.parameters` or `PDLPData.parameters.termination_criteria` fields. E.g.:
```julia
using BatchPDLP
PDLP_data = PDLPData(10000, 3, 3,
                    iteration_limit = 50000, # Max PDLP iterations per LP
                    termination_evaluation_frequency = 100, # PDLP iterations before checking termination criteria
                    abs_tol = 1E-6, # Absolute tolerance for convergence
                    rel_tol = 1E-6) # Relative tolerance for convergence

# Change primal/dual infeasibility tolerances from default (1E-8) to 1E-7
PDLPData.parameters.termination_criteria.eps_primal_infeasible = 1E-7
PDLPData.parameters.termination_criteria.eps_dual_infeasible = 1E-7
```

There are several functions built into `lower_level_subroutines.jl` and `kernels.jl` that can be used to add constraints to the necessary fields in `PDLPData`, although they are written specifically to interact with `SourceCodeMcCormick.jl` and the `ParBB` extension of EAGO. An example of how some of these functions are meant to be used is included in the `./benchmarks` folder, along with the script that was run to generate benchmark results for the paper. Generally, if a constraint is to be added to each LP, the following actions should be taken:

1) Ensure that the LP constraint to be added is in `>=` form, with only a constant on the right-hand side.
2) Add variable coefficients for this constraint to the `PDLPData.original_problem.constraint_matrix` field, in the line `PDLPData.dims.current_LP_length + 1` for each LP. Note that `PDLPData.dims.current_LP_length + 1` should not exceed `PDLPData.dims.total_LP_length`, or else data will spill into the constraint space designated for the following LP(s). 
3) Add the right-hand side constant term to `PDLPData.original_problem.right_hand_side`, at the index `PDLPData.dims.current_LP_length + 1` for each LP
4) Change `PDLPData.active_constraint` at the index `PDLPData.dims.current_LP_length + 1` to `true`, for each LP where a constraint was added.
5) Set `PDLPData.current_LP_length` equal to `PDLPData.current_LP_length + 1`.

Once all necessary constraints have been added, the LPs can be solved by calling the `PDLP` function with the `PDLPData` struct as an argument. Storage for the LP solutions (with a size of (`PDLPData.dims.n_LPs`, `PDLPData.dims.n_vars`)) and objective values (with a size of (`PDLPData.dims.n_LPs`)) are provided as arguments, though are not technically required. The global upper bound used in a global optimization routine may also be provided as an argument, which will cause `BatchPDLP.jl` to terminate individual LPs if dual feasible solutions are found with dual objective values greater than the global upper bound (with termination code `TERMINATION_REASON_GLOBAL_UPPER_BOUND_HIT`). The call to `PDLP` may look like the following:
```julia
using BatchPDLP, CUDA

# A generic PDLPData struct
PDLP_data = PDLPData(10000, 3, 3)

# Add constraints to each LP
# [Add constraints in some way]

# Pre-allocate storage for solutions and dual objective values
PDLP_solutions = CuArray{Float64}(undef, PDLP_data.dims.n_LPs, PDLP_data.dims.n_vars)
dual_objectives = CuArray{Float64}(undef, PDLP_data.dims.n_LPs)
upper_bound = 0.0

# Run PDLP
PDLP(PDLP_data,
     solutions=PDLP_solutions
     objectives=dual_objectives,
     global_upper_bound=upper_bound)
```

In this example, 10000 LPs are prepared to be solved, and results will be saved to `PDLP_solutions` and `dual_objectives`. Note that by default, `BatchPDLP` returns dual objective values rather than primal objective values so as to give a conservative lower bound for global optimization purposes. By the Strong Duality Theorem, the primal and dual objective values should be equal at a solution, though in practice PDLP terminates when the gap between these objective values is within the specified absolute and relative tolerances. Note also that the `global_upper_bound` is not required to be added, but in this case, it specifies that any LPs with dual feasible solutions having objective values above 0.0 should be terminated, regardless of the duality gap. 


## Citing BatchPDLP
Please cite the following paper when using `BatchPDLP.jl`. In plain text form this is:
```
Gottlieb, R. X., Alston, D., and Stuber, M. D. Re-Architected PDLP for Batch Parallelization of Linear Programs. Under Review.
```
A BibTeX entry is given below:
```bibtex
@Article{,
  author    = {Robert X. Gottlieb, Dimitri Alston, and Matthew D. Stuber},
  journal   = {Under Review},
  title     = {Re-Architected PDLP for Batch Parallelization of Linear Programs},
  year      = {2026},
  pages     = {},
  doi       = {},
  eprint    = {},
  publisher = {},
  url       = {},
}
```


## References
1. Applegate, D., Díaz, M., Hinder, O., Lu, H., Lubin, M., O’Donoghue, B., Schudy, W.: Practical large-scale linear programming using primal-dual hybrid gradient (2021) https://doi.org/10.48550/ARXIV.2106.04756 arXiv:2106.04756 [math.OC]
2. Besard, T., Foket, C., and De Sutter, B. Effective Extensible Programming: Unleashing Julia on GPUs. IEEE Transactions on Parallel and Distributed Systems (2018). https://doi.org/10.1109/TPDS.2018.2872064
3. Wilhelm, M.E., Stuber, M.D.: EAGO.jl: easy advanced global optimization in Julia. Optimization Methods and Software 37(2), 425–450 (2022) https://doi.org/10.1080/10556788.2020.1786566
4. Gottlieb, R.X., Stuber, M.D.: Automatic generation of GPU kernels for evaluators of McCormick-based relaxations and subgradients. Under Revision (2025)
5. Lu, H., Yang, J.: cuPDLP.jl: A GPU implementation of restarted primal-dual hybrid gradient for linear programming in Julia (2024) https://doi.org/10.48550/ARXIV.2311.12180 arXiv:2311.12180 [math.OC]