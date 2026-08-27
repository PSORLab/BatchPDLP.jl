# This package contains a version of the PDLP algorithm that is designed to 
# solve multiple LPs simultaneously and take CuArrays as inputs. This section 
# borrows extremely heavily from cuPDLPx:
# https://github.com/MIT-Lu-Lab/cuPDLPx/
#
# The version implemented here is meant to be the same as cuPDLPx, but with
# the lowest-level functions modified explicitly to solve multiple LPs. Additionally,
# the following assumptions can be made:
# 1) There are no equality constraints. Because the relaxations and subgradients
#    are relaxed constraints, we always have either less-than or greater-than
#    constraints. Equality constraints in the original problem formulation are
#    split into pairs of less-than and greater-than constraints.
# 2) cuPDLPx does bounds and objective vector rescaling in addition to Ruiz and Pock-Chambolle
#    rescaling methods. Here, we do not do bounds and objective rescaling.
#    
# 3) Since we're getting relaxations of the original optimization problem, it's
#    easiest to structure the LP using an epigraph reformulation. So, the objective
#    functions of the LPs being constructed will always be identical. 

module BatchPDLP
    # Import CUDA and frequently used structs, functions, and macros
    import CUDA
    import CUDA: CuArray, sync_threads, unsafe_load, threadIdx, blockIdx, blockDim, gridDim
    import CUDA: @cuda, @cuDynamicSharedMem, @cuStaticSharedMem
    
    # Export the main struct and the PDLP function itself
    export PDLPData, PDLP

    # Export ways of adding constraints
    export add_LP_objective_constraint, add_LP_constraint, add_LP_lower_bound,
            add_best_obj_LP_constraints, add_best_cons_LP_constraints, add_multiple_LP_lower_bound

    include(joinpath(@__DIR__, "structs.jl"))
    include(joinpath(@__DIR__, "kernels.jl"))
    include(joinpath(@__DIR__, "main_loop.jl"))
    include(joinpath(@__DIR__, "main_loop_rHalpern.jl"))
    include(joinpath(@__DIR__, "lower_level_subroutines.jl"))
    include(joinpath(@__DIR__, "primary_subroutines.jl"))
end