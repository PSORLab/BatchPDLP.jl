# This package contains a version of the PDLP algorithm that is designed to 
# solve multiple LPs simultaneously and take CuArrays as inputs. This section 
# borrows extremely heavily from cuPDLP.jl:
# https://github.com/jinwen-yang/cuPDLP.jl
#
# which utilizes FirstOrderLp.jl:
# https://github.com/google-research/FirstOrderLp.jl
#
# The version implemented here is meant to be the same as cuPDLP.jl, but with
# the lowest-level functions modified explicitly to solve multiple LPs. Additionally,
# the following assumptions can be made:
# 1) There are no equality constraints. Because the relaxations and subgradients
#    are relaxed constraints, we always have either less-than or greater-than
#    constraints. Equality constraints in the original problem formulation are
#    split into pairs of less-than and greater-than constraints.
# 2) Only the adaptive stepsize method will be used. The purpose isn't to recreate
#    cuPDLP for general use; it's only being used as a subsolver. Making only
#    one of the stepsize cases is simpler.
# 3) The original cuPDLP can solve quadratic programming problems, but we will only
#    ever need to solve LPs. So, we can remove the `objective_matrix` in the
#    cuPDLP.QuadraticProgrammingProblem and any lines of code that address that
#    field. 
# 4) Since we're getting relaxations of the original optimization problem, it's
#    easiest to structure the LP using an epigraph reformulation. So, the objective
#    functions of the LPs being constructed will always be identical. 

module BatchPDLP
    # Add necessary packages
    using CUDA
    
    # This constant checks to see how many streaming multiprocessors the user's GPU has,
    # and uses this number to determine block sizes for all CUDA kernels
    const GPU_blocks = Int32(CUDA.attribute(CUDA.device(), CUDA.DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT))

    # Export the main struct and the PDLP function itself
    export PDLPData, PDLP

    include(joinpath(@__DIR__, "structs.jl"))
    include(joinpath(@__DIR__, "kernels.jl"))
    include(joinpath(@__DIR__, "main_loop.jl"))
    include(joinpath(@__DIR__, "lower_level_subroutines.jl"))
    include(joinpath(@__DIR__, "primary_subroutines.jl"))
end