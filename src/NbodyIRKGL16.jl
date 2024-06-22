__precompile__()

module NbodyIRKGL16

using Reexport
@reexport using DiffEqBase

using Parameters
using OrdinaryDiffEq
using SIMD

export fbirkgl16_gen, fbirkgl16_simd
export IRKGLstep_fixed!, IRKGLstep_adap!
export IRKNGLstep_fixed!, IRKNGLstep_adap!
export IRKGLstep_SIMD_fixed!, IRKGLstep_SIMD_adap!
export IRKNGLstep_SIMD_fixed!, IRKNGLstep_SIMD_adap!
export PolInterp, PolInterp!
export Rdigits
export VecArray

include("./common/IRKGL_Coefficients.jl")
include("./common/aux_functions_SIMD.jl")
include("./common/aux_functions.jl")

include("./generic/IRKGL_Solver.jl")
include("./generic/IRKGL_Step_Functions.jl")
include("./generic/IRKNGL_Step_Functions.jl")

include("./simd/VecArray_def.jl")
include("./simd/IRKGL_SIMD_Solver.jl")
include("./simd/IRKGL_SIMD_Step_Functions.jl")
include("./simd/IRKNGL_SIMD_Step_Functions.jl")


end # module

