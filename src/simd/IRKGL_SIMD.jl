__precompile__()

module IRKGL_SIMD

using Reexport
@reexport using DiffEqBase

using Parameters
using OrdinaryDiffEq
using SIMD

export  IRKGL_simd, IRKNGL_simd
export IRKGLstep_SIMD_fixed!, IRKGLstep_SIMD_adap!
export IRKNGLstep_SIMD_fixed!, IRKNGLstep_SIMD_adap!
export PolInterp, PolInterp!
export Rdigits
export  VecArray

include("VecArray_def.jl")
include("../common/IRKGL_Coefficients.jl")
include("../common/aux_functions.jl")
include("../common/aux_functions_SIMD.jl")
include("IRKGL_SIMD_Solver.jl")
include("IRKGL_SIMD_Step_Functions.jl")
include("IRKNGL_SIMD_Solver.jl")
include("IRKNGL_SIMD_Step_Functions.jl")


end # module
