__precompile__()

module IRKGL_SIMD

using Reexport
@reexport using DiffEqBase

using Parameters
using OrdinaryDiffEq
using SIMD

export  IRKGL_simd, IRKNGL_simd
export  VecArray


include("../IRKGL_Coefficients.jl")
include("VecArray_def.jl")
include("IRKGL_SIMD_Solver.jl")
include("IRKGL_SIMD_Step_Functions.jl")
include("IRKNGL_SIMD_Solver.jl")
include("IRKNGL_SIMD_Step_Functions.jl")


end # module
