__precompile__()

module IRKGL_SIMD

using Reexport
@reexport using DiffEqBase

using LinearAlgebra
using Parameters
using OrdinaryDiffEq
using SIMD

export  IRKGL_simd, IRKNGL_simd, IRKAlgorithm
export  VecArray, floatType, IRKGLCoefficients_adap

include("IRKGL_SIMD_Coefficients.jl")
include("IRKGL_SIMD_Solver.jl")
include("IRKGL_SIMD_Step_Functions.jl")
include("IRKNGL_SIMD_Solver.jl")
include("IRKNGL_SIMD_Step_Functions.jl")



end # module
