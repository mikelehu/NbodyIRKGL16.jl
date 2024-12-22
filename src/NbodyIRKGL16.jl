__precompile__()

module NbodyIRKGL16

using Reexport
@reexport using DiffEqBase
@reexport using OrdinaryDiffEq

using Parameters
using SIMD

export nbirkgl16
export IRKGLstep_fixed!, IRKGLstep_adap!
export IRKNGLstep_fixed_simpl!, IRKNGLstep_adap_simpl!
export IRKGLstep_SIMD_fixed!, IRKGLstep_SIMD_adap!
export IRKNGLstep_SIMD_fixed_simpl!, IRKNGLstep_SIMD_adap_simpl!
export PolInterp, PolInterp!
export Rdigits
export VecArray

include("./simd/VecArray_def.jl")
include("./IRKGL_Solver.jl")

include("./common/IRKGL_Coefficients.jl")
include("./common/aux_functions_SIMD.jl")
include("./common/aux_functions.jl")

include("./generic/IRKGL_gen_Step_Functions.jl")
include("./generic/IRKNGL_gen_Step_Functions.jl")

include("./simd/IRKGL_SIMD_Step_Functions.jl")
include("./simd/IRKNGL_SIMD_Step_Functions.jl")


end # module

