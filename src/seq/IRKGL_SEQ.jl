__precompile__()

module IRKGL_SEQ

using Reexport
@reexport using DiffEqBase

using Parameters
using OrdinaryDiffEq

export IRKGL_seq, IRKNGL_seq
export IRKGLstep_fixed!, IRKGLstep_adap!
export IRKNGLstep_fixed!, IRKNGLstep_adap!
export PolInterp, PolInterp!
export Rdigits

include("../common/IRKGL_Coefficients.jl")
include("../common/aux_functions.jl")
include("IRKGL_Seq_Solver.jl")
include("IRKGL_Step_Functions.jl")
include("IRKNGL_Seq_Solver.jl")
include("IRKNGL_Step_Functions.jl")

end # module
