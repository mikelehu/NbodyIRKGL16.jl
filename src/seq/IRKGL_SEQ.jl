__precompile__()

module IRKGL_SEQ

using Reexport
@reexport using DiffEqBase

using Parameters
using OrdinaryDiffEq

export  IRKGL_seq, IRKNGL_seq

include("../IRKGL_Coefficients.jl")
include("IRKGL_Seq_Solver.jl")
include("IRKGL_Step_Functions.jl")
include("IRKNGL_Seq_Solver.jl")
include("IRKNGL_Step_Functions.jl")

end # module
