
module ConicBenchmarkUtilities

using GZip
using SparseArrays
using LinearAlgebra

include("cbf_io.jl")
export readcbfdata, writecbfdata

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

include("moi.jl")
export cbftomoi!#, moitocbf

end # module
