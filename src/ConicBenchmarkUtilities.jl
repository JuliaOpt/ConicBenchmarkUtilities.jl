
module ConicBenchmarkUtilities

using GZip

include("cbf_io.jl")
export readcbfdata, writecbfdata

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

include("moi.jl")
export cbftomoi!, moitocbf_solution#, moitocbf

end # module
