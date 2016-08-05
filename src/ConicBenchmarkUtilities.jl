module ConicBenchmarkUtilities

using GZip
import Compat: String

export readcbfdata, cbftompb

include("cbf_input.jl")
include("mpb.jl")

end # module
