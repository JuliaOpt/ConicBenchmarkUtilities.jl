module ConicBenchmarkUtilities

using GZip
import Compat: String

export readcbfdata, cbftompb, remove_zero_varcones, socrotated_to_soc

include("cbf_input.jl")
include("mpb.jl")
include("preprocess_mpb.jl")

end # module
