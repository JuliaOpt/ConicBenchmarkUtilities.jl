module ConicBenchmarkUtilities

using GZip
import Compat: String

export readcbfdata, cbftompb
export remove_zero_varcones, socrotated_to_soc, remove_ints_in_nonlinear_cones

include("cbf_input.jl")
include("mpb.jl")
include("preprocess_mpb.jl")

end # module
