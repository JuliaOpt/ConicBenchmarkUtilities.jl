module ConicBenchmarkUtilities

using GZip

using SparseArrays
using LinearAlgebra

# this is required because findnz does not support arrays by default in julia v1
function SparseArrays.findnz(A::AbstractMatrix)
    I = findall(!iszero, A)
    return (getindex.(I, 1), getindex.(I, 2), A[I])
end

export readcbfdata, cbftompb, mpbtocbf, writecbfdata
export remove_zero_varcones, socrotated_to_soc, remove_ints_in_nonlinear_cones, dualize

include("cbf_input.jl")
include("cbf_output.jl")
include("mpb.jl")
include("preprocess_mpb.jl")
include("convex_to_cbf.jl")
include("jump_to_cbf.jl")

end # module
