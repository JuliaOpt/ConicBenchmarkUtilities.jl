
module ConicBenchmarkUtilities

mutable struct CBFData
    sense::Symbol
    nvar::Int
    intlist::Vector{Int}
    var::Vector{Tuple{String, Int}}
    psdvar::Vector{Int}
    ncon::Int
    con::Vector{Tuple{String, Int}}
    psdcon::Vector{Int}
    objfcoord::Vector{Tuple{NTuple{3, Int}, Float64}}
    objacoord::Vector{Tuple{NTuple{1, Int}, Float64}}
    objbcoord::Float64
    fcoord::Vector{Tuple{NTuple{4, Int}, Float64}}
    acoord::Vector{Tuple{NTuple{2, Int}, Float64}}
    bcoord::Vector{Tuple{NTuple{1, Int}, Float64}}
    hcoord::Vector{Tuple{NTuple{4, Int}, Float64}}
    dcoord::Vector{Tuple{NTuple{3, Int}, Float64}}
end
CBFData() = CBFData(:xxx, 0, [], [], [], 0, [], [], [], [], 0.0, [], [], [], [], [])
export CBFData

using GZip
include("file_io.jl")
export readcbfdata, writecbfdata

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
include("MOI_io.jl")
export cbftomoi!, moitocbf_solution, moitocbf

end
