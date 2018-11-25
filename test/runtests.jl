
using ConicBenchmarkUtilities
using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
using GZip
using Test
# using SCS
#
# MOIU.@model(SCSModelData, (MOI.ZeroOne, MOI.Integer), (),
#     (MOI.Zeros, MOI.Reals, MOI.Nonnegatives, MOI.Nonpositives,
#         MOI.SecondOrderCone, MOI.RotatedSecondOrderCone, MOI.ExponentialCone,
#         MOI.DualExponentialCone, MOI.PowerCone, MOI.DualPowerCone,
#         MOI.PositiveSemidefiniteConeTriangle,),
#     (), (), (), (MOI.VectorOfVariables,), (MOI.VectorAffineFunction,)
#     )
#
# optimizer = MOIU.CachingOptimizer(SCSModelData{Float64}(), SCS.Optimizer(eps=1e-6, verbose=false))
using Hypatia
 MOIU.@model(HypatiaModelData,
    (),
    (
        MOI.EqualTo, MOI.GreaterThan, MOI.LessThan, MOI.Interval,
    ),
    (
        MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives,
        MOI.SecondOrderCone, MOI.RotatedSecondOrderCone,
        MOI.ExponentialCone, MOI.PowerCone, MOI.GeometricMeanCone,
        MOI.PositiveSemidefiniteConeTriangle,
        MOI.LogDetConeTriangle,
    ),
    (),
    (MOI.SingleVariable,),
    (MOI.ScalarAffineFunction,),
    (MOI.VectorOfVariables,),
    (MOI.VectorAffineFunction,),
    )

optimizer = MOIU.CachingOptimizer(HypatiaModelData{Float64}(), Hypatia.Optimizer(verbose=true))

examples = Dict{String, NamedTuple}(
    # "exA.cbf" => (objval = 0.705710, xval = [0.179894, 0.254408, 0.179894], Xval = [[0.217251 -0.25997 0.217251; -0.25997 0.31109 -0.25997; 0.217251 -0.25997 0.217251]]),
    # "exB.cbf" => (objval = 5.0, xval = [1.0, 1.0], Xval = [ones(2, 2)]),
    # "exC.cbf" => (objval = 5.098445, xval = [1.94819, 4.92228], Xval = []),
    # "exD.cbf" => (objval = 0.0, xval = [], Xval = [zeros(2, 2)]),
    "exE.cbf" => (objval = 0.0, xval = [], Xval = [zeros(2, 2)]),
)

getcbflines(file) = strip.(split(strip(read(file, String)), "\n"))

@testset "ConicBenchmarkUtilities tests" begin

@testset "example $exname" for (exname, r) in examples
    filename = joinpath(@__DIR__, "data", exname)

    # file IO
    dat = readcbfdata(filename)

    if endswith(filename, "cbf.gz")
        fd = GZip.gzopen(filename, "r")
    else
        fd = open(filename, "r")
    end
    cbflines = getcbflines(fd)
    close(fd)

    writecbfdata("example_out.cbf", dat, cbflines[1])
    @test getcbflines("example_out.cbf") == cbflines
    rm("example_out.cbf")

    # MOI IO
    MOI.empty!(optimizer)
    cbftomoi!(optimizer, dat)
    MOI.optimize!(optimizer)

    # @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.Success
    # @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FeasiblePoint
    # @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ r.objval atol=1e-4 rtol=1e-4

    @show MOI.get(optimizer, MOI.ObjectiveValue())
    (x, X) = moitocbf_solution(dat, optimizer)
    @show x
    @show X
    # @test isempty(x) == isempty(r.xval)
    # if !isempty(x)
    #     @test x ≈ r.xval atol=1e-4 rtol=1e-4
    # end
    # for j in eachindex(X)
    #     @test X[j] ≈ r.Xval[j] atol=1e-4 rtol=1e-4
    # end

    # dat2 = moitocbf(optimizer)
    #
    # MOI.empty!(optimizer)
    # cbftomoi!(optimizer, dat2)
    # MOI.optimize!(optimizer)
    #
    # @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.Success
    # @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FeasiblePoint
    # @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ r.objval atol=1e-4 rtol=1e-4
    #
    # (x2, X2) = moitocbf_solution(dat2, optimizer)
    # @test isempty(X2)
    # if !isempty(x2)
    #     @test x2[1:length(r.xval)] ≈ x atol=1e-4 rtol=1e-4
    # end
end

# @testset "integer problem" begin
#     filename = joinpath(@__DIR__, "data", "exAint.cbf")
#     dat = readcbfdata(filename)
#     cbflines = getcbflines(filename)
#     writecbfdata("example_out.cbf", dat, cbflines[1])
#     @test getcbflines("example_out.cbf") == cbflines
#     rm("example_out.cbf")
# end
#
# @testset "GZipped CBF file" begin
#     filename = joinpath(@__DIR__, "data", "exCzip.cbf.gz")
#     dat = readcbfdata(filename)
#     GZip.gzopen(filename, "r") do fd
#         cbflines = getcbflines(fd)
#         writecbfdata("example_out.cbf", dat, cbflines[1])
#         @test getcbflines("example_out.cbf") == cbflines
#         rm("example_out.cbf")
#     end
# end

end
