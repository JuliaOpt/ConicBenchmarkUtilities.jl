
using ConicBenchmarkUtilities
using MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
using Test

using SCS

MOIU.@model(SCSModelData, (MOI.ZeroOne, MOI.Integer), (),
    (MOI.Zeros, MOI.Reals, MOI.Nonnegatives, MOI.Nonpositives,
        MOI.SecondOrderCone, MOI.RotatedSecondOrderCone, MOI.ExponentialCone,
        MOI.DualExponentialCone, MOI.PowerCone, MOI.DualPowerCone,
        MOI.PositiveSemidefiniteConeTriangle,),
    (), (), (), (MOI.VectorOfVariables,), (MOI.VectorAffineFunction,)
    )

optimizer = MOIU.CachingOptimizer(SCSModelData{Float64}(), SCS.Optimizer(eps=1e-6, verbose=false))

# using Hypatia

# MOIU.@model(HypatiaModelData,
#     (),
#     (
#         MOI.EqualTo, MOI.GreaterThan, MOI.LessThan, MOI.Interval,
#     ),
#     (
#         MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives,
#         MOI.SecondOrderCone, MOI.RotatedSecondOrderCone,
#         MOI.ExponentialCone, MOI.PowerCone, MOI.GeometricMeanCone,
#         MOI.PositiveSemidefiniteConeTriangle,
#         MOI.LogDetConeTriangle,
#     ),
#     (),
#     (MOI.SingleVariable,),
#     (MOI.ScalarAffineFunction,),
#     (MOI.VectorOfVariables,),
#     (MOI.VectorAffineFunction,),
#     )
#
# optimizer = MOIU.CachingOptimizer(HypatiaModelData{Float64}(), Hypatia.Optimizer(verbose=false))


@testset "ConicBenchmarkUtilities tests" begin

# CBF data input/output tests
@testset "read/write CBF $filename" for filename in ("example1", "example3", "example4", "psdvaronly")
    dat = readcbfdata(filename * ".cbf")
    if startswith(filename, "example")
        comment = "# Example C.$(last(filename)) from the CBF documentation version 2"
    else
        comment = "# Generated by ConicBenchmarkUtilities.jl"
    end
    writecbfdata("example_out.cbf", dat, comment)
    @test strip(read(filename * ".cbf", String)) == strip(read("example_out.cbf", String))
    rm("example_out.cbf")

    # TODO delete
    @show filename

    MOI.empty!(optimizer)
    cbftomoi!(optimizer, dat)
    # @show optimizer
    # println(optimizer.optimizer.c)
    # println(optimizer.optimizer.A)
    # println(optimizer.optimizer.b)
    # println(optimizer.optimizer.G)
    # println(optimizer.optimizer.h)
    # println(optimizer.optimizer.cone)
    # println(optimizer.optimizer.x)

    MOI.optimize!(optimizer)
    (x, X) = moitocbf_solution(dat, optimizer)
    # @show x
    # @show X

    dat2 = moitocbf(optimizer)
    # @show dat.nvar
    # @show dat2.nvar
    # @show dat.intlist
    # @show dat2.intlist
    # @show dat.var
    # @show dat2.var
    # @show dat.psdvar
    # @show dat2.psdvar
    # @show dat.ncon
    # @show dat2.ncon
    # @show dat.con
    # @show dat2.con
    # @show dat.psdcon
    # @show dat2.psdcon
    #
    # @show dat.objacoord
    # @show dat.objfcoord
    # @show dat2.objacoord
    # @show dat.objbcoord
    # @show dat2.objbcoord
    #
    # @show dat.acoord
    # @show dat.fcoord
    # @show dat2.acoord
    # @show dat.bcoord
    # @show dat2.bcoord

    # println("\n\n")
end

# MathOptInterface conversion tests


# TODO use MOI conic tests: write then read CBF, convert to MOI, then solve and run the optimizer and result checks



end
