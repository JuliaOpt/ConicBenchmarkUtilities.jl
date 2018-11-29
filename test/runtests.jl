
using ConicBenchmarkUtilities
using MathOptInterface
MOI = MathOptInterface
MOIU = MathOptInterface.Utilities
MOIB = MOI.Bridges
using GZip
using Test
using SCS

MOIU.@model(ModelData,
    (MOI.Integer,),
    (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan),
    (MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives, MOI.SecondOrderCone, MOI.ExponentialCone, MOI.PositiveSemidefiniteConeTriangle),
    (),
    (MOI.SingleVariable,),
    (MOI.ScalarAffineFunction,),
    (MOI.VectorOfVariables,),
    (MOI.VectorAffineFunction,)
    )

optimizer = MOIB.RSOC{Float64}(MOIU.CachingOptimizer(ModelData{Float64}(), SCS.Optimizer(eps=1e-6, verbose=false)))

examples = Dict{String, NamedTuple}(
    "exA.cbf" => (objval = 0.705710, xval = [0.179894, 0.254408, 0.179894], Xval =
        [[0.217251 -0.25997 0.217251; -0.25997 0.31109 -0.25997; 0.217251 -0.25997 0.217251]]),
    "exB.cbf" => (objval = 5.0, xval = [1.0, 1.0], Xval = [ones(2, 2)]),
    "exC.cbf" => (objval = 5.098445, xval = [1.94819, 4.92228], Xval = []),
    "exD.cbf" => (objval = 0.0, xval = [], Xval = [zeros(2, 2)]),
    "exE.cbf" => (objval = -4.808370, xval = [-3.66096, 3.40549, 3.15003, 1.14741], Xval = []),
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

    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.Success
    @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FeasiblePoint
    @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ r.objval atol=1e-4 rtol=1e-4

    (x, X) = moitocbf_solution(dat, optimizer)
    @test isempty(x) == isempty(r.xval)
    if !isempty(x)
        @test x ≈ r.xval atol=1e-4 rtol=1e-4
    end
    for j in eachindex(X)
        @test X[j] ≈ r.Xval[j] atol=1e-4 rtol=1e-4
    end

    dat2 = moitocbf(optimizer)

    MOI.empty!(optimizer)
    cbftomoi!(optimizer, dat2)
    MOI.optimize!(optimizer)

    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.Success
    @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FeasiblePoint
    @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ r.objval atol=1e-4 rtol=1e-4

    (x2, X2) = moitocbf_solution(dat2, optimizer)
    @test isempty(X2)
    if !isempty(x2)
        @test x2[1:length(r.xval)] ≈ x atol=1e-4 rtol=1e-4
    end
end

@testset "integer problem" begin
    filename = joinpath(@__DIR__, "data", "exAint.cbf")

    # file IO
    dat = readcbfdata(filename)
    cbflines = getcbflines(filename)
    writecbfdata("example_out.cbf", dat, cbflines[1])
    @test getcbflines("example_out.cbf") == cbflines
    rm("example_out.cbf")

    # MOI IO
    MOI.empty!(optimizer)
    cbftomoi!(optimizer, dat)
    intcons = MOI.get(optimizer, MOI.ListOfConstraintIndices{MOI.SingleVariable, MOI.Integer}())
    ints = [var.variable.value for var in MOI.get(optimizer, MOI.ConstraintFunction(), intcons)]
    @test length(ints) == 2
    @test sort(ints) == sort(dat.intlist)
end

@testset "GZipped CBF file" begin
    filename = joinpath(@__DIR__, "data", "exCzip.cbf.gz")
    dat = readcbfdata(filename)
    GZip.gzopen(filename, "r") do fd
        cbflines = getcbflines(fd)
        writecbfdata("example_out.cbf", dat, cbflines[1])
        @test getcbflines("example_out.cbf") == cbflines
        rm("example_out.cbf")
    end
end

# tests taken from MOI.Test contlinear.jl and contconic.jl and intconic.jl
@testset "MathOptInterface linear and conic tests" begin
@testset "linear feasible" begin
    # min -3x - 2y - 4z
    # st    x +  y +  z == 3
    #            y +  z == 2
    #       x>=0 y>=0 z>=0
    # Opt obj = -11, soln x = 1, y = 0, z = 2
    MOI.empty!(optimizer)

    # set up MOI model
    v = MOI.add_variables(optimizer, 3)

    vc = MOI.add_constraint(optimizer, MOI.VectorOfVariables(v), MOI.Nonnegatives(3))
    c = MOI.add_constraint(optimizer, MOI.VectorAffineFunction(MOI.VectorAffineTerm.([1,1,1,2,2], MOI.ScalarAffineTerm.(1.0, [v;v[2];v[3]])), [-3.0,-2.0]), MOI.Zeros(2))

    MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([-3.0, -2.0, -4.0], v), 0.0))
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MinSense)

    # convert to CBF, write, read, and convert to MOI
    dat = moitocbf(optimizer)
    writecbfdata("example_out.cbf", dat, "# test")
    dat2 = readcbfdata("example_out.cbf")
    rm("example_out.cbf")
    MOI.empty!(optimizer)
    cbftomoi!(optimizer, dat2)

    # optimize and check solution
    MOI.optimize!(optimizer)

    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.Success
    @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FeasiblePoint
    @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ -11 atol=1e-3 rtol=1e-3

    (x, X) = moitocbf_solution(dat2, optimizer)
    @test x ≈ [1, 0, 2] atol=1e-3 rtol=1e-3
    @test isempty(X)
end

@testset "linear infeasible" begin
    # min x
    # s.t. 2x+y <= -1
    # x,y >= 0
    # infeasible
    MOI.empty!(optimizer)

    # set up MOI model
    x = MOI.add_variable(optimizer)
    y = MOI.add_variable(optimizer)

    c = MOI.add_constraint(optimizer, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([2.0,1.0], [x,y]), 0.0), MOI.LessThan(-1.0))
    bndx = MOI.add_constraint(optimizer, MOI.SingleVariable(x), MOI.GreaterThan(0.0))
    bndy = MOI.add_constraint(optimizer, MOI.SingleVariable(y), MOI.GreaterThan(0.0))

    MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, x)], 0.0))
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MinSense)

    # convert to CBF, write, read, and convert to MOI
    dat = moitocbf(optimizer)
    writecbfdata("example_out.cbf", dat, "# test")
    dat2 = readcbfdata("example_out.cbf")
    rm("example_out.cbf")
    MOI.empty!(optimizer)
    cbftomoi!(optimizer, dat2)

    # optimize and check solution
    MOI.optimize!(optimizer)

    @test MOI.get(optimizer, MOI.ResultCount()) >= 1
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.Success
    @test MOI.get(optimizer, MOI.DualStatus()) == MOI.InfeasibilityCertificate

    (x, X) = moitocbf_solution(dat2, optimizer)
    @test length(x) == 2
    @test isempty(X)
end

@testset "rotated SOC infeasible" begin
    # min 0
    # s.t.
    #      x ≤ 1
    #      y = 1/2
    #      z ≥ 2
    #      z^2 ≤ 2x*y
    # infeasible
    MOI.empty!(optimizer)

    # set up MOI model
    x = MOI.add_variables(optimizer, 3)

    vc1 = MOI.add_constraint(optimizer, MOI.SingleVariable(x[1]), MOI.LessThan(1.0))
    vc2 = MOI.add_constraint(optimizer, MOI.SingleVariable(x[2]), MOI.EqualTo(0.5))
    vc3 = MOI.add_constraint(optimizer, MOI.SingleVariable(x[3]), MOI.GreaterThan(2.0))
    rsoc = MOI.add_constraint(optimizer, MOI.VectorOfVariables(x), MOI.RotatedSecondOrderCone(3))

    MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(Float64[0, 0, 0], x), 0.0))
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MinSense)

    # convert to CBF, write, read, and convert to MOI
    dat = moitocbf(optimizer)
    writecbfdata("example_out.cbf", dat, "# test")
    dat2 = readcbfdata("example_out.cbf")
    rm("example_out.cbf")
    MOI.empty!(optimizer)
    cbftomoi!(optimizer, dat2)

    # optimize and check solution
    MOI.optimize!(optimizer)

    @test MOI.get(optimizer, MOI.ResultCount()) >= 1
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.Success
    @test MOI.get(optimizer, MOI.DualStatus()) == MOI.InfeasibilityCertificate

    (x, X) = moitocbf_solution(dat2, optimizer)
    @test length(x) == 3
    @test isempty(X)
end

@testset "primal exponential cone" begin
    # min x + y + z
    #  st  y e^(x/y) <= z, y > 0 (i.e (x, y, z) are in the exponential primal cone)
    #      x == 1
    #      y == 2
    MOI.empty!(optimizer)

    # set up MOI model
    v = MOI.add_variables(optimizer, 3)

    vc = MOI.add_constraint(optimizer, MOI.VectorAffineFunction{Float64}(MOI.VectorOfVariables(v)), MOI.ExponentialCone())
    cx = MOI.add_constraint(optimizer, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, v[1])], 0.), MOI.EqualTo(1.))
    cy = MOI.add_constraint(optimizer, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, v[2])], 0.), MOI.EqualTo(2.))

    MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, v), 0.0))
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MinSense)

    # convert to CBF, write, read, and convert to MOI
    dat = moitocbf(optimizer)
    writecbfdata("example_out.cbf", dat, "# test")
    dat2 = readcbfdata("example_out.cbf")
    rm("example_out.cbf")
    MOI.empty!(optimizer)
    cbftomoi!(optimizer, dat2)

    # optimize and check solution
    MOI.optimize!(optimizer)

    @test MOI.get(optimizer, MOI.ResultCount()) >= 1
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.Success
    @test MOI.get(optimizer, MOI.DualStatus()) == MOI.FeasiblePoint
    @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ 3 + 2exp(1/2) atol=1e-3 rtol=1e-3

    (x, X) = moitocbf_solution(dat2, optimizer)
    @test x ≈ [1, 2, 2exp(1/2)] atol=1e-3 rtol=1e-3
    @test isempty(X)
end

@testset "PSD cone" begin
    # min X[1,1] + X[2,2]    max y
    #     X[2,1] = 1         [0   y/2     [ 1  0
    #                         y/2 0    <=   0  1]
    #     X >= 0              y free
    # Optimal solution:
    #     ⎛ 1   1 ⎞
    # X = ⎜       ⎟           y = 2
    #     ⎝ 1   1 ⎠
    MOI.empty!(optimizer)

    # set up MOI model
    X = MOI.add_variables(optimizer, 3)

    cX = MOI.add_constraint(optimizer, MOI.VectorAffineFunction{Float64}(MOI.VectorOfVariables(X)), MOI.PositiveSemidefiniteConeTriangle(2))
    c = MOI.add_constraint(optimizer, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, X[2])], 0.0), MOI.EqualTo(1.0))

    MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(1.0, [X[1], X[end]]), 0.0))
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MinSense)

    # convert to CBF, write, read, and convert to MOI
    dat = moitocbf(optimizer)
    writecbfdata("example_out.cbf", dat, "# test")
    dat2 = readcbfdata("example_out.cbf")
    rm("example_out.cbf")
    MOI.empty!(optimizer)
    cbftomoi!(optimizer, dat2)

    # optimize and check solution
    MOI.optimize!(optimizer)

    @test MOI.get(optimizer, MOI.ResultCount()) >= 1
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.Success
    @test MOI.get(optimizer, MOI.DualStatus()) == MOI.FeasiblePoint
    @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ 2 atol=1e-3 rtol=1e-3

    (x, X) = moitocbf_solution(dat2, optimizer)
    @test x ≈ ones(3) atol=1e-3 rtol=1e-3
    @test isempty(X)
end

@testset "integer SOC (no solve)" begin
    # min 0x - 2y - 1z
    #  st  x            == 1
    #      x >= ||(y,z)||
    #      (y,z) integer
    MOI.empty!(optimizer)

    # set up MOI model
    (x,y,z) = MOI.add_variables(optimizer, 3)

    MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([-2.0,-1.0], [y,z]), 0.0))
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MinSense)

    ceq = MOI.add_constraint(optimizer, MOI.VectorAffineFunction([MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x))], [-1.0]), MOI.Zeros(1))
    csoc = MOI.add_constraint(optimizer, MOI.VectorOfVariables([x,y,z]), MOI.SecondOrderCone(3))

    bin1 = MOI.add_constraint(optimizer, MOI.SingleVariable(y), MOI.Integer())
    bin2 = MOI.add_constraint(optimizer, MOI.SingleVariable(z), MOI.Integer())

    # convert to CBF, write, read, and convert to MOI
    dat = moitocbf(optimizer)
    writecbfdata("example_out.cbf", dat, "# test")
    dat2 = readcbfdata("example_out.cbf")
    rm("example_out.cbf")
    MOI.empty!(optimizer)
    cbftomoi!(optimizer, dat2)

    @test MOI.get(optimizer, MOI.NumberOfConstraints{MOI.SingleVariable, MOI.Integer}()) == 2
end
end
end
