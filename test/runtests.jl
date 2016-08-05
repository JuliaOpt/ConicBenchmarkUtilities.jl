using ConicBenchmarkUtilities
using Base.Test
using ECOS, MathProgBase, JuMP


dat = readcbfdata("example.cbf")

c, A, b, con_cones, var_cones, vartypes, dat.sense, dat.objoffset = cbftompb(dat)

@test_approx_eq c [1.0, 0.64]
@test_approx_eq A [-50.0 -31; -3.0 2.0]
@test_approx_eq b [-250.0, 4.0]
@test vartypes == [:Cont, :Cont]
@test dat.sense == :Max
@test dat.objoffset == 0.0
@test con_cones == [(:NonPos,[1]),(:NonNeg,[2])]

m = MathProgBase.ConicModel(ECOSSolver(verbose=0))
MathProgBase.loadproblem!(m, -c, A, b, con_cones, var_cones)
MathProgBase.optimize!(m)

x_sol = MathProgBase.getsolution(m)
objval = MathProgBase.getobjval(m)

mj = Model(solver=ECOSSolver(verbose=0))
@variable(mj, x[1:2])
@objective(mj, Max, x[1] + 0.64x[2])
@constraint(mj, 50x[1] + 31x[2] <= 250)
@constraint(mj, 3x[1] - 2x[2] >= -4)
status = solve(mj)

@test_approx_eq_eps x_sol getvalue(x) 1e-6
@test_approx_eq_eps -objval getobjectivevalue(mj) 1e-6
