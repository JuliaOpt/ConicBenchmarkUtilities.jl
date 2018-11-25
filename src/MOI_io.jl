
# get CBF vector cone string from MOI vector cone object
moitocbf_cone(::MOI.EqualTo{Float64}) = "L="
moitocbf_cone(::MOI.GreaterThan{Float64}) = "L+"
moitocbf_cone(::MOI.LessThan{Float64}) = "L-"
moitocbf_cone(::MOI.Zeros) = "L="
moitocbf_cone(::MOI.Reals) = "F"
moitocbf_cone(::MOI.Nonnegatives) = "L+"
moitocbf_cone(::MOI.Nonpositives) = "L-"
moitocbf_cone(::MOI.SecondOrderCone) = "Q"
moitocbf_cone(::MOI.RotatedSecondOrderCone) = "QR"
moitocbf_cone(::MOI.ExponentialCone) = "EXP"
moitocbf_cone(::MOI.DualExponentialCone) = "EXP*"

# get vector length of triangle of symmetric n*n matrix
psd_len(n) = div(n*(n+1), 2)

# get vector index for (i,j) term in n*n matrix (lower triangle, row major)
function mattovecidx(i, j)
    if i < j
        (i,j) = (j,i)
    end
    return div((i-1)*i, 2) + j
end

# get (i,j) matrix indices (lower triangle, row major) for kth term in vector
function vectomatidx(k)
    i = div(1 + isqrt(8k-7), 2)
    j = k - div(i*(i-1), 2)
    return (i,j)
end

# fill a PSD matrix from a PSD vector (lower triangle, row major)
function vectomat!(m::Matrix{Float64}, v::Vector{Float64}, side::Int)
    k = 1
    for i in 1:side, j in 1:i
        if j == i
            m[i,j] = v[k]
        else
            m[i,j] = m[j,i] = v[k]
        end
        k += 1
    end
    return m
end

# convert a CBFData object to an MOI model by modifying an empty MOI model
function cbftomoi!(model::MOI.ModelLike, dat::CBFData)
    @assert dat.nvar == (isempty(dat.var) ? 0 : sum(c -> c[2], dat.var))
    @assert dat.ncon == (isempty(dat.con) ? 0 : sum(c -> c[2], dat.con))

    if !MOI.is_empty(model)
        error("MOI model object is not empty")
    end

    # objective sense
    if dat.sense == :Min
        MOI.set(model, MOI.ObjectiveSense(), MOI.MinSense)
    elseif dat.sense == :Max
        MOI.set(model, MOI.ObjectiveSense(), MOI.MaxSense)
    else
        error("objective sense $(dat.sense) not recognized")
    end

    # non-PSD variables
    x = MOI.add_variables(model, dat.nvar) # MOI non-PSD variables
    for j in dat.intlist # integer constraints
        MOI.add_constraint(model, MOI.SingleVariable(x[j]), MOI.Integer())
    end

    # PSD variables
    npsdvar = 0 # number of columns
    psdvaridxs = UnitRange[] # column ranges
    for i in eachindex(dat.psdvar)
        ilen = psd_len(dat.psdvar[i])
        push!(psdvaridxs, npsdvar .+ (1:ilen))
        npsdvar += ilen
    end
    X = MOI.add_variables(model, npsdvar) # MOI PSD variables

    # objective function
    objterms = [MOI.ScalarAffineTerm(v, x[a]) for ((a,), v) in dat.objacoord]
    for ((a, b, c), v) in dat.objfcoord
        if b != c
            v += v # scale off-diagonals
        end
        push!(objterms, MOI.ScalarAffineTerm(v, X[psdvaridxs[a][mattovecidx(b, c)]]))
    end
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(objterms, dat.objbcoord))

    # non-PSD variable cones
    k = 0
    for (cname, clen) in dat.var
        if cname == "F" # ignore free cones
            k += clen
            continue
        end

        if cname in ("EXP", "EXP*") # reverse order
            @assert clen == 3
            F = MOI.VectorOfVariables(x[k+3,k+2,k+1])
            S = (cname == "EXP") ? MOI.ExponentialCone() : MOI.DualExponentialCone()
        else
            F = MOI.VectorOfVariables(x[k .+ (1:clen)])

            if cname in ("POWER", "POWER*") # TODO power cones equivalent to MathOptInterface's 3D-PowerCone definition
                error("Power cones not currently supported by ConicBenchmarkUtilities")
                # if clen != 3 || length(params) != 2
                #     error("currently cannot handle power cones that aren't equivalent to MathOptInterface's 3D-PowerCone definition (or its dual cone)")
                # end
                # sigma = sum(params) # TODO get power cone parameters vector
                # exponent = params[1]/sigma
                # if cname == "POWER"
                #     S = MOI.PowerCone(exponent)
                # else
                #     S = MOI.DualPowerCone(exponent)
                # end
            elseif cname == "L="
                S = MOI.Zeros(clen)
            elseif cname == "L-"
                S = MOI.Nonpositives(clen)
            elseif cname == "L+"
                S = MOI.Nonnegatives(clen)
            elseif cname == "Q"
                @assert clen >= 3
                S = MOI.SecondOrderCone(clen)
            elseif cname == "QR"
                @assert clen >= 3
                S = MOI.RotatedSecondOrderCone(clen)
            else
                error("cone type $cname is not recognized")
            end
        end

        @assert MOI.output_dimension(F) == MOI.dimension(S)
        MOI.add_constraint(model, F, S)
        k += clen
    end
    @assert k == dat.nvar

    # PSD variable cones
    for i in eachindex(dat.psdvar)
        F = MOI.VectorOfVariables(X[psdvaridxs[i]])
        S = MOI.PositiveSemidefiniteConeTriangle(dat.psdvar[i])
        @assert MOI.output_dimension(F) == MOI.dimension(S)
        MOI.add_constraint(model, F, S)
    end

    # non-PSD constraints
    if length(dat.con) > 0
        conterms = [Vector{MOI.ScalarAffineTerm{Float64}}() for i in 1:dat.ncon]
        for ((a, b), v) in dat.acoord # non-PSD terms
            push!(conterms[a], MOI.ScalarAffineTerm(v, x[b]))
        end
        for ((a, b, c, d), v) in dat.fcoord # PSD terms
            if c != d
                v += v # scale off-diagonals
            end
            idx = psdvaridxs[b][mattovecidx(c, d)]
            push!(conterms[a], MOI.ScalarAffineTerm(v, X[idx]))
        end

        conoffs = zeros(dat.ncon)
        for ((a,), v) in dat.bcoord # constants
            conoffs[a] = v
        end

        k = 0
        for (cname, clen) in dat.con
            if cname == "F" # ignore free cones
                k += clen
                continue
            end

            if cname in ("EXP", "EXP*") # reverse order
                @assert clen == 3
                vats = [MOI.VectorAffineTerm(4-l, t) for l in 1:clen for t in conterms[k+l]]
                F = MOI.VectorAffineFunction(vats, conoffs[[k+3,k+2,k+1]])
                S = (cname == "EXP") ? MOI.ExponentialCone() : MOI.DualExponentialCone()
            else
                vats = [MOI.VectorAffineTerm(l, t) for l in 1:clen for t in conterms[k+l]]
                F = MOI.VectorAffineFunction(vats, conoffs[k .+ (1:clen)])

                if cname in ("POWER", "POWER*") # TODO power cones equivalent to MathOptInterface's 3D-PowerCone definition
                    error("Power cones not currently supported by ConicBenchmarkUtilities")
                    # if clen != 3 || length(params) != 2
                    #     error("currently cannot handle power cones that aren't equivalent to MathOptInterface's 3D-PowerCone definition (or its dual cone)")
                    # end
                    # sigma = sum(params) # TODO get power cone parameters vector
                    # exponent = params[1]/sigma
                    # if cname == "POWER"
                    #     S = MOI.PowerCone(exponent)
                    # else
                    #     S = MOI.DualPowerCone(exponent)
                    # end
                elseif cname == "L="
                    S = MOI.Zeros(clen)
                elseif cname == "L-"
                    S = MOI.Nonpositives(clen)
                elseif cname == "L+"
                    S = MOI.Nonnegatives(clen)
                elseif cname == "Q"
                    @assert clen >= 3
                    S = MOI.SecondOrderCone(clen)
                elseif cname == "QR"
                    @assert clen >= 3
                    S = MOI.RotatedSecondOrderCone(clen)
                else
                    error("cone type $cname is not recognized")
                end
            end

            @assert MOI.output_dimension(F) == MOI.dimension(S)
            MOI.add_constraint(model, F, S)
            k += clen
        end
        @assert k == dat.ncon
    end

    # PSD constraints
    if length(dat.psdcon) > 0
        npsdcon = 0 # number of rows
        psdconidxs = UnitRange[] # row ranges
        for i in eachindex(dat.psdcon)
            ilen = psd_len(dat.psdcon[i])
            push!(psdconidxs, npsdcon .+ (1:ilen))
            npsdcon += ilen
        end
        @assert npsdcon == last(last(psdconidxs))

        psdconterms = [Vector{MOI.ScalarAffineTerm{Float64}}() for i in 1:npsdcon]
        for ((a, b, c, d), v) in dat.hcoord # PSD terms
            idx = first(psdconidxs[a]) - 1 + mattovecidx(c, d)
            push!(psdconterms[idx], MOI.ScalarAffineTerm(v, x[b]))
        end

        psdconoffs = zeros(npsdcon)
        for ((a, b, c), v) in dat.dcoord # constants
            idx = first(psdconidxs[a]) - 1 + mattovecidx(b, c)
            psdconoffs[idx] = v
        end

        k = 0
        for i in eachindex(dat.psdcon)
            idxs = psdconidxs[i]
            vats = [MOI.VectorAffineTerm(l-k, t) for l in idxs for t in psdconterms[l]]
            F = MOI.VectorAffineFunction(vats, psdconoffs[idxs])
            S = MOI.PositiveSemidefiniteConeTriangle(dat.psdcon[i])
            @assert MOI.output_dimension(F) == MOI.dimension(S)
            MOI.add_constraint(model, F, S)
            k += length(idxs)
        end
        @assert k == npsdcon
    end

    return model
end

# convert an MOI model to a CBFData object
# NOTE does not use variable cones
# NOTE need to apply split interval bridge and PSD square to triangle bridge
function moitocbf(model)
    if MOI.is_empty(model)
        error("MOI model object is empty")
    end

    # helper functions for MOI
    getmodelcons(F, S) = MOI.get(model, MOI.ListOfConstraintIndices{F, S}())
    getconfun(conidx) = MOI.get(model, MOI.ConstraintFunction(), conidx)
    getconset(conidx) = MOI.get(model, MOI.ConstraintSet(), conidx)

    dat = CBFData() # CBF data object to be returned

    # objective sense
    objsense = MOI.get(model, MOI.ObjectiveSense())
    if objsense == MOI.MinSense
        dat.sense = :Min
    elseif objsense == MOI.MaxSense
        dat.sense = :Max
    else
        @assert objsense == MOI.FeasibilitySense
        println("Objective sense is feasibility sense")
        dat.sense = :Feas
    end

    # variables
    dat.nvar = MOI.get(model, MOI.NumberOfVariables())
    dat.intlist = [getconfun(ci).variable.value for ci in getmodelcons(
        MOI.SingleVariable, MOI.Integer)]
    dat.var = [("F", dat.nvar)] # not using variable cones: only free variables

    # objective function
    F = MOI.get(model, MOI.ObjectiveFunctionType())
    obj = MOI.get(model, MOI.ObjectiveFunction{F}())
    if F == MOI.SingleVariable
        dat.objacoord = [((F.variable.value,), 1.0)]
    else
        @assert F == MOI.ScalarAffineFunction{Float64}
        dat.objacoord = [((t.variable_index.value,), t.coefficient) for t in obj.terms]
        dat.objbcoord = obj.constant
    end

    # MOI scalar constraints
    for S in (
        MOI.EqualTo{Float64},
        MOI.GreaterThan{Float64},
        MOI.LessThan{Float64},
        )
        for ci in getmodelcons(MOI.SingleVariable, S)
            dat.ncon += 1
            push!(dat.con, (moitocbf_cone(si), 1))
            push!(dat.acoord, ((dat.ncon, getconfun(ci).variable.value), 1.0))
            push!(dat.bcoord, ((dat.ncon,), -getconset(ci).value))
        end

        for ci in getmodelcons(MOI.ScalarAffineFunction{Float64}, S)
            dat.ncon += 1
            push!(dat.con, (moitocbf_cone(si), 1))
            fi = getconfun(ci)
            for vt in fi.terms
                push!(dat.acoord, ((dat.ncon, vt.variable_index.value), vt.coefficient))
            end
            push!(dat.bcoord, ((dat.ncon,), fi.constant - getconset(ci).value))
        end
    end

    # MOI vector constraints
    for S in (
        MOI.Zeros,
        MOI.Reals,
        MOI.Nonnegatives,
        MOI.Nonpositives,
        MOI.SecondOrderCone,
        MOI.RotatedSecondOrderCone,
        MOI.ExponentialCone,
        MOI.DualExponentialCone,
        )
        for ci in getmodelcons(MOI.VectorOfVariables, S)
            vars = getconfun(ci).variables
            if S in (MOI.ExponentialCone, MOI.DualExponentialCone) # reverse order
                reverse!(vars)
            end
            for vj in vars
                dat.ncon += 1
                push!(dat.acoord, ((dat.ncon, vj.value), 1.0))
            end
            si = getconset(ci)
            push!(dat.con, (moitocbf_cone(si), MOI.dimension(si)))
        end

        for ci in getmodelcons(MOI.VectorAffineFunction{Float64}, S)
            fi = getconfun(ci)
            si = getconset(ci)
            dim = MOI.dimension(si)
            if S in (MOI.ExponentialCone, MOI.DualExponentialCone) # reverse order
                @assert dim == 3
                for vt in fi.terms
                    idx = dat.ncon + 4 - vt.output_index
                    push!(dat.acoord, ((idx, vt.scalar_term.variable_index.value),
                        vt.scalar_term.coefficient))
                end
                for row in [3,2,1]
                    dat.ncon += 1
                    push!(dat.bcoord, ((dat.ncon,), fi.constants[row]))
                end
            else
                for vt in fi.terms
                    idx = dat.ncon + vt.output_index
                    push!(dat.acoord, ((idx, vt.scalar_term.variable_index.value),
                        vt.scalar_term.coefficient))
                end
                for row in 1:dim
                    dat.ncon += 1
                    push!(dat.bcoord, ((dat.ncon,), fi.constants[row]))
                end
            end
            push!(dat.con, (moitocbf_cone(si), dim))
        end
    end

    # TODO power cone constraints

    # PSD constraints
    npsdcon = 0
    for ci in getmodelcons(MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeTriangle)
        npsdcon += 1
        side = getconset(ci).side_dimension
        vars = getconfun(ci).variables
        k = 0
        for i in 1:side, j in 1:i
            k += 1
            push!(dat.hcoord, ((npsdcon, vars[k].value, i, j), 1.0))
        end
        @assert k == length(vars)
        push!(dat.psdcon, side)
    end

    for ci in getmodelcons(MOI.VectorAffineFunction{Float64}, MOI.PositiveSemidefiniteConeTriangle)
        npsdcon += 1
        side = getconset(ci).side_dimension
        fi = getconfun(ci)
        for vt in fi.terms
            (i,j) = vectomatidx(vt.output_index)
            push!(dat.hcoord, ((npsdcon, vt.scalar_term.variable_index.value, i, j),
                vt.scalar_term.coefficient))
        end
        k = 0
        for i in 1:side, j in 1:i
            k += 1
            if !iszero(fi.constants[k])
                push!(dat.dcoord, ((npsdcon, i, j), fi.constants[k]))
            end
        end
        @assert k == MOI.output_dimension(fi)
        push!(dat.psdcon, side)
    end

    return dat
end

# convert an MOI solution to a (scalar, PSD) solution pair
function moitocbf_solution(dat::CBFData, model::MOI.ModelLike)
    vars = MOI.get(model, MOI.ListOfVariableIndices())
    soln = MOI.get(model, MOI.VariablePrimal(), vars)

    scalar_soln = soln[1:dat.nvar]

    matrix_soln = Vector{Matrix{Float64}}()
    k = dat.nvar
    for i in eachindex(dat.psdvar)
        side = dat.psdvar[i]
        veclen = psd_len(side)
        veci = soln[k .+ (1:veclen)]
        mati = Matrix{Float64}(undef, side, side)
        vectomat!(mati, veci, side)
        push!(matrix_soln, mati)
        k += veclen
    end

    return (scalar_soln, matrix_soln)
end
