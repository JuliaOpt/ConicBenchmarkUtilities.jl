
# get MOI cone object from CBF cone string and dimension
function cbftomoi_cones(cname, clen::Int)
    if cname == "L="
        # equal to
        return MOI.Zeros(clen)
    elseif cname == "F"
        # free
        return MOI.Reals(clen)
    elseif cname == "L-"
        # nonpositive
        return MOI.Nonpositives(clen)
    elseif cname == "L+"
        # nonnegative
        return MOI.Nonnegatives(clen)
    elseif cname == "Q"
        # second-order cone
        @assert clen >= 3
        return MOI.SecondOrderCone(clen)
    elseif cname == "QR"
        # rotated second-order cone
        @assert clen >= 3
        return MOI.RotatedSecondOrderCone(clen)
    elseif cname == "EXP"
        # exponential
        @assert clen == 3
        return MOI.ExponentialCone()
    elseif cname == "EXP*"
        # dual exponential
        @assert clen == 3
        return MOI.DualExponentialCone()
    # elseif cname in ("POWER", "POWER*")
    #     # power (parametrized)
    #     if clen != 3 || length(params) != 2
    #         error("currently cannot handle power cones that aren't equivalent to MathOptInterface's 3D-PowerCone definition (or its dual cone)")
    #     end
    #     sigma = sum(params)
    #     exponent = params[1]/sigma
    #     if cname == "POWER"
    #         return MOI.PowerCone(exponent)
    #     else
    #         return MOI.DualPowerCone(exponent)
    #     end
    else
        error("cone type $cname is not recognized")
    end
end

# get vector length of triangle of symmetric n*n matrix
psd_len(n) = div(n*(n+1), 2)

# get vector index (lower triangle, row major) for (i,j) term in n*n matrix
function mattovecidx(i, j)
    if i < j
        (i,j) = (j,i)
    end
    return div((i-1)*i, 2) + j
end

# convert an MOI solution to a CBF solution
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

# fill a PSD matrix from a PSD vector (lower triangle, row major)
function vectomat!(m::Matrix{Float64}, v::Vector{Float64}, side::Int)
    k = 1
    for i in 1:side, j in i:side
        if j == i
            m[i,j] = v[k]
        else
            m[i,j] = m[j,i] = v[k]
        end
        k += 1
    end
    return m
end

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
    append!(objterms, MOI.ScalarAffineTerm(v,
        X[psdvaridxs[a][mattovecidx(b, c)]]
        ) for ((a, b, c), v) in dat.objfcoord)
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(objterms, dat.objoffset))

    # non-PSD variable cones
    k = 0
    for (cname, clen) in dat.var
        S = cbftomoi_cones(cname, clen)
        if !(S isa MOI.Reals)
            F = MOI.VectorOfVariables(x[k .+ (1:clen)])
            @assert MOI.output_dimension(F) == MOI.dimension(S)
            MOI.add_constraint(model, F, S)
        end
        k += clen
    end
    @assert k == dat.nvar

    # PSD variable cones
    for i in eachindex(dat.psdvar)
        F = MOI.VectorOfVariables(X[psdvaridxs[i]])
        S = MOI.PositiveSemidefiniteConeTriangle(dat.psdvar[i])
        # println()
        # println(F)
        # println(S)
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
            idx = psdvaridxs[b][mattovecidx(c, d)]
            push!(conterms[a], MOI.ScalarAffineTerm(v, X[idx]))
        end

        conoffs = zeros(dat.ncon)
        for ((a,), v) in dat.bcoord # constants
            conoffs[a] = v
        end

        k = 0
        for (cname, clen) in dat.con
            vats = [MOI.VectorAffineTerm(l, t) for l in 1:clen for t in conterms[k+l]]
            F = MOI.VectorAffineFunction(vats, conoffs[k .+ (1:clen)])
            S = cbftomoi_cones(cname, clen)
            # println()
            # println(cname)
            # println(F)
            # println(S)
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
            # println()
            # println(F)
            # println(S)
            @assert MOI.output_dimension(F) == MOI.dimension(S)
            MOI.add_constraint(model, F, S)
            k += length(idxs)
        end
        @assert k == npsdcon
    end

    return model
end



# function mpbtocbf(name, c, A, b, con_cones, var_cones, vartypes, sense=:Min)
#     num_scalar_var = 0
#     for (cone, idx) in var_cones
#         if cone != :SDP
#             num_scalar_var += length(idx)
#         end
#     end
#     num_scalar_con = 0
#     for (cone, idx) in con_cones
#         if cone != :SDP
#             num_scalar_con += length(idx)
#         end
#     end
#
#     # need to shuffle rows and columns to put them in order
#     var_idx_old_to_new = zeros(Int, length(c))
#     con_idx_old_to_new = zeros(Int, length(b))
#     var_idx_new_to_old = zeros(Int, num_scalar_var)
#     con_idx_new_to_old = zeros(Int, num_scalar_con)
#
#     # CBF fields
#     var = Vector{Tuple{String, Int}}()
#     con = Vector{Tuple{String, Int}}()
#
#     i = 1
#     for (cone, idx) in var_cones
#         if cone == :ExpPrimal
#             @assert all(var_idx_old_to_new[idx] .== 0)
#             @assert length(idx) == 3
#             # MPB: (x,y,z) : y*exp(x/y) <= z
#             # CBF: (z,y,x) : y*exp(x/y) <= z
#             var_idx_old_to_new[idx] = i+2:-1:i
#             var_idx_new_to_old[i+2:-1:i] = idx
#             i += 3
#             push!(var, (conemap_rev[cone], length(idx)))
#         elseif cone != :SDP
#             for k in idx
#                 var_idx_old_to_new[k] = i
#                 var_idx_new_to_old[i] = k
#                 i += 1
#             end
#             push!(var, (conemap_rev[cone], length(idx)))
#         end
#     end
#     @assert i - 1 == num_scalar_var
#
#     i = 1
#     for (cone, idx) in con_cones
#         if cone == :ExpPrimal
#             @assert all(con_idx_old_to_new[idx] .== 0)
#             @assert length(idx) == 3
#             con_idx_old_to_new[idx] = i+2:-1:i
#             con_idx_new_to_old[i+2:-1:i] = idx
#             i += 3
#             push!(con, (conemap_rev[cone], length(idx)))
#         elseif cone != :SDP
#             for k in idx
#                 @assert con_idx_old_to_new[k] == 0
#                 con_idx_old_to_new[k] = i
#                 con_idx_new_to_old[i] = k
#                 i += 1
#             end
#             push!(con, (conemap_rev[cone], length(idx)))
#         end
#     end
#     @assert i - 1 == num_scalar_con
#
#     objacoord = collect(zip(findnz(sparse(c[var_idx_new_to_old]))...))
#     bcoord = collect(zip(findnz(sparse(b[con_idx_new_to_old]))...))
#     # MPB is b - Ax ∈ K, CBF is b + Ax ∈ K
#     Acbf = -A[con_idx_new_to_old,var_idx_new_to_old]
#
#     acoord = collect(zip(findnz(Acbf)...))::Vector{Tuple{Int,Int,Float64}}
#
#     intlist = Int[]
#     for i in 1:length(vartypes)
#         if var_idx_old_to_new[i] == 0 && vartypes[i] != :Cont
#             error("CBF format does not support integer restrictions on PSD variables")
#         end
#         if vartypes[i] == :Cont
#         elseif vartypes[i] == :Int
#             push!(intlist,var_idx_old_to_new[i])
#         elseif vartypes[i] == :Bin
#             # TODO: Check if we need to add variable bounds also
#             push!(intlist,var_idx_old_to_new[i])
#         else
#             error("Unrecognized variable category $(vartypes[i])")
#         end
#     end
#
#     psdvar = Int[]
#     psdcon = Int[]
#
#     # Map from MPB linear variable index to (psdvar,i,j)
#     psdvar_idx_old_to_new = fill((0,0,0), length(c))
#     # Map from MPB linear constraint index to (psdcon,i,j)
#     psdcon_idx_old_to_new = fill((0,0,0), length(b))
#
#     for (cone, idx) in var_cones
#         if cone == :SDP
#             y = length(idx)
#             conedim = round(Int, sqrt(0.25 + 2y) - 0.5)
#             push!(psdvar, conedim)
#             k = 1
#             for i in 1:conedim, j in i:conedim
#                 psdvar_idx_old_to_new[idx[k]] = (length(psdvar), i, j)
#                 k += 1
#             end
#             @assert length(idx) == k - 1
#         end
#     end
#
#     for (cone, idx) in con_cones
#         if cone == :SDP
#             y = length(idx)
#             conedim = round(Int, sqrt(0.25 + 2y) - 0.5)
#             push!(psdcon, conedim)
#             k = 1
#             for i in 1:conedim, j in i:conedim
#                 psdcon_idx_old_to_new[idx[k]] = (length(psdcon), i, j)
#                 k += 1
#             end
#             @assert length(idx) == k - 1
#         end
#     end
#
#     objfcoord = Vector{Tuple{Int,Int,Int,Float64}}()
#     for i in 1:length(c)
#         if c[i] != 0.0 && psdvar_idx_old_to_new[i] != (0,0,0)
#             varidx, vari, varj = psdvar_idx_old_to_new[i]
#             scale = (vari == varj) ? 1.0 : sqrt(2)
#             push!(objfcoord, (varidx, vari, varj, c[i]/scale))
#         end
#     end
#     dcoord = Vector{Tuple{Int,Int,Int,Float64}}()
#     for i in 1:length(b)
#         if b[i] != 0.0 && psdcon_idx_old_to_new[i] != (0,0,0)
#             conidx, coni, conj = psdcon_idx_old_to_new[i]
#             scale = (coni == conj) ? 1.0 : sqrt(2)
#             push!(dcoord, (conidx, coni, conj, b[i]/scale))
#         end
#     end
#
#     A_I, A_J, A_V = findnz(A)
#     fcoord = Vector{Tuple{Int,Int,Int,Int,Float64}}()
#     hcoord = Vector{Tuple{Int,Int,Int,Int,Float64}}()
#
#     for (i,j,v) in zip(A_I,A_J,A_V)
#         if psdvar_idx_old_to_new[j] != (0,0,0)
#             if psdcon_idx_old_to_new[i] != (0,0,0)
#                 error("CBF format does not allow PSD variables to appear in affine expressions defining PSD constraints")
#             end
#             newrow = con_idx_old_to_new[i]
#             @assert newrow != 0
#             varidx, vari, varj = psdvar_idx_old_to_new[j]
#             scale = (vari == varj) ? 1.0 : sqrt(2)
#             push!(fcoord, (newrow, varidx, vari, varj, -v/scale))
#         elseif psdcon_idx_old_to_new[i] != (0,0,0)
#             newcol = var_idx_old_to_new[j]
#             conidx, coni, conj = psdcon_idx_old_to_new[i]
#             scale = (coni == conj) ? 1.0 : sqrt(2)
#             push!(hcoord, (conidx, newcol, coni, conj, -v/scale))
#         end
#     end
#
#     return CBFData(name,sense,var,psdvar,con,psdcon,objacoord,objfcoord,0.0,fcoord,acoord,bcoord,hcoord,dcoord,intlist,num_scalar_var,num_scalar_con)
# end
