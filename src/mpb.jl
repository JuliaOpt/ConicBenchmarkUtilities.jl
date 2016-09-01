
const conemap = Dict("L=" => :Zero, "F" => :Free,
                     "L-" => :NonPos, "L+" => :NonNeg,
                     "Q" => :SOC, "QR" => :SOCRotated,
                     "EXP" => :ExpPrimal)
                     #"EXP*" => :ExpDual)
const conemap_rev = Dict(:Zero => "L=", :Free => "F",
                     :NonPos => "L-", :NonNeg => "L+",
                     :SOC => "Q", :SOCRotated => "QR",
                     :ExpPrimal => "EXP")
                     #, :ExpDual => "EXP*")

function cbfcones_to_mpbcones(c::Vector{Tuple{String,Int}},total)
    i = 1
    mpb_cones = Vector{Tuple{Symbol,Vector{Int}}}(0) 

    for (cname,count) in c
        conesymbol = conemap[cname]
        if conesymbol == :ExpPrimal
            @assert count == 3
            indices = i+2:-1:i
        else
            indices = i:(i+count-1)
        end
        push!(mpb_cones, (conesymbol, collect(indices)))
        i += count
    end
    @assert i == total + 1
    return mpb_cones
end

# https://github.com/JuliaLang/julia/issues/13942#issuecomment-217324812
function unzip{T<:Tuple}(A::Array{T})
    res = map(x -> x[], T.parameters)
    res_len = length(res)
    for t in A
        for i in 1:res_len
            push!(res[i], t[i])
        end
    end
    res
end

function cbftompb(dat::CBFData)

    @assert dat.nvar == sum(c->c[2],dat.var)
    @assert dat.nconstr == sum(c->c[2],dat.con)

    c = dat.objvec

    var_cones = cbfcones_to_mpbcones(dat.var, dat.nvar) 
    con_cones = cbfcones_to_mpbcones(dat.con, dat.nconstr)
    
    I, J, V = unzip(dat.acoord)
    A = sparse(I,J,-V,dat.nconstr,dat.nvar)
    b = dat.bcoord
    vartypes = fill(:Cont, dat.nvar)
    if !isempty(dat.isint)
        vartypes[dat.isint] = :Int
    end

    return c, A, b, con_cones, var_cones, vartypes, dat.sense, dat.objoffset

end

function mpbtocbf(name, c, A, b, con_cones, var_cones, vartypes, sense=:Min)

    # need to shuffle rows and columns to put them in order
    var_idx_old_to_new = zeros(Int,length(c))
    con_idx_old_to_new = zeros(Int,length(b))
    var_idx_new_to_old = zeros(Int,length(c))
    con_idx_new_to_old = zeros(Int,length(b))

    # CBF fields
    var = Vector{Tuple{String,Int}}(0)
    con = Vector{Tuple{String,Int}}(0)

    i = 1
    for (cone,idx) in var_cones
        if cone == :ExpPrimal
            @assert all(var_idx_old_to_new[idx] .== 0)
            @assert length(idx) == 3
            # MPB: (x,y,z) : y*exp(x/y) <= z
            # CBF: (z,y,x) : y*exp(x/y) <= z
            var_idx_old_to_new[idx] = i+2:-1:i
            var_idx_new_to_old[i+2:-1:i] = idx
            i += 3
        else
            for k in idx
                var_idx_old_to_new[k] = i
                var_idx_new_to_old[i] = k
                i += 1
            end
        end
        push!(var, (conemap_rev[cone],length(idx)))
    end
    @assert i - 1 == length(c)

    i = 1
    for (cone,idx) in con_cones
        if cone == :ExpPrimal
            @assert all(con_idx_old_to_new[idx] .== 0)
            @assert length(idx) == 3
            con_idx_old_to_new[idx] = i+2:-1:i
            con_idx_new_to_old[i+2:-1:i] = idx
            i += 3
        else
            for k in idx
                @assert con_idx_old_to_new[k] == 0
                con_idx_old_to_new[k] = i
                con_idx_new_to_old[i] = k
                i += 1
            end
        end
        push!(con, (conemap_rev[cone],length(idx)))
    end
    @assert i - 1 == length(b)

    objvec = c[var_idx_new_to_old]
    bcoord = b[con_idx_new_to_old]
    # MPB is b - Ax ∈ K, CBF is b + Ax ∈ K
    Acbf = -A[con_idx_new_to_old,var_idx_new_to_old]

    acoord = collect(zip(findnz(Acbf)...))::Vector{Tuple{Int,Int,Float64}}

    isint = Vector{Bool}(0)
    for i in 1:length(vartypes)
        if vartypes[i] == :Cont
            push!(isint,false)
        elseif vartypes[i] == :Int
            push!(isint,true)
        elseif vartypes[i] == :Bin
            # TODO: Check if we need to add variable bounds also
            push!(isint,true)
        else
            error("Unrecognized variable category $(vartypes[i])")
        end
    end

    return CBFData(name,sense,var,con,objvec,0.0,acoord,bcoord,isint,length(c),length(b))


end
