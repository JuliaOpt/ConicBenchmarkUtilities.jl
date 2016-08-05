
const conemap = Dict("L=" => :Zero, "F" => :Free,
                     "L-" => :NonPos, "L+" => :NonNeg,
                     "Q" => :SOC, "QR" => :SOCRotated)

function cbfcones_to_mpbcones(c::Vector{Tuple{String,Int}},total)
    i = 1
    mpb_cones = Vector{Tuple{Symbol,Vector{Int}}}(0) 

    for (cname,count) in c
        conesymbol = conemap[cname]
        indices = i:(i+count-1)
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
