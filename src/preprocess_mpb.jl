# Some MPB conic solvers don't support all equivalent representatons,
# although they should.

# A few utilities to normalize the MPB representation

function remove_zero_varcones(c, A, b, con_cones, var_cones, vartypes)
    old_to_new_idx = zeros(Int,length(c))
    last_idx = 1
    new_varcones = Vector{Tuple{Symbol,Vector{Int}}}(0) 
    for (cname, cidx) in var_cones
        if cname != :Zero
            for i in cidx
                old_to_new_idx[i] = last_idx
                last_idx += 1
            end
            push!(new_varcones,(cname, old_to_new_idx[cidx]))
        else
            # dropped
        end
    end
    keep_indices = find(old_to_new_idx)
    c = c[keep_indices]
    A = A[:,keep_indices]
    vartypes = vartypes[keep_indices]
    var_cones = new_varcones

    return c, A, b, con_cones, var_cones, vartypes
end


function socrotated_to_soc(c, A, b, con_cones, var_cones, vartypes)

    con_cones = copy(con_cones)
    var_cones = copy(var_cones)
    vartypes = copy(vartypes)
    c = copy(c)
    b = copy(b)
    I, J, V = findnz(A)
    nslack = 0
    # introduce slack variables and put them into SOCRotated cones
    for i in 1:length(con_cones)
        cname, cidx = con_cones[i]
        if cname == :SOCRotated
            for j in cidx
                nslack += 1
                push!(I, j)
                push!(J, nslack+length(c))
                push!(V, 1.0)
                push!(vartypes,:Cont)
            end
            con_cones[i] = (:Zero, cidx)
            push!(var_cones, (:SOCRotated, (length(c)+nslack-length(cidx)+1):(length(c)+nslack)))
        end
    end
    append!(c, zeros(nslack))
    A = sparse(I,J,V,size(A,1),size(A,2)+nslack)

    # new rows to add to constraint matrix
    I = Int[]
    J = Int[]
    V = Float64[]
    rowidx = 1
    for i in 1:length(var_cones)
        cname, cidx = var_cones[i]
        if cname == :SOCRotated
            var_cones[i] = (:Free,cidx)
            # (y,z,x) in RSOC <=> (y+z,y-z,sqrt(2)*x) in SOC
            push!(I, rowidx)
            push!(J, cidx[1])
            push!(V, -1.0)
            push!(I, rowidx)
            push!(J, cidx[2])
            push!(V, -1.0)
            
            push!(I, rowidx+1)
            push!(J, cidx[1])
            push!(V, -1.0)
            push!(I, rowidx+1)
            push!(J, cidx[2])
            push!(V, 1.0)
            
            append!(I, (rowidx+2):(rowidx+length(cidx)-1))
            append!(J, cidx[3:end])
            append!(V, fill(-sqrt(2), length(cidx)-2))
            push!(con_cones, (:SOC, (size(A,1)+rowidx):(size(A,1)+rowidx+length(cidx)-1)))
            rowidx += length(cidx)
            append!(b, zeros(length(cidx)))
        end
    end

    A = vcat(A, sparse(I,J,V, rowidx-1, size(A,2)))

    return c, A, b, con_cones, var_cones, vartypes
end


