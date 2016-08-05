type CBFData
    name::String
    sense::Symbol
    var::Vector{Tuple{String,Int}}
    con::Vector{Tuple{String,Int}}
    objvec::Vector{Float64}
    objoffset::Float64
    acoord::Vector{Tuple{Int,Int,Float64}} # linear coefficients
    bcoord::Vector{Float64} # linear offsets
    isint::Vector{Bool}
    nvar::Int
    nconstr::Int
end

CBFData() = CBFData("",:xxx,Vector{Tuple{String,Int}}(0),Vector{Tuple{String,Int}}(0),Vector{Float64}(),0.0,Vector{Tuple{Int,Int,Float64}}(0),Vector{Float64}(),Vector{Bool}(0),0,0)

function readcbfdata(filename)

    if endswith(filename,"cbf.gz")
        fd = gzopen(filename,"r")
    else
        @assert endswith(filename, "cbf")
        fd = open(filename,"r")
    end

    dat = CBFData()
    dat.name = split(basename(filename),".")[1]

    while !eof(fd)
        line = readline(fd)
        startswith(line,"#") && continue # comments
        length(line) == 1 && continue # blank lines
        
        # new block
        
        if startswith(line,"VER")
            nextline = readline(fd)
            @assert startswith(nextline,"1") || startswith(nextline,"2") 
            continue
        end

        if startswith(line,"OBJSENSE")
            nextline = readline(fd)
            if strip(nextline) == "MIN"
                dat.sense = :Min
            else
                dat.sense = :Max
            end
            continue
        end

        if startswith(line,"VAR")
            nextline = readline(fd)
            totalvars, lines = split(nextline)
            totalvars = parse(Int,strip(totalvars))
            lines = parse(Int,strip(lines))
            varcnt = 0

            for k in 1:lines
                nextline = readline(fd)
                cone, size = split(nextline)
                size = parse(Int,strip(size))
                push!(dat.var, (cone, size))
                varcnt += size
            end
            @assert totalvars == varcnt
            dat.nvar = varcnt
            continue
        end

        if startswith(line, "INT")
            nextline = readline(fd)
            intvar = parse(Int,strip(nextline))
            dat.isint = falses(dat.nvar)
            # ignore integer variables
            for k in 1:intvar
                nextline = readline(fd)
                idx = parse(Int,strip(nextline))
                dat.isint[idx+1] = true
            end
            continue
        end

        if startswith(line,"CON")
            nextline = readline(fd)
            totalconstr, lines = split(nextline)
            totalconstr = parse(Int,strip(totalconstr))
            lines = parse(Int,strip(lines))
            constrcnt = 0

            for k in 1:lines
                nextline = readline(fd)
                cone, size = split(nextline)
                size = parse(Int,strip(size))
                push!(dat.con, (cone, size))
                constrcnt += size
            end
            @assert totalconstr == constrcnt
            dat.nconstr = constrcnt
            continue
        end

        if startswith(line,"PSDVAR") || startswith(line,"PSDCON")
            error("Input problem is an SDP")
            break
        end

        if startswith(line,"OBJACOORD")
            nextline = readline(fd)
            nnz = parse(Int,strip(nextline))
            dat.objvec = zeros(dat.nvar)
            for k in 1:nnz
                nextline = readline(fd)
                i, val = split(strip(nextline))
                dat.objvec[parse(Int,i)+1] = float(val)
            end
            continue
        end

        if startswith(line,"OBJBCOORD")
            nextline = readline(fd)
            dat.objoffset = float(strip(nextline))
            warn("Instance has objective offset")
        end

        if startswith(line,"ACOORD")
            nextline = readline(fd)
            nnz = parse(Int,strip(nextline))
            for k in 1:nnz
                nextline = readline(fd)
                i, j, val = split(strip(nextline))
                push!(dat.acoord, (parse(Int,i)+1,parse(Int,j)+1,parse(Float64,val)))
            end
        end

        if startswith(line,"BCOORD")
            nextline = readline(fd)
            nnz = parse(Int,strip(nextline))
            dat.bcoord = zeros(dat.nconstr)
            for k in 1:nnz
                nextline = readline(fd)
                i, val = split(strip(nextline))
                dat.bcoord[parse(Int,i)+1] = float(val)
            end
        end


    end

    return dat
end


