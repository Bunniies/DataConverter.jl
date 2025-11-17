function read_GHeader(path::String)
    data = open(path, "r")
    g_header = zeros(Int32, 4)
    read!(data, g_header)
    close(data)
    a = GHeader(g_header)
    return a
end

function read_CHeader(path::String; legacy::Bool=false)
    gh = read_GHeader(path)
    
    data = open(path, "r")
    seek(data, gh.hsize)
    
    a = Vector{CHeader}(undef, gh.ncorr)
    if !legacy
        aux_f = zeros(Float64, 6)
        aux_i = zeros(Int32, 4)
        theta = zeros(Float64, 6)

        for k = 1:gh.ncorr
            read!(data, aux_f)
            read!(data, theta)

            qs1 = read(data, Int32)
            if qs1 != 0
                qn1 = read(data, Int32)
                qeps1 = read(data, Float64)
                q1 = Sm(qs1, qn1, qeps1, 1)
            else
                q1 = Sm(qs1, 1)
            end

            qs2 = read(data, Int32)
            if qs2 != 0
                qn2 = read(data, Int32)
                qeps2 = read(data, Float64)
                q2 = Sm(qs2, qn2, qeps2, 1)
            else
                q2 = Sm(qs2, 1)
            end


            gs1 = read(data, Int32)
            if gs1 != 0 && gs1 != 3 && gs1 != 4
                gn1 = read(data, Int32)
                geps1 = read(data, Float64)
                g1 = Sm(gs1, gn1, geps1, 2)
            elseif gs1 == 3 || gs1 == 4
                g1 = Sm(gs1, q1.niter, q1.eps, 2)
            else
                g1 = Sm(gs1, 2)
            end
            
            gs2 = read(data, Int32)
            if gs2 != 0 && gs2 != 3 && gs2 != 4
                gn2 = read(data, Int32)
                geps2 = read(data, Float64)
                g2 = Sm(gs2, gn2, geps2, 2)
            elseif gs1 == 3 || gs1 == 4
                g2 = Sm(gs2, q2.niter, q2.eps, 2)
            else
                g2 = Sm(gs2, 2)
            end
            
            read!(data, aux_i)
            a[k] = CHeader(aux_f, aux_i, theta, [q1, q2, g1, g2])
        end
    else
        aux_f = zeros(Float64, 4)
        aux_i = zeros(Int32, 4)
        for k = 1:gh.ncorr
            read!(data, aux_f)
            read!(data, aux_i)
            a[k] = CHeader(aux_f, aux_i)
        end
    end
    close(data)
    return a
end

@doc raw"""
    read_mesons(path::String, g1::Union{String, Nothing}=nothing, g2::Union{String, Nothing}=nothing; id::Union{String, Nothing}=nothing, legacy::Bool=false)

    read_mesons(path::Vector{String}, g1::Union{String, Nothing}=nothing, g2::Union{String, Nothing}=nothing; id::Union{String, Nothing}=nothing, legacy::Bool=false)

This function read a mesons dat file at a given path and returns a vector of `CData` structures for different masses and Dirac structures.
Dirac structures `g1` and/or `g2` can be passed as string arguments in order to filter correaltors.
ADerrors id can be specified as argument. If is not specified, the `id` is fixed according to the ensemble name (example: "H400"-> id = "H400")

*For the old version (without smearing, distance preconditioning and theta) set legacy=true.

Examples:
```@example
read_mesons(path)
read_mesons(path, "G5")
read_mesons(path, nothing, "G5")
read_mesons(path, "G5", "G5")
read_mesons(path, "G5", "G5", id="H100")
read_mesons(path, "G5_d2", "G5_d2", legacy=true)
```
"""
function read_mesons(path::String, g1::Union{String, Nothing}=nothing, g2::Union{String, Nothing}=nothing; id::Union{String, Nothing}=nothing, legacy::Bool=false,
    nnoise_trunc::Union{Int64, Nothing}=nothing)    
    t1 = isnothing(g1) ? nothing : findfirst(x-> x==g1, gamma_name) - 1
    t2 = isnothing(g2) ? nothing : findfirst(x-> x==g2, gamma_name) - 1
    if isnothing(id)
        bname = basename(path)
        m = findfirst(r"[A-Z][0-9]{3}r[0-9]{3}", bname)
        id = bname[m[1:4]]
        #id = parse(Int64, bname[m[2:4]])
    end

    data = open(path, "r")
    g_header = read_GHeader(path)
    c_header = read_CHeader(path, legacy=legacy)

    ncorr = g_header.ncorr
    tvals = g_header.tvals
    nnoise = g_header.nnoise

    nnoise_trunc = isnothing(nnoise_trunc) ? nnoise : min(nnoise, nnoise_trunc)

    fsize = filesize(path)

    datsize = 4 + sum(getfield.(c_header, :dsize)) * tvals * nnoise #data_size / ncnfg
    ncfg = div(fsize - g_header.hsize - sum(getfield.(c_header, :hsize)), datsize) #(total size - header_size) / data_size

    corr_match = findall(x-> (x.type1==t1 || isnothing(t1)) && (x.type2==t2 || isnothing(t2)), c_header)

        
    seek(data, g_header.hsize + sum(c.hsize for c in c_header))

    res = Array{CData}(undef, length(corr_match))

    data_re = Array{Float64}(undef, length(corr_match), ncfg, tvals)
    data_im = zeros(length(corr_match), ncfg, tvals)
    vcfg = Array{Int32}(undef, ncfg)

    for icfg = 1:ncfg
        vcfg[icfg] = read(data, Int32)
        c=1
        for k = 1:ncorr
            if k in corr_match 
                if c_header[k].is_real == 1
                    tmp = Array{Float64}(undef, tvals*nnoise)
                    read!(data, tmp)
                    tmp2 = reshape(tmp, (nnoise, tvals))
                    tmp2 = mean(tmp2[1:nnoise_trunc, :], dims=1)
                    data_re[c, icfg, :] = tmp2[1, :]
                elseif c_header[k].is_real == 0
                    tmp = Array{Float64}(undef, 2*tvals*nnoise)
                    read!(data, tmp)
                    tmp2 = reshape(tmp, (2, nnoise, tvals))
                    tmp2 = mean(tmp2[:, 1:nnoise_trunc, :], dims=2)
                    data_re[c, icfg, :] = tmp2[1, 1, :]
                    data_im[c, icfg, :] = tmp2[2, 1, :]

                end
                c += 1
            else
                seek(data, position(data)  + c_header[k].dsize*tvals*nnoise)
            end
        

        end
    end

    for c in eachindex(corr_match)
        res[c] = CData(c_header[corr_match[c]], vcfg, data_re[c, :, :], data_im[c, :, :], id)
    end
    close(data)

    return res
end

function read_mesons(path::Vector{String}, g1::Union{String, Nothing}=nothing, g2::Union{String, Nothing}=nothing; id::Union{String, Nothing}=nothing, legacy::Bool=false,
    nnoise_trunc::Union{Int64, Nothing}=nothing)
    res = read_mesons.(path, g1, g2, id=id, legacy=legacy, nnoise_trunc=nnoise_trunc)
    nrep = length(res)
    ncorr = length(res[1])

    cdata = Vector{Vector{CData}}(undef, ncorr)
    for icorr = 1:ncorr
        cdata[icorr] = Vector{CData}(undef, nrep) 
        for r = 1:nrep
            cdata[icorr][r] = res[r][icorr]
        end
    end
    return cdata
end

@doc"""
    read_mesons_multiple_files(path::Vector{String}, g1::Union{String, Nothing}=nothing, g2::Union{String, Nothing}=nothing)

    This function reads dat files produced by job arrays, where each config dat file is saved as an independent file.
    The function reads all the dat files, respecting configuration order, and returns a CData structure where all the measurements
    on each configuration are concatenated.

    Currently this only works for a single replica. 
"""
function read_mesons_multiple_files(path::Vector{String}, g1::Union{String, Nothing}=nothing, g2::Union{String, Nothing}=nothing)
    idx_perm  = sortperm([parse(Int, match(r"cnfg(\d+)", basename(f)).captures[1]) for f in path])
    path = path[idx_perm]

    cd_aux = read_mesons(path[1], g1, g2)

    for k in eachindex(path)
        if k == 1
            continue
        end

        cd_tmp = read_mesons(path[k], g1, g2)

        for j in eachindex(cd_aux)
            cd_aux[j].re_data = vcat(cd_aux[j].re_data, cd_tmp[j].re_data)
            cd_aux[j].im_data = vcat(cd_aux[j].im_data, cd_tmp[j].im_data)
        end
        push!(cd_aux[k].vcfg, cd_tmp[k].vcfg[1])
    end

    return cd_aux
end
