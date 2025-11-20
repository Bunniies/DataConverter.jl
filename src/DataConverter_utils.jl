@doc"""
    ave_cdata_src(a::Vector{CData}, k1::Float64, k2::Float64)
    
    This function takes an array of CData with k1, k2 and perform the average
    over  pair of source positions (t1,t2) such that (T - t2 - 1 = t1), followed by a second 
    average over the resulting noise sources len(src_pos)/2.
    Such a setup was introduced to spread multiple source positions in a span of 0.5 from the left and right
    side of the source at thalf. 

    It returns:
        - CData if a==Vector{CData} # single replica 
        - Vector{Cdata} if a==Vector{Vector{CData}} # multiple replica
    
    where the CData attributes are taken from ‘a‘, with x0 set to thalf and 
    the re_data and im_data set to the averaged results computed as described above.
"""
function ave_cdata_src(a::Vector{CData}, k1::Float64, k2::Float64)
    a = filter(x -> getfield(getfield(x, :header), :k1) == k1 && getfield(getfield(x, :header), :k2) == k2 , a)
    if length(a) == 0
        error("No CData structure found with k1=$(k1) and k2=$(k2)")
    end

    kappas = (getfield.(getfield.(a, :header), :k1), getfield.(getfield.(a, :header), :k2))
    if all(kappas[k] != unique(kappas)[1] for k in eachindex(kappas))
        error("Mismatch in kappa values.")
    end
    
    tsrcs = getfield.(getfield.(a, :header), :x0)
    NCFG, T = size(a[1].re_data)
    thalf = Int64(T/2)
    nnoise = Int64(length(tsrcs)/2) 

    re_data_noise = Array{Float64}(undef, (NCFG, T, nnoise))
    im_data_noise = Array{Float64}(undef, (NCFG, T, nnoise))
    # println("tsrc:        ", tsrcs )
    # println("len(nnoise): ", length(nnoise), " has to be len(tsrc)/2")
    idx_noise = 1
    for i in eachindex(tsrcs)
        for j in eachindex(tsrcs)
            if j>i
                if (T - tsrcs[j] -1) == tsrcs[i] 
                    # println(tsrcs[i], " ", tsrcs[j], " ", (T - tsrcs[j] -1) )
                    re_data_noise[:,:, idx_noise] = (circshift(a[i].re_data, (0,-tsrcs[i])) + circshift(a[j].re_data, (0,T-tsrcs[j]))) / 2
                    im_data_noise[:,:, idx_noise] = (circshift(a[i].im_data, (0,-tsrcs[i])) + circshift(a[j].im_data, (0,T-tsrcs[j]))) / 2
                    idx_noise +=1
                end
            end
        end
    end

    re_mean = dropdims(mean(re_data_noise, dims=3), dims=3)
    im_mean = dropdims(mean(im_data_noise, dims=3), dims=3)
    header = deepcopy(a[1].header)
    header.x0 = thalf

    resCData = CData(header, a[1].vcfg, re_mean, im_mean, a[1].id)
    return resCData
end
function ave_cdata_src(a::Vector{Vector{CData}}, k1::Float64, k2::Float64)
    b = hcat(a...)
    nrep, _ = size(b)
    store_cd = Vector{CData}(undef, nrep)
    
    for k in 1:nrep
        store_cd[k] = ave_cdata_src(b[k,:], k1, k2)
    end
    return store_cd
end

@doc"""
    get_unique_kappas(cd::Vector{CData})
    get_unique_kappas(cd::Vector{Vector{CData}})

    Given an vector of CData structures returns the unique kappa values k1 and k2 
    corresponding to the measured masses.
    For ensemble with multiple replicas a vector of vector of CData is passed instead

    k1, k2 = get_unique_kappas(cd)
"""
function get_unique_kappas(cd::Vector{CData})
    k1_tot = []
    k2_tot = []
    for k in eachindex(cd)
        k1_aux = getfield(getfield(cd[k], :header), :k1)
        k2_aux = getfield(getfield(cd[k], :header), :k2)
        push!(k1_tot, k1_aux)
        push!(k2_tot, k2_aux)
    end
    k1 = sort(unique(k1_tot), rev=true)
    k2 = sort(unique(k2_tot), rev=true)
    return k1, k2
end
function get_unique_kappas(cd::Vector{Vector{CData}})
    k1_tot = []
    k2_tot = []
    for k in eachindex(cd)
        k1_aux = getfield(getfield.(cd[k], :header)[1], :k1)
        k2_aux = getfield(getfield.(cd[k], :header)[1], :k2)
        push!(k1_tot, k1_aux)
        push!(k2_tot, k2_aux)
    end
    k1 = sort(unique(k1_tot), rev=true)
    k2 = sort(unique(k2_tot), rev=true)
    return k1, k2 
end

@doc"""
    get_sector(k::Float64, k_list; ens_deg=true)
    
    Given a value of kappa k and a list of unique kappa values k_list
    the function returns the character l or h for light and heavy sector respectively.
    If ens_deg=false it returns also the character s for strange sector, assuming it its corresponding kappa is the second element in k_list
"""
function get_sector(k::Float64, k_list; ens_deg::Bool=true)
    if k == k_list[1]
        return "l"
    elseif ens_deg==false && k == k_list[2]
        return "s"
    else
        idx = if ens_deg
            findall(x->x==k, k_list)[1]-1
        else
            findall(x->x==k, k_list)[1]-2
        end
        return string("h",idx)
    end
end

@doc"""
    get_fname(k1, k2, k1_list, k2_list; ens_def=true)

    Given k1, k2 values of kappas and two list of uniques kappa values k1_list and k2_list
    returns the string used for saving the txt files in appropriate folders. 
    Possible outcomes are:
    - mll
    - mls (if ens_deg=false only)
    - mlh
    - msh (if ens_deg=false only)
    - mhh
"""
function get_fname(k1, k2, k1_list, k2_list; ens_deg::Bool=true)
    sec1 = get_sector(k1, k1_list, ens_deg=ens_deg)
    sec2 = get_sector(k2, k2_list, ens_deg=ens_deg)
    return string("m",sec2, sec1)
end
