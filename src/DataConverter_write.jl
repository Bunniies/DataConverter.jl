function write_2pt_data(cd::CData, path::String, name::String, rep::String)

    open(joinpath(path,name), "w") do io
    # open(path, "w") do io
        writedlm(io, ["nb replicas: $(length(rep))"])
        writedlm(io, ["nb confs replica $(rep): $(length(cd.vcfg))"] )
        writedlm(io, " ")
        writedlm(io, ["# kappa 1: $(cd.header.k1)"])
        writedlm(io, ["# kappa 2: $(cd.header.k2)"])
        writedlm(io, ["# theta 1: $(cd.header.theta1)"])
        writedlm(io, ["# theta 2: $(cd.header.theta2)"])
        writedlm(io, " ")

        T = collect(Int64, 0:size(cd.re_data,2)-1)
        for cnfg in cd.vcfg
            writedlm(io, ["# rep $(rep) conf $(cnfg)"])
            writedlm(io, " ")
            for t in T
                writedlm(io,  [t '\t' cd.re_data[cnfg, t+1] '\t' cd.im_data[cnfg, t+1]])
                # writedlm(io,  zip(t, cd.re_data[cnfg, t+1],  cd.im_data[cnfg, t+1]) )
            end
            writedlm(io, " ")
        end

    end
end


function write_2pt_data(cd::Vector{CData}, path::String, name::String, rep::Vector{String})
    if length(cd) != length(rep)
        error("Mismatch in number of replicas!")
    end

    open(joinpath(path,name), "w") do io
        writedlm(io, ["nb replicas: $(length(rep))"])
        for k in eachindex(rep)
            writedlm(io, ["nb confs replica $(rep[k]): $(length(cd[k].vcfg))"] )
        end
        writedlm(io, " ")
        writedlm(io, ["# kappa 1: $(cd[1].header.k1)"])
        writedlm(io, ["# kappa 2: $(cd[1].header.k2)"])
        writedlm(io, ["# theta 1: $(cd[1].header.theta1)"])
        writedlm(io, ["# theta 2: $(cd[1].header.theta2)"])
        writedlm(io, " ")

        T = collect(Int64, 0:size(cd[1].re_data,2)-1)

        for k in eachindex(rep)
            for cnfg in cd[k].vcfg
                writedlm(io, ["# rep $(rep[k]) conf $(cnfg)"])
                writedlm(io, " ")
                for t in T
                    writedlm(io,  [t '\t' cd[k].re_data[cnfg, t+1] '\t' cd[k].im_data[cnfg, t+1]])
                end
                writedlm(io, " ")
            end
        end
    end
end