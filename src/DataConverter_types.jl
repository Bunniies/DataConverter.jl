const noise_name=["Z2", "GAUSS", "U1", "Z4"]
const gamma_name = ["G0", "G1", "G2", "G3",
"invalid", "G5", "1", "G0G1", "G0G2", "G0G3",
"G0G5", "G1G2", "G1G3", "G1G5", "G2G3", "G2G5", "G3G5", 
"G0_d1", "G1_d1", "G2_d1", "G3_d1",
"invalid", "G5_d1", "1_d1", "G0G1_d1", "G0G2_d1", "G0G3_d1",
"G0G5_d1", "G1G2_d1", "G1G3_d1", "G1G5_d1", "G2G3_d1", "G2G5_d1", "G3G5_d1",
"G0_d2", "G1_d2", "G2_d2", "G3_d2",
"invalid", "G5_d2", "1_d2", "G0G1_d2", "G0G2_d2", "G0G3_d2",
"G0G5_d2", "G1G2_d2", "G1G3_d2", "G1G5_d2", "G2G3_d2", "G2G5_d2", "G3G5_d2"]
const qs_name=["Local", "Wuppertal", "3D Gradient Flow", "Gradient Flow"]
const gs_name=["Local", "APE", "3D Wilson Flow", "Quark 3D Gradient Flow", "Quark Gradient Flow"]

mutable struct GHeader
    ncorr::Int32
    nnoise::Int32
    tvals::Int32
    noisetype::Int32

    hsize::Int32 #headersize
    GHeader(a, b, c, d) = new(a, b, c, d, 4*4)
    function GHeader(x::Vector{Int32})
        a = new(x[1], x[2], x[3], x[4], 4*4)
        return a
    end
end

mutable struct Sm
    type::Int32
    niter::Int32
    eps::Float64
    qg::Int32 #1->fermionic 2->gluonic
    Sm(t,q) = new(t, 0, 0.0, q)
    Sm(t, n, e, q) = new(t, n, e, q)
end

mutable struct CHeader
    k1::Float64
    k2::Float64
    mu1::Float64
    mu2::Float64
    dp1::Float64
    dp2::Float64
    theta1::Vector{Float64}
    theta2::Vector{Float64}
    q1::Sm
    q2::Sm
    g1::Sm
    g2::Sm
    type1::Int32
    type2::Int32
    x0::Int32
    is_real::Int32

    hsize::Int32 #headersize
    dsize::Int32 #datasize / (nnoise * T * ncfg)

    function CHeader(aux_f::Vector{Float64}, aux_i::Vector{Int32}, theta::Vector{Float64}, sm_par::Vector{Sm})
        a = new()
        a.k1 = aux_f[1]
        a.k2 = aux_f[2]
        a.mu1 = aux_f[3]
        a.mu2 = aux_f[4]
        a.dp1 = aux_f[5]
        a.dp2 = aux_f[6]
        a.type1 = aux_i[1]
        a.type2 = aux_i[2]
        a.x0 = aux_i[3]
        a.is_real = aux_i[4]

        a.theta1 = theta[1:3]
        a.theta2 = theta[4:6]
        a.q1 = sm_par[1]
        a.q2 = sm_par[2]
        a.g1 = sm_par[3]
        a.g2 = sm_par[4]
        a.hsize = 8*12 + 4*8 #without smearing parameters niter, neps
        if a.q1.type != 0
            a.hsize += 8+4
        end
        if a.q2.type != 0
            a.hsize += 8+4
        end
        if a.g1.type != 0 && a.g1.type != 3 && a.g1.type != 4
            a.hsize += 8+4
        end
        if a.g2.type != 0 && a.g2.type != 3 && a.g2.type != 4
            a.hsize += 8+4
        end
        a.dsize = 16 - 8* a.is_real 
        
        return a
    end
    function CHeader(aux_f::Vector{Float64}, aux_i::Vector{Int32})
        a = new()
        a.k1 = aux_f[1]
        a.k2 = aux_f[2]
        a.mu1 = aux_f[3]
        a.mu2 = aux_f[4]
        a.dp1 = 0.0
        a.dp2 = 0.0
        a.type1 = aux_i[1]
        a.type2 = aux_i[2]
        a.x0 = aux_i[3]
        a.is_real = aux_i[4]

        a.theta1 = zeros(3)
        a.theta2 = zeros(3)
        a.q1 = Sm(0, 1)
        a.q2 = Sm(0, 1)
        a.g1 = Sm(0, 2)
        a.g2 = Sm(0, 2)
        a.hsize = 8*4 + 4*4
        a.dsize = 16 - 8* a.is_real 

        return a
    end
end
function Base.:(==)(a::CHeader, b::CHeader)
    for s in [:k1, :k2, :mu1, :mu2, :dp1, :dp2, :type1, :type2, :x0, :is_real, :theta1, :theta2]
        if getfield(a, s) != getfield(b, s)
            return false
        end
    end
    return true
end
#Base.:(!=)(a::CHeader, b::CHeader) = !(a == b)

mutable struct CData
    header::CHeader
    vcfg::Array{Int32}
    re_data::Array{Float64}
    im_data::Array{Float64}
    id::String
    CData(a, b, c, d, e) = new(a, b, c, d, e)

end
Base.copy(a::CData) = CData(a.header, a.vcfg, a.re_data, a.im_data, a.id)

function Base.show(io::IO, a::CData)
    g1 = gamma_name[a.header.type1+1]
    g2 = gamma_name[a.header.type2+1]
    print(io, typeof(a), "(")
    print(io, "id: ", a.id )
    print(io, " channel: ", g1, "-", g2, " kappas: (", a.header.k1, ", ", a.header.k2,")", " source: ", a.header.x0, ")")
end


