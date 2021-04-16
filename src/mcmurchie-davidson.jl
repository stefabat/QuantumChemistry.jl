

"""
    hermite_expansion(imax::Int, jmax::Int, Kαβ::Real, PA::Real, PB::Real, p::Real)

Compute the Hermite expansion coefficients for a 1-dimensional Cartesian overlap distribution
using a two-term recursion relation.
"""
function hermite_expansion(imax::Int, jmax::Int, Kαβ::Real, PA::Real, PB::Real, p::Real)

    # initialize array of expansion coefficients
    # E = OffsetArray(zeros(imax+jmax+1,imax+1,jmax+1),-1,-1,-1)
    E = zeros(imax+jmax+1,imax+1,jmax+1)

    # enter recursion
    for j = 0:jmax
        for i = 0:imax
            for t = i+j:-1:0
                if t > 0
                    if i > 0
                        # E[t,i,j] += (1/(2*p*t)) * i * E[t-1, i-1, j]
                        E[t+1,i+1,j+1] += (1/(2*p*t)) * i * E[t, i, j+1]
                    end
                    if j > 0
                        # E[t,i,j] += (1/(2*p*t)) * j * E[t-1, i, j-1]
                        E[t+1,i+1,j+1] += (1/(2*p*t)) * j * E[t, i+1, j]
                    end
                else
                    if i == j == 0
                        # E[t,i,j] = Kαβ
                        E[t+1,i+1,j+1] = Kαβ
                    elseif j == 0
                        # E[t,i,j] = PA * E[0, i-1, j] + E[1, i-1, j]
                        E[t+1,i+1,j+1] = PA * E[1, i, j+1] + E[2, i, j+1]
                    else
                        # E[t,i,j] = PB * E[0, i, j-1] + E[1, i, j-1]
                        E[t+1,i+1,j+1] = PB * E[1, i+1, j] + E[2, i+1, j]
                    end
                end
            end
        end
    end

    return E
end



"""
    hermite_integral(tmax::Int, umax::Int, vmax::Int, p::Real, RPC::StaticVector{3})

Compute the integral of an Hermite Gaussian divided by the Coulomb operator.
"""
function hermite_integral(tmax::Int, umax::Int, vmax::Int, p::Real, RPC::AbstractVector)

    # initialize array
    R = OffsetArray(zeros(tmax+umax+vmax+1,tmax+1,umax+1,vmax+1),-1,-1,-1,-1)

    # compute auxiliary integrals Rⁿ₀₀₀
    for n = 0:tmax+umax+vmax
        R[n,0,0,0] = (-2.0*p)^n * boys(n,p*RPC⋅RPC)
    end

    # transfer angular momentum n to t
    for t=0:max-1
        for n=0:tmax+umax+vmax-t-1
            R[n,t+1,0,0] = RPC[1]*R[n+1, t, 0, 0]
            if t > 0
                R[n,t+1,0,0] += t*R[n+1, t-1, 0, 0]
            end
        end
    end

    # transfer angular momentum n to u
    for u=0:umax-1
        for t=0:tmax
            for n=0:tmax+umax+vmax-t-u-1
                R[n,t,u+1,0] = RPC[2]*R[n+1, t, u, 0]
                if u > 0
                    R[n,t,u+1,0] += u*R[n+1, t, u-1, 0]
                end
            end
        end
    end

    # transfer angular momentum n to v
    for v=0:vmax-1
        for u=0:umax
            for t=0:tmax
                for n=0:tmax+umax+vmax-t-u-v-1
                    R[n,t,u,v+1] = RPC[3]*R[n+1, t, u, v]
                    if v > 0
                    R[n,t,u,v+1] += v*R[n+1, t, u, v-1]
                    end
                end
            end
        end
    end

    return R

end





"""
    overlap(α::Real, la::Int, RA::AbstractVector{Real}, β::Real, lb::Int, RB::AbstractVector{Real})

Compute the overlap integral matrix between two primitive Cartesian shells centered on
`RA` and `RB`, with exponents `α` and `β` and angular momenta `la` and `lb`.
The centers `RA` and `RB` are expected to be 3-dimensional vectors.
"""
function overlap(α::Real, la::Int, RA::AbstractVector{T},
                 β::Real, lb::Int, RB::AbstractVector{T}) where {T<:Real}

    # precomputing all required quantities
    p = α + β
    μ = (α * β)/(α + β)
    RP = (RA.*α .+ RB.*β)./(α + β)
    RAB = RA .- RB; RPA = RP .- RA; RPB = RP .- RB
    Kαβ = exp.(-μ.*RAB.^2)

    # number of Cartesian primitive functions in shell
    Nla = (la+1)*(la+2)÷2
    Nlb = (lb+1)*(lb+2)÷2
    S = zeros(Nla,Nlb)

    for (b,(jx,jy,jz)) in enumerate(get_ijk(lb))
        for (a,(ix,iy,iz)) in enumerate(get_ijk(la))
            Ex = Etij(0, ix, jx, Kαβ[1], RPA[1], RPB[1], p)
            Ey = Etij(0, iy, jy, Kαβ[2], RPA[2], RPB[2], p)
            Ez = Etij(0, iz, jz, Kαβ[3], RPA[3], RPB[3], p)
            @inbounds S[a,b] = Ex * Ey * Ez
        end
    end

    return ((π / p)^1.5) .* S
end




"""Compute the kinetic energy integral -0.5*<Ga|∇^2|Gb> of two Cartesian PGFs."""
function kinetic(α::Real, ikm::NTuple{3,Int}, RA::NTuple{3,Float64},
                 β::Real, jln::NTuple{3,Int}, RB::NTuple{3,Float64})

    # extract the angular momentum elementwise
    (i,k,m) = ikm
    (j,l,n) = jln

    # precompute common terms
    dαSab = 2*α*overlap(α, (i,k,m), RA, β, (j,l,n), RB)
    fα2 = 4*α^2

    # Dij^2 * Skl * Smn
    Dij = i*(i-1) * overlap(α, (i-2,k,m), RA, β, (j,l,n), RB) -
          (2*i+1) * dαSab +
              fα2 * overlap(α, (i+2,k,m), RA, β, (j,l,n), RB)

    # Sij * Dkl^2 * Smn
    Dkl = k*(k-1) * overlap(α, (i,k-2,m), RA, β, (j,l,n), RB) -
          (2*k+1) * dαSab +
              fα2 * overlap(α, (i,k+2,m), RA, β, (j,l,n), RB)

    # Sij * Skl * Dmn^2
    Dmn = m*(m-1) * overlap(α, (i,k,m-2), RA, β, (j,l,n), RB) -
          (2*m+1) * dαSab +
              fα2 * overlap(α, (i,k,m+2), RA, β, (j,l,n), RB)

    return -0.5 * (Dij + Dkl + Dmn)
end


"""Compute the electron-nuclear attraction energy integral <Ga|1/r|Gb> of two PGFs."""
function attraction(α::Real, ikm::NTuple{3,Int}, RA::NTuple{3,Float64},
                    β::Real, jln::NTuple{3,Int}, RB::NTuple{3,Float64},
                                                 RC::NTuple{3,Float64})

    # precomputing all required quantities
    (i,k,m) = ikm; (j,l,n) = jln
    p = α + β
    RP = (α.*RA + β.*RB) ./ p
    RAB = RA - RB; RPA = RP - RA; RPB = RP - RB; RPC = RP - RC
    μ = (α * β)/(α + β)
    Kαβ = exp.(-μ .* RAB.^2)

    vne = 0.0
    for t = 0:i+j
        Eij = Etij(t, i, j, Kαβ[1], RPA[1], RPB[1], p)
        for u = 0:k+l
            Ekl = Etij(u, k, l, Kαβ[2], RPA[2], RPB[2], p)
            for v = 0:m+n
                Emn = Etij(v, m, n, Kαβ[3], RPA[3], RPB[3], p)
                vne +=  Eij * Ekl * Emn * Rtuv(t, u, v, 0, p, RPC)
            end
        end
    end

    return 2.0*π*vne/p
end


"""Compute the two-electron integral <Ga(r1)Gb(r1)|1/r12|Gc(r2)Gd(r2)>."""
function repulsion(α::Real, ikm1::NTuple{3,Int}, RA::NTuple{3,Float64},
                   β::Real, jln1::NTuple{3,Int}, RB::NTuple{3,Float64},
                   γ::Real, ikm2::NTuple{3,Int}, RC::NTuple{3,Float64},
                   δ::Real, jln2::NTuple{3,Int}, RD::NTuple{3,Float64})

    # precompute quantities for electron 1
    (i1,k1,m1) = ikm1; (j1,l1,n1) = jln1
    RAB = RA - RB
    p = α + β
    RP = (α.*RA + β.*RB) ./ p
    Kαβ = exp.(-α*β/p.*RAB.^2)
    RPA = RP - RA; RPB = RP - RB

    # precompute quantities for electron 2
    (i2,k2,m2) = ikm2; (j2,l2,n2) = jln2
    RCD = RC - RD
    q = γ + δ
    RQ = (γ.*RC + δ.*RD) ./ q
    Kγδ = exp.(-γ*δ/q.*RCD.^2)
    RQC = RQ - RC; RQD = RQ - RD

    # precompute quantities for auxiliary integral R
    RPQ = RP - RQ
    ξ = (p * q)/(p + q)

    vee = 0.0
    for t = 0:i1+j1
        Eij1 = Etij(t, i1, j1, Kαβ[1], RPA[1], RPB[1], p)
        for u = 0:k1+l1
            Ekl1 = Etij(u, k1, l1, Kαβ[2], RPA[2], RPB[2], p)
            for v = 0:m1+n1
                Emn1 = Etij(v, m1, n1, Kαβ[3], RPA[3], RPB[3], p)
                for τ = 0:i2+j2
                    Eij2 = Etij(τ, i2, j2, Kγδ[1], RQC[1], RQD[1], q)
                    for ν = 0:k2+l2
                        Ekl2 = Etij(ν, k2, l2, Kγδ[2], RQC[2], RQD[2], q)
                        for ϕ = 0:m2+n2
                            Emn2 = Etij(ϕ, m2, n2, Kγδ[3], RQC[3], RQD[3], q)
                            vee += Eij1 * Ekl1 * Emn1 * Eij2 * Ekl2 * Emn2 *
                                   Rtuv(t+τ, u+ν, v+ϕ, 0, ξ, RPQ) * (-1)^(τ+ν+ϕ)
                        end
                    end
                end
            end
        end
    end

    return (vee*2.0*π^2.5)/(p*q*sqrt(p+q))
end


# TODO: implement the boys function instead of relying on GSL
using GSL: sf_hyperg_1F1

"""Compute the Boys fucntion F_n(x)."""
function boys(n::Int, x::Real)
    return sf_hyperg_1F1(n+0.5, n+1.5, -x) / (2.0*n+1.0)
end

# This does not work I believe due to missing flag during standard
# compilation of GSL library from Arch Linux repo
# ERROR: could not load library "/usr/lib/libgsl.so"
# /usr/lib/libgsl.so: undefined symbol: cblas_ctrmv
#function _gsl_sf_hyperg_1F1(a::AbstractFloat, b::AbstractFloat, x::AbstractFloat)
#    return ccall( (:gsl_sf_hyperg_1F1, "/usr/lib/libgsl.so"),
#               Cdouble, (Cdouble, Cdouble, Cdouble), a, b, x)
#end


"""Compute the double factorial n!! of an integer number."""
function dfactorial(n::Int)
    if n == 0
        return 1.0
    elseif iseven(n) && n > 0
        return reduce(*,n:-2:2)[end]
    elseif isodd(n) && n > 0
        return reduce(*,n:-2:1)[end]
    elseif isodd(n) && n < 0
        return 1/reduce(*,n+2:2:1)[end]
    else
        error("n!! undefined for even negative n values")
    end
end


# using Memoize

"""
    Etij(t::Int ,i::Int, j::Int, Kαβx::Real, XPA::Real, XPB::Real, p::Real)

Compute the Hermite expansion coefficients for a 1-dimensional Cartesian overlap distribution
using a two-term recursion relation.
"""
function Etij(t::Int ,i::Int, j::Int, Kαβx::Real, XPA::Real, XPB::Real, p::Real)
# @memoize function Etij(t::Int ,i::Int, j::Int, Kαβx::Real, XPA::Real, XPB::Real, α::Real, β::Real)

    # enter recursion
    if t < 0 || t > i+j
        return 0.0
    elseif t == 0
        if i == j == 0
            return Kαβx
        elseif j == 0
            return XPA * Etij(0, i-1, j, Kαβx, XPA, XPB, p) +
                         Etij(1, i-1, j, Kαβx, XPA, XPB, p)
        else
            return XPB * Etij(0, i, j-1, Kαβx, XPA, XPB, p) +
                         Etij(1, i, j-1, Kαβx, XPA, XPB, p)
        end
    else
        return (1/(2*p*t)) * (i * Etij(t-1, i-1, j, Kαβx, XPA, XPB, p) +
                              j * Etij(t-1, i, j-1, Kαβx, XPA, XPB, p) )
    end
end




"""Compute the integral of an Hermite Gaussian divided by the Coulomb operator."""
function Rtuv(t::Int, u::Int, v::Int, n::Int, p::Real, RPC::NTuple{3,Float64})
    if t == u == v == 0
        return (-2.0*p)^n * boys(n,p*abs(sum(RPC.*RPC)))
    elseif u == v == 0
        if t > 1
            return  (t-1)*Rtuv(t-2, u, v, n+1, p, RPC) +
                   RPC[1]*Rtuv(t-1, u, v, n+1, p, RPC)
        else
            return RPC[1]*Rtuv(t-1, u, v, n+1, p, RPC)
        end
    elseif v == 0
        if u > 1
            return  (u-1)*Rtuv(t, u-2, v, n+1, p, RPC) +
                   RPC[2]*Rtuv(t, u-1, v, n+1, p, RPC)
        else
            return RPC[2]*Rtuv(t, u-1, v, n+1, p ,RPC)
        end
    else
        if v > 1
            return  (v-1)*Rtuv(t, u, v-2, n+1, p, RPC) +
                   RPC[3]*Rtuv(t, u, v-1, n+1, p, RPC)
        else
            return RPC[3]*Rtuv(t, u, v-1, n+1, p, RPC)
        end
    end
end


"""
    get_ijk(l::Int)

Return an array with the Cartesian quantum numbers summing to angular momentum `l`.
"""
function get_ijk(l::Int)
    # number of angular momentum projections
    Nijk = (l+1)*(l+2)÷2
    # preallocate vector of ml values
    ijk = Vector{NTuple{3,Int}}(undef, Nijk)
    # initialize counter
    it = 1
    for a = 1:l+1
        for b = 1:a
            i = l + 1 - a
            j = a - b
            k = b - 1
            ijk[it] = (i, j, k)
            it +=1
        end
    end
    return ijk
end

