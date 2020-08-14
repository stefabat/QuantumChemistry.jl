
# module McMurchieDavidson





"""Compute the overlap integral <Ga|Gb> of two Cartesian PGFs."""
function overlap(α::Real, ikm::NTuple{3,Int}, RA::NTuple{3,Float64},
                 β::Real, jln::NTuple{3,Int}, RB::NTuple{3,Float64})

    # precomputing all required quantities
    μ = (α * β)/(α + β)
    RP = (RA.*α .+ RB.*β)./(α + β)
    RAB = RA .- RB; RPA = RP .- RA; RPB = RP .- RB
    Kαβ = exp.(-μ.*RAB.^2)

    # calculate overlaps in the 3 Cartesian directions
    @inbounds Sx = Etij(0, ikm[1], jln[1], Kαβ[1], RPA[1], RPB[1], α, β)
    @inbounds Sy = Etij(0, ikm[2], jln[2], Kαβ[2], RPA[2], RPB[2], α, β)
    @inbounds Sz = Etij(0, ikm[3], jln[3], Kαβ[3], RPA[3], RPB[3], α, β)

    return Sx * Sy * Sz * (π / (α + β) )^1.5
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
        Eij = Etij(t, i, j, Kαβ[1], RPA[1], RPB[1], α, β)
        for u = 0:k+l
            Ekl = Etij(u, k, l, Kαβ[2], RPA[2], RPB[2], α, β)
            for v = 0:m+n
                Emn = Etij(v, m, n, Kαβ[3], RPA[3], RPB[3], α, β)
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
        Eij1 = Etij(t, i1, j1, Kαβ[1], RPA[1], RPB[1], α, β)
        for u = 0:k1+l1
            Ekl1 = Etij(u, k1, l1, Kαβ[2], RPA[2], RPB[2], α, β)
            for v = 0:m1+n1
                Emn1 = Etij(v, m1, n1, Kαβ[3], RPA[3], RPB[3], α, β)
                for τ = 0:i2+j2
                    Eij2 = Etij(τ, i2, j2, Kγδ[1], RQC[1], RQD[1], γ, δ)
                    for ν = 0:k2+l2
                        Ekl2 = Etij(ν, k2, l2, Kγδ[2], RQC[2], RQD[2], γ, δ)
                        for ϕ = 0:m2+n2
                            Emn2 = Etij(ϕ, m2, n2, Kγδ[3], RQC[3], RQD[3], γ, δ)
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
Compute the Hermite expansion coefficients for a 1D Cartesian overlap distribution
using a two-term recursion relation.
"""
# @memoize function Etij(t::Int ,i::Int, j::Int, Kαβx::Real, XPA::Real, XPB::Real, α::Real, β::Real)
function Etij(t::Int ,i::Int, j::Int, Kαβx::Real, XPA::Real, XPB::Real, α::Real, β::Real)

    # compute overlap exponent
    p = α + β

    # enter recursion
    if t < 0 || t > i+j
        return 0.0
    elseif t == 0
        if i == j == 0
            return Kαβx
        elseif j == 0
            return XPA * Etij(0, i-1, j, Kαβx, XPA, XPB, α, β) +
                         Etij(1, i-1, j, Kαβx, XPA, XPB, α, β)
        else
            return XPB * Etij(0, i, j-1, Kαβx, XPA, XPB, α, β) +
                         Etij(1, i, j-1, Kαβx, XPA, XPB, α, β)
        end
    else
        return (1/(2*p*t)) * (i * Etij(t-1, i-1, j, Kαβx, XPA, XPB, α, β) +
                              j * Etij(t-1, i, j-1, Kαβx, XPA, XPB, α, β) )
    end
end


# function Etij(i⃗::NTuple{3,Int}, j⃗::NTuple{3,Int}, Kαβ::NTuple{3,Real},
#               RPA::NTuple{3,Real}, RPB::NTuple{3,Real}, α::Real, β::Real)

    
#     for tx = 0:i⃗[1]+j⃗[1]
#         Ex = Etij(tx, i⃗[1], j⃗[1], Kαβ[1], RPA[1], RPB[1], α, β)
#         for ty = 0:i⃗[2]+j⃗[2]
#             Ey = Etij(tx, i⃗[2], j⃗[2], Kαβ[2], RPA[2], RPB[2], α, β)
#             for tz = 0:i⃗[3]+j⃗[3]
#                 Ez = Etij(tz, i⃗[3], j⃗[3], Kαβ[3], RPA[3], RPB[3], α, β)
#                 E += Ex * Ey * Ez
#             end
#         end
#     end

# end

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



# end # of Module
