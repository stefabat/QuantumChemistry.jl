
using StaticArrays
# module Integrals


#--------------------------------- Shell -----------------------------------#

"""
    overlap(Sa::Shell, Sb::Shell)

Compute the overlap integral between two contracted atomic orbitals shells.
"""
function overlap(Sa::Shell, Sb::Shell)

    # precomputing all required quantities
    # α = shell_a.exponents
    # β = shell_b.exponents
    la = Sa.angular_momentum
    lb = Sb.angular_momentum
    da = Sa.coefficients
    db = Sb.coefficients
    Ra = Sa.center
    Rb = Sb.center

    Nla = (la+1)*(la+2)÷2
    Nlb = (lb+1)*(lb+2)÷2
    S = zeros(Nla, Nlb)

    # loop over primitives
    for (b,β) in enumerate(Sb.exponents)
        for (a,α) in enumerate(Sa.exponents)
            S += da[a] * db[b] * overlap(α, la, Ra, β, lb, Rb)
        end
    end

    return S

    # Ta = cart2sph(la)
    # Tb = cart2sph(lb)
    # # transform to spherical harmonics
    # return Ta*S*Tb'

end



"""Returns the kinetic energy integral matrix T of two GTO shells."""
function kinetic(Ga::Shell, Gb::Shell)

    # extract info from GTOs
    α = exponent(Ga); β = exponent(Gb)
    RA = center(Ga); RB = center(Gb)
    la = ltot(Ga); lb = ltot(Gb)
    Nal = (la+1)*(la+2)÷2
    Nbl = (lb+1)*(lb+2)÷2
    T = zeros(Nal,Nbl)

    fα2 = 4*α^2

    for (a,(i,k,m)) in enumerate(angular[la])
        for (b,(j,l,n)) in enumerate(angular[lb])
            # precompute common terms
            dαSab = 2*α*overlap(α, (i,k,m), RA, β, (j,l,n), RB)

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

            T[a,b] = Dij + Dkl + Dmn
        end
    end

    return -0.5 .* T
end


"""Returns the nuclear attraction energy integral of the distribution Ωab from center C."""
function attraction(Ga::Shell, Gb::Shell, RC::NTuple{3,Float64})

    # precompute all required quantities
    α = exponent(Ga); β = exponent(Gb)
    RA = center(Ga); RB = center(Gb)
    p = α + β
    RP = (α.*RA .+ β.*RB) ./ p
    RAB = RA .- RB; RPA = RP .- RA; RPB = RP .- RB; RPC = RP .- RC
    μ = (α * β)/p
    Kαβ = exp.(-μ .* RAB.^2)

    la = ltot(Ga); lb = ltot(Gb)
    Nal = (la+1)*(la+2)÷2
    Nbl = (lb+1)*(lb+2)÷2
    V = zeros(Nal,Nbl)

    for (a,(i,k,m)) in enumerate(angular[la])
        for (b,(j,l,n)) in enumerate(angular[lb])
            for t = 0:i+j
                Eij = Etij(t, i, j, Kαβ[1], RPA[1], RPB[1], α, β)
                for u = 0:k+l
                    Ekl = Etij(u, k, l, Kαβ[2], RPA[2], RPB[2], α, β)
                    for v = 0:m+n
                        Emn = Etij(v, m, n, Kαβ[3], RPA[3], RPB[3], α, β)
                        V[a,b] += Eij * Ekl * Emn * Rtuv(t, u, v, 0, p, RPC)
                    end
                end
            end
        end
    end

    return 2.0*π/p .* V
end


"""Returns the two-electron integral gabcd = < Ga(r1) Gb(r1) | 1/r12 | Gc(r2) Gd(r2) >."""
function repulsion(Ga::Shell, Gb::Shell, Gc::Shell, Gd::Shell)

    # extract info from gaussians for electron 1
    α = exponent(Ga); β = exponent(Gb)
    # (i1,k1,m1) = ltot(Ga); (j1,l1,n1) = ltot(Gb)
    la = ltot(Ga); lb = ltot(Gb)
    Nal = (la+1)*(la+2)÷2
    Nbl = (lb+1)*(lb+2)÷2
    RA = center(Ga); RB = center(Gb); RAB = RA .- RB
    # precompute quantities for electron 1
    p = α + β
    RP = (α.*RA .+ β.*RB) ./ p
    Kαβ = exp.(-α*β / p .* RAB.^2)
    RPA = RP .- RA; RPB = RP .- RB

    # extract info from gaussians for electron 2
    γ = exponent(Gc); δ = exponent(Gd)
    # (i2,k2,m2) = ltot(Gc); (j2,l2,n2) = ltot(Gd)
    lc = ltot(Gc); ld = ltot(Gd)
    Ncl = (lc+1)*(lc+2)÷2
    Ndl = (ld+1)*(ld+2)÷2
    RC = center(Gc); RD = center(Gd); RCD = RC .- RD
    # precompute quantities for electron 2
    q = γ + δ
    RQ = (γ.*RC .+ δ.*RD) ./ q
    Kγδ = exp.(-γ*δ / q .* RCD.^2)
    RQC = RQ .- RC; RQD = RQ .- RD

    # precompute quantities for auxiliary integral R
    RPQ = RP .- RQ
    ξ = p * q / (p + q)

    Vee = zeros(Nal,Nbl,Ncl,Ndl)

    for (a,(i1,k1,m1)) in enumerate(angular[la])
        for (b,(j1,l1,n1)) in enumerate(angular[lb])
            for (c,(i2,k2,m2)) in enumerate(angular[la])
                for (d,(j2,l2,n2)) in enumerate(angular[lb])
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
                                            Vee[a,b,c,d] += Eij1 * Ekl1 * Emn1 * Eij2 * Ekl2 * Emn2 *
                                            Rtuv(t+τ, u+ν, v+ϕ, 0, ξ, RPQ) * (-1)^(τ+ν+ϕ)
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return (2.0*π^2.5)/(p*q*sqrt(p+q)) .* Vee
end



##### utils

# Cartesian to spherical index
# function cart2sph_index(i::Int, j::Int, k::Int)
#     #index = findall(x->x==(i,j,k),get_ijk(i+j+k))

#     index = indexin([(i,j,k)],get_ijk(i+j+k))

#     if length(index) == 1
#         return index[1]
#     else
#         error("Something went wrong")
#     end
# end


function NlmS(l::Int, m::Int)
    # factor in front of the sqrt
    fac = 1.0/(2^abs(m)*factorial(l))

    # numerator divided by half
    numhalf = factorial(l+abs(m))*factorial(l-abs(m))

    # check m for the denominator, multiply and return
    if m == 0
        return fac*sqrt(numhalf)
    else
        return fac*sqrt(2.0*numhalf)
    end
end

# note that `v` can take half integer values, thus is `Real`
function Ctuvlm(t::Int, u::Int, v::Real, l::Int, m::Int)
    # if statement to set v_m directly in the expression
    if m >= 0
        return (-1)^(t+v) * (0.25)^t * binomial(l,t) *
        binomial(l-t,abs(m)+t) * binomial(t,u) *
        binomial(abs(m),Int(2*v))
    else
        # note the -0.5 due to v_m
        return (-1)^(t+v-0.5) * (0.25)^t * binomial(l,t) *
        binomial(l-t,abs(m)+t) * binomial(t,u) *
        binomial(abs(m),Int(2*v))
    end
end

"""Compute cartesian to spherical harmonic transformation matrix for angular momentum l"""
function cart2sph(l::Int)

    # for the p shell, l=1, m=1 corresponds to px, l=1, m=-1 to py and l=1, m=0 to pz
    # d shell: dxy,  dyz, dz^2, dxz, dx^2-y^2
    # (2,m) : m=-2, m=-1,  m=0, m=1, m=2

    # Initialize matrix
    Nsph = 2*l + 1
    Ncart = (l+1)*(l+2)÷2
    T = zeros(Nsph,Ncart)

    # loop over magnetic quantum numbers
    for m = -l:1:l
        # compute v_m
        vm = 0.0
        if m < 0
            vm = 0.5
        end
        # since the loop over v runs with half integers
        # for m < 0, we run it instead over 2*v, where
        # vmin and vmax are beginning and end of the
        # loop over 2*v, such that this works for all
        # cases, i.e. also when v is an integer
        vmin = Int(2 * vm)
        vmax = Int(2 * (floor(Int,abs(m)/2-vm) + vm))

        # loop over auxiliary indices
        for t = 0:floor(Int,(l-abs(m))/2)
            for u = 0:t
                for v = vmin:2:vmax
                    # compute the Cartesian angular momenta
                    # to know the associated GTO to the
                    # transformation coefficient calculated here
                    i = Int(2*t+abs(m)-2*(u+v/2))
                    j = Int(2*(u+v/2))
                    k = Int(l-2*t-abs(m))

                    # compute index in matrix
                    # index = c2s_index(i,j,k)
                    index = indexin([(i,j,k)],get_ijk(i+j+k))[1]

                    # compute coeff and save it
                    T[m+l+1,index] = NlmS(l,m)*Ctuvlm(t,u,v/2,l,m)
                end
            end
        end
    end

    return SMatrix{Nsph,Ncart}(T)
end

# end # of Module
