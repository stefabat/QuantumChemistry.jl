

# module Integrals

include("mcmurchie-davidson.jl")
# using .McMurchieDavidson

include("basisfunctions.jl")
# using .BasisFunctions


angular = Dict(
0=>[(0,0,0)],
1=>[(1,0,0),(0,1,0),(0,0,1)],
2=>[(2,0,0),(1,1,0),(1,0,1),
(0,2,0),(0,1,1),(0,0,2)],
3=>[(3,0,0),(2,1,0),(2,0,1),
(1,2,0),(1,1,1),(1,0,2),
(0,3,0),(0,2,1),(0,1,2),
(0,0,3)],
4=>[(4,0,0),(3,1,0),(3,0,1),
(2,2,0),(2,1,1),(2,0,2),
(1,3,0),(1,2,1),(1,1,2),
(1,0,3),(0,4,0),(0,3,1),
(0,2,2),(0,1,3),(0,0,4)],
5=>[(5,0,0),(4,1,0),(4,0,1),
(3,2,0),(3,1,1),(3,0,2),
(2,3,0),(2,2,1),(2,1,2),
(2,0,3),(1,4,0),(1,3,1),
(1,2,2),(1,1,3),(1,0,4),
(0,5,0),(0,4,1),(0,3,2),
(0,2,3),(0,1,4),(0,0,5)]
)


function c2s_index(i,j,k)
    index = findall(x->x==(i,j,k),angular[i+j+k])
    if length(index) == 1
        return index[1]
    else
        error("Something went wrong")
    end
end


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
        # vvmin and vvmax are beginning and end of the
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
                    index = c2s_index(i,j,k)
                    
                    # compute coeff and save it
                    T[m+l+1,index] = NlmS(l,m)*Ctuvlm(t,u,v/2,l,m)
                end
            end
        end
    end
    return T
end


#----------------------------------- PGF --------------------------------------#

"""Returns the self-overlap integral <Ga|Ga> of a `PGF` Ga."""
function overlap(Ga::PGF)
    α = exponent(Ga)
    num = dfactorial(2*lx(Ga)-1)*dfactorial(2*ly(Ga)-1)*dfactorial(2*lz(Ga)-1)
    den = (4.0*α)^ltot(Ga)
    return ((π/(2.0*α))^1.5 * num)/den
end


"""Returns the normalization factor of a `PGF` Ga."""
function normalization(Ga::PGF)
    α = exponent(Ga)
    num = (2.0*α/π)^0.75 * (4.0*α)^(ltot(Ga)/2.0)
    den = sqrt(dfactorial(2*lx(Ga)-1)*dfactorial(2*ly(Ga)-1)*dfactorial(2*lz(Ga)-1))
    return num/den
end


"""Returns the overlap integral <Ga|Gb> of two Cartesian GTOs."""
function overlap(Ga::PGF, Gb::PGF)
    
    # extract info from PGFs
    α = exponent(Ga); β = exponent(Gb)
    RA = center(Ga); RB = center(Gb)
    ikm = shell(Ga); jln = shell(Gb)
    
    return overlap(α, ikm, RA, β, jln, RB)
end


"""Returns the kinetic energy integral -0.5*<Ga|∇^2|Gb> of two GTOs."""
function kinetic(Ga::PGF, Gb::PGF)
    
    # extract info from PGFs
    α = exponent(Ga); β = exponent(Gb)
    RA = center(Ga); RB = center(Gb)
    ikm = shell(Ga); jln = shell(Gb)
    
    return kinetic(α, ikm, RA, β, jln, RB)
end


"""Returns the nuclear attraction energy integral of the distribution Ωab from center C."""
function attraction(Ga::PGF, Gb::PGF, RC::NTuple{3,Float64})
    
    # precomputing all required quantities
    α = exponent(Ga); β = exponent(Gb)
    ikm = shell(Ga); jln = shell(Gb)
    RA = center(Ga); RB = center(Gb)
    
    return attraction(α, ikm, RA, β, jln, RB, RC)
end


"""Returns the two-electron integral gabcd = < Ga(r1) Gb(r1) | 1/r12 | Gc(r2) Gd(r2) >."""
function repulsion(Ga::PGF, Gb::PGF, Gc::PGF, Gd::PGF)
    
    # extract info from gaussians for electron 1
    α = exponent(Ga); β = exponent(Gb)
    ikm1 = shell(Ga); jln1 = shell(Gb)
    RA = center(Ga); RB = center(Gb)
    
    # extract info from gaussians for electron 2
    γ = exponent(Gc); δ = exponent(Gd)
    ikm2 = shell(Gc); jln2 = shell(Gd)
    RC = center(Gc); RD = center(Gd)
    
    return repulsion(α, ikm1, RA, β, jln1, RB, γ, ikm2, RC, δ, jln2, RD)
end


#----------------------------------- CGF -------------------------------------#

"""Returns the self-overlap integral of a contracted CGF."""
function overlap(μ::CGF)
    return overlap(μ,μ)
end


"""Returns the overlap integral between two contracted GTOs."""
function overlap(μ::CGF, ν::CGF)
    Gμ = primitives(μ); dμ = coefs(μ); Nμ = norms(μ)
    Gν = primitives(ν); dν = coefs(ν); Nν = norms(ν)
    S = 0.0
    for a in 1:nprimitives(μ)
        for b in 1:nprimitives(ν)
            S += Nμ[a] * Nν[b] * dμ[a] * dν[b] * overlap(Gμ[a],Gν[b])
        end
    end
    return S
end


"""Returns the kinetic energy integral between two contracted GTOs."""
function kinetic(μ::CGF, ν::CGF)
    Gμ = primitives(μ); dμ = coefs(μ); Nμ = norms(μ)
    Gν = primitives(ν); dν = coefs(ν); Nν = norms(ν)
    T = 0.0
    for a in 1:nprimitives(μ)
        for b in 1:nprimitives(ν)
            T += Nμ[a] * Nν[b] * dμ[a] * dν[b] * kinetic(Gμ[a],Gν[b])
        end
    end
    return T
end


"""
Returns the nuclear attraction integral between two contracted GTOs and
the nucleus centered at `C`.
"""
function attraction(μ::CGF, ν::CGF, C::NTuple{3,Float64})
    Gμ = primitives(μ); dμ = coefs(μ); Nμ = norms(μ)
    Gν = primitives(ν); dν = coefs(ν); Nν = norms(ν)
    V = 0.0
    for a in 1:nprimitives(μ)
        for b in 1:nprimitives(ν)
            V += Nμ[a] * Nν[b] * dμ[a] * dν[b] * attraction(Gμ[a],Gν[b],C)
        end
    end
    return V
end


"""Returns the two-electron repulsion integral over four contracted GTOs."""
function repulsion(μ::CGF, ν::CGF, λ::CGF, σ::CGF)
    Gμ = primitives(μ); dμ = coefs(μ); Nμ = norms(μ)
    Gν = primitives(ν); dν = coefs(ν); Nν = norms(ν)
    Gλ = primitives(λ); dλ = coefs(λ); Nλ = norms(λ)
    Gσ = primitives(σ); dσ = coefs(σ); Nσ = norms(σ)
    V = 0.0
    for a in 1:nprimitives(μ)
        for b in 1:nprimitives(ν)
            for c in 1:nprimitives(λ)
                for d in 1:nprimitives(σ)
                    V += Nμ[a] * Nν[b] * Nλ[c] * Nσ[d] *
                    dμ[a] * dν[b] * dλ[c] * dσ[d] *
                    repulsion(Gμ[a],Gν[b],Gλ[c],Gσ[d])
                end
            end
        end
    end
    return V
end


#--------------------------------- Shell -----------------------------------#

"""Returns the overlap integral matrix S of two GTO shells."""
function overlap(Ga::Shell,Gb::Shell)
    
    # precomputing all required quantities
    α = exponent(Ga); β = exponent(Gb)
    RA = center(Ga); RB = center(Gb)
    μ = (α * β)/(α + β)
    RP = (RA.*α .+ RB.*β)./(α + β)
    RAB = RA .- RB; RPA = RP .- RA; RPB = RP .- RB
    Kαβ = exp.(-μ.*RAB.^2)
    
    la = ltot(Ga); lb = ltot(Gb)
    Sx = zeros(la+1,lb+1)
    Sy = zeros(la+1,lb+1)
    Sz = zeros(la+1,lb+1)
    Nal = (la+1)*(la+2)÷2
    Nbl = (lb+1)*(lb+2)÷2
    S = zeros(Nal,Nbl)
    
    # loop over Cartesian quantum numbers
    for a = 0:la
        for b = 0:lb
            @inbounds Sx[a+1,b+1] = Etij(0, a, b, Kαβ[1], RPA[1], RPB[1], α, β)
            @inbounds Sy[a+1,b+1] = Etij(0, a, b, Kαβ[2], RPA[2], RPB[2], α, β)
            @inbounds Sz[a+1,b+1] = Etij(0, a, b, Kαβ[3], RPA[3], RPB[3], α, β)
        end
    end
    
    for (a,(i,k,m)) in enumerate(angular[la])
        for (b,(j,l,n)) in enumerate(angular[lb])
            S[a,b] = Sx[i+1,j+1] * Sy[k+1,l+1] * Sz[m+1,n+1]
        end
    end

    # S = ((π / (α + β) )^1.5) .* S

    return ((π / (α + β) )^1.5) .* S
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



# end # of Module
