

# module BasisFunctions

using LinearAlgebra: norm
using StaticArrays
import Base: exponent

# type export
# export PGF, CGF, Shell

# method export
# export center, exponent, shell
# export Rx, Ry, Rz, lx, ly, lz, ltot
# export primitives, nprimitives, exponents, coefs, norms


"""Abstract type representing a basis function."""
abstract type BasisFunction end


#----------------------------------- PGF --------------------------------------#

"""Concrete type defining a primitive Gaussian function."""
struct PGF <: BasisFunction
    R::SVector{3,Float64}    # center of the PGF
    α::Float64              # exponent
    l::SVector{3,Int}        # Cartesian quantum numbers
end

center(pgf::PGF)   = pgf.R
exponent(pgf::PGF) = pgf.α
shell(pgf::PGF)    = pgf.l

Rx(pgf::PGF) = pgf.R[1]
Ry(pgf::PGF) = pgf.R[2]
Rz(pgf::PGF) = pgf.R[3]

lx(pgf::PGF) = pgf.l[1]
ly(pgf::PGF) = pgf.l[2]
lz(pgf::PGF) = pgf.l[3]
ltot(pgf::PGF) = sum(pgf.l)
#------------------------------------------------------------------------------#


#----------------------------------- CGF --------------------------------------#

"""Concrete type defining a contracted Gaussian function."""
struct CGF
    funcs::Vector{PGF}          # primitive PGFs
    coefs::Vector{Float64}      # contraction coefs
    norms::Vector{Float64}      # normalization factors
end


"""Construct a `CGF` centered on `R` with angular momenutm `l`."""
function CGF(R::SVector{3,Float64}, l::SVector{3,Int},
             α::AbstractVector{T}, d::AbstractVector{T}) where {T<:Real}

    # number of coefs and exponents must match
    @assert(length(d) == length(α))
    n = length(d)

    # initialize data vector
    funcs = Vector{PGF}(undef,n)
    norms = Vector{Float64}(undef,n)

    for i in 1:n
        funcs[i] = PGF(R,α[i],l)
        norms[i] = normalization(funcs[i])
    end

    temp_cgto = CGF(funcs, d, norms)
    norms .*= 1.0/sqrt(overlap(temp_cgto))

    return CGF(funcs, d, norms)
end

primitives(cgf::CGF) = cgf.funcs
nprimitives(cgf::CGF) = length(cgf.funcs)
coefs(cgf::CGF) = cgf.coefs
norms(cgf::CGF) = cgf.norms

center(cgf::CGF) = cgf.funcs[1].R
exponents(cgf::CGF) = map(x->exponent(x),primitives(cgf))
shell(cgf::CGF) = shell(cgf.funcs[1])

Rx(cgf::CGF) = Rx(cgf.funcs[1])
Ry(cgf::CGF) = Ry(cgf.funcs[1])
Rz(cgf::CGF) = Rz(cgf.funcs[1])

lx(cgf::CGF) = lx(cgf.funcs[1])
ly(cgf::CGF) = ly(cgf.funcs[1])
lz(cgf::CGF) = lz(cgf.funcs[1])
ltot(cgf::CGF) = ltot(cgf.funcs[1])
#------------------------------------------------------------------------------#



################ copied over from integrals ##############


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




# end # of Module
