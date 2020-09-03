

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





# end # of Module
