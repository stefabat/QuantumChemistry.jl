
# A shell is a set of atomic orbitals with the same angular momentum and center.
# For instance, a p shell contains three atomic orbitals: px, py and pz.
# The atomic orbitals composing a shell are assumed to be contracted Gaussian
# functions, such that a shell also contains the list of exponents and the
# matrix of contraction coefficients

"""
Concrete type representing an atomic orbital shell.
`M` is the number of primitive functions.
`N` is the number of contracted functions.
"""
struct Shell{M,N}
    center::SVector{3,Float64}   # center of the Shell
    angular_momentum::Int                  # angular momentum
    coefficients::SMatrix{M,N,Float64}     # contraction coefficients matrix
    exponents::SVector{M,Float64}     # Gaussian exponents
end


"""Construct a `Shell` centered on `R` with angular momenutm `l`."""
function Shell(center::AbstractVector{T}, angular_momentum::Int,
            coefficients::AbstractArray{T}, exponents::AbstractVector{T}) where {T<:Number}

    # extract sizes
    M=size(coefficients,1)
    N=size(coefficients,2)

    # input sanity check
    @assert length(center) == 3
    @assert angular_momentum ≥ 0
    @assert length(exponents) == M
    @assert M ≥ N

    # promote to static arrays
    R = SVector{3,Float64}(center)
    l = angular_momentum
    d = SMatrix{M,N,Float64}(coefficients)
    α = SVector{M,Float64}(exponents)

    return Shell(R, l, d, α)
end

get_num_primitives(shell::Shell) = length(shell.exponents)
get_num_contractions(shell::Shell) = size(shell.coefficients,2)

coefficients(shell::Shell) = shell.coefficients
center(shell::Shell) = shell.center
exponents(shell::Shell) = shell.exponents
get_am(shell::Shell) = shell.angular_momentum

Rx(shell::Shell) = shell.center[1]
Ry(shell::Shell) = shell.center[2]
Rz(shell::Shell) = shell.center[3]


