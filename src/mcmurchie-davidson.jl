

"""
    hermite_expansion!(E::AbstractArray, imax::Int, jmax::Int,
                       Kαβ::Real, p::Real, PA::Real, PB::Real)

Compute the Hermite expansion coefficients for a 1-dimensional Cartesian overlap distribution
using a two-term recursion relation.
The array `E` must have dimensions at least `(imax+jmax+1,imax+1,jmax+1)`.
To increase performance there are no bounds check.
"""
function hermite_expansion!(E::AbstractArray, imax::Int, jmax::Int,
                            Kαβ::Real, p::Real, PA::Real, PB::Real)

    # Assuming E is zerod out!!!

    # note that Eₜⁱʲ = E[t+1,i+1,j+1] because in Julia we start
    # to count at 1, e.g. E₀¹² = E[1,2,3]
    for j = 0:jmax
        for i = 0:imax
            @inbounds for t = i+j:-1:0
                if t > 0
                    if i > 0
                        E[t+1,i+1,j+1] += (1/(2*p*t)) * i * E[t, i, j+1]
                    end
                    if j > 0
                        E[t+1,i+1,j+1] += (1/(2*p*t)) * j * E[t, i+1, j]
                    end
                else
                    if i == j == 0
                        E[t+1,i+1,j+1] = Kαβ
                    elseif j == 0
                        E[t+1,i+1,j+1] = PA * E[1, i, j+1] + E[2, i, j+1]
                    else
                        E[t+1,i+1,j+1] = PB * E[1, i+1, j] + E[2, i+1, j]
                    end
                end
            end
        end
    end

    return nothing
end


"""
    hermite_expansion(imax::Int, jmax::Int, Kαβ::Real, p::Real, PA::Real, PB::Real)

Compute the Hermite expansion coefficients for a 1-dimensional Cartesian overlap distribution
using a two-term recursion relation.
"""
function hermite_expansion(imax::Int, jmax::Int, Kαβ::Real, p::Real, PA::Real, PB::Real)

    # initialize array of expansion coefficients
    E = zeros(imax+jmax+1,imax+1,jmax+1)
    # compute expansion coefficients
    hermite_expansion!(E, imax, jmax, Kαβ, p, PA, PB)

    return E
end


"Compute the Hermite expansion coefficients for two primitive shells."
function hermite_expansion!(E::AbstractArray,
                            α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
                            β::Real, lb::Int, Bx::Real, By::Real, Bz::Real)

    # exponent of Gaussian product
    p  =  α + β
    # center of Gaussian product
    Px = (α*Ax + β*Bx)/p
    Py = (α*Ay + β*By)/p
    Pz = (α*Az + β*Bz)/p
    # distances between Gaussian centers
    PAx = Px - Ax; PAy = Py - Ay; PAz = Pz - Az
    PBx = Px - Bx; PBy = Py - By; PBz = Pz - Bz

    # exponential prefactor
    μ  = (α * β)/p
    Kαβx = exp(-μ * (Ax-Bx)^2)
    Kαβy = exp(-μ * (Ay-By)^2)
    Kαβz = exp(-μ * (Az-Bz)^2)

    E .= 0.0
    Ex = @view E[1:la+lb+1 , 1:la+1 , 1:lb+1 , 1]
    Ey = @view E[1:la+lb+1 , 1:la+1 , 1:lb+1 , 2]
    Ez = @view E[1:la+lb+1 , 1:la+1 , 1:lb+1 , 3]
    hermite_expansion!(Ex, la, lb, Kαβx, p, PAx, PBx)
    hermite_expansion!(Ey, la, lb, Kαβy, p, PAy, PBy)
    hermite_expansion!(Ez, la, lb, Kαβz, p, PAz, PBz)

    return nothing
end


"""
    hermite_integral!(R::AbstractArray, tmax::Int, umax::Int, vmax::Int, p::Real, PQx::Real, PQy::Real, PQz::Real)

Compute the integral of an Hermite Gaussian divided by the Coulomb operator.
The array `R` must have dimensions at least `(tmax+umax+vmax+1,tmax+1,umax+1,vmax+1)`.
To increase performance there are no bounds check.
"""
function hermite_integral!(R::AbstractArray, tmax::Int, umax::Int, vmax::Int,
    p::Real, PQx::Real, PQy::Real, PQz::Real)

    # note that Rⁿₜᵤᵥ = R[n+1,t+1,u+1,v+1], e.g. R⁰₂₀₁ = R[1,3,1,2]

    # compute auxiliary integrals Rⁿ₀₀₀
    @inbounds for n = 0:tmax+umax+vmax
        R[n+1, 1, 1, 1] = (-2.0*p)^n * boys(n,p * (PQx^2 + PQy^2 + PQz^2))
    end

    # transfer angular momentum n to t
    for t = 0:tmax-1
        @inbounds for n = 0:tmax+umax+vmax-t-1
            R[n+1, t+2, 1, 1] = PQx * R[n+2, t+1, 1, 1]
            if t > 0
                R[n+1, t+2, 1, 1] += t * R[n+2, t, 1, 1]
            end
        end
    end

    # transfer angular momentum n to u
    for u = 0:umax-1
        for t = 0:tmax
            @inbounds for n = 0:tmax+umax+vmax-t-u-1
                R[n+1, t+1, u+2, 1] = PQy * R[n+2, t+1, u+1, 1]
                if u > 0
                    R[n+1, t+1, u+2, 1] += u * R[n+2, t+1, u, 1]
                end
            end
        end
    end

    # transfer angular momentum n to v
    for v = 0:vmax-1
        for u = 0:umax
            for t = 0:tmax
                @inbounds for n = 0:tmax+umax+vmax-t-u-v-1
                    R[n+1, t+1, u+1, v+2] = PQz * R[n+2, t+1, u+1, v+1]
                    if v > 0
                        R[n+1, t+1, u+1, v+2] += v * R[n+2, t+1, u+1, v]
                    end
                end
            end
        end
    end

    return R
end


"Compute the Hermite integral for two primitive shells."
function hermite_integral!(R::AbstractArray,
                           α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
                           β::Real, lb::Int, Bx::Real, By::Real, Bz::Real,
                                             Cx::Real, Cy::Real, Cz::Real)

    # exponent of Gaussian product
    p  =  α + β
    # center of Gaussian product
    Px = (α*Ax + β*Bx)/p
    Py = (α*Ay + β*By)/p
    Pz = (α*Az + β*Bz)/p
    # distances between Gaussian centers
    # PAx = Px - Ax; PAy = Py - Ay; PAz = Pz - Az
    # PBx = Px - Bx; PBy = Py - By; PBz = Pz - Bz

    PCx = Px - Cx; PCy = Py - Cy; PCz = Pz - Cz
    # # exponential prefactor
    # μ  = (α * β)/p
    # Kαβx = exp(-μ * (Ax-Bx)^2)
    # Kαβy = exp(-μ * (Ay-By)^2)
    # Kαβz = exp(-μ * (Az-Bz)^2)

    tmax = umax = vmax = la+lb

    R .= 0.0
    Rtuv = @view R[1:tmax+umax+vmax+1 , 1:tmax , 1:umax , 1:vmax]
    hermite_integral!(Rtuv, tmax, umax, vmax, p, PCx, PCy, PCz)

    return nothing
end


"""
    hermite_integral(tmax::Int, umax::Int, vmax::Int, p::Real, PQx::Real, PQy::Real, PQz::Real)

Compute the integral of an Hermite Gaussian divided by the Coulomb operator.
"""
function hermite_integral(tmax::Int, umax::Int, vmax::Int, p::Real, PQx::Real, PQy::Real, PQz::Real)

    # initialize array
    R = zeros(tmax+umax+vmax+1,tmax+1,umax+1,vmax+1)

    return hermite_integral!(R, tmax, umax, vmax, p, PQx, PQy, PQz)
end


"""
    overlap!(S::AbstractArray, E::AbstractArray
             α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
             β::Real, lb::Int, Bx::Real, By::Real, Bz::Real)

Compute the overlap integral matrix between two primitive Cartesian shells centered on
`(Ax,Ay,Az)` and `(Bx,By,Bz)`, with exponents `α` and `β` and angular momenta `la` and `lb`.
The array `S` must have dimensions at least `((la+1)*(la+2)÷2,(lb+1)*(lb+2)÷2)`.
To increase performance there are no bounds check.
"""
function overlap!(S::AbstractArray, E::AbstractArray,
                  α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
                  β::Real, lb::Int, Bx::Real, By::Real, Bz::Real)

    # E = zeros(la+lb+1,la+1,lb+1,3)
    hermite_expansion!(E, α, la, Ax, Ay, Az, β, lb, Bx, By, Bz)

    î = get_ijk(la)
    ĵ = get_ijk(lb)

    for (b,(jx,jy,jz)) in enumerate(ĵ)
        @inbounds for (a,(ix,iy,iz)) in enumerate(î)
            S[a,b] = sqrt(π/(α+β))^3  *
                     E[1,ix+1,jx+1,1] *
                     E[1,iy+1,jy+1,2] *
                     E[1,iz+1,jz+1,3]
        end
    end

    return nothing
end

# function overlap!(S::AbstractArray,
#                   α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
#                   β::Real, lb::Int, Bx::Real, By::Real, Bz::Real)

#     E = zeros(la+lb+1,la+1,lb+1,3)
#     S = overlap!(S, E, α, la, Ax, Ay, Az, β, lb, Bx, By, Bz)

#     return nothing
# end

"""
    overlap(α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
            β::Real, lb::Int, Bx::Real, By::Real, Bz::Real)

Compute the overlap integral matrix between two primitive Cartesian shells centered on
`(Ax,Ay,Az)` and `(Bx,By,Bz)`, with exponents `α` and `β` and angular momenta `la` and `lb`.
"""
function overlap(α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
                 β::Real, lb::Int, Bx::Real, By::Real, Bz::Real)

    # number of Cartesian primitive functions in shell
    Nla = (la+1)*(la+2)÷2
    Nlb = (lb+1)*(lb+2)÷2
    S = zeros(Nla,Nlb)
    E = zeros(la+lb+1,la+1,lb+1,3)

    overlap!(S, E, α, la, Ax, Ay, Az, β, lb, Bx, By, Bz)

    return S
end


"""
    kinetic!(T::AbstractArray, E::AbstractArray,
            α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
            β::Real, lb::Int, Bx::Real, By::Real, Bz::Real)

Compute the kinetic energy integral matrix between two primitive
Cartesian shells centered on `(Ax,Ay,Az)` and `(Bx,By,Bz)`, with exponents `α` and
`β` and angular momenta `la` and `lb`.
"""
function kinetic!(T::AbstractArray, E::AbstractArray,
                  α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
                  β::Real, lb::Int, Bx::Real, By::Real, Bz::Real)

    hermite_expansion!(E, α, la+2, Ax, Ay, Az, β, lb, Bx, By, Bz)

    Ex = @view E[1:la+lb+3 , 1:la+3 , 1:lb+1 , 1]
    Ey = @view E[1:la+lb+3 , 1:la+3 , 1:lb+1 , 2]
    Ez = @view E[1:la+lb+3 , 1:la+3 , 1:lb+1 , 3]

    î = get_ijk(la)
    ĵ = get_ijk(lb)

    for (b,(jx,jy,jz)) in enumerate(ĵ)
        @inbounds for (a,(ix,iy,iz)) in enumerate(î)
            Tx = (-2*α^2*Ex[1,ix+3,jx+1]) + (α*(2*ix+1)*Ex[1,ix+1,jx+1])
            if ix > 1
                Tx -= 0.5*ix*(ix-1)*Ex[1,ix-1,jx+1]
            end

            Ty = (-2*α^2*Ey[1,iy+3,jy+1]) + (α*(2*iy+1)*Ey[1,iy+1,jy+1])
            if iy > 1
                Ty -= 0.5*iy*(iy-1)*Ey[1,iy-1,jy+1]
            end

            Tz = (-2*α^2*Ez[1,iz+3,jz+1]) + (α*(2*iz+1)*Ez[1,iz+1,jz+1])
            if iz > 1
                Tz -= 0.5*iz*(iz-1)*Ez[1,iz-1,jz+1]
            end

            T[a,b] = sqrt(π/(α+β))^3 * (
                     Tx * Ey[1,iy+1,jy+1] * Ez[1,iz+1,jz+1] +
                     Ex[1,ix+1,jx+1] * Ty * Ez[1,iz+1,jz+1] +
                     Ex[1,ix+1,jx+1] * Ey[1,iy+1,jy+1] * Tz )
        end
    end

    return nothing
end


"""
    kinetic(α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
            β::Real, lb::Int, Bx::Real, By::Real, Bz::Real)

Compute the kinetic energy integral matrix between two primitive
Cartesian shells centered on `(Ax,Ay,Az)` and `(Bx,By,Bz)`, with exponents `α` and
`β` and angular momenta `la` and `lb`.
"""
function kinetic(α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
                 β::Real, lb::Int, Bx::Real, By::Real, Bz::Real)

    # number of Cartesian primitive functions in shell
    Nla = (la+1)*(la+2)÷2
    Nlb = (lb+1)*(lb+2)÷2
    T = zeros(Nla,Nlb)
    E = zeros(la+lb+3,la+3,lb+1,3)

    kinetic!(T, E, α, la, Ax, Ay, Az, β, lb, Bx, By, Bz)

    return T
end


"""
    attraction(α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
               β::Real, lb::Int, Bx::Real, By::Real, Bz::Real,
                                 Cx::Real, Cy::Real, Cz::Real)

Compute the electron-nuclear attraction energy integral <Ga|1/r|Gb> between two
primitive Cartesian shells centered on `(Ax,Ay,Az)` and `(Bx,By,Bz)`, with
exponents `α` and `β` and angular momenta `la` and `lb`.
"""
function attraction!(V::AbstractArray, E::AbstractArray, R::AbstractArray,
                     α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
                     β::Real, lb::Int, Bx::Real, By::Real, Bz::Real,
                                       Cx::Real, Cy::Real, Cz::Real)

    # precomputing all required quantities
    # (i,k,m) = ikm; (j,l,n) = jln
    # p = α + β
    # RP = (α.*RA + β.*RB) ./ p
    # RAB = RA - RB; RPA = RP - RA; RPB = RP - RB; RPC = RP - RC
    # μ = (α * β)/(α + β)
    # Kαβ = exp.(-μ .* RAB.^2)

    hermite_expansion!(E, α, la, Ax, Ay, Az, β, lb, Bx, By, Bz)
    Ex = @view E[1:la+lb+1 , 1:la+1 , 1:lb+1 , 1]
    Ey = @view E[1:la+lb+1 , 1:la+1 , 1:lb+1 , 2]
    Ez = @view E[1:la+lb+1 , 1:la+1 , 1:lb+1 , 3]

    hermite_integral!(R, α, la, Ax, Ay, Az, β, lb, Bx, By, Bz, Cx, Cy, Cz)
    R0 = @view R[1, 1:la+lb+1, 1:la+lb+1, 1:la+lb+1]

    î = get_ijk(la)
    ĵ = get_ijk(lb)

    vne = 0.0
    for (b,(jx,jy,jz)) in enumerate(ĵ)
        for (a,(ix,iy,iz)) in enumerate(î)
            t1 = Ez[:, iz, jz] .* R0[t, u, :]
            t2 = Ey[:, iy, jy] .* t1

    for t = 0:ix+jx
        for u = 0:iy+jz
            Ekl = Ey[u, k, l]
            for v = 0:iz+jz
                Emn = Ez[v, m, n]
                vne +=  E[t,ix,jx] * E[u,iy,jy] * E[v,iz,jz] * R[1,t, u, v]
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


"""
    get_ijk!(ijk::AbstractArray, l::Int)

Compute an array with the Cartesian quantum numbers summing to angular momentum `l`.
"""
function get_ijk!(ijk::AbstractArray, l::Int)
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

    return nothing
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
    # get the vector
    get_ijk!(ijk, l)

    return ijk
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

## Works only up to f-type orbitals!!
"""Compute cartesian to spherical harmonic transformation matrix for angular momentum l"""
function cart2sph(l::Int)

    # for the p shell, l=1, m=1 corresponds to px, l=1, m=-1 to py and l=1, m=0 to pz
    # d shell: dxy,  dyz, dz^2, dxz, dx^2-y^2
    # (2,m) : m=-2, m=-1,  m=0, m=1, m=2
    # and so on

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

    # return SMatrix{Nsph,Ncart}(T)
    return T
end


# using Memoize

"""
    Etij(t::Int ,i::Int, j::Int, Kαβx::Real, p::Real, XPA::Real, XPB::Real)

Compute the Hermite expansion coefficients for a 1-dimensional Cartesian overlap distribution
using a two-term recursion relation.
"""
function Etij(t::Int ,i::Int, j::Int, Kαβx::Real, p::Real, XPA::Real, XPB::Real)
# @memoize function Etij(t::Int ,i::Int, j::Int, Kαβx::Real, XPA::Real, XPB::Real, α::Real, β::Real)

    # enter recursion
    if t < 0 || t > i+j
        return 0.0
    elseif t == 0
        if i == j == 0
            return Kαβx
        elseif j == 0
            return XPA * Etij(0, i-1, j, Kαβx, p, XPA, XPB) +
                         Etij(1, i-1, j, Kαβx, p, XPA, XPB)
        else
            return XPB * Etij(0, i, j-1, Kαβx, p, XPA, XPB) +
                         Etij(1, i, j-1, Kαβx, p, XPA, XPB)
        end
    else
        return (1/(2*p*t)) * (i * Etij(t-1, i-1, j, Kαβx, p, XPA, XPB) +
                              j * Etij(t-1, i, j-1, Kαβx, p, XPA, XPB) )
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
Inefficient overlap
"""
function overlap2(α::Real, la::Int, Ax::Real, Ay::Real, Az::Real,
                  β::Real, lb::Int, Bx::Real, By::Real, Bz::Real)

    # precomputing all required quantities
    p  =  α + β
    μ  = (α * β)/(α + β)

    Px = (α*Ax + β*Bx)/p
    Py = (α*Ay + β*By)/p
    Pz = (α*Az + β*Bz)/p

    PAx = Px - Ax; PAy = Py - Ay; PAz = Pz - Az
    PBx = Px - Bx; PBy = Py - By; PBz = Pz - Bz

    Kαβx = exp(-μ * (Ax-Bx)^2)
    Kαβy = exp(-μ * (Ay-By)^2)
    Kαβz = exp(-μ * (Az-Bz)^2)

    # number of Cartesian primitive functions in shell
    Nla = (la+1)*(la+2)÷2
    Nlb = (lb+1)*(lb+2)÷2
    S = zeros(Nla,Nlb)

    for (b,(jx,jy,jz)) in enumerate(get_ijk(lb))
        for (a,(ix,iy,iz)) in enumerate(get_ijk(la))
            Ex = Etij(0, ix, jx, Kαβx, p, PAx, PBx)
            Ey = Etij(0, iy, jy, Kαβy, p, PAy, PBy)
            Ez = Etij(0, iz, jz, Kαβz, p, PAz, PBz)
            @inbounds S[a,b] = Ex * Ey * Ez
        end
    end

    return ((π / p)^1.5) .* S
end

