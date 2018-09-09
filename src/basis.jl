

include("integrals.jl")

#######

using JSON.Parser.parsefile
sto3g = parsefile("basis/sto-3g.json");

# STO-3G for Hydrogen
α1s = convert(Vector{Float64},sto3g["H"][1]["prim"])
d1s = convert(Vector{Float64},sto3g["H"][1]["cont"][1])
CA  = [0.0, 0.0, 0.0]; ZA = 1
GA1s = CGTO(CA,(0,0,0),α1s,d1s)

# STO-3G for Helium
# α1s = convert(Vector{Float64},sto3g["He"][1]["prim"])
# d1s = convert(Vector{Float64},sto3g["He"][1]["cont"][1])
CB  = [0.0, 0.0, 1.4632]; ZB = 2
GB1s = CGTO(CB,(0,0,0),[2.22766,0.405771,0.109818]*(2.0925^2),d1s) # to reproduce Szabo

######

function atom_basis(α0::Real, γ::Real, ns::Int, np::Int, nd::Int)

    # atomic center
    A = zeros(3)

    nα = ns + 3*np + 6*nd
    gto = Vector{Integrals.GTO}(nα)

    # create the even-tempered series of s-type Gaussians
    for k = 1:ns
        gto[k]    = Integrals.GTO(A, α0*γ^k, (0,0,0))
    end

    # create the even-tempered series of p-type Gaussians
    for k = 0:np-1
        gto[ns+3*k+1] = Integrals.GTO(A, α0*γ^k, (1,0,0))
        gto[ns+3*k+2] = Integrals.GTO(A, α0*γ^k, (0,1,0))
        gto[ns+3*k+3] = Integrals.GTO(A, α0*γ^k, (0,0,1))
    end

    # create the even-tempered series of d-type Gaussians
    for k = 0:nd-1
        gto[ns+3*np+6*k+1] = Integrals.GTO(A, α0*γ^k, (2,0,0))
        gto[ns+3*np+6*k+2] = Integrals.GTO(A, α0*γ^k, (0,2,0))
        gto[ns+3*np+6*k+3] = Integrals.GTO(A, α0*γ^k, (0,0,2))
        gto[ns+3*np+6*k+4] = Integrals.GTO(A, α0*γ^k, (1,1,0))
        gto[ns+3*np+6*k+5] = Integrals.GTO(A, α0*γ^k, (1,0,1))
        gto[ns+3*np+6*k+6] = Integrals.GTO(A, α0*γ^k, (0,1,1))
    end

    return gto
end


function compute_ERI(basis, Tint::Real = 1e-16)

    N = length(basis)
    eri = zeros(N*N,N*N)
    discarded = 0

    for (i,α) in enumerate(basis)
        for (j,β) in enumerate(basis)
            for (k,γ) in enumerate(basis)
                for (l,δ) in enumerate(basis)
                    Vijkl = Integrals.repulsion(α, β, γ, δ)
                    if abs(Vijkl) > Tint
                        p = (i-1)*N + j
                        q = (k-1)*N + l
                        @inbounds eri[p,q] = Vijkl
                    else
                        discarded += 1
                    end
                end
            end
        end
        println(100*i/N," %")
    end

    println("Number of integrals: ",N^4)
    println("Discarded integrals: ",discarded)
    println(" Ratio of discarded: ",round(100*discarded/(N^4),2)," %")

    return eri

end
