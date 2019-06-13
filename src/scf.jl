

include("molecule.jl")
include("basis.jl")

using Printf
using LinearAlgebra


function rhf(mol::Molecule, basis::AOBasis, maxiter::Int = 3, Etol::Float64 = 1e-6)

    N = electrons(mol)
    occ = Int(N/2)
    M = dimension(basis)

    # initialize matries
    D = zeros(M,M)
    F = zeros(M,M)

    # compute integrals
    S = overlap(basis)
    T = kinetic(basis)
    Vne = attraction(basis, mol)
    Vee = repulsion_debug(basis)
    # Vee = repulsion(basis)
    Vnn = nuclear(mol)

    # symmetric orthogonalization
    s,U = eigen(Symmetric(S))
    X = U*diagm(0 => 1.0./sqrt.(s))*U'

    Eold = 0.0

    @printf("%12s%12s%12s\n","Iteration","Energy","ΔE")
    # enter SCF loop
    for iter = 1:maxiter

        # Fock matrix construction
        for i = 1:M
            for j = 1:M
                F[i,j] = T[i,j] + Vne[i,j]
                for k = 1:M
                    for l = 1:M
                        # F[i,j] += D[k,l] * (Vee[index(i,j),index(k,l)] -
                        #               0.5 * Vee[index(i,k),index(j,l)])
                        F[i,j] += D[k,l] * (Vee[i,j,k,l] - 0.5 * Vee[i,k,j,l])
                    end
                end
                # F[j,i] = F[i,j]
            end
        end

        # G = F - T - Vne
        # println("G = ",G)
        # println("F = ",F)
        # println("F' = ",X'*F*X)

        ## WHEN I TRANSFORM MANUALLY IT DOESN'T WORK
        ## I need to make the matrix Symmetric first
        ## This sort the eigenvalues!
        ## The problem were the eigenvalues not sorted
        # solve the Roothaan-Hall equations FC = SCE
        Efock,Cfock = eigen(Symmetric(X'*F*X))#;permute=true,scale=false)
        # println("C' = ",Cfock)

        # transform back coefs
        C = X*Cfock
        E = Efock

        # THIS WORKS PERFECTLY
        # E,C = eigen(F,S)
        # println("ε = ",E)
        # println("C = ",C)
        
        # compute density matrix
        for μ = 1:M
            for ν = 1:M
                D[μ,ν] = 2.0 * dot(C[μ,1:occ],C[ν,1:occ])
            end
        end
        # println("D = ",D)


        # compute new estimate of electronic energy
        Enew = 0.0
        for i = 1:M
            for j = 1:M
                Enew += 0.5*D[i,j]*(T[i,j] + Vne[i,j] + F[i,j])
            end
        end

        ΔE = Enew - Eold

        @printf("%12i%12.6f%12.6f\n",iter,Enew+Vnn,ΔE)

        if abs(ΔE) < Etol
            println("SCF converged!")
            return Enew+Vnn
        end
        Eold = Enew

    end # End of SCF loop


end
