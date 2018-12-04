

include("molecule.jl")
include("basis.jl")

using Printf
using LinearAlgebra


function rhf(mol::Molecule, basis::AOBasis)

    Etol = 10.0^(-6)
    maxiter = 15

    N = electrons(mol)
    occ = Int(N/2)
    M = dimension(basis)

    # compute integrals
    S = overlap(basis)
    D = zeros(M,M)
    F = zeros(M,M)
    T = kinetic(basis)
    Vne = attraction(basis, mol)
    Vee = repulsion(basis)
    Vnn = nuclear(mol)

    # symmetric orthogonalization
    s,U = eigen(S;permute=false,scale=false)
    X = U*diagm(0 => 1.0./sqrt.(s))*U'

    # form Hcore and diagonalize
    Hcore = T + Vne
    Ecore,Ccore = eigen(X'*Hcore*X;permute=false,scale=false)

    ### in alternative, I can avoid to compute X and directly do
    # Ecore,Ccore = eig(Hcore,S)

    # transform back Ccore in the old basis
    Ccore = X*Ccore

    # sort energies and orbitals
    # idx = sortperm(Ecore)
    # Ccore = Ccore[:,idx]
    # Ecore = Ecore[idx]

    # compute density matrix
    for μ = 1:M
        for ν = 1:M
            D[μ,ν] = 2.0 * dot(Ccore[μ,1:occ],Ccore[ν,1:occ])
        end
    end

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
                        F[i,j] += D[k,l] * (Vee[i,j,k,l] - 0.5*Vee[i,k,j,l])
                    end
                end
            end
        end

        # solve the Roothaan-Hall equations FC = SCE
        Efock,Cfock = eigen(X'*F*X;permute=false,scale=false)

        # order evals and evecs
        # idx = sortperm(Efock)
        # Cfock = Cfock[:,idx]
        # Efock = Efock[idx]

        # transform back coefs
        C = X*Cfock
        # C = Cfock
        # construct density matrix
        # D = 2.0 * (C[:,1:occ] * C[:,1:occ]')
        # compute density matrix
        for μ = 1:M
            for ν = 1:M
                D[μ,ν] = 2.0 * dot(C[μ,1:occ],C[ν,1:occ])
            end
        end
    

        Enew = 0.0
        for i = 1:M
            for j = 1:M
                Enew += 0.5*D[i,j]*(Hcore[i,j] + F[i,j])
            end
        end

        ΔE = Enew - Eold

        @printf("%12i%12.6f%12.6f\n",iter,Enew+Vnn,ΔE)

        if abs(ΔE) < Etol
            println("SCF converged!")
            break
        end
        Eold = Enew

    end # End of SCF loop

end
