

function rhf(molecule::Molecule)

Etol = 10.0^(-6)
maxiter = 15

N = electrons(molecule)
occ = Int(N/2)
M = dimension(basis(molecule))
AOs = aos(basis(molecule))

# Initialize matrices
S = zeros(M,M)
D = zeros(M,M)
F = zeros(M,M)
T = zeros(M,M)
Vne = zeros(M,M)
Vee = zeros(M,M,M,M)
Vnn = 0.0

# Compute AO integrals
for (i,χi) in enumerate(AOs)
    for (j,χj) in enumerate(AOs)
        T[i,j] = kinetic(χi,χj)
        S[i,j] = overlap(χi,χj)
        for (k,χk) in enumerate(AOs)
            for (l,χl) in enumerate(AOs)
                Vee[i,j,k,l] = repulsion(χi,χj,χk,χl)
            end
        end
    end
end

# nulcear-nuclear repulsion
for (I,Inuc) in enumerate(atoms(molecule))
    for (J,Jnuc) in enumerate(atoms(molecule))
        if I<J
            Vnn += Z(Inuc)*Z(Jnuc)/norm(xyz(Inuc)-xyz(Jnuc))
        end
    end
end

# nuclear-electron attraction
for (i,χi) in enumerate(AOs)
    for (j,χj) in enumerate(AOs)
        for I in atoms(molecule)
            Vne[i,j] += nuclear(χi,χj,xyz(I))*(-Z(I))
        end
    end
end

# symmetric orthogonalization
s,U = eig(S)
X = U*diagm(1.0./sqrt.(s))*U'

# form Hcore and diagonalize
Hcore = T + Vne
Ecore,Ccore = eig(X'*Hcore*X)

### in alternative, I can avoid to compute X and directly do
# Ecore,Ccore = eig(Hcore,S)

# sort energies and orbitals
idx = sortperm(Ecore)
Ccore = Ccore[:,idx]
Ecore = Ecore[idx]

# transform coefs to old basis
Ccore = X*Ccore

# compute density matrix
# D = 2.0 * (Ccore[:,1:occ] * Ccore[:,1:occ]')

Eold = 0.0

@printf("%12s%12s%12s\n","Iteration","Energy","ΔE")
# enter SCF loop
for iter = 1:maxiter

    # Fock matrix construction
    for i = 1:M
        for j = 1:M
            Vmf = 0.0
            for k = 1:M
                for l = 1:M
                    Vmf += D[k,l] * (Vee[i,j,k,l] - 0.5*Vee[i,k,j,l])
                end
            end
            F[i,j] = T[i,j] + Vne[i,j] + Vmf
        end
    end

    # solve the Roothaan-Hall equations FC = SCE
    Efock,Cfock = eig(X'*F*X)

    # order evals and evecs
    idx = sortperm(Efock)
    Cfock = Cfock[:,idx]
    Efock = Efock[idx]

    # transform back coefs
    C = X*Cfock
    # construct density matrix
    D = 2.0 * (C[:,1:occ] * C[:,1:occ]')

    Enew = 0.0
    for i = 1:M
        for j = 1:M
            Enew += D[j,i]*(Hcore[i,j] + F[i,j])
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
