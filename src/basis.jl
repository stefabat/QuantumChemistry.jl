
"""Type representing an atomic orbital basis."""
struct Basis
    name::String
    shells::Vector{Shell}
end

### some getter functions ###
get_name(basis::Basis) = basis.name
get_shells(basis::Basis) = basis.shells
get_num_shells(basis::Basis) = length(basis.shells)

function get_num_primitives(basis::Basis)
    return sum([get_num_primitives(shell) for shell in basis.shells])
end

function get_num_contractions(basis::Basis)
    return sum([get_num_contractions(shell) for shell in basis.shells])
end

function dimension(basis::Basis)
    m = 0
    n = 0
    for shell in get_shells(basis)
        l = get_am(shell)
        m += get_num_primitives(shell)   * (2*l+1)
        n += get_num_contractions(shell) * (2*l+1)
    end
    return m,n
end
###


# """Load a basis set file stored in JSON format."""
# loadbasis(filename::String) = parsefile(filename)


"""Construct an atomic orbital basis for a given molecule."""
function Basis(name::String, molecule::Molecule)

    # parse basis set data from JSON file
    # basisdata = loadbasis("src/basis/$name.json")
    basis = parsefile("src/basis/$name.1.json")

    shells = Vector{Shell}()

    for atom in molecule
        Z = get_atomic_number(atom)
        R = get_center(atom)
        for shell in basis["elements"]["$Z"]["electron_shells"]
            # the angular momentum is in a 1-element array, so we always access the
            # first element, which is already returned as an integer
            l = shell["angular_momentum"][1]
            # the exponents in the JSON file are encoded as an arrays of strings, so we need
            # to parse each element of the array into a float
            α = parse.(Float64,shell["exponents"])
            # the coefficients matrix in the JSON file is encoded as an array of arrays of strings,
            # so, first we "splat it" to a matrix of strings and then we parse each element to float
            d = parse.(Float64,hcat(shell["coefficients"]...))

            # add the shell to the
            push!(shells,Shell(R,l,d,α))
        end
    end

    return Basis(name,shells)
end


# """Compute the overlap matrix S of a given `AOBasis`."""
# function overlap(basis::AOBasis)
#     M = dimension(basis)
#     AOs = contractions(basis)
#     S = zeros(M,M)
#     # permutational symmetry Sij = Sji
#     for i = 1:M
#         for j = 1:i
#             S[i,j] = S[j,i] = overlap(AOs[i],AOs[j])
#         end
#     end
#     return S
# end


# """Compute the kinetic energy matrix T of a given `AOBasis`."""
# function kinetic(basis::AOBasis)
#     M = dimension(basis)
#     AOs = contractions(basis)
#     T = zeros(M,M)
#     # permutational symmetry Tij = Tji
#     for i = 1:M
#         for j = 1:i
#             T[i,j] = T[j,i] = kinetic(AOs[i],AOs[j])
#         end
#     end
#     return T
# end


# """Compute the nuclei-electron attraction potential matrix Vne of a given `AOBasis`."""
# function attraction(basis::AOBasis, mol::Molecule)
#     M = dimension(basis)
#     AOs = contractions(basis)
#     Vne = zeros(M,M)
#     # permutational symmetry Vij = Vji
#     for i = 1:M
#         for j = 1:i
#             for I in atoms(mol)
#                 Vne[i,j] += attraction(AOs[i],AOs[j],xyz(I))*(-Z(I))
#             end
#             Vne[j,i] = Vne[i,j]
#         end
#     end
#     return Vne
# end


# """Compute the electron-electron repulsion tensor Vee of a given `AOBasis`."""
# function repulsion(basis::AOBasis)
#     M = dimension(basis)
#     AOs = contractions(basis)
#     Vee = zeros((M*M-M)÷2+M,(M*M-M)÷2+M)
#     for i = 1:M
#         for j = 1:i
#             ij = index(i,j)
#             for k = 1:i
#                 for l = 1:k
#                     kl = index(k,l)
#                     if ij >= kl
#                         # println(i,j,k,l)
#                         Vee[ij,kl] = Vee[kl,ij] = repulsion(AOs[i],AOs[j],AOs[k],AOs[l])
#                     end
#                 end
#             end
#         end
#     end
#     return Vee
# end


# """Return compound index of the ERI tensor."""
# index(i::Int, j::Int) = i >= j ? (i*i - i) ÷ 2 + j : (j*j - j) ÷ 2 + i


# """Compute the nuclear-nuclear repulsion for a given `Molecule`."""
# function nuclear(mol::Molecule)
#     Vnn = 0.0
#     for (I,Inuc) in enumerate(atoms(mol))
#         for (J,Jnuc) in enumerate(atoms(mol))
#             if I<J
#                 Vnn += Z(Inuc)*Z(Jnuc)/norm(xyz(Inuc)-xyz(Jnuc))
#             end
#         end
#     end

#     return Vnn
# end


# """
# Compute the electron-electron repulsion tensor Vee of a given `AOBasis`.
# Use only for test purposes on small systems, EXTREMELY SLOW!
# """
# function repulsion_debug(basis::AOBasis)
#     M = dimension(basis)
#     AOs = contractions(basis)
#     Vee = zeros(M,M,M,M)
#     for i = 1:M
#         for j = 1:M
#             for k = 1:M
#                 for l = 1:M
#                     Vee[i,j,k,l] = repulsion(AOs[i],AOs[j],AOs[k],AOs[l])
#                 end
#             end
#         end
#     end
#     return Vee
# end

