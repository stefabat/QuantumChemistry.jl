
include("integrals.jl")


"""Supertype of all types of basis sets."""
abstract type AbstractBasis end


"""Type representing an AO basis."""
struct AOBasis <: AbstractBasis
    name::String
    cgto::Vector{CGTO}
end

### some getter functions ###
name(basis::AOBasis) = basis.name
dimension(basis::AOBasis) = length(basis.cgto)
contractions(basis::AOBasis) = basis.cgto
###

using JSON: JSON.Parser.parsefile

"""Load a basis set file stored in JSON format."""
loadbasis(filename::String) = parsefile(filename)


# TODO: generalize it to higher angular momenta
# TODO: group functions by angular momentum instead of atoms?
"""Construct the AO basis for a given molecule."""
function AOBasis(basisname::String, mol::Molecule)

    # load basis data from disk to memory
    basisdata = loadbasis("src/basis/$basisname.json")

    # TODO: find a way to preallocate this?
    basis = Vector{CGTO}()

    for atom in atoms(mol)
        for fn in basisdata[name(atom)]
            α = convert(Vector{Float64},fn["prim"])
            d = convert(Vector{Float64},fn["cont"][1])
            if fn["angular"] == "s"
                push!(basis,CGTO(xyz(atom),(0,0,0),α,d))
            elseif fn["angular"] == "p"
                push!(basis,CGTO(xyz(atom),(1,0,0),α,d))
                push!(basis,CGTO(xyz(atom),(0,1,0),α,d))
                push!(basis,CGTO(xyz(atom),(0,0,1),α,d))
            elseif fn["angular"] == "d"
                push!(basis,CGTO(xyz(atom),(2,0,0),α,d))
                push!(basis,CGTO(xyz(atom),(0,2,0),α,d))
                push!(basis,CGTO(xyz(atom),(0,0,2),α,d))
                push!(basis,CGTO(xyz(atom),(1,1,0),α,d))
                push!(basis,CGTO(xyz(atom),(1,0,1),α,d))
                push!(basis,CGTO(xyz(atom),(0,1,1),α,d))
            elseif fn["angular"] == "f"
                push!(basis,CGTO(xyz(atom),(3,0,0),α,d))
                push!(basis,CGTO(xyz(atom),(0,3,0),α,d))
                push!(basis,CGTO(xyz(atom),(0,0,3),α,d))
                push!(basis,CGTO(xyz(atom),(2,1,0),α,d))
                push!(basis,CGTO(xyz(atom),(2,0,1),α,d))
                push!(basis,CGTO(xyz(atom),(1,2,0),α,d))
                push!(basis,CGTO(xyz(atom),(0,2,1),α,d))
                push!(basis,CGTO(xyz(atom),(1,0,2),α,d))
                push!(basis,CGTO(xyz(atom),(0,1,2),α,d))
                push!(basis,CGTO(xyz(atom),(1,1,1),α,d))
            end
        end
    end

    return AOBasis(basisname,basis)
end


function overlap(basis::AOBasis)
    M = dimension(basis)
    AOs = contractions(basis)
    S = zeros(M,M)
    for (i,χi) in enumerate(AOs)
        for (j,χj) in enumerate(AOs)
            S[i,j] = overlap(χi,χj)
        end
    end
    return S
end

function kinetic(basis::AOBasis)
    M = dimension(basis)
    AOs = contractions(basis)
    T = zeros(M,M)
    for (i,χi) in enumerate(AOs)
        for (j,χj) in enumerate(AOs)
            T[i,j] = kinetic(χi,χj)
        end
    end
    return T
end

# nuclear-electron attraction
function potential(basis::AOBasis, mol::Molecule)
    M = dimension(basis)
    AOs = contractions(basis)
    Vne = zeros(M,M)
    for (i,χi) in enumerate(AOs)
        for (j,χj) in enumerate(AOs)
            for I in atoms(mol)
                Vne[i,j] += nuclear(χi,χj,xyz(I))*(-Z(I))
            end
        end
    end
    return Vne
end

# Compute ERI integrals
function repulsion(basis::AOBasis)
    M = dimension(basis)
    AOs = contractions(basis)
    Vee = zeros(M,M,M,M)
    for (i,χi) in enumerate(AOs)
        for (j,χj) in enumerate(AOs)
            for (k,χk) in enumerate(AOs)
                for (l,χl) in enumerate(AOs)
                    Vee[i,j,k,l] = repulsion(χi,χj,χk,χl)
                end
            end
        end
    end
    return Vee
end

