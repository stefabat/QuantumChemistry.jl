
# include("basis.jl")


"""Type representing an atom."""
struct Atom
    name::String
    xyz ::Vector{Float64}
    Z   ::Integer
    mass::Float64
end


"""Explicit constructor for `Atom` type."""
function Atom(name::String, coords::Vector{T}) where {T<:Real}
    if haskey(ptable, name)
        Z  = ptable[name]["Z"]      # Atomic number
        Ar = ptable[name]["Ar"]     # Relative atomic weight
        return Atom(name, coords, Z, Ar)
    else
        error("element $name not found")
    end
end


### some getter functions ###
name(atom::Atom) = atom.name
xyz(atom::Atom) = atom.xyz
Z(atom::Atom) = atom.Z
mass(atom::Atom) = atom.mass
###


"""Type representing a molecule."""
struct Molecule
    atoms::Vector{Atom}
    charge::Integer
end


### some getter functions ###
atoms(mol::Molecule) = mol.atoms
natoms(mol::Molecule) = length(atoms(mol))
# basis(mol::Molecule) = mol.basis
charge(mol::Molecule) = mol.charge
electrons(mol::Molecule) = sum(map(x->Z(x),atoms(mol))) - charge(mol)
"""Returns the coordinates of `mol` in a Nx3 matrix."""
function xyz(mol::Molecule)
    coords = Matrix{Float64}(undef,natoms(mol),3)
    for i = 1:natoms(mol)
        coords[i,:] = xyz(atoms(mol)[i])
    end
    return coords
end
###


"""Simple constructor for `Molecule` type."""
function Molecule(coords::Matrix, charge::Int = 0)
    if size(coords,2) != 4
        error("expected a matrix with format ['Atom' x y z;...] ")
    end

    natoms = size(coords,1)
    atoms = Vector{Atom}(undef,natoms)

    for i = 1:natoms
        atoms[i] = Atom(coords[i,1],convert(Vector{Float64},coords[i,2:4]))
    end

    return Molecule(atoms, charge)
end


"""Constructor of type `Molecule` from an xyz file."""
function Molecule(xyzfile::String, charge::Int = 0)

    types,coords = readxyz(xyzfile)

    natoms = size(coords,1)
    atoms = Vector{Atom}(undef,natoms)

    for i = 1:natoms
        atoms[i] = Atom(types[i],coords[i,:])
    end

    return Molecule(atoms, charge)
end


"""Parse xyz file."""
function readxyz(xyzfile::String)
    f = open(xyzfile,"r")
    natoms = parse(readline(f))
    title = readline(f)
    atoms = read(f)
    close(f)

    # check that the xyz file is valid
    if natoms != size(atoms,1)
        error("number of atoms doesn't match the number declared")
    end

    println("Reading xyzfile '$xyzfile'")

    types  = convert(Array{String,1}, atoms[:,1])
    coords = convert(Array{Float64,2}, atoms[:,2:end])

    return types,coords
end



# periodic table
ptable = Dict("H"  => Dict("Z"=>1, "Ar"=>  1.008),
              "He" => Dict("Z"=>2, "Ar"=>  4.003),
              "Li" => Dict("Z"=>3, "Ar"=>  6.940),
              "Be" => Dict("Z"=>4, "Ar"=>  9.012),
              "B"  => Dict("Z"=>5, "Ar"=> 10.810),
              "C"  => Dict("Z"=>6, "Ar"=> 12.010),
              "N"  => Dict("Z"=>7, "Ar"=> 14.010),
              "O"  => Dict("Z"=>8, "Ar"=> 16.000),
              "F"  => Dict("Z"=>9, "Ar"=> 19.000),
              "Ne" => Dict("Z"=>10,"Ar"=> 20.180),
              "Br" => Dict("Z"=>35,"Ar"=> 79.904),
              "I"  => Dict("Z"=>53,"Ar"=>126.904))

