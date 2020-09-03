

# module Geometry

using StaticArrays
using DelimitedFiles
using JSON.Parser: parsefile


"""Type representing an element."""
struct Element
    name::String
    symbol::String
    number::Int
    mass::Float64
end

get_element_number(element::Element) = element.number

"""Explicit constructor for `Element` type"""
function Element(symbol::String)
    table = parsefile("src/data/PeriodicTable.json")["elements"]
    for element in table
        if isequal(element["symbol"],symbol)
            name = element["name"]
            number = element["number"]
            mass = element["atomic_mass"]
            return Element(name, symbol, number, mass)
        else
            error("Element not found!")
        end
    end
end


"""Type representing an atom."""
struct Atom
    element::Element
    center::SVector{3,Float64}
end

get_center(atom::Atom) = atom.center
get_atomic_number(atom::Atom) = get_element_number(atom.element)

"""Explicit constructor for `Atom` type."""
function Atom(symbol::String, xyz::Vector{T}) where {T<:Real}
    return Atom(Element(symbol), xyz)
end


"""Type representing a molecule."""
struct Molecule
    atoms::Vector{Atom}
    charge::Int
end

function Base.getindex(mol::Molecule, index)
    return mol.atoms[index]
end
  
function Base.length(mol::Molecule)
    return length(mol.atoms)
end
  
function Base.iterate(mol::Molecule)
    return iterate(mol.atoms)
end
  
function Base.iterate(mol::Molecule, state)
    return iterate(mol.atoms, state)
end
  
function Base.push!(mol::Molecule, atom::Atom)
    return push!(mol.atoms, atom)
end

"""Returns the coordinates of `mol` in a Nx3 matrix."""
function coordinates(mol::Molecule)
    xyz = @SMatrix zeros(length(mol),3)
    for (i,atom) in enumerate(mol)
        xyz[i,:] = atom.center
    end
    return xyz
end


"""
Simple constructor for `Molecule` type.
Geometry expected in Angstrom.    
"""
function Molecule(coords::Matrix, charge::Int = 0)

    # input sanity check
    if size(coords,2) != 4
        error("expected a matrix with format ['Atom' x y z;...] ")
    end

    natoms = size(coords,1)
    atoms = Vector{Atom}(undef,natoms)

    for i = 1:natoms
        atoms[i] = Atom(coords[i,1],tobohr*coords[i,2:4])
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
    natoms = parse(Int,readline(f))
    atoms = readdlm(f,skipstart=1)
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
# ptable = Dict("H"  => Dict("Z"=>1, "Ar"=>  1.008),
#               "He" => Dict("Z"=>2, "Ar"=>  4.003),
#               "Li" => Dict("Z"=>3, "Ar"=>  6.940),
#               "Be" => Dict("Z"=>4, "Ar"=>  9.012),
#               "B"  => Dict("Z"=>5, "Ar"=> 10.810),
#               "C"  => Dict("Z"=>6, "Ar"=> 12.010),
#               "N"  => Dict("Z"=>7, "Ar"=> 14.010),
#               "O"  => Dict("Z"=>8, "Ar"=> 16.000),
#               "F"  => Dict("Z"=>9, "Ar"=> 19.000),
#               "Ne" => Dict("Z"=>10,"Ar"=> 20.180),
#               "Mg" => Dict("Z"=>12,"Ar"=> 24.305),
#               "Ti" => Dict("Z"=>22,"Ar"=> 47.867),
#               "Br" => Dict("Z"=>35,"Ar"=> 79.904),
#               "I"  => Dict("Z"=>53,"Ar"=>126.904))

# using Psi4 physical constants
toang = 0.52917720859
tobohr = 1.8897261328856432 # == 1.0/toang

# end
