

"""Type representing an atom."""
struct Atom
    name::String
    xyz ::Vector{Float64}
    Z   ::Integer
    mass::Float64
end


"""Explicit constructor for `Atom` type."""
function Atom(name::String, coords::Vector{T}) where {T<:Real}
    Z = ptable[name]["Z"]
    M = ptable[name]["M"]
    return Atom(name, coords, Z, M)
end


### some getter functions ###
name(atom::Atom) = atom.name
xyz(atom::Atom) = atom.xyz
Z(atom::Atom) = atom.Z
mass(atom::Atom) = atom.mass
###


"""Type representing an AO basis."""
struct AOBasis
    name::String
    cgto::Vector{CGTO}
end


### some getter functions ###
name(basis::AOBasis) = basis.name
dimension(basis::AOBasis) = length(basis.cgto)
contractions(basis::AOBasis) = basis.cgto
###


"""Type representing a molecule."""
struct Molecule
    atoms::Vector{Atom}
    basis::AOBasis
    charge::Integer
end


### some getter functions ###
atoms(mol::Molecule) = mol.atoms
basis(mol::Molecule) = mol.basis
charge(mol::Molecule) = mol.charge
electrons(mol::Molecule) = sum(map(x->Z(x),atoms(mol))) - charge(mol)
###


using JSON.Parser.parsefile


"""Load a basis set file stored in JSON format."""
loadbasis(filename::String) = parsefile(filename)


"""Parse xyz file"""
function readxyz(filename::String)
    f = open(filename,"r")
    natoms = parse(readline(f))
    title = readline(f)
    atoms = readdlm(f)
    close(f)

    # check that the xyz file is valid
    if natoms != size(atoms,1)
        error("number of atoms doesn't match the number declared")
    end

    println("Reading xyzfile '$filename'")

    types  = convert(Array{String,1}, atoms[:,1])
    coords = convert(Array{Float64,2}, atoms[:,2:end])

    return types,coords
end


"""Specialized constructor of the Molecule type."""
function Molecule(xyzfile::String, basisname::String, charge::Int = 0)

    atm,coords = readxyz(xyzfile)
    natm = length(atm)

    atoms = Vector{Atom}(natm)

    # generate list of `Atom`
    for i = 1:natm
        atoms[i] = Atom(atm[i], coords[i,:])
    end

    # load basis from file
    basisdata = loadbasis("src/basis/$basisname.json")
    
    basis = Vector{CGTO}()
    
    for atom in atoms
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

    return Molecule(atoms, AOBasis(basisname,basis), charge)
end


# periodic table
ptable = Dict("H"  => Dict("Z"=>1, "M"=>  1.008),
              "He" => Dict("Z"=>2, "M"=>  4.003),
              "Li" => Dict("Z"=>3, "M"=>  6.940),
              "Be" => Dict("Z"=>4, "M"=>  9.012),
              "B"  => Dict("Z"=>5, "M"=> 10.810),
              "C"  => Dict("Z"=>6, "M"=> 12.010),
              "N"  => Dict("Z"=>7, "M"=> 14.010),
              "O"  => Dict("Z"=>8, "M"=> 16.000),
              "F"  => Dict("Z"=>9, "M"=> 19.000),
              "Ne" => Dict("Z"=>10,"M"=> 20.180),
              "Br" => Dict("Z"=>35,"M"=> 79.904),
              "I"  => Dict("Z"=>53,"M"=>126.904))

