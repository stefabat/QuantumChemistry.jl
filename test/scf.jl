
include("../src/scf.jl")

mol = Molecule(["H" 0 0 0;"H" 0 0 1.889726132886])
basis = AOBasis("sto-3g",mol)

rhf(mol, basis)
