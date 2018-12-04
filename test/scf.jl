
include("../src/scf.jl")

mol = Molecule(["H" 0 0 0;"H" 0 0 1.8897162])
basis = AOBasis("sto-3g",H2)

rhf(mol, basis)
