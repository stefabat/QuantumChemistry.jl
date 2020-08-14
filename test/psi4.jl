using PyCall

psi4 = pyimport("psi4")

mol = psi4.geometry("""
                    symmetry c1
                    C
                    C 1 1.421
                    """)

#psi4.set_options(Dict("scf_type"=>"pk"))
#psi4.set_options(Dict("basis"=>"sto-3g"))
psi4.set_options(Dict("puream"=>false))
#psi4.set_options(Dict("guess"=>"core"))
#psi4.set_options(Dict("diis"=>false))

# en,wfn = psi4.energy("scf/sto-3g", return_wfn=true)
wfn = psi4.core.Wavefunction.build(mol,"sto-3g")

mints = psi4.core.MintsHelper(wfn.basisset())

mints.ao_overlap().to_array()

