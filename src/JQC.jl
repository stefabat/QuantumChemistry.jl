
module JQC

    include("molecule.jl")
    export Atom,name,xyz,Z,mass
    export Molecule,atoms,natoms,charge,electrons,xyz

    include("integrals.jl")
    export BasisFunction,GTO,CGTO
    export center,exponent,shell,Ax,Ay,Az,lx,ly,lz,ltot
    export overlap,normalization,kinetic,attraction,repulsion
    export primitives,nprimitives,coefs,norms,exponents

    include("basis.jl")
    export AbstractBasis,AOBasis
    export name,dimension,contractions,index
    export overlap,kinetic,attraction,repulsion
    export repulsion_debug

    include("scf.jl")
    export rhf

end
