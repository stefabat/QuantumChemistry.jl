#
# Integrals are tested against values obtained with Psi4
# using the libint library on an Arch Linux x86_64 machine,
# kernal version 4.19.4-arch1-1-ARCH and CPU Intel i5-7200U
#
# Conda version:  4.5.11
# Python version: 3.6.7         (build h0371630_0)
# Psi4 version:   1.2.1+406f4de (build py36hf231b52_0)
# Libint version: 1.2.1         (build h87b9b30_4)
# MKL version:    2018.0.3
# Numpy version:  1.15.4        (build py36h1d66e8a_0)
#

using Test

include("../src/molecule.jl")
include("../src/basis.jl")

@testset "integrals" begin
    # H2/sto-3g at 1.0 Ang distance
    mol = Molecule(["H" 0 0 0;"H" 0 0 1.0])
    basis = AOBasis("sto-3g",mol)
    tol = 1e-12

    # test overlap
    S = overlap(basis)
    @test S[1,1] ≈ 0.99999999999999970 atol=tol
    @test S[1,2] ≈ 0.49648468902916987 atol=tol

    # test kinetic energy
    T = kinetic(basis)
    @test T[1,1] ≈ 0.76003188356660870 atol=tol
    @test T[1,2] ≈ 0.11131040276424069 atol=tol

    # test nuclear-electron attraction
    V = attraction(basis,mol)
    @test V[1,1] ≈ -1.7395282585031457 atol=tol
    @test V[1,2] ≈ -0.7941753316682432 atol=tol

    # test electron-electron repulsion
    W = repulsion(basis)
    @test W[index(1,1),index(1,1)] ≈ 0.77460594391989710 atol=tol
    @test W[index(2,1),index(1,1)] ≈ 0.30930896514686950 atol=tol
    @test W[index(2,2),index(1,1)] ≈ 0.47804137172137210 atol=tol
    @test W[index(1,2),index(2,1)] ≈ 0.15786577605541174 atol=tol

    # test nuclear-nuclear repulsion
    @test nuclear(mol) ≈ 0.5291772085899999 atol=tol

    # C2/cc-pvdz
    # mol = Molecule([
    #     "C" 0.00000 0.00000 0.00000;
    #     "C" 0.00000 0.00000 1.42100])
    # basis = AOBasis("cc-pvdz",mol)
    # tol = 5e-8

end
