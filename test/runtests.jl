
using QuantumChemistry
using Test

tests = ["mcmurchie-davidson"]

for test in tests
    include("$test.jl")
end
