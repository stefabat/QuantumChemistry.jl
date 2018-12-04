
using JQC
using Test

tests = ["integrals"]

for test in tests
    include("$test.jl")
end
