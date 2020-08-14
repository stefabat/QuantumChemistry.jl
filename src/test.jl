
using BenchmarkTools

includet("src/integrals.jl")

Ga100 = CGF(zeros(3),0.5,(1,0,0))
Ga010 = CGF(zeros(3),0.5,(0,1,0))
Ga001 = CGF(zeros(3),0.5,(0,0,1))

Gb100 = CGF(ones(3),0.5,(1,0,0))
Gb010 = CGF(ones(3),0.5,(0,1,0))
Gb001 = CGF(ones(3),0.5,(0,0,1))

Gabasis = [Ga100,Ga010,Ga001]
Gbbasis = [Gb100,Gb010,Gb001]

Gabas = [CGF(zeros(3),0.5,(0,0,2)),
         CGF(zeros(3),0.5,(0,1,1)),
         CGF(zeros(3),0.5,(0,2,0)),
         CGF(zeros(3),0.5,(1,0,1)),
         CGF(zeros(3),0.5,(1,1,0)),
         CGF(zeros(3),0.5,(2,0,0))]

Gbbas = [CGF(ones(3),0.5,(0,0,2)),
         CGF(ones(3),0.5,(0,1,1)),
         CGF(ones(3),0.5,(0,2,0)),
         CGF(ones(3),0.5,(1,0,1)),
         CGF(ones(3),0.5,(1,1,0)),
         CGF(ones(3),0.5,(2,0,0))]

mom0=[(0,0,0)]
mom1=[(1,0,0),(0,1,0),(0,0,1)]
mom2=[(0,0,2),(0,1,1),(0,2,0),(1,0,1),(1,1,0),(2,0,0)]

function overlap_basis(Ga,Gb)
    na = length(Ga)
    nb = length(Gb)
    S = zeros(na,nb)
    for i = 1:na
        for j = 1:nb
            S[i,j] = overlap(Ga[i],Gb[j])
        end
    end
    return S
end

function kinetic_basis(Ga,Gb)
    na = length(Ga)
    nb = length(Gb)
    T = zeros(na,nb)
    for i = 1:na
        for j = 1:nb
            T[i,j] = kinetic(Ga[i],Gb[j])
        end
    end
    return S
end

function attraction_basis(Ga,Gb,C)
    Na = length(Ga)
    Nb = length(Gb)
    V = zeros(Na,Nb)
    for i = 1:Na
        for j = 1:Nb
            V[i,j] = attraction(Ga[i],Gb[j],C)
        end
    end
    return V
end

Ga = GTOShell(zeros(3),0.5,2)
Gb = GTOShell(ones(3),0.5,2)

@btime overlap_basis($Gabas,$Gbbas)

@btime overlap($Ga,$Gb)
