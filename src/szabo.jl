

module Integrals

export CGF
export center,ltot,alpha
export sz_overlap,sz_kinetic,sz_nuclear

struct CGF
    A::Vector{Float64}
    α::Float64
    ijk::Tuple{Int,Int,Int}
end

center(gto::CGF) = gto.A
ltot(gto::CGF) = sum(gto.ijk)
alpha(gto::CGF) = gto.α

function sz_overlap(Ga::CGF, Gb::CGF)
    assert(ltot(Ga) == ltot(Gb) == 0)
    α  = alpha(Ga);  β  = alpha(Gb)
    RA = center(Ga); RB = center(Gb)
    return (π/(α+β))^1.5 * exp(-α*β*norm(RA-RB)^2/(α+β))
end

function sz_kinetic(Ga::CGF, Gb::CGF)
    assert(ltot(Ga) == ltot(Gb) == 0)
    α  = alpha(Ga);  β  = alpha(Gb)
    RA = center(Ga); RB = center(Gb)
    μ = α*β/(α+β)
    p = α + β
    RAB = norm(RA-RB)
    return μ*(3-2*μ*RAB^2)*(π/p)^1.5*exp(-μ*RAB^2)
end

function sz_nuclear(Ga::CGF, Gb::CGF, RC::Vector{Float64})
    assert(ltot(Ga) == ltot(Gb) == 0)
    α  = alpha(Ga);  β  = alpha(Gb)
    RA = center(Ga); RB = center(Gb)
    μ = α*β/(α+β)
    p = α + β
    RP = (α*RA + β*RB)/(α + β)
    RAB = norm(RA-RB)
    RPC = norm(RP-RC)
    if p*RPC^2 < 1e-14
        return 2*π/p *exp(-μ*RAB^2)
    else
        return 2*π/p *exp(-μ*RAB^2)*F0(p*RPC^2)
    end
end

using SpecialFunctions

function F0(t)
    return 0.5*sqrt(π/t)*erf(sqrt(t))
end

end # enf of module
