
"""
Returns the Hermite expansion coefficients for a 1D Cartesian overlap distribution
using a two-term recursion relation.
"""
function Etij(t::Int ,i::Int, j::Int, Kαβx::Real, XPA::Real, XPB::Real, α::Real, β::Real)

    # compute overlap exponent
    p = α + β

    # enter recursion
    if t < 0 || t > i+j
        return 0.0
    elseif t == 0
        if i == j == 0
            return Kαβx
        elseif j == 0
            return XPA * Etij(0, i-1, j, Kαβx, XPA, XPB, α, β) +
                         Etij(1, i-1, j, Kαβx, XPA, XPB, α, β)
        else
            return XPB * Etij(0, i, j-1, Kαβx, XPA, XPB, α, β) +
                         Etij(1, i, j-1, Kαβx, XPA, XPB, α, β)
        end
    else
        return (1/(2*p*t)) * (i * Etij(t-1, i-1, j, Kαβx, XPA, XPB, α, β) +
                              j * Etij(t-1, i, j-1, Kαβx, XPA, XPB, α, β) )
    end
end


"""
Returns the Hermite expansion coefficients for a 1D Cartesian overlap distribution
using a two-term recursion relation.
"""
function Etij_exp(t::Int ,i::Int, j::Int, Kαβx, XPA, XPB, α, β)

    # compute overlap exponent
    p = :($α + $β)

    # enter recursion
    if t < 0 || t > i+j
        return 0.0
    elseif t == 0
        if i == j == 0
            return :($Kαβx)
        elseif j == 0
            return :($XPA * $(Etij_exp(0, i-1, j, Kαβx, XPA, XPB, α, β)) +
                            $(Etij_exp(1, i-1, j, Kαβx, XPA, XPB, α, β)))
        else
            return :($XPB * $(Etij_exp(0, i, j-1, Kαβx, XPA, XPB, α, β)) +
                            $(Etij_exp(1, i, j-1, Kαβx, XPA, XPB, α, β)))
        end
    else
        return :((1/(2*$p*$t)) * ($i * $(Etij_exp(t-1, i-1, j, Kαβx, XPA, XPB, α, β)) +
                                  $j * $(Etij_exp(t-1, i, j-1, Kαβx, XPA, XPB, α, β))) )
    end
end

@generated function Etij_gen(::Type{Val{t}}, ::Type{Val{i}}, ::Type{Val{j}},
    Kαβx, XPA, XPB, α, β) where {t,i,j}
    return Etij_exp(t, i, j, :Kαβx, :XPA, :XPB, :α, :β)
end
