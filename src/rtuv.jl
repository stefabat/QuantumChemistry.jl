# DEVELOPMENT FILE

include("utils.jl")

# using Memoize
# @memoize function Rtuv(t::Int, u::Int, v::Int, n::Int, p::Real, RPC::Vector{Float64})::Float64

# """Returns the integral of an Hermite Gaussian divided by the Coulomb operator."""
# function Rtuv{T<:Real}(t::Int, u::Int, v::Int, n::Int, p::Real, RPC::Vector{T})

### Note that for t+u+v odd -> Rtuv is zero!!
## Note also that there are symmetries in t, u and v that can be used

"""Returns the integral of an Hermite Gaussian divided by the Coulomb operator."""
function Rtuv(t::Int, u::Int, v::Int, n::Int, p::Real, RPC::Vector{T}) where {T<:Real}
    if t == u == v == 0
        return (-2.0*p)^n * boys(n,p*abs(sum(RPC.*RPC)))
    elseif u == v == 0
        if t > 1
            return  (t-1)*Rtuv(t-2, u, v, n+1, p, RPC) +
                   RPC[1]*Rtuv(t-1, u, v, n+1, p, RPC)
        else
            return RPC[1]*Rtuv(t-1, u, v, n+1, p, RPC)
        end
    elseif v == 0
        if u > 1
            return  (u-1)*Rtuv(t, u-2, v, n+1, p, RPC) +
                   RPC[2]*Rtuv(t, u-1, v, n+1, p, RPC)
        else
            return RPC[2]*Rtuv(t, u-1, v, n+1, p ,RPC)
        end
    else
        if v > 1
            return  (v-1)*Rtuv(t, u, v-2, n+1, p, RPC) +
                   RPC[3]*Rtuv(t, u, v-1, n+1, p, RPC)
        else
            return RPC[3]*Rtuv(t, u, v-1, n+1, p, RPC)
        end
    end
end


function Rtuv_exp(t::Int, u::Int, v::Int, n, p, RPC)
    
    nn = :($n+1)

    if t == u == v == 0
        return :((-2.0 * $p)^$n * boys($n, $p * sum($RPC.^2)))
    elseif u == v == 0
        if t > 1
            return :($(t-1) * $(Rtuv_exp(t-2, u, v, nn, p, RPC)) +
                    $RPC[1] * $(Rtuv_exp(t-1, u, v, nn, p, RPC)))
        else
            return :($RPC[1] * $(Rtuv_exp(t-1, u, v, nn, p, RPC)))
        end
    elseif v == 0
        if u > 1
            return :($(u-1) * $(Rtuv_exp(t, u-2, v, nn, p, RPC)) +
                    $RPC[2] * $(Rtuv_exp(t, u-1, v, nn, p, RPC)))
        else
            return :($RPC[2] * $(Rtuv_exp(t, u-1, v, nn, p ,RPC)))
        end
    else
        if v > 1
            return :($(v-1) * $(Rtuv_exp(t, u, v-2, nn, p, RPC)) +
                    $RPC[3] * $(Rtuv_exp(t, u, v-1, nn, p, RPC)))
        else
            return :($RPC[3] * $(Rtuv_exp(t, u, v-1, nn, p, RPC)))
        end
    end
end


@generated function Rtuv_gen(::Type{Val{t}}, ::Type{Val{u}}, ::Type{Val{v}},
                            n::Int, p::Real, RPC::Vector{Float64}) where {t,u,v}
    return Rtuv_exp(t, u, v, :n, :p, :RPC)
end
