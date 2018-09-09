
include("utils.jl")

# using Memoize
# @memoize function Rtuv(t::Int, u::Int, v::Int, n::Int, p::Real, RPC::Vector{Float64})::Float64

# """Returns the integral of an Hermite Gaussian divided by the Coulomb operator."""
# function Rtuv{T<:Real}(t::Int, u::Int, v::Int, n::Int, p::Real, RPC::Vector{T})

### Note that for t+u+v odd -> Rtuv is zero!!
## Note also that there are symmetries in t, u and v that can be used

function Rtuv(t::Int, u::Int, v::Int, n::Int, p::Real, RPC::Vector{Float64})
    
    if t < 0 || u < 0 || v < 0
        return 0.0
    elseif t == u == v == 0
        return (-2.0*p)^n * boys(n,p*norm(RPC)^2)
        # return p
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


function Rtuv_exp(t::Int, u::Int, v::Int, n::Int, p, RPC)

    if t < 0 || u < 0 || v < 0
        return 0.0
    elseif t == u == v == 0
        # return :($p)
        return :((-2.0 * $p)^$n * boys($n, $p * norm($RPC)^2))
    elseif u == v == 0
        if t > 1
            return :($(t-1) * $(Rtuv_exp(t-2, u, v, n+1, p, RPC)) +
                    $RPC[1] * $(Rtuv_exp(t-1, u, v, n+1, p, RPC)))
        else
            return :($RPC[1] * $(Rtuv_exp(t-1, u, v, n+1, p, RPC)))
        end
    elseif v == 0
        if u > 1
            return :($(u-1) * $(Rtuv_exp(t, u-2, v, n+1, p, RPC)) +
                    $RPC[2] * $(Rtuv_exp(t, u-1, v, n+1, p, RPC)))
        else
            return :($RPC[2] * $(Rtuv_exp(t, u-1, v, n+1, p ,RPC)))
        end
    else
        if v > 1
            return :($(v-1) * $(Rtuv_exp(t, u, v-2, n+1, p, RPC)) +
                    $RPC[3] * $(Rtuv_exp(t, u, v-1, n+1, p, RPC)))
        else
            return :($RPC[3] * $(Rtuv_exp(t, u, v-1, n+1, p, RPC)))
        end
    end
end


@generated function Rtuv_gen(::Type{Val{t}}, ::Type{Val{u}}, ::Type{Val{v}},
                             ::Type{Val{n}}, p::Real, RPC::Vector{Float64}) where {t,u,v,n}
    return Rtuv_exp(t, u, v, n, :p, :RPC)
end

@generated function Rtuv_tot(::Type{Val{t}}, ::Type{Val{u}}, ::Type{Val{v}},
                            n::Int, p::Real, RPC::Vector{Float64}) where {t,u,v}
    if t == u == v == 0
        return :((-2.0*p)^n * boys(n,p*norm(RPC)^2))
    elseif u == v == 0
        if t > 1
            return :((t-1)*Rtuv_tot(Val{t-2}, Val{u}, Val{v}, n+1, p, RPC) +
                    RPC[1]*Rtuv_tot(Val{t-1}, Val{u}, Val{v}, n+1, p, RPC))
        else
            return :(RPC[1]*Rtuv_tot(Val{t-1}, Val{u}, Val{v}, n+1, p, RPC))
        end
    elseif v == 0
        if u > 1
            return :((u-1)*Rtuv_tot(Val{t}, Val{u-2}, Val{v}, n+1, p, RPC) +
                    RPC[2]*Rtuv_tot(Val{t}, Val{u-1}, Val{v}, n+1, p, RPC))
        else
            return :(RPC[2]*Rtuv_tot(Val{t}, Val{u-1}, Val{v}, n+1, p ,RPC))
        end
    else
        if v > 1
            return :((v-1)*Rtuv_tot(Val{t}, Val{u}, Val{v-2}, n+1, p, RPC) +
                    RPC[3]*Rtuv_tot(Val{t}, Val{u}, Val{v-1}, n+1, p, RPC))
        else
            return :(RPC[3]*Rtuv_tot(Val{t}, Val{u}, Val{v-1}, n+1, p, RPC))
        end
    end
end

