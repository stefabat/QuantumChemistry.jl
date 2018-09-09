
using SpecialFunctions


function boys(n,x,maxIter = 100)
    
    acc = 0.0
    for k = 0:maxIter
        acc += (-x)^k/(reduce(*,1:k)*(2*n + 2*k + 1))
    end

    return acc
end

function gsl_sf_hyperg_1F1(a::AbstractFloat, b::AbstractFloat, x::AbstractFloat)
    return ccall( (:gsl_sf_hyperg_1F1, "/usr/lib/libgsl.so"), Cdouble, (Cdouble, Cdouble, Cdouble), a, b, x)
end

function gsl_boys(n,x)
    return gsl_sf_hyperg_1F1(n+0.5,n+1.5,-x)/(2*n+1)
end


function gsl_sf_gamma_inc(a::AbstractFloat, x::AbstractFloat)
    return ccall( (:gsl_sf_gamma_inc_P, "/usr/lib/libgsl.so"), Cdouble, (Cdouble, Cdouble), a, x)
end


function gser(a::AbstractFloat, x::AbstractFloat, ϵ::AbstractFloat = 1e-16)
    
    ap = a
    lnΓ = lgamma(a)
    acc = 1.0/a
    frc = 1.0/a

    while(frc > acc*ϵ)
        ap  += 1
        frc *= x/ap
        acc += frc
    end

    return acc * exp(-x+a*log(x)-lnΓ)
end


function gcf(a::AbstractFloat, x::AbstractFloat,
             ϵ::AbstractFloat = 1e-16,
            ϵ2::AbstractFloat = 1e-16)

    lnΓ = lgamma(a)
    b = x + 1.0 - a
    c = 1.0/ϵ2
    d = 1.0/b
    h = d
    i = 0
    while true
        i += 1
        an = -i * (i-a)
        b += 2.0
        d = an * d + b
        if abs(d) < ϵ2
            d = ϵ2
        end
        c = b + an/c
        if abs(c) < ϵ2
            c = ϵ2
        end
        d = 1.0/d
        del = d*c
        h *= del
        if abs(del-1.0) <= ϵ
            break
        end
    end
    return 1.0 - exp(-x + a*log(x) - lnΓ) * h
end

