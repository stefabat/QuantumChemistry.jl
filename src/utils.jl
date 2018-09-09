
"""Returns the Boys fucntion F_n(x)."""
function boys(n::Int, x::Real)
    return _gsl_sf_hyperg_1F1(n+0.5, n+1.5, -x) / (2.0*n+1.0)
end

function _gsl_sf_hyperg_1F1(a::AbstractFloat, b::AbstractFloat, x::AbstractFloat)
    return ccall( (:gsl_sf_hyperg_1F1, "/usr/lib/libgsl.so"),
               Cdouble, (Cdouble, Cdouble, Cdouble), a, b, x)
end
