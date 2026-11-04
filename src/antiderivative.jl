#
# antiderivative.jl --
#
# Antiderivatives of piecewise polynomial interpolation kernels.
#
#-------------------------------------------------------------------------------------------
#
# This file is part of the InterpolationKernels package licensed under the MIT "Expat"
# License.
#

"""
    InterpolationKernels.Antiderivative{T,S,K}(ker)

Callable object representing the antiderivative of the kernel `ker` of support size `S` and
floating-point type `T`, normalized to vanish at the left end of the support. 

When called, `A(x)` is the integral of the kernel `ker` from the
left end of its support to `x`. `A(x)` is zero for `x` below the support and, for `x` above
the support, equals the integral of `ker` over its whole support (one for normalized
kernels).

The kernel must be piecewise polynomial of degree at most 5 between unit-spaced breakpoints
starting at `-S/2` with `S = length(ker)`. This is the case of B-splines, generic cubic
splines, their derivatives, and the Keys, Catmull-Rom and Mitchell & Netravali kernels.

The antiderivative of `ker` is exact, up to rounding errors, for the supported kernels.
"""
struct Antiderivative{T <: AbstractFloat, S, K <: Kernel{T, S}} <: Function
    kernel::K
    function Antiderivative(ker::K) where {T, S, K <: Kernel{T, S}}
        ker isa PiecewisePolynomialKernel || throw(DomainError(ker, "Antiderivative is only implemented for piecewise polynomial kernels of degree at most 5."))
        return new{T, S, K}(ker)
    end
end

const PiecewisePolynomialKernel{T} = Union{
    BSpline{<:Any, T}, BSplinePrime{<:Any, T},
    AbstractCubicSpline{T}, AbstractCubicSplinePrime{T},
}


Base.length(::Antiderivative{T, S}) where {T, S} = S
Base.eltype(::Type{<:Antiderivative{T}}) where {T} = T
Base.eltype(A::Antiderivative) = eltype(typeof(A))

# 3-point Gauss-Legendre rule on [-1,1], exact for polynomials up to degree 5.
function gauss_legendre_3(::Type{T}) where {T <: AbstractFloat}
    x = sqrt(T(3) / 5)
    return (-x, zero(T), x), (T(5) / 9, T(8) / 9, T(5) / 9)
end

# n-point Gauss-Legendre rule on [-1,1], exact for polynomials up to degree 5.
@inline function gauss_legendre(::Type{T}, n::Integer) where {T <: AbstractFloat}
    if n == 1
        return (zero(T),), (T(2),)
    elseif n == 2
        x = T(1 / √3)
        return (-x, x), (T(1), T(1))
    elseif n == 3
        x = T(sqrt(3 / 5))
        return (-x, zero(T), x), (T(5 / 9), T(8 / 9), T(5 / 9))
    end
end


# Integral of the kernel over `[a,b]` which must be inside a single polynomial piece.
function integrate_piece(ker::Kernel{T}, a::T, b::T) where {T}
    nodes, weights = gauss_legendre(T, length(ker) ÷ 2 + 1)
    mid = (a + b) / 2
    half = (b - a) / 2
    return half * sum(w * ker(mid + half * x) for (x, w) in zip(T.(nodes), T.(weights)))
end

(A::Antiderivative{T, S})(x::Real) where {T, S} = A(T(x)::T)

function (A::Antiderivative{T, S})(x::T) where {T <: AbstractFloat, S}
    ker = A.kernel
    hw = T(S) / 2
    isnan(x) && return x
    s = clamp(x, -hw, hw) + hw # position along the support
    c = min(floor(Int, s), S - 1) # index of the piece containing `x`
    result = zero(T)
    for k in 0:(c - 1)
        result += integrate_piece(ker, k - hw, k + 1 - hw)
    end
    return result + integrate_piece(ker, c - hw, s - hw)
end

Base.show(io::IO, A::Antiderivative) = print(io, "antiderivative(", A.kernel, ")")
