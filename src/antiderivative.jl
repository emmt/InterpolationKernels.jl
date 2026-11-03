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
floating-point type `T`, normalized to vanish at the left end of the support. It is
usually built by [`antiderivative(ker)`](@ref antiderivative).

"""
struct Antiderivative{T<:AbstractFloat,S,K<:Kernel{T,S}} <: Function
    kernel::K
end

Antiderivative(ker::K) where {T,S,K<:Kernel{T,S}} = Antiderivative{T,S,K}(ker)

"""
    antiderivative(ker) -> A

Return a callable object `A` such that `A(x)` is the integral of the kernel `ker` from the
left end of its support to `x`. `A(x)` is zero for `x` below the support and, for `x` above
the support, equals the integral of `ker` over its whole support (one for normalized
kernels).

The kernel must be piecewise polynomial of degree at most 5 between unit-spaced breakpoints
starting at `-S/2` with `S = length(ker)`. This is the case of B-splines, generic cubic
splines, their derivatives, and the Keys, Catmull-Rom and Mitchell & Netravali kernels.
Other kernels (e.g. Lanczos) throw a `MethodError`.

The antiderivative of `ker` is exact, up to rounding errors, for the supported kernels.

""" antiderivative

const PiecewisePolynomialKernel{T} = Union{
    BSpline{<:Any,T}, BSplinePrime{<:Any,T},
    AbstractCubicSpline{T}, AbstractCubicSplinePrime{T}}

antiderivative(ker::PiecewisePolynomialKernel) = Antiderivative(ker)

Base.length(::Antiderivative{T,S}) where {T,S} = S
Base.eltype(::Type{<:Antiderivative{T}}) where {T} = T
Base.eltype(A::Antiderivative) = eltype(typeof(A))

# 3-point Gauss-Legendre rule on [-1,1], exact for polynomials up to degree 5.
function gauss_legendre_3(::Type{T}) where {T<:AbstractFloat}
    x = sqrt(T(3)/5)
    return (-x, zero(T), x), (T(5)/9, T(8)/9, T(5)/9)
end

# Integral of the kernel over `[a,b]` which must be inside a single polynomial piece.
function integrate_piece(ker::Kernel{T}, a::T, b::T) where {T}
    nodes, weights = gauss_legendre_3(T)
    mid = (a + b)/2
    half = (b - a)/2
    return half*sum(w*ker(mid + half*x) for (x, w) in zip(nodes, weights))
end

(A::Antiderivative{T,S})(x::Real) where {T,S} = A(T(x)::T)

function (A::Antiderivative{T,S})(x::T) where {T<:AbstractFloat,S}
    ker = A.kernel
    hw = T(S)/2
    isnan(x) && return x
    s = clamp(x, -hw, hw) + hw # position along the support
    c = min(floor(Int, s), S - 1) # index of the piece containing `x`
    result = zero(T)
    for k in 0:c-1
        result += integrate_piece(ker, k - hw, k + 1 - hw)
    end
    return result + integrate_piece(ker, c - hw, s - hw)
end

Base.show(io::IO, A::Antiderivative) = print(io, "antiderivative(", A.kernel, ")")
