#
# InterpolationKernels.jl --
#
# Kernel functions used for linear filtering, windowing or linear
# interpolation.
#
#------------------------------------------------------------------------------
#
# This file is part of the InterpolationKernels package licensed under the MIT
# "Expat" License.
#
# Copyright (C) 2016-2021, Éric Thiébaut.
#

module InterpolationKernels

export
    BSpline,
    BSplinePrime,
    CardinalCubicSpline,
    CardinalCubicSplinePrime,
    CatmullRomSpline,
    CatmullRomSplinePrime,
    CubicSpline,
    CubicSplinePrime,
    Kernel,
    LanczosKernel,
    LanczosKernelPrime,
    MitchellNetravaliSpline,
    MitchellNetravaliSplinePrime,
    iscardinal,
    isnormalized

import Base: convert, adjoint, values

# List and union of floating-point types which are specially handled.
const FLOATS = (:BigFloat, :Float16, :Float32, :Float64)
const Floats = @eval $(Expr(:curly, :Union, FLOATS...))

#------------------------------------------------------------------------------
# SUPPORTS

abstract type Bound end
struct Open <: Bound end
struct Closed <: Bound end

"""

`InterpolationKernels.Support{T,S,L,R}` is the abstract type for the support of
an interpolation kernel parameterized by the floating-point type `T` used for
computations, the integer size `S` of the support and the types `L` and `R` of
the left and right bound which can be `InterpolationKernels.Open` or
`InterpolationKernels.Closed`.

"""
abstract type Support{T<:AbstractFloat,S,L<:Bound,R<:Bound} end

Base.length(sup::Support) = length(typeof(sup))
Base.length(::Type{<:Support{T,S}}) where {T,S} = S

Base.eltype(sup::Support) = eltype(typeof(sup))
Base.eltype(::Type{<:Support{T,S}}) where {T,S} = T

"""
    InterpolationKernels.SymmetricSupport{T,S,L,R}() -> sup

yields an instance of a symmetric support parameterized by the floating-point
type `T`, the integer size `S` of the support and the types `L` and `R` of the
left and right bound which can be `InterpolationKernels.Open` or
`InterpolationKernels.Closed`.

Depending on `L` and `R`, the support is:

    (L,R) = (Open,Open) -------> sup = (-S/2,S/2)
    (L,R) = (Closed,Open) -----> sup = [-S/2,S/2)
    (L,R) = (Open,Closed) -----> sup = (-S/2,S/2]
    (L,R) = (Closed,Closed) ---> sup = [-S/2,S/2]

"""
struct SymmetricSupport{T<:AbstractFloat,S,L,R} <: Support{T,S,L,R} end

"""
    InterpolationKernels.LeftAnchoredSupport{T,S,L,R}(a)

yields an instance of a support with lower bound `a` and parameterized by the
floating-point type `T`, the integer size `S` of the support and the types `L`
and `R` of the left and right bound which can be `InterpolationKernels.Open` or
`InterpolationKernels.Closed`.

Depending on `L` and `R`, the support is:

    (L,R) = (Open,Open) -------> sup = (a,a+S)
    (L,R) = (Closed,Open) -----> sup = [a,a+S)
    (L,R) = (Open,Closed) -----> sup = (a,a+S]
    (L,R) = (Closed,Closed) ---> sup = [a,a+S]

"""
struct LeftAnchoredSupport{T,S,L,R} <: Support{T,S,L,R}
    a::T
end

"""
    InterpolationKernels.RightAnchoredSupport{T,S,L,R}(b)

yields an instance of a support with upper bound `b` and parameterized by the
floating-point type `T`, the integer size `S` of the support and the types `L`
and `R` of the left and right bound which can be `InterpolationKernels.Open` or
`InterpolationKernels.Closed`.

Depending on `L` and `R`, the support is:

    (L,R) = (Open,Open) -------> sup = (b-S,b)
    (L,R) = (Closed,Open) -----> sup = [b-S,b)
    (L,R) = (Open,Closed) -----> sup = (b-S,b]
    (L,R) = (Closed,Closed) ---> sup = [b-S,b]

"""
struct RightAnchoredSupport{T,S,L,R} <: Support{T,S,L,R}
    b::T
end

"""
    InterpolationKernels.supremum(sup) -> b

yields the least upper bound `b` of the kernel support `sup`.

"""
supremum(::SymmetricSupport{T,S}) where {T,S} = T(S)/2
supremum(sup::LeftAnchoredSupport{T,S}) where {T,S} = infimum(sup) + T(S)
supremum(sup::RightAnchoredSupport) = sup.b

"""
    InterpolationKernels.infimum(sup) -> a

yields the greatest lower bound `a` of the kernel support `sup`.

"""
infimum(::SymmetricSupport{T,S}) where {T,S} = T(-S)/2
infimum(sup::LeftAnchoredSupport) = sup.a
infimum(sup::RightAnchoredSupport{T,S}) where {T,S} = supremum(sup) - T(S)

#------------------------------------------------------------------------------
# INTERPOLATION KERNELS

"""
# Interpolation Kernels

An interpolation kernel `Kernel{T,S}` is parametrized by the floating-point
type `T` of its coefficients and by the integer size `S` of its support.  For
efficiency reasons, only kernels with (small) finite size supports are
implemented.

A kernel `ker` is a callable object which may be used as a function with a real
argument:

    ker(x::Real)

yields kernel value at offset `x`.  Whatever the type of `x`, `ker(x)` is
always of type `T = eltype(ker)` the floating-point type associated with `ker`.
All kernel supports are symmetric; that is `ker(x)` is zero if `abs(x) > S/2`
with `S = length(ker)` the size of the kernel support.


## Kernel floating-point type conversion

Called as a function with a real argument, a given kernel returns a value of
its associated floating-point type.  This has been chosen to have fast
interpolation methods.  Converting a kernel `ker` to use floating-point type
`T` is simply done by one of:

    T(ker)
    Kernel{T}(ker)
    convert(Kernel{T}, ker)

Beware that changing the floating-point type may lead to a loss of precision if
the kernel has numerical parameters.


## Common methods

A few common methods are specialized for any interpolation kernel `ker`:

    eltype(ker) -> T

yields the floating-point type for calculations,

    length(ker) -> S

yield the size the support of `ker` which is also the number of neighbors
involved in an interpolation by this kernel,

    values(ker)

yields a tuple of the parameters of `ker` such that an identical instance can
be built by:

    typeof(ker)(values(ker)...)

finally:

    compute_offset_and_weights(ker, x) -> off, (w1, w2, ..., wS)

yields the offset `off` and an `S`-tuple of interpolation weights to
interpolate an array at coordinate `x` (in fractional index units).

"""
abstract type Kernel{T<:AbstractFloat,S} <: Function end

Base.eltype(ker::Kernel) = eltype(typeof(ker))
Base.eltype(::Type{<:Kernel{T,S}}) where {T,S} = T

Base.length(ker::Kernel) = length(typeof(ker))
Base.length(::Type{<:Kernel{T,S}}) where {T,S} = S

"""
    isnormalized(ker)

yields whether the kernel `ker` has the partition of unity property.  That is,
the sum of the values computed by the kernel `ker` on a unit spaced grid is
equal to one.

""" isnormalized

"""
    iscardinal(ker)

yields whether the kernel `ker` is zero for non-zero integer arguments.
Cardinal kernels are directly suitable for interpolation.

""" iscardinal

"""
    InterpolationKernels.compute_offset_and_weights(ker, x) -> off, wgt

yields the index offset `off` and the weights `wgt` to interpolate with kernel
`ker` at position `x` in fractional index units.  The offset is a scalar and
the weights are an `n`-tuple with `n = length(ker)` the size of the support of
the kernel, all returned values have the same floating point type `eltype(ker)`
as the kernel.

Not taking into account boundary conditions, interpolating a vector `A` at
position `x` would then write:

    off, wgt = InterpolationKernels.compute_offset_and_weights(ker, x)
    k = Int(off) # here boundary conditions should be imposed
    result = wgt[1]*A[k+1] + ... + wgt[n]*A[k+n]

Note that 1-based indexing is assumed by `compute_offset_and_weights` to
interpret the position `x` and compute the offset `off`.  If this is not the
case, the code should be:

    j1 = first(axes(A,1)) # first index in A
    off, wgt = InterpolationKernels.compute_offset_and_weights(ker, x - (j1 - 1))
    k = Int(off) + (j1 - 1) # here boundary conditions should be imposed
    result = wgt[1]*A[k+1] + ... + wgt[n]*A[k+n]

where expression `x - (j1 - 1)` is assuming that the position `x` is in
fractional index for `A`, that is `x = j1` at the first entry of `A`.

!!! note
    For fast computations, this method should be specialized for specific
    kernel types.  For kernels with symmetric support, the method
    [`InterpolationKernels.compute_weights`](@ref) is called by
    `compute_offset_and_weights` to calculate the interpolation weights; for
    such kernels it is sufficient to specialize `compute_weights` instead of
    `compute_offset_and_weights`.

""" compute_offset_and_weights

@inline compute_offset_and_weights(ker::Kernel{T}, x::T) where {T} =
    compute_offset_and_weights(support(ker), ker, x)

# The following version is specialized for symmetric supports and rely on an
# optimized version of `compute_weights`.
@generated function compute_offset_and_weights(sup::SymmetricSupport{T,S},
                                               ker::Kernel{T,S},
                                               x::T) where {T,S}
    if isodd(S)
        quote
            $(Expr(:meta, :inline))
            round_x = round(x)
            wgts = compute_weights(ker, x - round_x)
            off = round_x - $((S + 1) >> 1)
            return off, wgts
        end
    else
        quote
            $(Expr(:meta, :inline))
            floor_x = floor(x)
            wgts = compute_weights(ker, x - floor_x)
            off = floor_x - $(S >> 1)
            return off, wgts
        end
    end
end

# This version is to call the generic version by default.
@inline function compute_offset_and_weights(sup::Support{T,S},
                                            ker::Kernel{T}, x::T) where {T,S}
    generic_compute_offset_and_weights(sup, ker, x)
end

# The generic version is only called if no optimized version exists.  It can
# also be used to check code implementing an optimized version.
@inline generic_compute_offset_and_weights(ker::Kernel{T}, x::T) where {T} =
    generic_compute_offset_and_weights(support(ker), ker, x)

@generated function generic_compute_offset_and_weights(sup::Support{T,S},
                                                       ker::Kernel{T,S},
                                                       x::T) where {T,S}
    wgt = ntuple(j -> Symbol("w_", j), Val(S))
    exprs = [:($(wgt[j]) = ker(u - $j)) for j in 1:S]
    quote
        $(Expr(:meta, :inline))
        off = offset(sup, x)
        u = x - off
        $(exprs...)
        return off, $(Expr(:tuple, wgt...))
    end
end

"""
    InterpolationKernels.offset(sup, x) -> off

yields the coordinate offset for computing the interpolation weights
at coordinate `x` for a kernel whose support is `sup`:

    off = floor(x - b)       if `sup` is right-open
          ceil(x - b) - 1    if `sup` is left-open

with `b = supremum(sup)` the upper bound of the support.

""" offset

offset(sup::Support{T,S,<:Bound,Open}, x::T) where {T,S} =
    floor(x - supremum(sup))

offset(sup::Support{T,S,Open,Closed}, x::T) where {T,S} =
    ceil(x - (supremum(sup) - 1))

"""
    InterpolationKernels.compute_weights(ker, t) -> wgt

computes the interpolation weights returned by
[`InterpolationKernels.compute_offset_and_weights`](@ref) for kernel `ker` with
symmetric support.  Assuming interpolation is performed at at position `x`,
argument `t` is given by:

     t = x - floor(x)     if length(ker) is even, hence t ∈ [0,1)
     t = x - round(x)     if length(ker) is odd,  hence t ∈ [-1/2,+1/2]

The returned weights are then:

     wgt = ntuple(i -> ker(t + k - i), length(ker))

where `k = (length(ker) + 1) >> 1` (i.e., integer division of `length(ker)+1`
by 2).  These conventions have been adopted so that, by specializing the
`compute_weights` method, computing the `length(ker)` weights at the same time
may be done in much fewer operations than calling `ker` as a function for each
weight.

""" compute_weights

# Call generic version by default.  The generic version has a different name so
# that it can be used to check optimized versions.  Must be in-lined.
@inline compute_weights(ker::Kernel{T}, t::T) where {T} =
    generic_compute_weights(ker, t)

@generated function generic_compute_weights(ker::Kernel{T,S}, t::T) where {T,S}
    k = ((S + 1) >> 1)
    exprs = [(j < 0 ? :(ker(t + $(-j))) :
              j > 0 ? :(ker(t - $j)) : :(ker(t))) for j in 1-k:S-k]
    return quote
        $(Expr(:meta, :inline))
        return $(Expr(:tuple, exprs...))
    end
end

"""
    values(ker::InterpolationKernels.Kernel)

yields a tuple of the parameters of the interpolation kernel `ker` such that an
identical instance can be built by:

    typeof(ker)(values(ker)...)

"""
values(::Kernel) = ()

"""
    InterpolationKernels.support(ker) -> sup

yields the support of the interpolation kernel `ker`.

""" support

#------------------------------------------------------------------------------
# B-SPLINES

"""
    BSpline{S,T}()

yields a B-spline (short for *basis spline*) of order `S` that is a piecewise
polynomial function of degree `S - 1` on a support of length `S`.  The
parameter `T` is the floating-point type for computations, `T = Float64` is
assuled if this parameter is not specified.

Fr now, not all B-spline are implemented in `InterpolationKernels`, `S` must be: `1`
(for a **rectangular** B-spline), `2` (for a **linear** B-spline), `3` (for a
**quadratic** B-spline), or `4` (for a **cubic** B-spline).

If `ker` is a B-spline, then `ker'` is its derivative which can also be
directly constructed by calling [`BSplinePrime`](@ref).

!!! warning
    The derivative of B-spline of order `S ≤ 2` is not defined everywhere.  It
    is allowed to take their derivative but it (arbitrarily) yields zero where
    not defined.  Returning `NaN` would have been more correct but it has been
    considered that it would do more harm than good in practice.

""" BSpline

struct BSpline{S,T} <: Kernel{T,S} end

"""
    BSplinePrime{S,T}()

yields the derivative of a B-spline of order `S` for floating-point `T`.

See the caveats in [`BSpline`](@ref) about taking the derivative of B-splines of
order `S ≤ 2`.

""" BSplinePrime

struct BSplinePrime{S,T} <: Kernel{T,S} end

# Outer Constructors.
BSpline{S}() where {S} = BSpline{S,Float64}()
BSplinePrime{S}() where {S} = BSplinePrime{S,Float64}()

# Support is right-open for 1st order splines and open for all others.
support(::BSpline{1,T}) where {T} = SymmetricSupport{T,1,Closed,Open}()
support(::BSplinePrime{1,T}) where {T} = SymmetricSupport{T,1,Closed,Open}()
support(::BSpline{S,T}) where {S,T} = SymmetricSupport{T,S,Open,Open}()
support(::BSplinePrime{S,T}) where {S,T} = SymmetricSupport{T,S,Open,Open}()

# Only the 2 first B-splines are interpolating.
iscardinal(::Union{K,Type{K}}) where {K<:BSpline{1}} = true
iscardinal(::Union{K,Type{K}}) where {K<:BSpline{2}} = true
iscardinal(::Union{K,Type{K}}) where {K<:BSpline} = false
iscardinal(::Union{K,Type{K}}) where {K<:BSplinePrime} = false

# All B-splines are normalized, not their derivative.
isnormalized(::Union{K,Type{K}}) where {K<:BSpline} = true
isnormalized(::Union{K,Type{K}}) where {K<:BSplinePrime} = false

#
# Rectangular B-spline
# --------------------
#
@inline (::BSpline{1,T})(x::T) where {T} =
    ifelse((x < frac(T,-1,2))|(x ≥ frac(T,1,2)), zero(T), one(T))
compute_offset_and_weights(::BSpline{1,T}, x::T) where {T} =
    (round(x) - 1, (one(T),))
compute_weights(::BSpline{1,T}, t::T) where {T} = (one(T),)

(::BSplinePrime{1,T})(x::T) where {T} = zero(T)
compute_offset_and_weights(::BSplinePrime{1,T}, x::T) where {T} =
    (round(x) - 1, (zero(T),))
compute_weights(::BSplinePrime{1,T}, t::T) where {T} = (zero(T),)

#
# Linear B-splines
# ----------------
#
@inline function (::BSpline{2,T})(x::T) where {T}
    abs_x = abs(x)
    return ifelse(abs_x < 1, 1 - abs_x, zero(T))
end

compute_weights(ker::BSpline{2,T}, t::T) where {T} = (1 - t, t)

# The linear B-spline is not C¹ continuous, its left and right derivatives do
# not match at x ∈ (-1,0,1).  Setting h'(x) = NaN for x ∈ (-1,0,1) is correct
# but not useful in practice.  Taking the mean of the left and right
# derivatives yiedls h'(±1) = ∓1/2 and h'(0) = 0 but, then the support would no
# longer be semi-open interval of size 2.  For now, h'(x) = 0 for x ∈ (-1,0,1).
@inline function (::BSplinePrime{2,T})(x::T) where {T}
    if (-1 < x)&(x < 0)
        return one(T)
    elseif (0 < x)&(x < 1)
        return -one(T)
    else
        return zero(T)
    end
end

compute_weights(ker::BSplinePrime{2,T}, t::T) where {T} =
    ifelse(t > 0, (-one(T), one(T)), (zero(T), zero(T)))

#
# Quadratic B-spline
# ------------------
#
@inline function (::BSpline{3,T})(x::T) where {T}
    abs_x = abs(x)
    if abs_x ≥ frac(T,3,2)
        return zero(T)
    elseif abs_x ≤ frac(T,1,2)
        return frac(T,3,4) - abs_x*abs_x
    else
        return square(abs_x - frac(T,3,2))*frac(T,1,2)
    end
end

@inline function compute_weights(ker::BSpline{3,T}, t::T) where {T}
    #
    # Given `t = x - round(x)`, the weights are:
    #
    #     w1 = (1/8)*(1 - 2*t)^2 = (1/2)*(1/sqrt(2) - t)^2
    #     w2 = (3/4) - t^2
    #     w3 = (1/8)*(1 + 2*t)^2 = (1/2)*(1/sqrt(2) + t)^2
    #
    # which can be computed in 9 operations using the following factorized
    # expressions:
    #
    #     w1 = (1/2)*(1/sqrt(2) - t)^2
    #     w2 = (3/4) - t^2
    #     w3 = (1/2)*(1/sqrt(2) + t)^2
    #
    # or in 6 operations:
    #
    ht = t/2
    w = frac(T,1,8) + ht*t
    w1 = w - ht
    w2 = frac(T,3,4) - t*t
    w3 = w + ht
    return (w1, w2, w3)
end

@inline function (::BSplinePrime{3,T})(x::T) where {T}
    if (x ≤ frac(T,-3,2))|(x ≥ frac(T,3,2))
        return zero(T)
    elseif x < frac(T,-1,2)
        return x + frac(T,3,2)
    elseif x ≤ frac(T,1,2)
        return -2x
    else
        return x - frac(T,3,2)
    end
end

@inline function compute_weights(ker::BSplinePrime{3,T}, t::T) where {T}
    # 3 operations
    h = frac(T,1,2)
    return (t - h, -2t, t + h)
end

#
# Cubic B-spline
# --------------
#
@inline function (::BSpline{4,T})(x::T) where {T}
    abs_x = abs(x)
    if abs_x ≥ 2
        return zero(T)
    elseif abs_x ≥ 1
        return cube(2 - abs_x)*frac(T,1,6)
    else
        return (frac(T,1,2)*abs_x - 1)*abs_x*abs_x + frac(T,2,3)
    end
end

@inline function compute_weights(ker::BSpline{4,T}, t::T) where {T}
    #
    # With `t = x - floor(x)`, the weights are given by:
    #
    #     w1 = 1/6 - t/2 + t^2/2 - t^3/6
    #        = 1/6 + (t^2 - t)/2 - t^3/6
    #        = (1 - t)^3/6
    #     w2 = 2/3 - t^2 + t^3/2
    #        = 2/3 + (t/2 - 1)*t^2
    #     w3 = 1/6 + t/2 + t^2/2 - t^3/2
    #        = 1/6 + (t + t^2 - t^3)/2
    #        = 1/6 - ((t - 1)*t - 1)*t/2
    #        = 4/6 - (1 - t)^2*(t + 1)/2
    #     w4 = t^3/6
    #
    # Horner's scheme takes 6 operations per cubic polynomial, 24 operations
    # for the 4 weights.  Precomputing the powers of t, t^2 and t^3, takes 2
    # operations, then 6 operations per cubic polynomial are needed.
    #
    # Using factorizations, I manage to only use 15 operations:
    #
    h = frac(T,1,2)
    p = frac(T,2,3)
    q = frac(T,1,6)
    u = 1 - t
    u2 = u*u
    t2 = t*t
    w1 = q*u2*u
    w2 = p + (h*t - 1)*t2
    w3 = p - (h*u2)*(t + 1)
    w4 = q*t2*t
    return (w1, w2, w3, w4)
end

@inline function (::BSplinePrime{4,T})(x::T) where {T}
    if (x ≤ -2)|(x ≥ 2)
        return zero(T)
    elseif x < -1
        return frac(T,1,2)*square(x + 2)
    elseif x <  0
        return frac(T,-3,2)*x*(x + frac(T,4,3))
    elseif x <  1
        return frac(T,+3,2)*x*(x - frac(T,4,3))
    else
        return frac(T,-1,2)*square(x - 2)
    end
end

@inline function compute_weights(ker::BSplinePrime{4,T}, t::T) where {T}
    w1 = frac(T,-1,2)*(t - 1)^2
    w2 = frac(T, 3,2)*(t - frac(T,4,3))*t
    w3 = frac(T, 1,2) + (1 - frac(T,3,2)*t)*t
    w4 = frac(T, 1,2)*t^2
    return (w1, w2, w3, w4)
end

#------------------------------------------------------------------------------
# GENERIC CUBIC SPLINE

# Abstract types to implement common traits of cubic splines.
abstract type AbstractCubicSpline{T} <: Kernel{T,4} end
abstract type AbstractCubicSplinePrime{T} <: Kernel{T,4} end

support(::AbstractCubicSpline{T}) where {T} =
    SymmetricSupport{T,4,Open,Open}()
support(::AbstractCubicSplinePrime{T}) where {T} =
    SymmetricSupport{T,4,Open,Open}()

"""
    CubicSpline{T}(a, b) -> ker

yields an instance of a generic cubic spline for floating-point type `T` and
parameters `a = ker'(1)` and `b = ker(1)` the slope and the value of the
function `ker(x)` at `x = 1`.

A cubic spline kernel is at least C¹ continuous, the expression `ker'` yields a
kernel instance implementing the 1st derivative of the generic cubic spline
`ker` (see [`CubicSplinePrime`](@ref) to directly build a derivative).

Depending on the values of the parameters `a` and `b`, more specific cubic
spline kernels can be emulated:

* `CubicSpline{T}(-1/2,1/6)` yields a cubic B-spline as built by
  `BSpline{4,T}()`.

* `CubicSpline{T}(a,0)` yields a cardinal cubic spline as built by
  `CardinalCubicSpline{T}(a)`.

* `CubicSpline{T}(-1/2,0)` yields a Catmull-Rom kernel as built by
  `CatmullRomSpline{T}()`.

* `CubicSpline{T}(-b/2-c,b/6)` yields Mitchell & Netravali cubic spline as
  built by `MitchellNetravaliSpline{T}(b,c)`.

Instances of `CubicSpline` are very well optimized and, in practice, they may
be as fast or even faster than these more specialized counterparts.

""" CubicSpline

"""
    CubicSplinePrime{T}(a, b)

yields a kernel instance that is the 1st derivative of the generic cubic spline
of parameters `a` and `b` (see [`CubicSpline`](@ref)) for floating-point
type `T` (`Float64` by default).

""" CubicSplinePrime

struct CubicSpline{T} <: AbstractCubicSpline{T}
    a::T
    b::T
    c0::T
    c1::T
    c2::T
    c3::T
    c4::T
    c5::T
    c6::T
    function CubicSpline{T}(_a::Real, _b::Real) where {T}
        # NOTE: Paramaters must be signed in the expressions of the kernel
        # constants.  Since there is no loss of precision in computing these
        # constants with floating-point arithmetic even though parameters are
        # integers, we convert the parameters to type T.
        a, b = T(_a), T(_b)
        new{T}(a, b,
               1 - 2b,     # c0
               9b - a - 3, # c1
               2 + a - 6b, # c2
               -a - b,     # c3
               a + 2b,     # c4
               3b - 1,     # c5
               a + 3b)     # c6
    end
end

@inline function (ker::CubicSpline{T})(x::T) where {T}
    abs_x = abs(x)
    if abs_x ≤ 1
        return (ker.c2*abs_x + ker.c1)*abs_x^2 + ker.c0
    elseif abs_x < 2
        return (ker.c3 + ker.c4*abs_x)*(abs_x - 2)^2
    else
        return zero(T)
    end
end

@inline function compute_weights(ker::CubicSpline{T}, t::T) where {T}
    #
    # With `t = x - floor(x)`, the weights are given by:
    #
    #     w1 = (1 - t)^2 (b + (a + 2b) t)
    #        = u^2 (b + c4*t)
    #     w2 = (1 - 2b) + t^2 ((-3 - a + 9b) + (2 + a - 6b) t)
    #        = c0 + t^2 (c1 + c2*t)
    #     w3 = (1 - 2b) + (1 - t)^2 ((3b - 1) - (2 + a - 6b) t)
    #        = c0 + u^2 (c5 - c2*t)
    #     w4 = t^2 ((a + 3b) - (a + 2b) t)
    #        = t^2 (c6 - c4 t)
    #
    # where:
    #
    #     u = 1 - t
    #     c5 = 3b - 1
    #     c6 = a + 3b
    #
    # can be done in 15 operations:
    #
    u = 1 - t
    c4t = ker.c4*t
    c2t = ker.c2*t
    u2 = u*u
    t2 = t*t
    w1 = (ker.b  + c4t)*u2
    w2 = (ker.c1 + c2t)*t2 + ker.c0
    w3 = (ker.c5 - c2t)*u2 + ker.c0
    w4 = (ker.c6 - c4t)*t2
    return (w1, w2, w3, w4)
end

struct CubicSplinePrime{T} <: AbstractCubicSplinePrime{T}
    a::T
    b::T
    c1::T
    c2::T
    c3::T
    c4::T
    c5::T
    function CubicSplinePrime{T}(_a::Real, _b::Real) where {T}
        a, b = T(_a), T(_b)
        new{T}(a, b,
               -6 - 2a + 18b, # c1
               6 + 3a - 18b,  # c2
               -4a - 6b,      # c3
               3a + 6b,       # c4
               2a + 6b)       # c5
    end
end

@inline function (ker::CubicSplinePrime{T})(x::T) where {T}
    sign_x, abs_x = signabs(x)
    if abs_x ≤ 1
        return (ker.c1 + ker.c2*abs_x)*x
    elseif abs_x < 2
        return sign_x*(abs_x - 2)*(ker.c3 + ker.c4*abs_x)
    else
        return zero(T)
    end
end

@inline function compute_weights(ker::CubicSplinePrime{T}, t::T) where {T}
    #
    # With `t = x - floor(x)`, the weights are given by:
    #
    #     w1 = (a - (6b + 3a) t) (1 - t)
    #        = (a - c4 t) (1 - t)
    #     w2 = ((-6 - 2a + 18b) + (6 + 3a - 18b) t) t
    #        = (c1 + c2 t) t
    #     w3 = ((6 + 3a - 18b) t - a) (1 - t)
    #        = (c2 t - a) (1 - t)
    #     w4 = ((2a + 6b) - (3a + 6b) t) t
    #        = (c5 - c4 t) t
    #
    # can be computed in 11 operations in total:
    #
    u = 1 - t
    c4t = ker.c4*t
    c2t = ker.c2*t
    w1 = (ker.a - c4t)*u
    w2 = (ker.c1 + c2t)*t
    w3 = (c2t - ker.a)*u
    w4 = (ker.c5 - c4t)*t
    return (w1, w2, w3, w4)
end

# Outer constructors, traits and standard methods.
for K in (:CubicSpline, :CubicSplinePrime)
    @eval begin
        $K(a::Real, b::Real) = $K(promote(a, b)...)
        $K(a::Integer, b::Integer) = $K{Float64}(a, b)
        $K(a::T, b::T) where {T<:AbstractFloat} = $K{T}(a, b)
        values(ker::$K) = (ker.a, ker.b)
        iscardinal(::Union{K,Type{K}}) where {K<:$K} = false
        isnormalized(::Union{K,Type{K}}) where {K<:$K} =
            $(K === :CubicSpline)
    end
end

#------------------------------------------------------------------------------
# KEYS' CARDINAL CUBIC SPLINE

# Abstract types to implement common traits of cardinal cubic splines.
abstract type AbstractCardinalCubicSpline{T} <: AbstractCubicSpline{T} end
abstract type AbstractCardinalCubicSplinePrime{T} <: AbstractCubicSplinePrime{T} end

"""
    CardinalCubicSpline{T}(a)

yields an instance of the Keys family of cardinal cubic splines for
floating-point type `T` and parameter `a = ker'(1)` the slope of the function
`ker(x)` at `x = 1`.

These kernels are C¹ continuous piecewise normalized cardinal cubic spline
which depend on one parameter `a` and defined by:

    ker(x) = 1 + ((2 + a)*|x| - (3 + a))*x^2    if |x| ≤ 1
             a*(|x| - 1)*(|x| - 2)^2            if 1 ≤ |x| ≤ 2
             0                                  if |x| ≥ 2

The expression `ker'` yields a kernel instance which is the 1st derivative of
the Keys kernel `ker` (also see the constructor
[`CardinalCubicSplinePrime`](@ref)).

Reference:

* Keys, Robert, G., "Cubic Convolution Interpolation for Digital Image
  Processing", IEEE Trans. Acoustics, Speech, and Signal Processing,
  Vol. ASSP-29, No. 6, December 1981, pp. 1153-1160.

""" CardinalCubicSpline

struct CardinalCubicSpline{T} <: AbstractCardinalCubicSpline{T}
    a::T   # a
    ap2::T # a + 2
    ap3::T # a + 3
    am2::T # a - 2
    function CardinalCubicSpline{T}(_a::Real) where {T}
        a = T(_a)
        new{T}(a, a + 2, a + 3, a - 2)
    end
end

@inline function (ker::CardinalCubicSpline{T})(x::T) where {T}
    abs_x = abs(x)
    if abs_x ≥ 2
        return zero(T)
    elseif abs_x ≤ 1
        return 1 + (ker.ap2*abs_x - ker.ap3)*abs_x^2
    else
        return ker.a*(abs_x - 1)*(abs_x - 2)^2
    end
end

@inline function compute_weights(ker::CardinalCubicSpline{T}, t::T) where {T}
    #
    # Given `t = x - floor(x)`, compute:
    #
    #     w1 = a*(1 - t)^2*t
    #     w2 = (1 + t - (2 + a) t^2)*(1 - t)
    #     w3 = -t (a - (2a + 3)*t + (a + 2)*t^2)
    #     w4 = a*(1 - t)*t^2
    #
    # Noting that, with `u = 1 - t` and `v = a*t*u`(1 - t)`, the following
    # hold:
    #
    #     w1 = v*u
    #     w4 = v*t
    #     w1 + w4 = v
    #     w2 + w3 = 1 - v
    #     w2 = 1 - (1 + (2 + a)*u)*t^2
    #
    # can be done in in 11 operations:
    #
    u = 1 - t
    v = ker.a*t*u
    w = (1 + ker.ap2*u)*(t*t) # ap2 = a + 2
    w1 = v*u
    w2 = 1 - w
    w3 = w - v
    w4 = v*t
    return (w1, w2, w3, w4)
end

"""
    CardinalCubicSplinePrime{T}(a)

yields a kernel instance that is the 1st derivative of the Keys cardinal cubic
spline (see [`CardinalCubicSpline`](@ref)) for floating-point type `T` and
parameter `a`.  This derivative is given by:

    ker′(x) = (3(a + 2)*|x| - 2(a + 3))*x           if |x| ≤ 1
              (3a)*(|x| - 2)*(|x| - 4/3)*sign(x)    if 1 ≤ |x| ≤ 2
              0                                     if |x| ≥ 2

""" CardinalCubicSplinePrime

struct CardinalCubicSplinePrime{T} <: AbstractCardinalCubicSplinePrime{T}
    a ::T
    c1::T # 3(a + 2)
    c2::T # 2(a + 3)
    c3::T # 3a
    CardinalCubicSplinePrime{T}(a::Real) where {T} =
        new{T}(a, 3(a + 2), 2(a + 3), 3a)
end

@inline function (ker::CardinalCubicSplinePrime{T})(x::T) where {T}
    sign_x, abs_x = signabs(x)
    if abs_x ≥ 2
        return zero(T)
    elseif abs_x ≤ 1
        return (ker.c1*abs_x - ker.c2)*x
    else
        return (sign_x*ker.c3)*(abs_x - 2)*(abs_x - frac(T,4,3))
    end
end

# FIXME: optimize compute_weights for CardinalCubicSplinePrime

# Traits and outer constructors.
for K in (:CardinalCubicSpline, :CardinalCubicSplinePrime)
    @eval begin
        $K(a::Real) = $K{Float64}(a)
        $K(a::T) where {T<:AbstractFloat} = $K{T}(a)
        values(ker::$K) = (ker.a,)
        iscardinal(::Union{K,Type{K}}) where {K<:$K} =
            $(K === :CardinalCubicSpline)
        isnormalized(::Union{K,Type{K}}) where {K<:$K} =
            $(K === :CardinalCubicSpline)
    end
end

#------------------------------------------------------------------------------
# CATMULL-ROM INTERPOLATION KERNEL

"""
    CatmullRomSpline{T}()

yields an instance of the Catmull-Rom interpolation kernel for floating-point
type `T` which is assumed to be `Float64` if omitted.

Catmull-Rom interpolation kernel is a cardinal piecewise cubic spline defined
by:

    ker(x) = ((3/2)*|x| - (5/2))*x^2 + 1             if |x| ≤ 1
             (((5/2) - (1/2)*|x|)*|x| - 4)*|x| + 2   if 1 ≤ |x| ≤ 2
             0                                       if |x| ≥ 2

The expression `ker'` yields a kernel instance which is the 1st derivative of
the Catmull-Rom interpolation kernel `ker` (also see the constructor
[`CatmullRomSplinePrime`](@ref)).

""" CatmullRomSpline

struct CatmullRomSpline{T} <: AbstractCardinalCubicSpline{T} end

@inline function (::CatmullRomSpline{T})(x::T) where {T}
    abs_x = abs(x)
    if abs_x ≥ 2
        return zero(T)
    elseif abs_x ≤ 1
        # (((3/2)*x - (5/2))*x^2 + 1)
        return ((T(3)/2)*abs_x - T(5)/2)*abs_x*abs_x + 1
    else
        # (((5/2) - (1/2)*x)*x - 4)*x + 2
        # = (2 - x)^2 (1 - x)/2
        return (2 - abs_x)^2*(1 - abs_x)/2
    end
end

@inline function compute_weights(::CatmullRomSpline{T}, t::T) where {T}
    # 10 operations:
    u = 1 - t
    v = frac(T,-1,2)*t*u
    w1 = v*u
    w4 = v*t
    w = w4 - w1
    w2 = u - w1 + w
    w3 = t - w4 - w
    return (w1, w2, w3, w4)
end

"""
    CatmullRomSplinePrime{T}()

yields a kernel instance that is the 1st derivative of the Catmull-Rom
interpolation kernel (see [`CatmullRomSpline`](@ref)) for floating-point type
`T` which is assumed to be `Float64` if omitted.

The 1st derivative of the Catmull-Rom interpolation kernel is given by:

    ker′(x) = ((9/2)*|x| - 5)*x                      if a = |x| ≤ 1
              (5 - (3/2)*|x|)*x - 4*sign(x)          if 1 ≤ |x| ≤ 2
              0                                      if |x| ≥ 2

""" CatmullRomSplinePrime

struct CatmullRomSplinePrime{T} <: AbstractCardinalCubicSplinePrime{T} end

@inline function (::CatmullRomSplinePrime{T})(x::T) where {T}
    if x < 0
        if x ≤ -2
            return zero(T)
        elseif x ≥ -1
            return (frac(T,-9,2)*x - 5)*x
        else
            return (5 + frac(T,3,2)*x)*x + 4
        end
    else
        if x ≥ 2
            return zero(T)
        elseif x ≤ 1
            return (frac(T,9,2)*x - 5)*x
        else
            return (5 - frac(T,3,2)*x)*x - 4
        end
    end
end

@inline function compute_weights(::CatmullRomSplinePrime{T}, t::T) where {T}
    #
    # Weights are given by (18 operations):
    #
    #   w1 = -(3/2)*t^2 + 2*t - (1/2);
    #   w2 =  (9/2)*t^2 - 5*t;
    #   w3 = -(9/2)*t^2 + 4*t + (1/2);
    #   w4 =  (3/2)*t^2 - t;
    #
    #   w4 = ((3/2)*t - 1)*t;
    #   w1 = t - w4 - (1/2);
    #   w2 = 3*w4 - 2*t;
    #   w3 = t - 3*w4 + (1/2);
    #
    # can be computed in 11 operations:
    #
    w4 = frac(T,3,2)*t*t - t;
    w1 = t - w4 - frac(T,1,2);
    w2 = T(3)*w4 - 2*t;
    w3 = t - T(3)*w4 + frac(T,1,2);
    return (w1, w2, w3, w4)
end

# Traits and outer constructors.
for K in (:CatmullRomSpline, :CatmullRomSplinePrime)
    @eval begin
        $K() = $K{Float64}()
        iscardinal(::Union{K,Type{K}}) where {K<:$K} =
            $(K === :CatmullRomSpline)
        isnormalized(::Union{K,Type{K}}) where {K<:$K} =
            $(K === :CatmullRomSpline)
    end
end

#------------------------------------------------------------------------------
# MITCHELL & NETRAVALI KERNELS

"""
    MitchellNetravaliSpline{T}(b=1/3, c=1/3)

yields an instance of the Mitchell & Netravali family of kernels for
floating-point type `T` and parameters `(b,c)`.

These kernels are cubic splines which depends on 2 parameters, `b` and `c`.
Whatever the values of `(b,c)`, Mitchell & Netravali kernels are normalized,
even and C¹ continuous functions (these kernels and their first derivatives are
continuous).

Taking `b = 0` yields the family of cardinal cubic splines (see
[`CardinalCubicSpline`](@ref)) and is a sufficient and necessary condition to
have Mitchell & Netravali kernels be cardinal functions.

Using the constraint: `b + 2c = 1` yields a cubic filter with, at least,
quadratic order approximation.

Some specific values of `(b,c)` yield other well known kernels:

    (b,c) = (1,0)      --> cubic B-spline
    (b,c) = (0,-a)     --> Keys's cardinal cubic spline CardinalCubicSpline(a)
    (b,c) = (0,1/2)    --> Catmull-Rom kernel CatmullRomSpline()
    (b,c) = (b,0)      --> Duff's tensioned B-spline
    (b,c) = (6β,-α-3β) --> generic cubic spline CubicSpline(α,β)
    (b,c) = (1/3,1/3)  --> recommended by Mitchell-Netravali

The expression `ker'` yields a kernel instance which is the 1st derivative of
the Mitchell & Netravali kernel `ker` (also see the constructor
[`MitchellNetravaliSplinePrime`](@ref)).

Mitchell & Netravali family of kernels are currently instances of
[`CubicSpline`](@ref).

Reference:

* [Mitchell & Netravali, "*Reconstruction Filters in Computer Graphics*", in
  Computer Graphics, Vol. 22, Num. 4
  (1988)](http://www.cs.utexas.edu/users/fussell/courses/cs384g/lectures/mitchell/Mitchell.pdf).

""" MitchellNetravaliSpline

abstract type MitchellNetravaliSpline{T} <: AbstractCubicSpline{T} end

"""
    MitchellNetravaliSplinePrime([T=Float64,] [b=1/3, c=1/3,] B=Flat)

yields a kernel instance that is the 1st derivative of the Mitchell & Netravali
kernel (see [`MitchellNetravaliSpline`](@ref)) for floating-point type `T`,
parameters `b` and `c` and boundary conditions `B`.

""" MitchellNetravaliSplinePrime

abstract type MitchellNetravaliSplinePrime{T} <: AbstractCubicSplinePrime{T} end

for K in (:MitchellNetravaliSpline, :MitchellNetravaliSplinePrime)
    @eval begin
        $K() = $K{Float64}()
        $K{T}() where {T<:AbstractFloat} = $K{T}(one(T)/3, one(T)/3)
        $K(b::Real, c::Real) = $K{floating_point_type(b, c)}(b, c)
        function $K{T}(_b::Real, _c::Real) where {T<:AbstractFloat}
            b, c = T(_b), T(_c)
            return CubicSpline{T}(-(b/2 + c), b/6)
        end
    end
end

floating_point_type(x::Real) = floating_point_type(typeof(x))
floating_point_type(args::Real...) =
    floating_point_type(map(typeof, args)...)

floating_point_type(T::Type{<:AbstractFloat}) = T
floating_point_type(T::Type{<:Real}) = Float64
floating_point_type(types::Type{<:Real}...) =
    floating_point_type(promote_type(types...))

#------------------------------------------------------------------------------
# LANCZOS RESAMPLING KERNEL

"""
    LanczosKernel{S,T}()

yields an instance of a Lanczos re-sampling kernel of support size `S` (which
must be even) and for floating-point type `T`.

The Lanczos re-sampling kernels are even cardinal functions which tend to be
normalized for large support size.  They are defined by:

    ker(x) = S/(2*(π*x)^2)*sin(π*x)*sin(2*π*x/S)     if |x| ≤ S/2
             0                                       if |x| ≥ S/2

The expression `ker'` yields the first derivative of a Lanczos re-sampling
kernel `ker` (also see the constructor [`LanczosKernelPrime`](@ref)).

""" LanczosKernel

"""
    LanczosKernelPrime{S,T}()

yields a kernel instance that is the 1st derivative of the Lanczos re-sampling
kernel (see [`LanczosKernel`](@ref)) of support size `S` and for floating-point
type `T`.

""" LanczosKernelPrime

struct LanczosKernel{S,T} <: Kernel{T,S}
    a::T   # 1/2 support
    b::T   # a/pi^2
    c::T   # pi/a
    function LanczosKernel{S,T}() where {S,T}
        (typeof(S) == Int && S > 0 && iseven(S)) || bad_lanczos_kernel_size()
        a = T(S>>1)
        new{S,T}(a, a/T(π)^2, T(π)/a)
    end
end

struct LanczosKernelPrime{S,T} <: Kernel{T,S}
    a::T   # 1/2 support
    c::T   # pi/a
    function LanczosKernelPrime{S,T}() where {S,T}
        (typeof(S) == Int && S > 0 && iseven(S)) || bad_lanczos_kernel_size()
        a = T(S>>1)
        new{S,T}(a, T(π)/a)
    end
end

@noinline bad_lanczos_kernel_size() =
    throw(ArgumentError("Lanczos kernel size must be an even integer"))

# Outer constructors.
LanczosKernel{S}(T::Type{<:AbstractFloat} = Float64) where {S} =
    LanczosKernel{S,T}()
LanczosKernelPrime{S}(T::Type{<:AbstractFloat} = Float64) where {S} =
    LanczosKernelPrime{S,T}()

iscardinal(::Union{K,Type{K}}) where {K<:LanczosKernel} = true
iscardinal(::Union{K,Type{K}}) where {K<:LanczosKernelPrime} = false

isnormalized(::Union{K,Type{K}}) where {K<:LanczosKernel} = false
isnormalized(::Union{K,Type{K}}) where {K<:LanczosKernelPrime} = false

support(::LanczosKernel{S,T}) where {S,T} = SymmetricSupport{T,S,Open,Open}()
support(::LanczosKernelPrime{S,T}) where {S,T} = SymmetricSupport{T,S,Open,Open}()

# Expression for non-zero argument in the range (-S/2,S/2).
@inline _p(ker::LanczosKernel{S,T}, x::T) where {S,T} =
    ker.b*sin(π*x)*sin(ker.c*x)/(x*x)

# Expression for non-zero argument in the range (-S/2,S/2).
@inline function _p(ker::LanczosKernelPrime{S,T}, x::T) where {S,T}
    x1 = π*x
    s1, c1 = sin(x1), cos(x1)
    r1 = s1/x1
    x2 = ker.c*x # π*x/a
    s2, c2 = sin(x2), cos(x2)
    r2 = s2/x2
    return (c1*r2 + c2*r1 - 2*r1*r2)/x
end

(ker::LanczosKernel{S,T})(x::T) where {S,T} =
    (abs(x) ≥ ker.a ? zero(T) : x == 0 ? one(T) : _p(ker, x))

(ker::LanczosKernelPrime{S,T})(x::T) where {S,T} =
    (abs(x) ≥ ker.a ? zero(T) : x == 0 ? one(T) : _p(ker, x))

@generated function compute_weights(ker::LanczosKernel{S,T}, t::T) where {S,T}
    c = (S >> 1) # central index
    W = [Symbol("w",i) for i in 1:S] # all weights
    Expr(:block,
         Expr(:meta, :inline),
         Expr(:local, [:($w::T) for w in W]...),
         Expr(:if, :(t == zero(T)),
              Expr(:block, [:($(W[i]) = $(i == c ? 1 : 0)) for i in 1:S]...),
              Expr(:block,
                   [:($(W[i]) = _p(ker, t + $(c - i))) for i in 1:c-1]...,
                   :($(W[c]) = _p(ker, t)),
                   [:($(W[i]) = _p(ker, t - $(i - c))) for i in c+1:S]...)),
         Expr(:return, Expr(:tuple, W...)))
end

@generated function compute_weights(ker::LanczosKernelPrime{S,T},
                                    t::T) where {S,T}
    c = (S >> 1) # central index
    W = [Symbol("w",i) for i in 1:S] # all weights
    Expr(:block,
         #Expr(:meta, :inline),
         Expr(:local, [:($w::T) for w in W]...),
         [:($(W[i]) = ker(t + $(c - i))) for i in 1:c-1]...,
         :($(W[c]) = ker(t)),
         [:($(W[i]) = ker(t - $(i - c))) for i in c+1:S]...,
         Expr(:return, Expr(:tuple, W...)))
end

#------------------------------------------------------------------------------

@inline frac(::Type{T}, a::Integer, b::Integer) where {T<:AbstractFloat} = (T(a)/T(b))
@inline square(x) = x*x
@inline cube(x) = x*x*x

"""
    signabs(x) -> sign(x), abs(x)

yields the sign and absolute value of `x` both with the same type as `x`.

"""
@inline signabs(x::Unsigned) = (one(x), x)
@inline signabs(x::Real) = ifelse(x > 0, (one(x), x),
                                  ifelse(x < 0, (oftype(one(x),-1), -x),
                                         (zero(x), zero(x))))

#------------------------------------------------------------------------------

# Manage to call the short version of `show` for MIME output.
Base.show(io::IO, ::MIME"text/plain", ker::Kernel) = show(io, ker)

function Base.show(io::IO, ker::Kernel)
    summary(io, ker)
    print(io, '(')
    join(io, values(ker), ',')
    print(io, ')')
end

"""
    InterpolationKernels.with_eltype(T, ker)

yields an instance (resp. a type) of interpolation kernel instance (resp. type)
`ker` but with floating-point type `T`.

"""
with_eltype(::Type{T}, ::Type{K}) where {T<:AbstractFloat,K<:Kernel{T}} = K
with_eltype(::Type{T}, ker::Kernel{T}) where {T<:AbstractFloat} = ker

# The first instance below is needed to automatically do nothing when
# converting to a kernel of the same type (remmeber that convert is called to
# instanciate structure fields).
convert(::Type{K}, ker::K) where {K<:Kernel} = ker
convert(::Type{K}, ker::Kernel) where {K<:Kernel} = K(ker)

Kernel(ker::Kernel) = ker
Kernel{T}(ker::Kernel) where {T<:AbstractFloat} = with_eltype(T, ker)
Kernel{T,S}(ker::Kernel{<:Any,S}) where {T<:AbstractFloat,S} = with_eltype(T, ker)

for T in FLOATS
    @eval Base.$T(ker::Kernel) = with_eltype($T, ker)
end

"""
    InterpolationKernels.brief(ker)

yields a brief description of the kernel type or instance `ker`.

"""
brief(ker::Kernel) = brief(typeof(ker))

brief(::Type{<:BSpline{1}}) = "rectangular B-spline"
brief(::Type{<:BSplinePrime{1}}) = "derivative of rectangular B-spline"

brief(::Type{<:BSpline{2}}) = "linear B-spline"
brief(::Type{<:BSplinePrime{2}}) = "derivative of linear B-spline"

brief(::Type{<:BSpline{3}}) = "quadratic B-spline"
brief(::Type{<:BSplinePrime{3}}) = "derivative of quadratic B-spline"

brief(::Type{<:BSpline{4}}) = "cubic B-spline"
brief(::Type{<:BSplinePrime{4}}) = "derivative of cubic B-spline"

@noinline brief(::Type{<:BSpline{S}}) where {S} = "B-spline of order $S"
@noinline brief(::Type{<:BSplinePrime{S}}) where {S} =
    "derivative of B-spline of order $S"

brief(::Type{<:CubicSpline}) = "generic cubic spline"
brief(::Type{<:CubicSplinePrime}) = "derivative of generic cubic spline"

brief(::Type{<:CardinalCubicSpline}) = "cardinal cubic spline"
brief(::Type{<:CardinalCubicSplinePrime}) =
    "derivative of cardinal cubic spline"

brief(::Type{<:CatmullRomSpline}) = "Catmull & Rom cubic spline"
brief(::Type{<:CatmullRomSplinePrime}) =
    "derivative of Catmull & Rom cubic spline"

@noinline brief(::Type{<:LanczosKernel{S}}) where {S} =
    "Lanczos resampling kernel of size $S"
@noinline brief(::Type{<:LanczosKernelPrime{S}}) where {S} =
    "derivative of Lanczos resampling kernel of size $S"

# Manage to yield the derivative of (some) kernels when the notation `ker'` is
# used.
adjoint(ker::BSpline{S,T}) where {S,T} = BSplinePrime{S,T}()
adjoint(ker::CubicSpline{T}) where {T} = CubicSplinePrime{T}(ker.a, ker.b)
adjoint(ker::CardinalCubicSpline{T}) where {T} = CardinalCubicSplinePrime{T}(ker.a)
adjoint(ker::CatmullRomSpline{T}) where {T} = CatmullRomSplinePrime{T}()
adjoint(ker::LanczosKernel{S,T}) where {S,T} = LanczosKernelPrime{S,T}()

# Provide methods for all kernels.
for K in (:BSpline,             :BSplinePrime,
          :CubicSpline,         :CubicSplinePrime,
          :CardinalCubicSpline, :CardinalCubicSplinePrime,
          :CatmullRomSpline,    :CatmullRomSplinePrime,
          :LanczosKernel,       :LanczosKernelPrime)

    has_size = (K === :BSpline || K === :BSplinePrime ||
        K === :LanczosKernel || K === :LanczosKernelPrime)

    # We want that calling the kernel on a different type of real argument than
    # the floting-point type of the kernel convert the argument.
    # Unfortunately, defining:
    #
    #     (ker::$K{T})(x::Real) where {T<:AbstractFloat} = ker(convert(T,x))
    #
    # leads to ambiguities, the following is ugly but works...
    for T in FLOATS, R in (FLOATS..., :Integer)
        if R != T
            if has_size
                @eval @inline (ker::$K{S,$T})(x::$R) where {S} = ker($T(x))
            else
                @eval @inline (ker::$K{$T})(x::$R) = ker($T(x))
            end
        end
    end

    # Calling the kernel on an array.  FIXME: should be deprecated!
    if has_size
        @eval (ker::$K{S,T})(A::AbstractArray) where {S,T<:AbstractFloat} =
            map(ker, A)
    else
        @eval (ker::$K{T})(A::AbstractArray) where {T<:AbstractFloat} =
            map(ker, A)
    end

    # Change floating-point type.
    @eval $K{T}(ker::$K) where {T<:AbstractFloat} = with_eltype(T, ker)
    if has_size
        @eval begin
            $K{S,T}(ker::$K{S}) where {S,T<:AbstractFloat} = with_eltype(T, ker)
            with_eltype(::Type{T}, ::Type{<:$K{S}}) where {S,T<:AbstractFloat} = $K{S,T}
            with_eltype(::Type{T}, ker::$K{S}) where {S,T<:AbstractFloat} = $K{S,T}(values(ker)...)
            with_eltype(::Type{T}, ker::$K{T,S}) where {S,T<:AbstractFloat} = ker
        end
    else
        @eval begin
            with_eltype(::Type{T}, ::Type{<:$K}) where {T<:AbstractFloat} = $K{T}
            with_eltype(::Type{T}, ker::$K{T}) where {T<:AbstractFloat} = ker
            with_eltype(::Type{T}, ker::$K) where {T<:AbstractFloat} = $K{T}(values(ker)...)
        end
    end

    # Method `summary` prints the type of an object.
    if has_size
        @eval Base.summary(io::IO, ::$K{S,T}) where {S,T} =
            print(io, $(string(K, "{")), S, ',', T, '}')
    else
        @eval Base.summary(io::IO, ::$K{T}) where {T} =
            print(io, $(string(K, "{")), T, '}')
    end
end

end # module
