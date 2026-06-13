abstract type Bound end
struct Open <: Bound end
struct Closed <: Bound end

"""
    InterpolationKernels.Support{T,S,L,R}

Abstract super-type for the support of an interpolation kernel parameterized by the
floating-point type `T` used for computations, the integer size `S` of the support and the
types `L` and `R` of the left and right bounds which can be `InterpolationKernels.Open` or
`InterpolationKernels.Closed`.

"""
abstract type Support{T<:AbstractFloat,S,L<:Bound,R<:Bound} end

Base.length(sup::Support) = length(typeof(sup))
Base.length(::Type{<:Support{T,S}}) where {T,S} = S

Base.eltype(sup::Support) = eltype(typeof(sup))
Base.eltype(::Type{<:Support{T,S}}) where {T,S} = T

"""
    InterpolationKernels.SymmetricSupport{T,S,L,R}() -> sup

Return an instance of a symmetric support parameterized by the floating-point type `T`, the
integer size `S` of the support and the types `L` and `R` of the left and right bounds which
can be `InterpolationKernels.Open` or `InterpolationKernels.Closed`.

Depending on `L` and `R`, the support is:

    (L,R) = (Open,Open) ───────→ sup = (-S/2,S/2)
    (L,R) = (Closed,Open) ─────→ sup = [-S/2,S/2)
    (L,R) = (Open,Closed) ─────→ sup = (-S/2,S/2]
    (L,R) = (Closed,Closed) ───→ sup = [-S/2,S/2]

"""
struct SymmetricSupport{T<:AbstractFloat,S,L,R} <: Support{T,S,L,R} end

"""
    InterpolationKernels.LeftAnchoredSupport{T,S,L,R}(a)

Return an instance of a support with lower bound `a` and parameterized by the floating-point
type `T`, the integer size `S` of the support and the types `L` and `R` of the left and
right bounds which can be `InterpolationKernels.Open` or `InterpolationKernels.Closed`.

Depending on `L` and `R`, the support is:

    (L,R) = (Open,Open) ───────→ sup = (a,a+S)
    (L,R) = (Closed,Open) ─────→ sup = [a,a+S)
    (L,R) = (Open,Closed) ─────→ sup = (a,a+S]
    (L,R) = (Closed,Closed) ───→ sup = [a,a+S]

"""
struct LeftAnchoredSupport{T,S,L,R} <: Support{T,S,L,R}
    a::T
end

"""
    InterpolationKernels.RightAnchoredSupport{T,S,L,R}(b)

Return an instance of a support with upper bound `b` and parameterized by the floating-point
type `T`, the integer size `S` of the support and the types `L` and `R` of the left and
right bounds which can be `InterpolationKernels.Open` or `InterpolationKernels.Closed`.

Depending on `L` and `R`, the support is:

    (L,R) = (Open,Open) ───────→ sup = (b-S,b)
    (L,R) = (Closed,Open) ─────→ sup = [b-S,b)
    (L,R) = (Open,Closed) ─────→ sup = (b-S,b]
    (L,R) = (Closed,Closed) ───→ sup = [b-S,b]

"""
struct RightAnchoredSupport{T,S,L,R} <: Support{T,S,L,R}
    b::T
end

"""
    InterpolationKernels.supremum(sup) -> b

Return the least upper bound `b` of the kernel support `sup`.

"""
supremum(::SymmetricSupport{T,S}) where {T,S} = T(S)/2
supremum(sup::LeftAnchoredSupport{T,S}) where {T,S} = infimum(sup) + T(S)
supremum(sup::RightAnchoredSupport) = sup.b

"""
    InterpolationKernels.infimum(sup) -> a

Return the greatest lower bound `a` of the kernel support `sup`.

"""
infimum(::SymmetricSupport{T,S}) where {T,S} = T(-S)/2
infimum(sup::LeftAnchoredSupport) = sup.a
infimum(sup::RightAnchoredSupport{T,S}) where {T,S} = supremum(sup) - T(S)

Base.in(t, I::Support{T}) where {T} = in(convert(T, t), I)
Base.in(t::T, I::Support{T}) where {T} = above_left_bound(t, I) & below_right_bound(t, I)

above_left_bound(t::T, I::Support{T,<:Any,Open,<:Any}) where {T} = (infimum(I) < t)
above_left_bound(t::T, I::Support{T,<:Any,Closed,<:Any}) where {T} = (infimum(I) ≤ t)
below_right_bound(t::T, I::Support{T,<:Any,<:Any,Open}) where {T} = (t < supremum(I) < t)
below_right_bound(t::T, I::Support{T,<:Any,<:Any,Closed}) where {T} = (t ≤ supremumI)

TypeUtils.get_precision(I::Support{T}) where {T} = get_precision(T)
function TypeUtils.adapt_precision(::Type{T},
                                   I::SymmetricSupport{E,S,L,R}) where {T<:Precision,E,S,L,R}
    return SymmetricSupport{adapt_precision(T,E),S,L,R}()
end
