# Basic usage

An interpolation kernel `Kernel{T,S,B}` is parametrized by the floating-point
type `T` of its coefficients, by the size `S` of its support and by the
boundary conditions `B` applied for extrapolation (see [Boundary
conditions](boundaries.md)).  For efficiency reasons, only kernels with (small)
finite size supports are implemented.  To create a kernel instance, call its
constructor; for example:

```julia
ker = LanczosKernel(6)
```

yields a Lanczos re-sampling kernel of size 6 (see [Interpolation
kernels](kernels.md) for an exhaustive list of kernels implemented in
`InterpolationKernels`).

Any interpolation kernel `ker` is a callable object which may be used as a
function with a real argument:

```julia
    ker(x::Real)
```

yields kernel value at offset `x`.  All kernel supports are symmetric; that is
`ker(x)` is zero if `abs(x) > S/2` with `S` the size of the kernel support.


## Simple methods

Some simple methods are available for any interpolation kernel `ker`:

- `eltype(ker)` yields the floating-point type `T` for calculations;

- `length(ker)` yields the number `S` of samples in the support of `ker` which
  is also the number of neighbors involved in an interpolation by this kernel;

- [`boundaries(ker)`](@ref) yields `B` the type of the boundary conditions
  applied for extrapolation.

Since the floating-point type `T`, the support size `S` and the boundary
conditions `B` are parameters of the interpolation kernel type, the above
methods can also be applied to the type of an interpolation kernel.


## Traits

Methods [`iscardinal`](@ref) and [`isnormalized`](@ref) can be used to query
whether a kernel is a cardinal function (that is a function which yields zero
for all non-zero integers) and whether a kernel has the partition of unity
property.  For some parametric kernels, these traits depend on the specific
values of the parameters so these methods take a kerenl insta,ce (not a type)
as argument.


## Derivative

The expression `ker'` yields a kernel instance which is the 1st derivative of
the kernel `ker`.  Having the derivative if a kernel is useful because, thanks
to linearity of the interpolation, interpolating an array with the derivative
`ker'` is equivalent of taking the 1st derivative of the continuous function
given by interpolating the same array with the kernel `ker`.


## Interpolation weights

Efficient computation of interpolation weights is implemented by the
[`getweights`](@ref).  The principles of interpolation are detailed in [another
section](interpolation.md).


## Kernel conversion

A kernel object `ker` can be *converted* to change its floating-point type
and/or its boundary conditions.  For example, assuming `ker` is a kernel
instance, one can do:

```julia
convert(Kernel{Float32}, ker)
convert(Kernel{eltype(ker),length(ker),SafeFlat}, ker)
```

to use `Float32` floating-point arithmetic or to choose [`SafeFlat`](@ref)
boundary conditions but keeping the same floating-point type and (of course)
the same support size.

Such expressions are a bit tedious to type, so shortcuts are provided.  It is
sufficient to call a kernel with the new floating-point type `T` and/or boundary
conditions `B` to perform the conversion:

```julia
ker([T::Type{<:AbstractFloat} = eltype(ker),] B::Type{<:Boundaries}=boundaries(ker))
```

where any of `T` or `B` can be omitted to keep the current kernel setting and
their order is irrelevant.  Beware that changing the floating-point type may
lead to a loss of precision if the new floating-point type has more digits.  It
is possible to change the floating-point type of a kernel or its boundary
conditions by something like:

```julia
Float32(ker)    # change floating-point type of kernel `ker`
SafeFlat(ker)   # change boundary conditions of kernel `ker`
```

The above calls do not follow the usuall conventions that a constructor may be
called to convert its argument(s) to an instance of its type, this is however
practical.
