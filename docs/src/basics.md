# Basic usage

An interpolation kernel `Kernel{T,S}` is parameterized by the floating-point
type `T` of its coefficients and by the (integer) size `S` of its support.  For
efficiency reasons, only kernels with (small) finite size supports are
implemented.  To create a kernel instance, call its constructor; for example:

```julia
ker = LanczosKernel{6}()
```

yields a Lanczos re-sampling kernel of size 6 (see [Interpolation
kernels](kernels.md) for an exhaustive list of kernels implemented in
`InterpolationKernels`) and using the default floating-point type (that is,
`Float64`) for computations.  The same kind of kernel with a specific
floating-point type, say `T`, can be created by:

```julia
ker = LanczosKernel{6,T}()
```

Any interpolation kernel `ker` is a callable object which may be used as a
function with a real argument:

```julia
ker(x::Real)
```

yields kernel value at `x`.


## Basic methods

Some simple methods are available for any interpolation kernel `ker`:

- `eltype(ker)` yields the floating-point type `T` for calculations;

- `length(ker)` yields the number `S` of samples in the support of `ker` which
  is also the number of neighbors involved in an interpolation by this kernel.

Since the floating-point type `T` and the support size `S` are parameters of
the interpolation kernel type, the above methods can also be applied to the
type of an interpolation kernel.

Some interpolation kernels have numerical parameters, these parameters can be
retrieved by:

```julia
values(ker)
```

which yields a tuple of parameters, possibly empty.  A kernel instance
identical to `ker` can be built as follows:

```julia
typeof(ker)(values(ker)...)
```


## Kernel floating-point type

Calling a kernel `ker` with a real argument `x`, as `ker(x)`, always yield a
floating-point of type `T = eltype(ker)`.  This property is imposed for
efficiency reasons when interpolating arrays.  Calling a kernel `ker` with an
argument `x` that has a different floating-point type is therefore less
efficient as it involves converting the value of the real `x`.  It is however
very easy to change the floating-point type used by a kernel.

A kernel object `ker` can be *converted* to use a given floating-point type.
For example, assuming `ker` is a kernel instance, one can do either of:

```julia
convert(Kernel{Float32}, ker)
Kernel{Float32}(ker)
Float32(ker)
```

to use `Float32` floating-point arithmetic.


## Kernel support

The call:

```julia
support(ker)
```

yields the support of the kernel `ker` (and instance of
[`InterpolationKernels.Support`](@ref).  All kernels implemented in
`InterpolationKernels` have symmetric supports; that is, `ker(x)` is zero if
`abs(x) > S/2` with `S` the size of the kernel support.  `InterpolationKernels`
however provides the framework for any kind of support (symmetric,
left-/right-anchored, open, closed, semi-open, ...).  Methods [`infimum`](@ref)
and [`supremum`](@ref) respectively yield the lower and upper bounds of a
kernel support.


## Traits

Methods [`iscardinal`](@ref) and [`isnormalized`](@ref) respectively yield
whether a kernel is a cardinal function (that is, a function which yields zero
for all non-zero integers) and whether a kernel has the partition of unity
property.  For some parametric kernels, these traits depend on the specific
values of the parameters, so these methods only take a kernel instance (not a
type) as argument.


## Derivative

The expression `ker'` yields a kernel instance which is the 1st derivative of
the kernel `ker`.  Having the derivative of a kernel is useful in a number of
practical cases.  For instance, thanks to linearity of the [interpolation
procedure](interpolation.md), interpolating an array `A` with the derivative
`ker'` is equivalent to taking the 1st derivative of the continuous function
modeled by interpolating the array `A` with the kernel `ker`.


## Interpolation weights

Efficient computation of interpolation weights is implemented by the
[`InterpolationKernels.compute_offset_and_weights`](@ref)
[`InterpolationKernels.compute_weights`](@ref) methods.  These methods are not
exported because they are only required to implement array interpolation, as
done by the [FineShift](https://github.com/emmt/FineShift.jl) or
[LinearInterpolators](https://github.com/emmt/LinearInterpolators.jl) packages.
The principles of interpolation are detailed in [another
section](interpolation.md).
