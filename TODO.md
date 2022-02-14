- Add prefiltering by Thévenaz and Unser.

- Kernel parameter `T` can be `AbstractFloat` to have a kernel whose result
  type depends on that of the argument (as many numerical floating-point
  functions).  This may be very simple to implement but can only work
  efficiently for parameter-less kernels.  Kernel evaluation could be written
  as:

  ```.jul
  (ker::SomeKernelType{S,T}(x::Real) where {S,T} = call(ker, convert(T, x))
  function call(ker::SomeKernelType{S}(x::T) where {S,T<:AbstractFloat}
      # Here assume that `T` is the concrete flating-point type of the result.
      ...
  end
  ```

  Note that kernel floating-point type is ignored in the `call`, only that of
  the argument `x` matters.  This exploits the fact that calling
  `convert(AbstractFloat,x)` is a no-op if `x` has concrete floating-point type
  and yield `x` converted to `Float64` otherwise (e.g. irrational).  Also note
  that `one(AbstractFloat)` and `zero(AbstractFloat)` both yield a `Float64`.

- `iscardinal` and `isnormalized` may give different result for parametric
  kernels depending on whether the argument is a kernel type or a kernel
  instance?

- Use SIMD.jl for `compute_weights`.
