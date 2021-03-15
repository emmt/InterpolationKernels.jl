
## Version 0.2.0

* Boundary conditions are no longer part of the type definition of an
  interpolation kernel.

* Kernel type constructors should be type-stable.

* The `values` method can be applied to an interpolation kernel to retrieve its
  parameters.

* Possible floating-point types are limited to the concrete types in
  `InterpolationKernels.FLOATS` or the union `InterpolationKernels.Floats`.

* The `summary` method now (correctly) prints the type of a kernel and takes an
  instance of `IO` as first parameter.  The floating-point type is part of the
  textual representation written by the `show` method.

* New method `interpcoefs` that eventually calls `getweights` to compute the
  interpolation weights in a carefully optimized way.

* The order of type parameters for Lanczos resampling kernels have been
  exchanged.  This is to be more consistent with B-splines types and to have
  the size parameter `S` playing a more prominent role than the floating-point
  type parameter `T` which defaults to `Float64` and can thus be omitted.

* `RectangularSpline`, `LinearSpline`, `QuadraticSpline` and `CubicSpline` (and
  their derivatives `RectangularSplinePrime`, `LinearSplinePrime`,
  `QuadraticSplinePrime` and `CubicSplinePrime`) have been replaced by
  `BSpline{S}` (and derivative with `BSplinePrime{S}`) with `S` the **order**
  of the B-spline: `S = 1` for a **rectangular** B-spline, `S = 2` for a
  **linear** B-spline, `S = 3` for a **quadratic** B-spline, and `S = 4` for a
  **cubic** B-spline.

* `CubicSpline` (and derivative `CubicSplinePrime`) is a npw a new generic
  cubic spline which is C¹ continuous and parameterized by `a = ker'(1)` and `b
  = ker(1)` the slope and value of the function `ker(x)` at `x = 1` with `ker =
  CubicSpline()`.  This new kernel type is better optimized than the former
  Mitchell & Netravali kernel.  In the 0.1.x versions of the package,
  instances of `CubicSpline` were in fact cubic B-splines.

* Mitchell & Netravali kernels built by `MitchellNetravaliSpline` are now
  emulated by instances of `CubicSpline` and benefit from a speed up by a
  factor of about 1.5 thanks to this new family of kernels.  Computing the
  interpolation offset and weights now takes 18 operations (14 for the
  derivative) instead of 28 operations (20 for the derivative) in the former
  implementation.

* There is only one generic family of cardinal cubic splines implemented by
  type `CardinalCubicSpline` (and derivative `CardinalCubicSplinePrime`).  This
  family is parameterized by the slope `a = ker'(1)` the function `ker(x)` at
  `x = 1` with `ker = CardinalCubicSpline(a)` as were Keys' cardinal cubic
  splines in the previous version of the package (formely `KeysSpline` and
  `KeysSplinePrime`).  In the versions 0.1.x of the package,
  `CardinalCubicSpline` were cardinal cubic splines parameterized by the
  so-called *tension* defined by `c = 2*ker'(1) + 1`.  These splines can be
  emalated by the new cardinal cubic splines with parameter `a = (c - 1)/2`.

* Tests are done with reference implementations of the kernel functions.


## Version 0.1.1

* Documentation.
