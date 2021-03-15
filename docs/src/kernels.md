# Kernels

`InterpolationKernels` provide the following kernels (`T` is an optional
floating-point type):

| Types                        | Support | Cardinal | Parameters | Description                                                     |
|:-----------------------------|:-------:|:--------:|:----------:|:----------------------------------------------------------------|
| `BSpline{1,T}`               | 1       | yes      |            | [rectangular B-spline](#Rectangular-interpolation-kernel)       |
| `BSpline{2,T}`               | 2       | yes      |            | [linear B-spline](#Triangular-interpolation-kernel)             |
| `BSpline{3,T}`               | 3       | no       |            | [quadratic B-spline](#Quadratic-B-spline)                       |
| `BSpline{4,T}`               | 4       | no       |            | [cubic B-spline](#Cubic-B-spline)                               |
| `CubicSpline{T}`             | 4       | yes/no   | `a`, `b`   | [cubic spline](#Cubic-splines)                                  |
| `CardinalCubicSpline{T}`     | 4       | yes      | `a`        | [cardinal cubic spline](#Cardinal-cubic-spline)                 |
| `CatmullRomSpline{T}`        | 4       | yes      |            | [Catmull-Rom kernel](#Catmull-Rom-kernel)                       |
| `MitchellNetravaliSpline{T}` | 4       | yes/no   | `b`, `c`   | [Mitchell & Netravali kernels](#Mitchell-and-Netravali-kernels) |
| `LanczosKernel{S,T}`         | `S`     | yes      | `S`        | [Lanczos re-sampling kernels](#Lanczos-re-sampling-kernels)     |

In this table, the *Cardinal* column indicates whether a kernel `ker(x)` is a
cardinal function, that is `ker(k) = 0` for all integers `k` except that
`ker(0) = 1` which makes such a kernel directly suitable for interpolation.

A *spline* is defined as a piecewise polynomial of given degree.  A *B-spline*
of order `S` is a piecesiwe polynomial of degree `S - 1` which is everywhere
nonnegative, an even and normalized function (its integral is equal to 1).


## Rectangular interpolation kernel

The **rectangular interpolation kernel** (also known as **box kernel** or
**Fourier window** or **Dirichlet window**) is the 1st order B-spline equals to
`1` on `[-1/2,+1/2)`, and `0` elsewhere.  An instance of a rectangular
interpolation kernel is created by:

```julia
ker = BSpline{1,T}()
```

where the floating-point type `T` is assumed to be `Float64` if omitted.

The expression `ker'` yields the first derivative of the kernel `ker`.  An
instance of a kernel implementing the first derivative of a rectangular
interpolation kernel may be directly created by:

```julia
 BSplinePrime{1,T}()
```

where, again, the floating-point type parameter `T` my be omitted.


## Triangular interpolation kernel

The **triangular interpolation kernel** linear spline (also known as **triangle
kernel** or **Bartlett window** or **Fejér window**) is the 2nd order B-spline
defined by:

```
ker(x) = 1 - |x|       if |x| ≤ 1
         0             if |x| ≥ 1
```

An instance of a triangular interpolation kernel is created by:

```julia
ker = BSpline{2,T}()
```

where the floating-point type `T` is assumed to be `Float64` if omitted.

The expression `ker'` yields the first derivative of the kernel `ker`.  An
instance of a kernel implementing the first derivative of a triangular
interpolation kernel may be directly created by:

```julia
 BSplinePrime{2,T}()
```

where, again, the floating-point type parameter `T` my be omitted.


## Quadratic B-spline

The **quadratic B-spline** is the 3rd order B-spline defined by:

```
ker(x) = 3/4 - x^2                   if |x| ≤ 1/2
         (1/2)*(|x| - 3/2)^2         if |x| ≤ 3/2
         0                           if |x| ≥ 3/2
```

An instance of a quadratic B-spline is created by:

```julia
BSpline{3,T}()
```

where the floating-point type `T` is assumed to be `Float64` if omitted.

The expression `ker'` yields the first derivative of the kernel `ker`.  An
instance of a kernel implementing the first derivative of a quadratic B-spline
may be directly created by:

```julia
 BSplinePrime{3,T}()
```

where, again, the floating-point type parameter `T` my be omitted.


## Cubic B-spline

The 4th order (cubic) B-spline kernel (also known as **Parzen window** or as
**de la Vallée Poussin window**) is defined by:

```
ker(x) = (|x|/2 - 1)*|x|^2 + 2/3     if |x| ≤ 1
         (1/6)*(2 - |x|)^3           if 1 ≤ |x| ≤ 2
         0                           if |x| ≥ 2
```

An of a cubic B-spline is created by:

```julia
BSpline{4,T}()
```

where the floating-point type `T` is assumed to be `Float64` if omitted.

The cubic B-spline is a C² continuous function (its derivatives up to the
second one are everywhere continuous).  The expression `ker'` yields the first
derivative of the kernel `ker`.  An instance of a kernel implementing the first
derivative of a cubic B-spline may be directly created by:

```julia
 BSplinePrime{4,T}()
```

where, again, the floating-point type parameter `T` my be omitted.


## Cubic splines

An instance of the familily of cubic spline kernels is created by:

```julia
ker = CubicSpline{T}(a, b)
```

and is defined by:

```
ker(x) = ((2 + a - 6b)*|x| +  (9b - a - 3))*x^2 + (1 - 2b)  if |x| ≤ 1
         ((a + 2b)*|x| - (a + b))*(|x| - 2)^2               if |x| ≤ 2
         0                                                  if |x| ≥ 2
```

The parameters `a = ker'(1)` and `b = ker(1)` are the slope and the value of
the function `ker(x)` at `x = 1`.  The type parameter `T` is the floating-point
type for computations `T`, it may be omitted in which case it is guessed from
the types of `a` and `b` as `T = float(promote(typeof(a), typeof(b)))`.

A cubic spline kernel is at least C¹ continuous, the expression `ker'` yields a
kernel instance implementing the 1st derivative of the generic cubic spline
`ker`.  Such a derivative may be directly built by calling:

```julia
CubicSplinePrime{T}(a, b)
```

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

Calling [`BSpline`](@ref), [`CatmullRomSpline`](@ref),
[`CardinalCubicSpline`](@ref), or [`MitchellNetravaliSpline`](@ref) to build
these more specialized cubic splines may yield more efficient kernels as
computations involve more simple expressions.  Instances of `CubicSpline` are
however very well optimized and, in practice, they may be as fast or even
faster than their more specialized counterparts.  If ultimate performances
matter, the [`BenchmarkTools`](https://github.com/JuliaCI/BenchmarkTools.jl)
package may helps you to decide which kernel to choose for a given machine.


## Keys cardinal cubic kernels

Keys kernels are parametric cardinal cubic splines defined by:

```
ker(x) = 1 - (a + 3)*x^2 + (a + 2)*|x|^3      if |x| ≤ 1
         -4a + 8a*|x| - 5a*x^2 + a*|x|^3      if 1 ≤ |x| ≤ 2
         0                                    if |x| ≥ 2
```

These kernels are C¹ continuous piecewise normalized cardinal cubic spline
which depend on one parameter `a = ker'(1)` the slope of the function `ker(x)`
at `x = 1`.

To create an instance of a cardinal cubic spline with parameter `a`, call:

```julia
ker = CardinalCubicSpline{T}(a)
```

where `T` is the floating-point type for computations.  If omitted, `T` is
`typeof(a)` if it is `Float16`, `Float32`, `Float64`, or `BigFloat` and
`Float64` otherwise.

The expression `ker'` yields the first derivative of the cardinal cubic spline
`ker`.  An instance of the first derivative of such a kernel can also be
directly created by:

```julia
KeysSplinePrime{T}(a)
```

with the same defalut for `T` if this parameter is omitted.

### References

* Keys, Robert, G., "*Cubic Convolution Interpolation for Digital Image
  Processing*", IEEE Trans. Acoustics, Speech, and Signal Processing,
  Vol. ASSP-29, No. 6, December 1981, pp. 1153-1160.


## Catmull-Rom kernel

The Catmull-Rom kernel is a cardinal piecewise cubic spline defined by:

```
ker(x) = ((3/2)*|x| - (5/2))*x^2 + 1             if |x| ≤ 1
         (((5/2) - (1/2)*|x|)*|x| - 4)*|x| + 2   if 1 ≤ |x| ≤ 2
         0                                       if |x| ≥ 2
```

Being a cardinal function, the Catmull-Rom kernel is suitable for
interpolation.  Its derivative is given by:

```
ker′(x) = ((9/2)*|x| - 5)*x                      if a = |x| ≤ 1
          (5 - (3/2)*|x|)*x - 4*sign(x)          if 1 ≤ |x| ≤ 2
          0                                      if |x| ≥ 2
```

To create an instance of a Catmull-Rom interpolation kernel, call

```julia
ker = CatmullRomSpline{T}()
```

where the floating-point type `T` is assumed to be `Float64` if omitted.  To
create an instance of the derivative of a Catmull-Rom interpolation kernel,
call one of:

```julia
kerp = ker'
kerp = CatmullRomSplinePrime{T}()
```


## Mitchell & Netravali kernels

Mitchell & Netravali kernels are piecewise cubic splines defined by:

```
ker(x) = (1/6)*(((12 - 9b - 6c)*|x| - 18 + 12b + 6c)*x^2 + (6 - 2b))   if |x| ≤ 1
         (1/6)*(2b + 6c - (b + 6c)*|x|)*(2 - |x|)^2                    if 1 ≤ |x| ≤ 2
         0                                                             if |x| ≥ 2
```

Instances of a Mitchell & Netravali kernel and of its derivative are
respectively created by:

```julia
MitchellNetravaliSpline{T}(b, c)
MitchellNetravaliSplinePrime{T}(b, c)
```

with `T` the floating-point type for computations.  If `T` is omitted but `b`
and `c` are specified, `T` is deduced from the floating-point type of `b` and
`c`.  If `b` and `c` are omitted, `(b,c) = (1/3,1/3)` is assumed as recommended
by Mitchell & Netravali.  If `T` is omitted and `b` and `c` are omitted or both
integers, `T = Float64` is assumed.

Whatever the values of the parameters `b` and `c`, Mitchell & Netravali kernels
are normalized and even functions of class C¹ (these kernels and their first
derivatives are continuous).  The expression `ker'` yields the first derivative
of a Mitchell & Netravali kernel `ker`.

Taking `b = 0` yields Keys's family of kernels and is a sufficient and
necessary condition to have Mitchell & Netravali kernels be cardinal functions.

Using the constraint: `b + 2c = 1` yields a cubic filter with, at least,
quadratic order approximation.

Some specific values of `(b,c)` yield other well known kernels:

```
(b,c) = (1,0)      --> cubic B-spline
(b,c) = (0,-a)     --> Keys's cardinal cubic spline CardinalCubicSpline(a)
(b,c) = (0,1/2)    --> Catmull-Rom kernel CatmullRomSpline()
(b,c) = (b,0)      --> Duff's tensioned B-spline
(b,c) = (6β,-α-3β) --> generic cubic spline CubicSpline(α,β)
(b,c) = (1/3,1/3)  --> recommended by Mitchell-Netravali
```

This family of kernels and their derivatives is implemented as instances of
`CubicSpline` and `CubicSplinePrime` using the following property:

```julia
MitchellNetravaliSpline{T}(b, c) = CubicSpline{T}(-b/2 - c, b/6)
```


### References

* Mitchell & Netravali, "*Reconstruction Filters in Computer Graphics*",
  in Computer Graphics, Vol. 22, Num. 4 (1988)
  [pdf](http://www.cs.utexas.edu/users/fussell/courses/cs384g/lectures/mitchell/Mitchell.pdf).


## Lanczos re-sampling kernels

The Lanczos re-sampling kernels are defined by:

```
ker(x) = (S/2π²)*sin(π*x)*sin(2π*x/S)/x^2   if |x| ≤ S/2
         0                                  if |x| ≥ S/2
```

The Lanczos re-sampling kernels are even cardinal functions which tend to be
normalized for large support size (see [this Wikipedia
page](https://en.wikipedia.org/wiki/Lanczos_resampling)).

To create an instance of a Lanczos re-sampling kernel of support size `S`
(which must be even), call:

```julia
LanczosKernel{S,T}()
```

where the floating-point type `T` is assumed to be `Float64` if omitted.

The expression `ker'` yields the first derivative of a Lanczos re-sampling
kernel `ker`.  An instance of a kernel implementing the first derivative of a
Lanczos re-sampling kernel can also be directly created by:

```julia
LanczosKernelPrime{S,T}()
```

where, again, it is possible to omit the floating-point type parameter `T`.
