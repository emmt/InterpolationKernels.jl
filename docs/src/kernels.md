# Kernels

`InterpolationKernels` provide the following kernels:

| Name                      | Support | Cardinal | Parameters | Description                  |
|:--------------------------|:-------:|:--------:|:----------:|:-----------------------------|
| `RectangularSpline`       | 1       | yes      |            |                              |
| `LinearSpline`            | 2       | yes      |            | linear B-spline              |
| `QuadraticSpline`         | 3       | no       |            | quadratic B-spline           |
| `CubicSpline`             | 4       | no       |            | cubic B-spline               |
| `CatmullRomSpline`        | 4       | yes      |            | Catmull-Rom kernels          |
| `CardinalCubicSpline`     | 4       | yes      | `a`        |                              |
| `KeysSpline`              | 4       | yes      | `a`        | Keys cardinal kernels        |
| `MitchellNetravaliSpline` | 4       | yes/no   | `b`, `c`   | Mitchell & Netravali kernels |
| `LanczosKernel`           | `S`     | yes      | `S`        | Lanczos re-sampling kernels  |

In this table, the *Cardinal* column indicates whether a kernel `ker(x)` is a
cardinal function, that is `ker(k) = 0` for all integers `k` except that
`ker(0) = 1` which makes such a kernel suitable for interpolation.

A *spline* is defined as a piecewise polynomial of given degree.  A *B-spline*
is a spline which is everywhere nonnegative, an even and normalized function
(its integral is equal to 1).



## Rectangular Spline

The rectangular spline (also known as box kernel or Fourier window or Dirichlet
window) is the 1st order (constant) B-spline equals to `1` on `[-1/2,+1/2)`,
and `0` elsewhere.

A rectangular spline instance is created by:

```julia
ker = RectangularSpline([T=Float64,] B=Flat)
```

Its derivative is created by:

```julia
dker = RectangularSplinePrime([T=Float64,] B=Flat)
```

or by:

```julia
dker = ker'
```

## Linear B-spline

The linear spline (also known as triangle kernel or Bartlett window or Fejér
window) is the 2nd order (linear) B-spline.

A linear spline instance is created by:

```julia
ker = LinearSpline([T=Float64,] B=Flat)
```

Its derivative is created by:

```julia
dker = LinearSplinePrime([T=Float64,] B=Flat)
```
or by:

```julia
dker = ker'
```

## Quadratic B-spline

The quadratic spline is the 3rd order (quadratic) B-spline.

A quadratic spline instance is created by:

```julia
QuadraticSpline([T=Float64,] B=Flat)
```

Its derivative is created by:

```julia
QuadraticSplinePrime([T=Float64,] B=Flat)
```

## Cubic B-spline

The 4th order (cubic) B-spline kernel is defined by:

```
ker(x) = (|x|/2 - 1)*|x|^2 + 2/3     if |x| ≤ 1
         (1/6)*(2 - |x|)^3           if 1 ≤ |x| ≤ 2
         0                           if |x| ≥ 2
```

The cubic B-spline is also known as *Parzen window* or as *de la Vallée Poussin
window*.

The cubic B-spline is a function of class C² (its derivatives up to the second
one are everywhere continuous), its first derivative is given by:

```
ker'(x) = ...
```

To create an instance of a cubic B-spline, call:

```julia
ker = CubicSpline([T=Float64,] B=Flat)
```

where optional arguments are the floating-point type `T<:AbstractFloat` and the
boundary conditions `B< Boundaries` (their order is irrelevant).

To create an instance of the derivative of a cubic B-spline, call one of:


```julia
kerp = ker'
kerp = CubicSplinePrime([T=Float64,] B=Flat)
```


## Catmull-Rom kernel

The Catmull-Rom kernel is a special case of Mitchell & Netravali kernel.  It is
a piecewise cardinal cubic spline defined by:

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
ker = CatmullRomSpline([T=Float64,] B=Flat)
```

where optional arguments are the floating-point type `T` and the boundary
conditions `B`.  To create an instance of the derivative of a Catmull-Rom
interpolation kernel, call one of:

```julia
kerp = ker'
kerp = CatmullRomSplinePrime([T=Float64,] B=Flat)
```


## Cardinal cubic spline

```julia
CardinalCubicSpline([T=Float64,] c, B=Flat) -> ker
```

yields a cardinal cubic spline interpolation kernel for floating-point type `T`
tension parameter `c` and boundary conditions `B`.  The slope at `x = ±1` is
`±(c - 1)/2`.  Usually `c ≤ 1`, choosing `c = 0` yields a Catmull-Rom spline,
`c = 1` yields all zero tangents, `c = -1` yields a truncated approximation of
a cardinal sine, `c = -1/2` yields an interpolating cubic spline with
continuous second derivatives (inside its support).

Its derivative is given by:

```julia
CardinalCubicSplinePrime([T=Float64,] c, B=Flat) -> ker
```

```julia
ker(x) = |x| < 1 ? ((r*|x| - 1)*|x| - 1)*(|x| - 1) :
         |x| < 2 ? p*(|x| - 1)*(2 - |x|)^2            : 0
```

with `p = ...` and `r = ...`.


## Keys cardinal kernels

Keys kernels are parametric cardinal cubic splines defined by:

```
ker(x) = p(|x|)   if |x| ≤ 1
         q(|x|)   if 1 ≤ |x| ≤ 2
         0        if |x| ≥ 2
```

with:

```
p(x) = 1 - (a + 3)*x^2 + (a + 2)*x^3
q(x) = -4a + 8a*x - 5a*x^2 + a*x^3
```

These kernels are piecewise normalized cardinal cubic spline which depend on
one parameter `a` which is the slope of the spline at abscissa 1.

The first derivative of a Keys spline is defined by:

```
ker′(x) = p′(|x|)*sign(x)   if |x| ≤ 1
          q′(|x|)*sign(x)   if 1 ≤ |x| ≤ 2
          0                 if |x| ≥ 2
```

with:

```
p'(x) = -2*(a + 3)*x + 3*(a + 2)*x^2
q'(x) = 8a - 10a*x + 3a*x^2
```

To create an instance of a Keys spline with parameter `a`, call:

```julia
ker = KeysSpline([T=Float64,] a, B=Flat)
```

where optional parameters are the floating-point type `T` and the boundary
conditions `B`.

The expression `ker'` yields the first derivative of a Keys kernel `ker`.  An
instance of the first derivative of a Keys kernel can also be directly created
by:

```julia
KeysSplinePrime([T=Float64,] a, B=Flat)
```

Reference:

* Keys, Robert, G., "Cubic Convolution Interpolation for Digital Image
  Processing", IEEE Trans. Acoustics, Speech, and Signal Processing,
  Vol. ASSP-29, No. 6, December 1981, pp. 1153-1160.


## Mitchell & Netravali kernels

Mitchell & Netravali kernels are cubic splines which depends on 2 parameters,
`b` and `c`.  Whatever the values of `(b,c)`, all these kernels are normalized
and even functions of class C¹ (these kernels and their first derivatives are
continuous).

Taking `b = 0` yields Keys's family of kernels and is a sufficient and
necessary condition to have Mitchell & Netravali kernels be cardinal functions.

Using the constraint: `b + 2c = 1` yields a cubic filter with, at least,
quadratic order approximation.

Some specific values of `(b,c)` yield other well known kernels:

    (b,c) = (1,0)     ==> cubic B-spline
    (b,c) = (0,-a)    ==> Keys's cardinal cubics
    (b,c) = (0,1/2)   ==> Catmull-Rom cubics
    (b,c) = (b,0)     ==> Duff's tensioned B-spline
    (b,c) = (1/3,1/3) ==> recommended by Mitchell-Netravali

Mitchell & Netravali kernels are defined by:

```
ker(x) = p(|x|)      if |x| ≤ 1
         q(|x|)      if 1 ≤ |x| ≤ 2
         0           if |x| ≥ 2
```

with:

```
p(x) = (p3*x + p2)*x*x + p0
q(x) = ((q3*x + q2)*x + q1)*x + q0
```

and:

```
p0 =  1 - b/3
p2 = -3 + 2b + c
p3 =  2 - 3b/2 - c
q0 =  4b/3 +4*c
q1 = -2b - 8c
q2 =  b + 5c
q3 =  -b/6 - c
```

The first derivative of a Mitchell & Netravali kernel is given by:

```
ker'(x) = sign(x)*p'(|x|)      if |x| ≤ 1
          sign(x)*q'(|x|)      if 1 ≤ |x| ≤ 2
          0                    if |x| ≥ 2
```

with:

```
p'(x) = (dp2*x + dp1)*x
q'(x) = ((dq3*x + dq2)*x + dq1)*x
```

```
dp1 = 2*p2 = -6 + 4b + 2c
dp2 = 3*p3 =  6 - 9b/2 - 3c
dq0 =   q1 = -2b - 8c
dq1 = 2*q2 = 2b + 10c
dq2 = 3*q3 = -b/2 - 3c
```

To create an instance of a Mitchell & Netravali kernel, call:


```julia
MitchellNetravaliSpline([T=Float64,] [b=1/3, c=1/3,] B=Flat)
```

where optional parameters are the floating-point type `T`, the family
parameters `b` and `c` and the boundary conditions `B`.

The expression `ker'` yields the first derivative of a Mitchell & Netravali
kernel `ker`.  An instance of the first derivative of a Mitchell & Netravali
kernel can also be directly created by:

```julia
MitchellNetravaliSplinePrime([T=Float64,] [b=1/3, c=1/3,] B=Flat) -> ker
```

Reference:

* Mitchell & Netravali, "*Reconstruction Filters in Computer Graphics*",
  in Computer Graphics, Vol. 22, Num. 4 (1988)
  [pdf](http://www.cs.utexas.edu/users/fussell/courses/cs384g/lectures/mitchell/Mitchell.pdf).


## Lanczos re-sampling kernels

The Lanczos re-sampling kernels are defined by:

```
ker(x) = (S/2π²)*sin(π*x)*sin(2π*x/S)/x^2
```

The Lanczos re-sampling kernels are even cardinal functions which tend to be
normalized for large support size (see [this Wikipedia
page](https://en.wikipedia.org/wiki/Lanczos_resampling)).

To create an instance of a Lanczos re-sampling kernel of support size `S`
(which must be even), call:

```julia
LanczosKernel([T=Float64,] S, B=Flat)
```

where optional parameters are the floating-point type `T` and the boundary
conditions `B`.

The expression `ker'` yields the first derivative of a Lanczos re-sampling
kernel `ker`.  An instance of the first derivative of a Lanczos re-sampling
kernel can also be directly created by:

```julia
LanczosKernelPrime([T=Float64,] S, B=Flat)
```
