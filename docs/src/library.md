# Reference

The following provides detailled documentation about types and methods provided
by the `InterpolationKernels` package.  This information is also available from
the REPL by typing `?` followed by the name of a method or a type.

# Methods

```@docs
iscardinal
isnormalized
brief
getweights
```

## Boundary conditions

```@docs
Boundaries
Flat
SafeFlat
boundaries
```

## Kernels

```@docs
Kernel
```

### Rectangular B-spline

```@docs
RectangularSpline
RectangularSplinePrime
```

### Linear B-spline

```@docs
LinearSpline
LinearSplinePrime
```

### Quadratic B-spline

```@docs
QuadraticSpline
QuadraticSplinePrime
```

### Cubic B-spline

```@docs
CubicSpline
CubicSplinePrime
```

### Catmull-Rom interpolation kernel

```@docs
CatmullRomSpline
CatmullRomSplinePrime
```

### Cardinal cubic splines

```@docs
CardinalCubicSpline
CardinalCubicSplinePrime
```

### Keys cardinal kernels

```@docs
KeysSpline
KeysSplinePrime
```

### Mitchell-Netravali kernels

```@docs
MitchellNetravaliSpline
MitchellNetravaliSplinePrime
```

### Lanczos re-sampling kernels

```@docs
LanczosKernel
LanczosKernelPrime
```
