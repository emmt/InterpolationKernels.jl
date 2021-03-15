# Reference

The following provides detailled documentation about types and methods provided
by the `InterpolationKernels` package.  This information is also available from
the REPL by typing `?` followed by the name of a method or a type.

# Methods

```@docs
iscardinal
isnormalized
InterpolationKernels.brief
InterpolationKernels.getweights
InterpolationKernels.getcoefs
```

## Kernels

```@docs
Kernel
```

### B-splines

```@docs
BSpline
BSplinePrime
```

### Generic cubic spline

```@docs
CubicSpline
CubicSplinePrime
```

### Cardinal cubic splines

```@docs
CardinalCubicSpline
CardinalCubicSplinePrime
```

### Catmull-Rom interpolation kernel

```@docs
CatmullRomSpline
CatmullRomSplinePrime
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
