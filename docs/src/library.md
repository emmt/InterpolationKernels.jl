# Reference

The following provides detailed documentation about types and methods provided
by the `InterpolationKernels` package. This information is also available from
the REPL by typing `?` followed by the name of a method or a type.

## Kernel supports

```@docs
InterpolationKernels.Support
InterpolationKernels.SymmetricSupport
InterpolationKernels.LeftAnchoredSupport
InterpolationKernels.RightAnchoredSupport
InterpolationKernels.infimum
InterpolationKernels.supremum
```

## Interpolation kernels

```@docs
Kernel
InterpolationKernels.iscardinal
InterpolationKernels.isnormalized
InterpolationKernels.support
InterpolationKernels.values
InterpolationKernels.brief
InterpolationKernels.compute_weights
InterpolationKernels.compute_offset_and_weights
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
