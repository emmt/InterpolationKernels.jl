# Boundary conditions

In `InterpolationKernels`, the boundary conditions indicate how missing values
should be extrapolated.

The boundary conditions are specified for convenience and to provide a unified
behavior.  The boundary conditions are not directly used by the implemented
kernels but are part of the kernel types and should be taken into account by
interpolation methods using the kernel instances.

Calling the [`boundaries`](@ref) method on a kernel (or on a kernel type) yields
the type of boundary conditions associated to the kernel.
