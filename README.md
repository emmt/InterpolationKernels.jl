# A library of interpolation kernels for Julia

[![Doc](https://img.shields.io/badge/docs-dev-blue.svg)](https://emmt.github.io/InterpolationKernels.jl/dev)
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](./LICENSE.md)
[![Build Status](https://github.com/emmt/InterpolationKernels.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/emmt/InterpolationKernels.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/InterpolationKernels.jl?branch=master)](https://ci.appveyor.com/project/emmt/InterpolationKernels-jl/branch/master)
[![Coverage](https://codecov.io/github/emmt/InterpolationKernels.jl/graph/badge.svg?token=ygyMUzFwXb)](https://codecov.io/github/emmt/InterpolationKernels.jl)
[![deps](https://juliahub.com/docs/General/InterpolationKernels/stable/deps.svg)](https://juliahub.com/ui/Packages/General/InterpolationKernels?t=2)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

The `InterpolationKernels` package provides a library of interpolation kernels for
[Julia](https://julialang.org/). As suggested by their name, interpolations kernels are
mostly designed for interpolation, they can be thought as smart functions that know their
derivative and that implement optimized computation of interpolation weights.

You may have a look at some sections of the [the latest
documentation](https://emmt.github.io/InterpolationKernels.jl/dev/):

- [List of provided kernels](https://emmt.github.io/InterpolationKernels.jl/dev/kernels/).
- [Principles of interpolation](https://emmt.github.io/InterpolationKernels.jl/dev/interpolation/).

The `InterpolationKernels` package is used by:

- [FineShift](https://github.com/emmt/FineShift.jl), a Julia package for fast sub-sample
  shifting of arrays.

- [LinearInterpolators](https://github.com/emmt/LinearInterpolators.jl) which implements
  various interpolation methods for Julia as linear mappings.
