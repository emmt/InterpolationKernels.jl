# A library of interpolation kernels for Julia

| **Documentation**               | **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][doc-dev-img]][doc-dev-url] | [![][license-img]][license-url] | [![][github-ci-img]][github-ci-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |

The `InterpolationKernels` package provides a library of interpolation kernels
for [Julia](https://julialang.org/).  As suggested by their name,
interpolations kernels are mostly designed for interpolation, they can be
thought as smart functions that know their derivative and that implement
optimized computation of interpolation weights.

You may have a look at some sections of the [the latest
documentation](https://emmt.github.io/InterpolationKernels.jl/dev/):

- [List of provided kernels](https://emmt.github.io/InterpolationKernels.jl/dev/kernels/).
- [Principles of interpolation](https://emmt.github.io/InterpolationKernels.jl/dev/interpolation/).
- [Instructions for installation](https://emmt.github.io/InterpolationKernels.jl/dev/install/).

The `InterpolationKernels` package is used by:

- [FineShift](https://github.com/emmt/FineShift.jl), a Julia package for
  fast sub-sample shifting of arrays.

- [LinearInterpolators](https://github.com/emmt/LinearInterpolators.jl) which implements
  various interpolation methods for Julia as linear mappings.

[doc-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[doc-stable-url]: https://emmt.github.io/InterpolationKernels.jl/stable

[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://emmt.github.io/InterpolationKernels.jl/dev

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat

[github-ci-img]: https://github.com/emmt/InterpolationKernels.jl/actions/workflows/CI.yml/badge.svg?branch=master
[github-ci-url]: https://github.com/emmt/InterpolationKernels.jl/actions/workflows/CI.yml?query=branch%3Amaster

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/emmt/InterpolationKernels.jl?branch=master
[appveyor-url]: https://ci.appveyor.com/project/emmt/InterpolationKernels-jl/branch/master

[coveralls-img]: https://coveralls.io/repos/emmt/InterpolationKernels.jl/badge.svg?branch=master&service=github
[coveralls-url]: https://coveralls.io/github/emmt/InterpolationKernels.jl?branch=master

[codecov-img]: http://codecov.io/github/emmt/InterpolationKernels.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/emmt/InterpolationKernels.jl?branch=master
