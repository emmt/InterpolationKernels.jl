# Installation

`InterpolationKernels` is not yet an [offical Julia
package](https://pkg.julialang.org/) but its installation can be as easy as:

```julia
… pkg> add https://github.com/emmt/InterpolationKernels.jl
```

where `… pkg>` represents the package manager prompt (the ellipsis `…` denote
your current environment).  To start Julia's package manager, launch Julia and,
at the [REPL of
Julia](https://docs.julialang.org/en/stable/manual/interacting-with-julia/),
hit the `]` key; you should get the above `… pkg>` prompt.  To revert to
Julia's REPL, just hit the `Backspace` key at the `… pkg>` prompt.

More detailed explanations are given below.


## Using the package manager

`InterpolationKernels` can be installed by Julia's package manager using
`https` protocol:

```julia
… pkg> add https://github.com/emmt/InterpolationKernels.jl
```

or `ssh` protocol:

```julia
… pkg> add git@github.com:emmt/InterpolationKernels.jl
```

To check whether the `InterpolationKernels` package works correctly, type:

```julia
… pkg> test InterpolationKernels
```

Later, to update to the last version (and run tests), you can type:

```julia
… pkg> update InterpolationKernels
… pkg> test InterpolationKernels
```

If something goes wrong, it may be because you already have an old version of
`InterpolationKernels`.  Uninstall `InterpolationKernels` as follows:

```julia
… pkg> rm InterpolationKernels
… pkg> gc
… pkg> add https://github.com/emmt/InterpolationKernels.jl
```

before re-installing.


## Installation in scripts

To install `InterpolationKernels` in a Julia script, write:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/emmt/InterpolationKernels.jl",
                    rev="master"));
```

or with `url="git@github.com:emmt/InterpolationKernels.jl"` if you want to use `ssh`.

This also works from the Julia REPL.
