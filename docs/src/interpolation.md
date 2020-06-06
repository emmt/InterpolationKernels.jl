# Interpolation

Interpolation kernels, as their name suggest, are designed for interpolating
arrays.  The `InterpolationKernels` package provides the [`getweights`](@ref)
method to compute interpolation weights for a given kernel and offset.


## Interpolation principles

To explain how linear interpolation works, let us assume that we want to
interpolate a source vector `src` using the kernel `ker: ℝ → ℝ` to produce a
continuous model function `mdl(x)` with `x ∈ ℝ` the continuous coordinate.
Such a model writes:

```
mdl(x) = sum_{j ∈ J} ker(x - j)*src[j]
```

where `sum_{j ∈ J}` denotes a sum over all indices `j` of `src`, that is `J =
axes(src,1)`.  The above equation is a **convolution** of the source `src` by
the kernel `ker`.  The model function is `mdl: ℝ → 𝕂` where `𝕂` is the type of
the product in the above sum.

If the kernel `ker` has a finite support of size `S ∈ ℕ`, there are at most `S`
nonzero `ker(…)` terms in the sum.  Moreover, if `S` is much smaller than the
length of `src` and to spare computations, it is worth rewriting the result of
the interpolation as:

```
mdl(x) = sum_{k ∈ 1:S} ker(x + p - k)*src[k - p]
```

for some well chosen `p ∈ ℤ` so that the sum over `k` accounts for all nonzero
`ker(…)` terms.  This means that if `[kmin,kmax]` is the support of `ker(x)`
(with `kmax = kmin + S`), then the following inequalities must hold:

```
kmin ≤ x + p - maximum(1:S) ≤ x + p - minimum(1:S) ≤ kmax
```

where `minimum(1:S) ≡ 1` and `maximum(1:S) ≡ S` are the minimum and the maximum
indices `k` in the sum `sum_{k ∈ 1:S}`.  Re-arranging terms yields that `pmin ≤
p ≤ pmax` must hold with:

```
pmin = kmin - x + S
pmax = kmax - x + 1
```

Note that the width of the interval `[pmin,pmax]` is `pmax - pmin = 1` since
`kmax = kmin + S`, so the interval `[pmin,pmax]` always contains at least one
integer and contains at most two integers when `pmin` (and `pmax`) is integer.
This leaves 2 possibilities for choosing `p ∈ [pmin,pmax] ∩ ℤ`:

```
p = ceil(pmax) - 1 =  ceil(kmax - x)
p = floor(pmax)    = floor(kmax - x + 1)
```

which are equivalent for most but not all values of `x` (they differ when
`kmax-x` happens to be integer).  If `ker(kmin) = ker(kmax) = 0`, which should
be the case if `ker(x)` is everywhere continuous, the choice between these two
possibilities is irrelevant for the result of the interpolation.

Now introducing the so-called **interpolation weights** given by:

```
w[k] = ker(x + p - k)
```

for `k ∈ 1:S`, the interpolation can be rewritten as

```
dst[i] = sum_{k ∈ 1:S} w[k]*src[k - p]
```

which is a **correlation** of `w` and `src`.


## Interpolation weights

In `InterpolationKernels` all kernels have a **symmetric support**, hence `kmin
= -kmax` and `kmax = S/2`.  The interpolation weights are then given by:

```
w[k] = ker(x + ceil(S/2 - x) - k)
```

for `k ∈ 1:S` and assuming that `p = ceil(kmax - x)` has been chosen among the
two possibilities.

Instead of computing the nonzero interpolation weights one by one by `S` calls
to the kernel function, it may be beneficial to compute all weights in a row
exploiting common sub-expressions.  This facility is provided by the kernels in
`InterpolationKernels` and the interpolation weights for kernel `ker` can be
computed by:

```julia
w = getweights(ker, t)
```

which yields an `S`-tuple `w` of weights with `S` the size of the support of
the kernel `ker`.  To interpolate around position `x` (in fractional index
units), the offset `t` is given by:

- if `S` is even:

  ```julia
  t = x - floor(x)
  ```

- if `S` is odd:

  ```julia
  t = x - floor(x + 1/2) = x - round(x)
  ```

The offset `t` is therefore in the range `[0,1]` for `S` even and in the range
`[-1/2,+1/2]` for `S` odd.  The bounds being inclusive or not is irrelevant for
the result (so the rounding direction does not matter).

One of the reasons of using these conventions is that the integer index:

```julia
j = (iseven(S) ? floor(Int, x) : round(Int, x))
```

such that `t = x - j` plays a central role in determining the indices of the
entries in the interpolated array involved in the interpolation formula.

For developers who would like to implement other kernels than those provided by
`InterpolationKernels`, it is necessary to derive the general formula for the
weights as functions of the offset `t`.  Recalling that `w[k] = ker(v - k)`
with `v = x + ceil(S/2 - x)`, the identity `ceil(-u) = floor(u)` for any real
`u` can be used to rewrite `v` as:

```julia
v = x + ceil(S/2 - x)
  = x - floor(x - S/2)
```

now 2 different cases must be considered depending on the parity of `S`:

- if `S` is even, then `S = 2c` with `c` integer, therefore `S/2 = c` and:

  ```julia
  v = x - floor(x - c) = x - floor(x) + c = t + c
  ```

- if `S` is odd, then `S = 2c - 1` with `c` integer, therefore `S/2 = c - 1/2`
  and:

  ```julia
  v = x - floor(x - c + 1/2) = x - floor(x + 1/2) + c = t + c
  ```

To summarize, the interpolations weights are simply given by:

```julia
w[k] = ker(t + c - k)
```

where:

```julia
c = (S + 1) ÷ 2
```

with `÷` the integer division in Julia.


## Example: fine shifting

Let us now assume that we want to compute:

```
dst[i] ≈ src[i - r]
```

for all indices `i ∈ I` of the destination vector `dst` and some non-integer
offset `r` where `≈` denotes the approximation by the interpolation model
`mdl(x)` described above.  Hence `dst` is the result of performing a sub-sample
shift of `src` by offset `r`.

Combining equations (that is just replace `x` by `i-r`) yields:

```
dst[i] = mdl(i - r)
       = sum_{j ∈ J} ker(i - r - j)*src[j]
```

In order to spare computations, we take `j = i - q + k` for some well chosen
integer `q` and rewrite the interpolation as:

```
dst[i] = sum_{k ∈ 1:R} ker(q - r - k)*src[i - q + k]
```

we have to chose `q` such that the following inequalities hold:

```
kmin ≤ q - r - S ≤ q - r - 1 ≤ kmax
```

with `[kmin,kmax]` the support of the kernel `ker` and `kmax - kmin = S` as
before.  These inequalities are equivalent to:

```
(qmin = kmin + r + S) ≤ q ≤ (qmax = kmax + r + 1)
```

as before `qmax - qmin = 1` and there are 2 possibilities for choosing `q`:

```
q = ceil(qmax) - 1 = ceil(kmax + r)
q = floor(qmax) = floor(kmax + 1 + r)
```

Taking the second choice for `q` and assuming a symmetric support (hence `kmax
= S/2`) yields:

```
q = floor(S/2 + 1 + r)
  = floor(c + r)
```

where `c = 1 + S/2`.  The interpolation weights are now given by:

```
w[k] = ker(v - k)
```

with:

```
v = q - r
  = floor(S/2 + 1 + r) - r
```

Finally the result of fine shifting writes:

```
dst[i] = sum_{k ∈ 1:R} w[k]*src[i - q + k]
```

where the weights `w` and offset `q` do not depend on `i` and can thus be
pre-computed resulting in very fast computations.  This is exploited by the
[FineShift](https://github.com/emmt/FineShift.jl) package.

If *flat* boundary conditions hold, we are assuming that `src[j] = src[1]` if
`j ≤ 1` and `src[j] = src[n]` if `j ≥ n` with `n = length(src)` the number of
samples in `src`.  Then, if the shift `r` is too large (in magnitude) all `j =
i + k - q` for `k = 1, ..., S` are below or above the limits in the
interpolation formula.  This occurs when:

```
i + k - q ≤ 1  (∀ i,k)  <=> i + S - 1 ≤ q  (∀ i)  <=> S + m - 1 ≤ q
i + k - q ≥ n  (∀ i,k)  <=>         i ≥ q  (∀ i)  <=> 1 ≥ q
```
