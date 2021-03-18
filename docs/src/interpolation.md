# Interpolation

Interpolation kernels, as their name suggest, are designed for interpolating
arrays.  The `InterpolationKernels` package provides the method
[`InterpolationKernels.compute_offset_and_weights`](@ref) to efficiently
compute interpolation weights for a given kernel and interpolating position.


## Interpolation principles

To explain how linear interpolation works, let us assume that we want to
interpolate a source vector `A ∈ 𝕂ⁿ` using the kernel `h: ℝ → ℝ` to produce a
continuous model function `f(x)` with `x ∈ ℝ` the continuous coordinate.
Such a model writes:

```
∀ x ∈ ℝ:    f(x) = sum_{k ∈ ℤ} h(x - k)*ℬ(A,k)
```

where `sum_{k ∈ ℤ}` denotes a sum over all integers and `ℬ: 𝕂ⁿ×ℤ → 𝕂`
implements the boundary conditions with `Sup(A)` the set of valid indices for
`A`.  In Julia code, `Sup(A) = axes(A,1)` and `n = length(A)`.  The above
equation is a **discrete convolution** of `A` by the kernel `h`.  The model
function is `f: ℝ → 𝕂` where `𝕂` is the field of the values taken by the
product in the above sum.  We make no restrictions for the boundary conditions
except that the following always holds:

```
∀ k ∈ Sup(A):    ℬ(A,k) = A[k]
```

For instance, **flat boundaries** are implemented by:

```
          / A[min(Sup(A))]  if k ≤ min(Sup(A))
ℬ(A,k) = |  A[max(Sup(A))]  if k ≥ max(Sup(A))
          \ A[k]            else
```

Provided the kernel `h` has a finite size support, we can assume (without
loss of generality) that the kernel support `Sup(h)` is enclosed in a
right-open interval of nonzero integer width `s ∈ ℕ \ {0}`.  That is, there
exists `a ∈ ℝ` such that `Sup(h) ⊂ [a,b)` with `b = a + s`.  In other words:

```
x < a   or   x ≥ b    ⟹     h(x) = 0
```

The indices `k ∈ ℤ` such that `h(x-k)` is not certainly zero are:

```
k ∈ ℤ, x-k ∈ [a,b)
⟺    k ∈ ℤ ∩ (x-b, x-a]
⟺    k ∈ ⟦⌊x-b+1⌋, ⌊x-a⌋⟧
```

with `⌊…⌋` the floor function.  Since the objective is to have the smallest
interval to restrict the number of terms in the discrete convolution, it may be
better to enclose the kernel support in a left-open interval, that is `Sup(h) ⊂
(a,b]`, the set of indices `k` such that `h(x-k)` is not certainly zero is then
given by:

```
k ∈ ℤ, x-k ∈ (a,b]
⟺    k ∈ ℤ ∩ [x-b, x-a)
⟺    k ∈ ⟦⌈x-b⌉, ⌈x-a-1⌉⟧
```

with `⌈…⌉` the ceil function.  To summarize, the discrete convolution in the
model `f(x)` can be limited to indices `k ∈ ⟦kmin(x),kmax(x)⟧` where:

```
kmin(x) = ⌊x-b+1⌋    if Sup(h) ⊂ [a,b)
          ⌈x-b⌉      if Sup(h) ⊂ (a,b]
```

and `kmax(x) = kmin(x) + s - 1` (either of these above expressions can be
chosen if `Sup(h) ⊂ (a,b)`).  The interpolation formula can finally be
rewritten as:

```
∀ x ∈ ℝ:    f(x) = sum_{j ∈ jmin:jmax} w_j(x)*ℬ(A, l(x) + j)
```

where `k = l(x) + j` and `w_j(x) = h(x-k)` for `j ∈ ⟦jmin,jmax⟧` account for
all the non-zero terms of the sum in the original formula.  Hence `jmax =
jmin+s-1` and:

```
  l(x) = kmin(x) - jmin
w_j(x) = h(x - l(x) - j)
```

Note that `jmin ∈ ℤ` can be chosen as is the most convenient.  For instance,
assuming Julia or Fortran indexing, one would choose `jmin = 1` and thus `l(x)
= ⌊x-b⌋` if `Sup(h) ⊂ [a,b)`.


## Methods and structures

### General methods

The above developments suggest that, for any interpolation kernel `h` and
interpolation coordinate `x`, we mostly need a function that yields the
**offset** `l(x)` and the **interpolation weights** `w_j(x)` for all `j ∈
⟦jmin,jmax⟧`.  This is exactly what is done by the method
`compute_offset_and_weights` provided by the `InterpolationKernels` package.
This method is called as:

```julia
off, wgt = compute_offset_and_weights(h, x)
```

with `h` the interpolation kernel, `x` the coordinate, `off = l(x)` the offset
and `wgt` the `s`-tuple of interpolation weights given by `wgt[j] = w_j(x)`.
In Julia, the first index of tuples is `1`, so we have `jmin = 1` and thus:

```
off = kmin(x) - 1 = ⌊x-b⌋          if Sup(h) ⊂ [a,b)
                    ⌈x-b⌉ - 1      if Sup(h) ⊂ (a,b]
wgt[j] = h(x - off - j)            for j = 1, 2, ..., s
```

In order to compute the offset, the coordinate `x` but also the support of the
kernel `h` must be known.  These are the reasons to have interpolation kernels
in `InterpolationKernels` be defined as smart functions that know their
support.

Assuming the following methods are available (they are indeed defined in
`InterpolationKernels` but not exported):

```julia
support(ker) ---> sup # the kernel support `Sup(ker)`
length(sup) ----> s   # the integer size of the kernel support
infimum(sup) ---> a   # the least upper bound of the support `sup`
supremum(sup) --> b   # the greatest lower bound of the support `sup`
```

the offset `off` can be computed by (there are two cases to consider):

```
offset(sup::Support{T,S,<:Bound,Open}, x::T) where {T,S} =
    floor(x - supremum(sup))       #  if Sup(ker) is ⊂ [a,b)
offset(sup::Support{T,S,Open,Closed}, x::T) where {T,S} =
    ceil(x - (supremum(sup) - 1))  #  if Sup(ker) is ⊂ (a,b]
```

and a generic implementation of `compute_offset_and_weights` is given by:

```julia
compute_offset_and_weights(ker::Kernel{T}, x::T) where {T} =
    compute_offset_and_weights(support(ker), ker, x)

function compute_offset_and_weights(sup::Support{T,S},
                                    ker::Kernel{T,S}, x::T) where {T,S}
    off = offset(sup, x)
    return off, ntuple(j -> ker(x - off - j), Val(length(ker)))
end
```

Of course `compute_offset_and_weights` is optimized for the most popular
kernels in order to reduce the number of operations.  But the generic code
above gives you the ideas.


### Symmetric supports

For kernels with symmetric support, the following hold `a = -b` and `b = s/2`.
The interpolation weights can be then expressed as:

```
wgt[j] = h(t - j)        for j ∈ 1-m:s-m
```

with `m = (s+1)>>1` (that is the integer division of `s+1` by `2`) and:

```
t = x - floor(x)      if s is even
    x - round(x)      if s is odd
```

and the offset is given by:

```
off = floor(x) - m      if s is even
      round(x) - m      if s is odd
```

These expressions are used by `compute_offset_and_weights` for kernels with
symmetric support because they may help reducing the number of operations (at
least for splines) for computing interpolation weights.  In that case, the
weights are computed by calling `compute_weights` with the kernel and the
computed value of `t` (not `x`):

```julia
wgt = compute_weights(k, t)
```

When a new kernel is implemented, `compute_offset_and_weights` or, if the
kernel has a symmetric support, `compute_weights` may be specialized to
optimize computations.



## Example: fine shifting

Let us now assume that we want to compute:

```
C[i] ≈ A[i - r]
```

for all indices `i ∈ Sub(C)` of the destination vector `C` and some non-integer
offset `r` where `≈` denotes the approximation by the interpolation model
`f(x)` described above.  Hence `C` is the result of performing a sub-sample
shift of `A` by offset `r`.

Combining equations (that is just replace `x` by `i-r` in the interpolation
formula) yields:

```
C[i] = f(i - r)
     = sum_{j ∈ 1:s} h(i - r - l(i - r) - j)*ℬ(A, l(i - r) + j)
```

From the definition of the offset `l(x)` it is obvious that:

```
∀ (r,i) ∈ ℝ×ℤ:    l(i - r) = i + l(-r) = i - k
```

with:

```
k = -l(-r)
```

Hence the interpolation writes:

```
C[i] = sum_{j ∈ 1:s} h(k - r - j)*ℬ(A, i - k + j)
```

But as `v = k - r` is a constant that does not depend on `i`, the weights:

```
W[j] = h(k - r - j) = h(v - j)
```

can be computed once for all indices `i` in the destination `B` to perform fine
shifting by the following formula:

```
C[i] = sum_{j ∈ 1:s} W[j]*ℬ(A, i - k + j)
```

which is a discrete correlation of `W` and `A`.  The weights `W` and the
integer offset `k` do not depend on `i`, resulting in very fast computations.
This is exploited by the [FineShift](https://github.com/emmt/FineShift.jl)
package.
