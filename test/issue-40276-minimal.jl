#
# inline-issue-minimal.jl -
#
# Demonstrate a reduction of performances with Julia-1.6.0 compared to 1.5.4.
# This version is similar to `inline-issue.jl` but does not use
# `MayOptimize.jl`.
#
module InlineIssue

using BenchmarkTools

# Catmull-Rom spline:
function spline(x::T) where {T<:AbstractFloat}
    abs_x = abs(x)
    if abs_x ≥ 2
        zero(T)
    elseif abs_x ≤ 1 # (((3/2)*x - (5/2))*x^2 + 1)
        ((T(3)/2)*abs_x - (T(5)/2))*abs_x^2 + 1
    else # (((5/2) - (1/2)*x)*x - 4)*x + 2 = (2 - x)^2 (1 - x)/2
        (2 - abs_x)^2*(1 - abs_x)/2
    end
end

# Compute all weights (10 operations).
function compute_weights(::typeof(spline),
                         t::T) where {T<:AbstractFloat}
    u = 1 - t
    v = (T(-1)/2)*t*u
    w1 = v*u
    w4 = v*t
    w = w4 - w1
    w2 = u - w1 + w
    w3 = t - w4 - w
    return (w1, w2, w3, w4)
end

function inbounds_map!(f::Function,
                       dst::Array{T,N},
                       src::Array{T,N}) where {T<:AbstractFloat,N}
    @inbounds for i in eachindex(dst, src)
        dst[i] = f(src[i])
    end
    return dst
end

function simd_map!(f::Function,
                   dst::Array{T,N},
                   src::Array{T,N}) where {T<:AbstractFloat,N}
    @inbounds @simd for i in eachindex(dst, src)
        dst[i] = f(src[i])
    end
    return dst
end

function compute_weights!(f::Function,
                          dst::Array{T,2},
                          src::Array{T,1}) where {T<:AbstractFloat}
    @assert size(dst) == (4,length(src))
    @inbounds @simd for i in eachindex(src)
        w1, w2, w3, w4 = compute_weights(f, src[i])
        dst[1,i] = w1
        dst[2,i] = w2
        dst[3,i] = w3
        dst[4,i] = w4
    end
    return dst
end

@inline function inlined_compute_weights!(f::Function,
                                          dst::Array{T,2},
                                          src::Array{T,1}) where {T<:AbstractFloat}
    @assert size(dst) == (4,length(src))
    @inbounds @simd for i in eachindex(src)
        w1, w2, w3, w4 = compute_weights(f, src[i])
        dst[1,i] = w1
        dst[2,i] = w2
        dst[3,i] = w3
        dst[4,i] = w4
    end
    return dst
end

function runtests(; T::Type{<:AbstractFloat}=Float32, n::Int=1_000)
    x = rand(T, n)
    y = Array{T}(undef, n)
    t = rand(T, n) .- one(T)/2
    z = Array{T}(undef, 4, n)
    print("Tests with Julia-", VERSION, ", T=", T, ", n=", n, "\n")
    print(" ├─ Call spline ($n times):\n");
    print(" │   ├─ map! ───────────────────────"); @btime $(map!)($spline, $y, $x)
    print(" │   ├─ inbounds_map! ──────────────"); @btime $(inbounds_map!)($spline, $y, $x)
    print(" │   └─ simd_map! ──────────────────"); @btime $(simd_map!)($spline, $y, $x)
    print(" └─ Computation of weights with spline ($n times):\n");
    print("     ├─ compute_weights! ───────────"); @btime $(compute_weights!)($spline, $z, $t)
    print("     └─ inlined_compute_weights! ───"); @btime $(inlined_compute_weights!)($spline, $z, $t)

end

runtests(T=Float32)
#println()
#runtests(T=Float64)

end # module
