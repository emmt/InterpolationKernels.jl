#
# inline-issue-solo.jl -
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

# Catmull-Rom spline (inlined):
@inline function inlined_spline(x::T) where {T<:AbstractFloat}
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

# Idem but force in-lining.
@inline function compute_weights(::typeof(inlined_spline),
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

function simple_map!(f::Function,
                     dst::Array{T,N},
                     src::Array{T,N}) where {T<:AbstractFloat,N}
    for i in eachindex(dst, src)
        dst[i] = f(src[i])
    end
    return dst
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

function simple_compute_weights!(f::Function,
                                 dst::Array{T,2},
                                 src::Array{T,1}) where {T<:AbstractFloat}
    @assert size(dst) == (4,length(src))
    for i in eachindex(src)
        w1, w2, w3, w4 = compute_weights(f, src[i])
        dst[1,i] = w1
        dst[2,i] = w2
        dst[3,i] = w3
        dst[4,i] = w4
    end
    return dst
end

@inline function inlined_simple_compute_weights!(f::Function,
                                                 dst::Array{T,2},
                                                 src::Array{T,1}) where {T<:AbstractFloat}
    @assert size(dst) == (4,length(src))
    for i in eachindex(src)
        w1, w2, w3, w4 = compute_weights(f, src[i])
        dst[1,i] = w1
        dst[2,i] = w2
        dst[3,i] = w3
        dst[4,i] = w4
    end
    return dst
end

function inbounds_compute_weights!(f::Function,
                                 dst::Array{T,2},
                                 src::Array{T,1}) where {T<:AbstractFloat}
    @assert size(dst) == (4,length(src))
    @inbounds for i in eachindex(src)
        w1, w2, w3, w4 = compute_weights(f, src[i])
        dst[1,i] = w1
        dst[2,i] = w2
        dst[3,i] = w3
        dst[4,i] = w4
    end
    return dst
end

@inline function inlined_inbounds_compute_weights!(f::Function,
                                                   dst::Array{T,2},
                                                   src::Array{T,1}) where {T<:AbstractFloat}
    @assert size(dst) == (4,length(src))
    @inbounds for i in eachindex(src)
        w1, w2, w3, w4 = compute_weights(f, src[i])
        dst[1,i] = w1
        dst[2,i] = w2
        dst[3,i] = w3
        dst[4,i] = w4
    end
    return dst
end

function simd_compute_weights!(f::Function,
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

@inline function inlined_simd_compute_weights!(f::Function,
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
    print(" │   ├─ simple:  "); @btime $(simple_map!)($spline, $y, $x)
    print(" │   ├─ inbounds:"); @btime $(inbounds_map!)($spline, $y, $x)
    print(" │   └─ simd:    "); @btime $(simd_map!)($spline, $y, $x)
    print(" ├─ Call inlined_spline ($n times):\n");
    print(" │   ├─ simple:  "); @btime $(simple_map!)($inlined_spline, $y, $x)
    print(" │   ├─ inbounds:"); @btime $(inbounds_map!)($inlined_spline, $y, $x)
    print(" │   └─ simd:    "); @btime $(simd_map!)($inlined_spline, $y, $x)
    print(" │\n");
    print(" ├─ Computation of weights with spline ($n times):\n");
    print(" │   ├─ simple:  "); @btime $(simple_compute_weights!)($spline, $z, $t)
    print(" │   ├─ inbounds:"); @btime $(inbounds_compute_weights!)($spline, $z, $t)
    print(" │   └─ simd:    "); @btime $(simd_compute_weights!)($spline, $z, $t)
    print(" ├─ Computation of weights with inlined_spline ($n times):\n");
    print(" │   ├─ simple:  "); @btime $(simple_compute_weights!)($inlined_spline, $z, $t)
    print(" │   ├─ inbounds:"); @btime $(inbounds_compute_weights!)($inlined_spline, $z, $t)
    print(" │   └─ simd:    "); @btime $(simd_compute_weights!)($inlined_spline, $z, $t)
    print(" ├─ Inlined computation of weights with spline ($n times):\n");
    print(" │   ├─ simple:  "); @btime $(inlined_simple_compute_weights!)($spline, $z, $t)
    print(" │   ├─ inbounds:"); @btime $(inlined_inbounds_compute_weights!)($spline, $z, $t)
    print(" │   └─ simd:    "); @btime $(inlined_simd_compute_weights!)($spline, $z, $t)
    print(" └─ Inlined computation of weights with inlined_spline ($n times):\n");
    print("     ├─ simple:  "); @btime $(inlined_simple_compute_weights!)($inlined_spline, $z, $t)
    print("     ├─ inbounds:"); @btime $(inlined_inbounds_compute_weights!)($inlined_spline, $z, $t)
    print("     └─ simd:    "); @btime $(inlined_simd_compute_weights!)($inlined_spline, $z, $t)

end

runtests(T=Float32)
println()
runtests(T=Float64)

end # module
