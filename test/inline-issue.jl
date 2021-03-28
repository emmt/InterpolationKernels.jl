module InlineIssue

using BenchmarkTools, MayOptimize

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

function Base.map!(opt::Type{<:OptimLevel},
                   f::Function,
                   dst::Array{T,N},
                   src::Array{T,N}) where {T<:AbstractFloat,N}
    @maybe_vectorized opt for i in eachindex(dst, src)
        dst[i] = f(src[i])
    end
    return dst
end

function compute_weights!(opt::Type{<:OptimLevel},
                          f::Function,
                          dst::Array{T,2},
                          src::Array{T,1}) where {T<:AbstractFloat}
    @assert size(dst) == (4,length(src))
    @maybe_vectorized opt for i in eachindex(src)
        wgts = compute_weights(f, src[i])
        dst[1,i] = wgts[1]
        dst[2,i] = wgts[2]
        dst[3,i] = wgts[3]
        dst[4,i] = wgts[4]
    end
    return dst
end

@inline function inlined_compute_weights!(opt::Type{<:OptimLevel},
                                          f::Function,
                                          dst::Array{T,2},
                                          src::Array{T,1}) where {T<:AbstractFloat}
    @assert size(dst) == (4,length(src))
    @maybe_vectorized opt for i in eachindex(src)
        wgts = compute_weights(f, src[i])
        dst[1,i] = wgts[1]
        dst[2,i] = wgts[2]
        dst[3,i] = wgts[3]
        dst[4,i] = wgts[4]
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
    print(" │   ├─ Debug:    "); @btime $(map!)($(Debug), $spline, $y, $x)
    print(" │   ├─ InBounds: "); @btime $(map!)($(InBounds), $spline, $y, $x)
    print(" │   └─ Vectorize:"); @btime $(map!)($(Vectorize), $spline, $y, $x)
    print(" ├─ Call inlined_spline ($n times):\n");
    print(" │   ├─ Debug:    "); @btime $(map!)($(Debug), $inlined_spline, $y, $x)
    print(" │   ├─ InBounds: "); @btime $(map!)($(InBounds), $inlined_spline, $y, $x)
    print(" │   └─ Vectorize:"); @btime $(map!)($(Vectorize), $inlined_spline, $y, $x)
    print(" ├─ Computation of weights with spline ($n times):\n");
    print(" │   ├─ Debug:    "); @btime $(compute_weights!)($(Debug), $spline, $z, $t)
    print(" │   ├─ InBounds: "); @btime $(compute_weights!)($(InBounds), $spline, $z, $t)
    print(" │   └─ Vectorize:"); @btime $(compute_weights!)($(Vectorize), $spline, $z, $t)
    print(" ├─ Computation of weights with inlined_spline ($n times):\n");
    print(" │   ├─ Debug:    "); @btime $(compute_weights!)($(Debug), $inlined_spline, $z, $t)
    print(" │   ├─ InBounds: "); @btime $(compute_weights!)($(InBounds), $inlined_spline, $z, $t)
    print(" │   └─ Vectorize:"); @btime $(compute_weights!)($(Vectorize), $inlined_spline, $z, $t)
    print(" ├─ Inlined computation of weights with spline ($n times):\n");
    print(" │   ├─ Debug:    "); @btime $(inlined_compute_weights!)($(Debug), $spline, $z, $t)
    print(" │   ├─ InBounds: "); @btime $(inlined_compute_weights!)($(InBounds), $spline, $z, $t)
    print(" │   └─ Vectorize:"); @btime $(inlined_compute_weights!)($(Vectorize), $spline, $z, $t)
    print(" └─ Inlined computation of weights with inlined_spline ($n times):\n");
    print("     ├─ Debug:    "); @btime $(inlined_compute_weights!)($(Debug), $inlined_spline, $z, $t)
    print("     ├─ InBounds: "); @btime $(inlined_compute_weights!)($(InBounds), $inlined_spline, $z, $t)
    print("     └─ Vectorize:"); @btime $(inlined_compute_weights!)($(Vectorize), $inlined_spline, $z, $t)

end

runtests(T=Float32)
runtests(T=Float64)

end # module
