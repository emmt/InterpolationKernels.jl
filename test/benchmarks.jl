module BenchmarkingInterpolationKernels

using BenchmarkTools, Statistics, Printf
using InterpolationKernels
using InterpolationKernels: compute_weights, compute_offset_and_weights, brief

function compute_offset_and_weights!(y::Array{T,2},
                   ker::Kernel{T,1},
                   x::Array{T,1}) where {T}
    @assert size(y) == (2,length(x))
    @inbounds @simd for i in eachindex(x)
        off, wgts = compute_offset_and_weights(ker, x[i])
        y[1,i] = off
        y[2,i] = wgts[1]
    end
    nothing
end

function compute_offset_and_weights!(y::Array{T,2},
                   ker::Kernel{T,2},
                   x::Array{T,1}) where {T}
    @assert size(y) == (3,length(x))
    @inbounds @simd for i in eachindex(x)
        off, wgts = compute_offset_and_weights(ker, x[i])
        y[1,i] = off
        y[2,i] = wgts[1]
        y[3,i] = wgts[2]
    end
    nothing
end

function compute_offset_and_weights!(y::Array{T,2},
                   ker::Kernel{T,3},
                   x::Array{T,1}) where {T}
    @assert size(y) == (4,length(x))
    @inbounds @simd for i in eachindex(x)
        off, wgts = compute_offset_and_weights(ker, x[i])
        y[1,i] = off
        y[2,i] = wgts[1]
        y[3,i] = wgts[2]
        y[4,i] = wgts[3]
    end
    nothing
end

function compute_offset_and_weights!(y::Array{T,2},
                   ker::Kernel{T,4},
                   x::Array{T,1}) where {T}
    @assert size(y) == (5,length(x))
    @inbounds @simd for i in eachindex(x)
        off, wgts = compute_offset_and_weights(ker, x[i])
        y[1,i] = off
        y[2,i] = wgts[1]
        y[3,i] = wgts[2]
        y[4,i] = wgts[3]
        y[5,i] = wgts[4]
    end
    nothing
end

gflops(nops::Real, ns::Real) = (ns, nops/ns)

const n = 1000


for T in (Float32, Float64)
    x = 200*rand(T, n) .- 100;
    for (ker, nops) in ((BSpline{1,T}(), 2*n),
                        (BSpline{2,T}(), 3*n),
                        (BSpline{3,T}(), 9*n),
                        (BSpline{4,T}(), 17*n),
                        (CatmullRomSpline{T}(), 13*n),
                        (CardinalCubicSpline{T}(-1/2), 14*n),
                        (MitchellNetravaliSpline{T}(1/3,1/3), 18*n),)
        print(stdout, "\nTests for $(summary(ker)): ")
        flush(stdout)
        y = Array{T}(undef, 1 + length(ker), n)
        b = @benchmark $(compute_offset_and_weights!)($y, $ker, $x)
        t = b.times
        @printf "%.3f Gflops\n" nops/minimum(t)
        @printf "  - memory: %d allocations, %d bytes\n" b.allocs b.memory
        @printf "  - minimum time: %7.0f ns, %.3f Gflops\n" gflops(nops, minimum(t))...
        @printf "  - median time:  %7.0f ns, %.3f Gflops\n" gflops(nops, median(t))...
        @printf "  - mean time:    %7.0f ns, %.3f Gflops\n" gflops(nops, mean(t))...
        @printf "  - maximum time: %7.0f ns, %.3f Gflops\n" gflops(nops, maximum(t))...
    end
end

end # module
