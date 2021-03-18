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

function runtests(n::Int = 1000; verbose::Bool=false, io::IO=stdout)
    for T in (Float32, Float64)
        x = 200*rand(T, n) .- 100;
        for (ker, nops) in ((BSpline{1,T}(), 2*n),
                            (BSpline{2,T}(), 3*n),
                            (BSpline{3,T}(), 9*n),
                            (BSpline{4,T}(), 17*n),
                            (CatmullRomSpline{T}(), 13*n),
                            (CardinalCubicSpline{T}(-1/2), 14*n),
                            (MitchellNetravaliSpline{T}(1/3,1/3), 18*n),)
            print(io, "Tests for $(summary(ker)): ")
            flush(io)
            y = Array{T}(undef, 1 + length(ker), n)
            b = @benchmark $(compute_offset_and_weights!)($y, $ker, $x)
            t = b.times
            @printf(io, "%.3f Gflops\n", nops/minimum(t))
            if verbose
                @printf(io, "  - memory: %d allocations, %d bytes\n", b.allocs, b.memory)
                @printf(io, "  - minimum time: %7.0f ns, %.3f Gflops\n", gflops(nops, minimum(t))...)
                @printf(io, "  - median time:  %7.0f ns, %.3f Gflops\n", gflops(nops, median(t))...)
                @printf(io, "  - mean time:    %7.0f ns, %.3f Gflops\n", gflops(nops, mean(t))...)
                @printf(io, "  - maximum time: %7.0f ns, %.3f Gflops\n", gflops(nops, maximum(t))...)
                println(io)
            end
        end
    end
end

end # module

if !isinteractive()
    let verbose = false,
        n = 1000
        i = 0
        while i < length(ARGS)
            i += 1
            arg = ARGS[i]
            if arg == "-h" || arg == "--help"
                println("Usage:  julia -O3 $(@__FILE__) [-h|--help] [-v] [--] [NUMBER]")
                println("Options:  julia -O3 SCRIPT [-h|--help] [-v] [--] [NUMBER]")
                println("    -h|--help      Print this help.")
                println("    -v|--verbose   Print detailed results.")
                println("    NUMBER         Vector size (default: $n")
                exit()
            elseif arg == "-v" || arg == "--verbose"
                verbose = true
            elseif arg == "--"
                i += 1
                break
            end
        end
        if i < length(ARGS)
            n = parse(Int, ARGS[i+1],base=10)
        end
        BenchmarkingInterpolationKernels.runtests(n, verbose=verbose)
    end
end
nothing
