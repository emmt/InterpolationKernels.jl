module BenchmarkingInterpolationKernels

using BenchmarkTools, Statistics, Printf
using InterpolationKernels
using InterpolationKernels:
    compute_weights, compute_offset_and_weights,
    brief, support, infimum, supremum

function compute_weights!(y::Array{T,2},
                          ker::Kernel{T,1},
                          x::Array{T,1}) where {T}
    @assert size(y) == (1,length(x))
    @inbounds @simd for i in eachindex(x)
        wgts = compute_weights(ker, x[i])
        y[1,i] = wgts[1]
    end
    nothing
end

function compute_weights!(y::Array{T,2},
                          ker::Kernel{T,2},
                          x::Array{T,1}) where {T}
    @assert size(y) == (2,length(x))
    @inbounds @simd for i in eachindex(x)
        wgts = compute_weights(ker, x[i])
        y[1,i] = wgts[1]
        y[2,i] = wgts[2]
    end
    nothing
end

function compute_weights!(y::Array{T,2},
                          ker::Kernel{T,3},
                          x::Array{T,1}) where {T}
    @assert size(y) == (3,length(x))
    @inbounds @simd for i in eachindex(x)
        wgts = compute_weights(ker, x[i])
        y[1,i] = wgts[1]
        y[2,i] = wgts[2]
        y[3,i] = wgts[3]
    end
    nothing
end

function compute_weights!(y::Array{T,2},
                          ker::Kernel{T,4},
                          x::Array{T,1}) where {T}
    @assert size(y) == (4,length(x))
    @inbounds @simd for i in eachindex(x)
        wgts = compute_weights(ker, x[i])
        y[1,i] = wgts[1]
        y[2,i] = wgts[2]
        y[3,i] = wgts[3]
        y[4,i] = wgts[4]
    end
    nothing
end

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

function simdmap!(dst::AbstractArray{<:Any,N},
                  f,
                  src::AbstractArray{<:Any,N}) where {N}
    @inbounds @simd for i in eachindex(dst, src)
        dst[i] = f(src[i])
    end
    return dst
end

number_of_operations(f, ker::Kernel) = number_of_operations(f, typeof(ker))
number_of_operations(ker::Kernel) = number_of_operations(ker, ker)

# For calling the function on a single value, the number of operations is the
# maximum possible.  In case a constant is returned, the number of operations is
# assumed to be one (not zero).  Function compute_offset_and_weights usually adds 3
# operations compared to compute_weights.
for (K, n1, n2, n3) in ((BSpline{1},               2,  2,  2),
                        (BSplinePrime{1},          1,  1,  2), # n1 should be 0
                        (BSpline{2},               2,  1,  4),
                        (BSplinePrime{2},          2,  1,  4),
                        (BSpline{3},               3,  6,  9),
                        (BSplinePrime{3},          2,  3,  6),
                        (BSpline{4},               6, 15, 18),
                        (BSplinePrime{4},          5, 13, 16),
                        (CubicSpline,              6, 15, 18),
                        (CubicSplinePrime,         6, 11, 14),
                        (CardinalCubicSpline,      6, 11, 14),
                        (CardinalCubicSplinePrime, 6, 24, 27), # weights not optimized
                        (CatmullRomSpline,         6, 10, 13),
                        (CatmullRomSplinePrime,    5, 11, 14),)
    @eval begin
        number_of_operations(::$K, ::Type{<:$K}) = $n1
        number_of_operations(::typeof(compute_weights),::Type{<:$K}) = $n2
        number_of_operations(::typeof(compute_offset_and_weights),::Type{<:$K}) = $n3
    end
end

gflops(nops::Real, ns::Real) = (ns, nops/ns)

function runtests(n::Int = 1000; table::Bool=false, verbose::Bool=false, io::IO=stdout)
    # Benchmark calling kernels.
    if table
        print(io, "Benchmark calling kernels ",
              "(minimum times in nanoseconds for $n calls):\n")
        print(io, "| Julia version                | $VERSION |\n")
        print(io, "|:-----------------------------|------:|\n")
    end
    for T in (Float32, Float64)
        for ker in (BSpline{1,T}(),
                    BSpline{2,T}(),
                    BSpline{3,T}(),
                    BSpline{4,T}(),
                    CubicSpline{T}(-1/2,1/6),
                    CatmullRomSpline{T}(),
                    CardinalCubicSpline{T}(-1/2),
                    MitchellNetravaliSpline{T}(1/3,1/3),)
            xmin = infimum(support(ker))
            xmax = supremum(support(ker))
            x = (xmax - xmin)*rand(T, n) .+ xmin;
            @assert eltype(x) === T
            y = Array{T}(undef, n)
            nops = number_of_operations(ker)*n
            if !table
                print(io, "Testing calling ", summary(ker), " ", n, " times: ")
                flush(io)
            end
            b = @benchmark $(simdmap!)($y, $ker, $x)
            t = b.times
            if table
                @printf(io, "| %-28s | %5.0f |\n", summary(ker), minimum(t))
            else
                @printf(io, "%.3f Gflops\n", nops/minimum(t))
                if verbose
                    @printf(io, "  - memory: %d allocations, %d bytes\n",
                            b.allocs, b.memory)
                    @printf(io, "  - minimum time: %7.0f ns, %.3f Gflops\n",
                            gflops(nops, minimum(t))...)
                    @printf(io, "  - median time:  %7.0f ns, %.3f Gflops\n",
                            gflops(nops, median(t))...)
                    @printf(io, "  - mean time:    %7.0f ns, %.3f Gflops\n",
                            gflops(nops, mean(t))...)
                    @printf(io, "  - maximum time: %7.0f ns, %.3f Gflops\n",
                            gflops(nops, maximum(t))...)
                    println(io)
                end
            end
        end
    end

    # Benchmark compute_weights.
    if table
        print(io, "\nBenchmarking `compute_weights` ",
              "(minimum times in nanoseconds for $n calls):\n")
        print(io, "| Julia version                | $VERSION |\n")
        print(io, "|:-----------------------------|------:|\n")
    end
    for T in (Float32, Float64)
        for ker in (BSpline{1,T}(),
                    BSpline{2,T}(),
                    BSpline{3,T}(),
                    BSpline{4,T}(),
                    CubicSpline{T}(-1/2,1/6),
                    CatmullRomSpline{T}(),
                    CardinalCubicSpline{T}(-1/2),
                    MitchellNetravaliSpline{T}(1/3,1/3),)
            if isodd(length(ker))
                x = rand(T, n) .- one(T)/2
            else
                x = rand(T, n)
            end
            @assert eltype(x) === T
            z = Array{T}(undef, length(ker), n)
            nops = number_of_operations(compute_weights, ker)*n
            if !table
                print(io, "Testing `compute_weights` for ",
                      summary(ker), " ", n, " times: ")
                flush(io)
            end
            b = @benchmark $(compute_weights!)($z, $ker, $x)
            t = b.times
            if table
                @printf(io, "| %-28s | %5.0f |\n", summary(ker), minimum(t))
            else
                @printf(io, "%.3f Gflops\n", nops/minimum(t))
                if verbose
                    @printf(io, "  - memory: %d allocations, %d bytes\n",
                            b.allocs, b.memory)
                    @printf(io, "  - minimum time: %7.0f ns, %.3f Gflops\n",
                            gflops(nops, minimum(t))...)
                    @printf(io, "  - median time:  %7.0f ns, %.3f Gflops\n",
                            gflops(nops, median(t))...)
                    @printf(io, "  - mean time:    %7.0f ns, %.3f Gflops\n",
                            gflops(nops, mean(t))...)
                    @printf(io, "  - maximum time: %7.0f ns, %.3f Gflops\n",
                            gflops(nops, maximum(t))...)
                    println(io)
                end
            end
        end
    end

    # Benchmark compute_offset_and_weights.
    if table
        print(io, "\nBenchmarking `compute_offset_and_weights` ",
              "(minimum times in nanoseconds for $n calls):\n")
        print(io, "| Julia version                | $VERSION |\n")
        print(io, "|:-----------------------------|------:|\n")
    end
    for T in (Float32, Float64)
        x = 200*rand(T, n) .- 100;
        for ker in (BSpline{1,T}(),
                    BSpline{2,T}(),
                    BSpline{3,T}(),
                    BSpline{4,T}(),
                    CubicSpline{T}(-1/2,1/6),
                    CatmullRomSpline{T}(),
                    CardinalCubicSpline{T}(-1/2),
                    MitchellNetravaliSpline{T}(1/3,1/3),)
            nops = number_of_operations(compute_offset_and_weights, ker)*n
            if !table
                print(io, "Testing `compute_offset_and_weights` for ",
                      summary(ker), " ", n, " times: ")
                flush(io)
            end
            z = Array{T}(undef, 1 + length(ker), n)
            b = @benchmark $(compute_offset_and_weights!)($z, $ker, $x)
            t = b.times
            if table
                @printf(io, "| %-28s | %5.0f |\n", summary(ker), minimum(t))
            else
                @printf(io, "%.3f Gflops\n", nops/minimum(t))
                if verbose
                    @printf(io, "  - memory: %d allocations, %d bytes\n",
                            b.allocs, b.memory)
                    @printf(io, "  - minimum time: %7.0f ns, %.3f Gflops\n",
                            gflops(nops, minimum(t))...)
                    @printf(io, "  - median time:  %7.0f ns, %.3f Gflops\n",
                            gflops(nops, median(t))...)
                    @printf(io, "  - mean time:    %7.0f ns, %.3f Gflops\n",
                            gflops(nops, mean(t))...)
                    @printf(io, "  - maximum time: %7.0f ns, %.3f Gflops\n",
                            gflops(nops, maximum(t))...)
                    println(io)
                end
            end
        end
    end
end

end # module

if !isinteractive()
    let verbose = false,
        table=false,
        n = 1000

        option = true
        for i in 1:length(ARGS)
            arg = ARGS[i]
            if option && (arg == "-h" || arg == "--help")
                println("Usage:  julia -O3 $(@__FILE__) [-h|--help] [-v] [--] [NUMBER]")
                println("Options:  julia -O3 SCRIPT [-h|--help] [-v] [--] [NUMBER]")
                println("    -h|--help      Print this help.")
                println("    -v|--verbose   Print detailed results.")
                println("    -t|--table     Print results as a table.")
                println("    NUMBER         Vector size (default: $n")
                exit(0)
            elseif option && (arg == "-v" || arg == "--verbose")
                verbose = true
                continue
            elseif option && (arg == "-t" || arg == "--table")
                table = true
                continue
            elseif option && (arg == "--")
                option = false
                continue
            end
            if i == length(ARGS)
                n = tryparse(Int, ARGS[i], base=10)
                if n === nothing
                    println(stderr, "invalid number of elements")
                    exit(1)
                end
            else
                println(stderr, "too many arguments")
                exit(1)
            end
        end
        BenchmarkingInterpolationKernels.runtests(n, verbose=verbose, table=table)
    end
end
nothing
