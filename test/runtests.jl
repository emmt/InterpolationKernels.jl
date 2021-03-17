#isdefined(Main, :InterpolationKernels) || include("../src/InterpolationKernels.jl")

module InterpolationKernelsTests

#include("../src/interp.jl")
using InterpolationKernels
using InterpolationKernels:
    infimum, supremum, support, brief,
    compute_weights, generic_compute_weights,
    compute_offset_and_weights, generic_compute_offset_and_weights

using Test

shortname(::Nothing) = ""
shortname(m::RegexMatch) = m.captures[1]
shortname(::Type{T}) where {T} = shortname(string(T))
shortname(str::AbstractString) =
    shortname(match(r"([_A-Za-z][_A-Za-z0-9]*)([({]|$)", str))

const FLOATS = (Float32, Float64, BigFloat)

#------------------------------------------------------------------------------
# Non-specialized, non-optimized spline-based kernels to serve as a reference.

abstract type ReferenceSpline{T,S} <: Kernel{T,S} end

struct Spline1{T,M,N1} <: ReferenceSpline{T,1}
    vals::NTuple{M,T}
    c1::NTuple{N1,T}
end

struct Spline1Prime{T,M,N1} <: ReferenceSpline{T,1}
    vals::NTuple{M,T}
    c1::NTuple{N1,T}
end

struct Spline2{T,M,N1} <: ReferenceSpline{T,2}
    vals::NTuple{M,T}
    c1::NTuple{N1,T}
end

struct Spline2Prime{T,M,N1} <: ReferenceSpline{T,2}
    vals::NTuple{M,T}
    c1::NTuple{N1,T}
end

struct Spline3{T,M,N1,N2} <: ReferenceSpline{T,3}
    vals::NTuple{M,T}
    c1::NTuple{N1,T}
    c2::NTuple{N2,T}
end

struct Spline3Prime{T,M,N1,N2} <: ReferenceSpline{T,3}
    vals::NTuple{M,T}
    c1::NTuple{N1,T}
    c2::NTuple{N2,T}
end

struct Spline4{T,M,N1,N2} <: ReferenceSpline{T,4}
    vals::NTuple{M,T}
    c1::NTuple{N1,T}
    c2::NTuple{N2,T}
end

struct Spline4Prime{T,M,N1,N2} <: ReferenceSpline{T,4}
    vals::NTuple{M,T}
    c1::NTuple{N1,T}
    c2::NTuple{N2,T}
end

#Base.convert(::Type{K}, ker::K) where {K<:ReferenceSpline} = ker

Base.values(ker::ReferenceSpline) = ker.vals

@static if VERSION < v"1.4"
    # Method `evalpoly` appeared in Julia 1.4.
    evalpoly(x::Real, c::NTuple{1,Real}) =
        convert(promote_type(typeof(x), typeof(c[1])), c[1])
    evalpoly(x::Real, c::NTuple{2,Real}) = c[1] + x*c[2]
    evalpoly(x::Real, c::NTuple{3,Real}) = c[1] + x*(c[2] + x*c[3])
    evalpoly(x::Real, c::NTuple{4,Real}) = c[1] + x*(c[2] + x*(c[3] + x*c[4]))
end

function (ker::Spline1{T})(_x::Real) where {T}
    x = T(_x)
    half = T(1)/2
    if -half ≤ x < half
        return evalpoly(abs(x), ker.c1)
    else
        return zero(T)
    end
end

function (ker::Spline1Prime{T})(_x::Real) where {T}
    x = T(_x)
    half = T(1)/2
    if -half ≤ x < half
        return sign(x)*evalpoly(abs(x), ker.c1)
    else
        return zero(T)
    end
end

function (ker::Spline2{T})(_x::Real) where {T}
    x = T(_x)
    abs_x = abs(x)
    if abs_x ≥ 1
        return zero(T)
    else
        return evalpoly(abs_x, ker.c1)
    end
end

function (ker::Spline2Prime{T})(_x::Real) where {T}
    x = T(_x)
    abs_x = abs(x)
    if abs_x ≥ 1
        return zero(T)
    else
        return sign(x)*evalpoly(abs_x, ker.c1)
    end
end

function (ker::Spline3{T})(_x::Real) where {T}
    x = T(_x)
    abs_x = abs(x)
    if abs_x ≥ T(3)/2
        return zero(T)
    elseif abs_x ≤ T(1)/2
        return evalpoly(abs_x, ker.c1)
    else
        return evalpoly(abs_x, ker.c2)
    end
end

function (ker::Spline3Prime{T})(_x::Real) where {T}
    x = T(_x)
    abs_x = abs(x)
    if abs_x ≥ T(3)/2
        return zero(T)
    elseif abs_x ≤ T(1)/2
        return sign(x)*evalpoly(abs_x, ker.c1)
    else
        return sign(x)*evalpoly(abs_x, ker.c2)
    end
end

function (ker::Spline4{T})(_x::Real) where {T}
    x = T(_x)
    abs_x = abs(x)
    if abs_x ≥ 2
        return zero(T)
    elseif abs_x ≤ 1
        return evalpoly(abs_x, ker.c1)
    else
        return evalpoly(abs_x, ker.c2)
    end
end

function (ker::Spline4Prime{T})(_x::Real) where {T}
    x = T(_x)
    abs_x = abs(x)
    if abs_x ≥ 2
        return zero(T)
    elseif abs_x ≤ 1
        return sign(x)*evalpoly(abs_x, ker.c1)
    else
        return sign(x)*evalpoly(abs_x, ker.c2)
    end
end

Base.adjoint(ker::Spline1) =
    Spline1Prime(values(ker), derivepoly(ker.c1))

Base.adjoint(ker::Spline2) =
    Spline2Prime(values(ker), derivepoly(ker.c1))

Base.adjoint(ker::Spline3) =
    Spline3Prime(values(ker), derivepoly(ker.c1),  derivepoly(ker.c2))

Base.adjoint(ker::Spline4) =
    Spline4Prime(values(ker), derivepoly(ker.c1),  derivepoly(ker.c2))

# FIXME: can be @generated
derivepoly(c::NTuple{1,T}) where {T<:AbstractFloat} = (zero(T),)
derivepoly(c::NTuple{2,T}) where {T<:AbstractFloat} = (c[2],)
derivepoly(c::NTuple{3,T}) where {T<:AbstractFloat} = (c[2], 2c[3])
derivepoly(c::NTuple{4,T}) where {T<:AbstractFloat} = (c[2], 2c[3], 3c[4])

bspline1(T::Type{<:AbstractFloat} = Float64) =
    Spline1((), (one(T),))

bspline2(T::Type{<:AbstractFloat} = Float64) =
    Spline2((), (one(T), -one(T)))

bspline3(T::Type{<:AbstractFloat} = Float64) =
    Spline3((),
            (T(3//4), T(0), T(-1)),
            (T(9//8), T(-3//2), T(1//2)))

bspline4(T::Type{<:AbstractFloat} = Float64) =
    Spline4((),
            (T(2//3), T(0), T(-1), T(1//2)),
            (T(4//3), T(-2), T( 1), T(-1//6)))

# A family of cubic splines parameterized by a = ker'(1) and b = ker(1).
generic_cubic_spline(a::Real, b::Real) = genericcubic_spline(Float64, a, b)
generic_cubic_spline(T::Type{<:AbstractFloat}, _a::Real, _b::Real) = begin
    a, b = T(_a), T(_b)
    # p1(x) = c0 + c1*x + c2*x^2
    # p2(x) = (c3 + c4*x)*(x - 2)^2
    # c0 = 1 - 2b
    # c1 = 9b - a - 3
    # c2 = 2 + a - 6b
    # c3 = -a - b
    # c4 = a + 2b
    # Spline4((a, b),
    #         (c0, zero(T), c1, c2),
    #         (4c3, 4(c4 - c3), c3 - 4c4, c4))
    Spline4((a, b),
            (1 - 2b, zero(T), 9b - a - 3, 2 + a - 6b),
            (-4a - 4b, 8a + 12b,  -5a - 9b,  a + 2b))
end

# A family of cardinal cubic splines.
keys_spline(a::Real) = keys_spline(Float64, a)
keys_spline(T::Type{<:AbstractFloat}, _a::Real) = begin
    a = T(_a)
    Spline4((a,),
            (one(T), zero(T), -a - 3, a + 2),
            (-4a, 8a, -5a, a))
end

catmull_rom_spline(T::Type{<:AbstractFloat} = Float64) = keys_spline(T, T(-1//2))

mitchell_netravali_spline(b::Real, c::Real) =
    mitchell_netravali_spline(Float64, b, c)

mitchell_netravali_spline(T::Type{<:AbstractFloat}, _b::Real, _c::Real) = begin
    b, c = T(_b), T(_c)
    Spline4((b,c),
            ((6 - 2b)/6, zero(T), (-18 + 12b + 6c)/6, (12 - 9b - 6c)/6),
            ((8b + 24c)/6, (-12b - 48c)/6, (6b + 30c)/6, (-b - 6c)/6))
end

kernel_model(ker::BSpline{1}) = bspline1(eltype(ker))
kernel_model(ker::BSplinePrime{1}) = bspline1(eltype(ker))'

kernel_model(ker::BSpline{2}) = bspline2(eltype(ker))
kernel_model(ker::BSplinePrime{2}) = bspline2(eltype(ker))'

kernel_model(ker::BSpline{3}) = bspline3(eltype(ker))
kernel_model(ker::BSplinePrime{3}) = bspline3(eltype(ker))'

kernel_model(ker::BSpline{4}) = bspline4(eltype(ker))
kernel_model(ker::BSplinePrime{4}) = bspline4(eltype(ker))'

kernel_model(ker::CatmullRomSpline) = catmull_rom_spline(eltype(ker))
kernel_model(ker::CatmullRomSplinePrime) = catmull_rom_spline(eltype(ker))'

kernel_model(ker::CubicSpline) = generic_cubic_spline(eltype(ker), values(ker)...)
#kernel_model(ker::CubicSplinePrime) = cubic_spline(eltype(ker), values(ker)...)'

kernel_model(ker::CardinalCubicSpline) = keys_spline(eltype(ker), values(ker)...)
kernel_model(ker::CardinalCubicSplinePrime) = keys_spline(eltype(ker), values(ker)...)'

kernel_model(::Kernel) = nothing

#------------------------------------------------------------------------------

raw_type(ker::Kernel) = raw_type(typeof(ker))
for K in (:BSpline,              :BSplinePrime,
          :CubicSpline,   :CubicSplinePrime,
          :CardinalCubicSpline,  :CardinalCubicSplinePrime,
          :CatmullRomSpline,     :CatmullRomSplinePrime,
          :LanczosKernel,        :LanczosKernelPrime)
    @eval raw_type(::Type{<:$K}) = $K
end
#raw_type(::Type{BSpline{S,T}}) where {S,T} = BSpline{S}
#raw_type(::Type{BSplinePrime{S,T}}) where {S,T} = BSplinePrime{S}

has_size_parameter(ker::Kernel) = has_size_parameter(typeof(ker))
has_size_parameter(::Type{<:Kernel}) = false
has_size_parameter(::Type{<:BSpline}) = true
has_size_parameter(::Type{<:BSplinePrime}) = true
has_size_parameter(::Type{<:LanczosKernel}) = true
has_size_parameter(::Type{<:LanczosKernelPrime}) = true

for (K, Kp) in ((:BSpline,             :BSplinePrime),
                (:CubicSpline,  :CubicSplinePrime),
                (:CardinalCubicSpline, :CardinalCubicSplinePrime),
                (:CatmullRomSpline,    :CatmullRomSplinePrime),
                (:LanczosKernel,       :LanczosKernelPrime))
    if K === :BSpline || K === :LanczosKernel
        @eval begin
            derivative_type(::Type{$K}) = $Kp
            derivative_type(::Type{$K{S}}) where {S} = $Kp{S}
            derivative_type(::Type{$K{S,T}}) where {S,T} = $Kp{S,T}
            primitive_type(::Type{$Kp}) = $K
            primitive_type(::Type{$Kp{S}}) where {S} = $K{S}
            primitive_type(::Type{$Kp{S,T}}) where {S,T} = $K{S,T}
        end
    else
        @eval begin
            derivative_type(::Type{$K}) = $Kp
            derivative_type(::Type{$K{T}}) where {T} = $Kp{T}
            primitive_type(::Type{$Kp}) = $K
            primitive_type(::Type{$Kp{T}}) where {T} = $K{T}
        end
    end
end

derivative_type(ker::Kernel) = derivative_type(typeof(ker))
primitive_type(ker::Kernel) = primitive_type(typeof(ker))

derivative_type(::Type{<:Kernel}) = Any
primitive_type(::Type{<:Kernel}) = Any

function simple_derivative(f::Kernel, x::Real, dx::Real)
    x1 = x - dx
    x2 = x + dx
    return (f(x1) - f(x2))/(x1 - x2)
end

function simple_derivative(f::Kernel{T},
                           x::AbstractArray{T},
                           dx::Real = sqrt(eps(T))) where {T}
    df = similar(x)
    for i in eachindex(x, df)
        df[i] = simple_derivative(f, x[i], dx)
    end
    return df
end

@testset "Interpolation kernels" begin
    @testset "Constructors" begin
        for S in 1:4
            @test eltype(BSpline{S}()) === Float64
            @test eltype(BSpline{S,Float32}()) === Float32
            @test eltype(BSplinePrime{S}()) === Float64
            @test eltype(BSplinePrime{S,Float32}()) === Float32
        end
        @test eltype(CubicSpline(1,0)) === Float64
        @test eltype(CubicSpline(1,0f0)) === Float32
        @test eltype(CubicSpline(1,0.0)) === Float64
        @test eltype(CubicSpline{Float32}(1,0.0)) === Float32

        @test eltype(CubicSplinePrime(1,0)) === Float64
        @test eltype(CubicSplinePrime(1,0f0)) === Float32
        @test eltype(CubicSplinePrime(1,0.0)) === Float64
        @test eltype(CubicSplinePrime{Float32}(1,0.0)) === Float32

        @test eltype(CardinalCubicSpline(0)) === Float64
        @test eltype(CardinalCubicSpline(0f0)) === Float32
        @test eltype(CardinalCubicSpline(0.0)) === Float64
        @test eltype(CardinalCubicSpline{Float32}(0.0)) === Float32

        @test eltype(CardinalCubicSplinePrime(0)) === Float64
        @test eltype(CardinalCubicSplinePrime(0f0)) === Float32
        @test eltype(CardinalCubicSplinePrime(0.0)) === Float64
        @test eltype(CardinalCubicSplinePrime{Float32}(0.0)) === Float32

        @test eltype(CatmullRomSpline()) === Float64
        @test eltype(CatmullRomSpline{Float32}()) === Float32
        @test eltype(CatmullRomSplinePrime()) === Float64
        @test eltype(CatmullRomSplinePrime{Float32}()) === Float32

        for S in (4,6)
            @test eltype(LanczosKernel{S}()) === Float64
            @test eltype(LanczosKernel{S,Float32}()) === Float32
            @test eltype(LanczosKernelPrime{S}()) === Float64
            @test eltype(LanczosKernelPrime{S,Float32}()) === Float32
        end
    end

    @testset "Values" begin
        for ker in (BSpline{1}(),
                    BSplinePrime{1}(),
                    BSpline{2}(),
                    #BSplinePrime{2}(), # FIXME:
                    BSpline{3}(),
                    BSplinePrime{3}(),
                    BSpline{4}(),
                    BSplinePrime{4}(),
                    CubicSpline(-1/3,1/5),
                    CubicSplinePrime(-2/3,1/7),
                    CardinalCubicSpline(-1/3),
                    CardinalCubicSplinePrime(-2/3),
                    CatmullRomSpline(),
                    CatmullRomSplinePrime(),
                    LanczosKernel{4}(),
                    LanczosKernelPrime{4}())
            mdl = kernel_model(ker)
            if mdl === nothing
                continue
            end
            s = 0.123 # step chosen so as to avoid integral values of x
            sup = support(ker)
            x = infimum(sup)-s:s:supremum(sup)+s
            @test ker.(x) ≈ mdl.(x)
        end
    end

    # Check the ability of the generic cubic spline to emulate other kernels.
    @testset "Equivalences" begin
        s = 0.1
        # Emulation of the cubic B-spline.
        let ker = BSpline{4}(),
            mdl = CubicSpline(-1/2,1/6)
            sup = support(ker)
            x = infimum(sup)-s:s:supremum(sup)+s
            @test mdl.(x) ≈ ker.(x)
            @test mdl'.(x) ≈ ker'.(x)
        end
        # Emulation of the family of cardinal cubic splines.
        for a in (-1/3, 0.1, -1/2)
            ker = CardinalCubicSpline(a)
            mdl = CubicSpline(a, 0)
            sup = support(ker)
            x = infimum(sup)-s:s:supremum(sup)+s
            @test mdl.(x) ≈ ker.(x)
            @test mdl'.(x) ≈ ker'.(x)
        end
        # Emulation of the Catmull-Rom kernel.
        let ker = CatmullRomSpline(),
            mdl = CubicSpline(-1/2, 0)
            sup = support(ker)
            x = infimum(sup)-s:s:supremum(sup)+s
            @test mdl.(x) ≈ ker.(x)
            @test mdl'.(x) ≈ ker'.(x)
        end
        # Emulation of the family of Mitchell & Netravali kernels.
        for (b,c) in ((1/3,1/3), (0.1,0.7), (-0.2,0.9))
            let ker = MitchellNetravaliSpline(b, c),
                mdl = mitchell_netravali_spline(b, c)
                sup = support(ker)
                x = infimum(sup)-s:s:supremum(sup)+s
                @test mdl.(x) ≈ ker.(x)
                @test mdl'.(x) ≈ ker'.(x)
            end
        end
    end

    @testset "Conversion" begin
        for ker in (BSpline{1}(),
                    BSplinePrime{1}(),
                    BSpline{2}(),
                    BSplinePrime{2}(),
                    BSpline{3}(),
                    BSplinePrime{3}(),
                    BSpline{4}(),
                    BSplinePrime{4}(),
                    CatmullRomSpline(),
                    CatmullRomSplinePrime(),
                    CubicSpline(-1/3,1/5),
                    CubicSplinePrime(-2/3,1/7),
                    CardinalCubicSpline(-1/3),
                    CardinalCubicSplinePrime(-2/3),
                    LanczosKernel{4}(),
                    LanczosKernelPrime{4}())
            K = raw_type(ker)
            @test eltype(ker) === Float64
            @test Kernel(ker) === ker
            @test Kernel{eltype(ker)}(ker) === ker
            if has_size_parameter(ker)
                @test Kernel{eltype(ker),length(ker)}(ker) === ker
            end
            for T in FLOATS
                if has_size_parameter(ker)
                    S = length(ker)
                    @test typeof(T(ker)) <: K{S,T}
                    if T !== BigFloat
                        @test T(ker) === K{S,T}(ker)
                        @test T(ker) === K{S,T}(values(ker)...)
                    end
                else
                    @test typeof(T(ker)) <: K{T}
                    if T !== BigFloat
                        @test T(ker) === K{T}(ker)
                        @test T(ker) === K{T}(values(ker)...)
                    end
                end
            end
        end
    end

    @testset "Methods" begin
        tol = 1e-14
        for (ker, len, nrml, card) in (
            (BSpline{1}(),                   1, true,  true),
            (BSplinePrime{1}(),              1, false, false),
            (BSpline{2}(),                   2, true,  true),
            (BSplinePrime{2}(),              2, false, false),
            (BSpline{3}(),                   3, true,  false),
            (BSplinePrime{3}(),              3, false, false),
            (BSpline{4}(),                   4, true,  false),
            (BSplinePrime{4}(),              4, false, false),
            (CatmullRomSpline(),             4, true,  true),
            (CatmullRomSplinePrime(),        4, false, false),
            (CubicSpline(-1/3,1/9),          4, true,  false),
            (CubicSplinePrime(-1/3,1/9),     4, false, false),
            (CardinalCubicSpline(-1),        4, true,  true),
            (CardinalCubicSplinePrime(-1/2), 4, false, false),
            (LanczosKernel{2}(),             2, false, true),
            (LanczosKernel{4}(),             4, false, true),
            (LanczosKernel{6}(),             6, false, true),
            (LanczosKernel{8}(),             8, false, true),
            (LanczosKernelPrime{4}(),        4, false, false),
            (LanczosKernelPrime{6}(),        6, false, false))
            name = brief(ker)
            @test isa(name, String)
            @testset "$name" begin
                @test isa(summary(ker), String)
                @test isnormalized(ker) == nrml
                @test iscardinal(ker) == card
                @test length(ker) == len
                @test length(typeof(ker)) == len
                @test shortname(summary(ker)) == shortname(typeof(ker))
                if iscardinal(ker)
                    @test ker(0) == 1
                    @test maximum(abs.(ker([-3,-2,-1,1,2,3]))) ≤ tol
                end
                sup = support(ker)
                @test length(sup) == len
                @test supremum(sup) == infimum(sup) + length(sup)

                # Check compute_weights.  Compare values computed by optimal,
                # generic and reference method.
                k = (length(ker) + 1) >> 1
                err_gen, err_opt = 0.0, 0.0
                for t in (0.0, 0.1, 0.2, 0.3, 0.4)
                    wgt_ref = ntuple(i -> ker(t + k - i), length(ker))
                    wgt_gen = generic_compute_weights(ker, t)
                    wgt_opt = compute_weights(ker, t)
                    err_gen = max(err_gen, maximum(abs.(wgt_gen .- wgt_ref)))
                    err_opt = max(err_opt, maximum(abs.(wgt_opt .- wgt_ref)))
                end
                @test err_gen ≤ tol
                @test err_opt ≤ tol

                # Check compute_offset_and_weights.  Compare values computed by optimal, generic
                # and reference method.
                k = (length(ker) + 1) >> 1
                off_gen_err, wgt_gen_err = 0.0, 0.0
                off_opt_err, wgt_opt_err = 0.0, 0.0
                for x in (50.0, -30.9, 8.2, 11.3, -27.6)
                    # FIXME: These rules are only valid for symmetric kernels.
                    if isodd(length(ker))
                        round_x = round(x)
                        t = x - round_x
                        off_ref = round_x - (length(ker) + 1)/2
                    else
                        floor_x = floor(x)
                        t = x - floor_x
                        off_ref = floor_x - length(ker)/2
                    end
                    wgt_ref = ntuple(i -> ker(t + k - i), length(ker))
                    off_gen, wgt_gen = generic_compute_offset_and_weights(ker, x)
                    off_gen_err = max(off_gen_err, abs(off_gen - off_ref))
                    wgt_gen_err = max(wgt_gen_err, maximum(abs.(wgt_gen .- wgt_ref)))
                    off_opt, wgt_opt = compute_offset_and_weights(ker, x)
                    off_opt_err = max(off_opt_err, abs(off_opt - off_ref))
                    wgt_opt_err = max(wgt_opt_err, maximum(abs.(wgt_opt .- wgt_ref)))
                end
                @test off_gen_err == 0
                @test wgt_gen_err ≤ tol
                @test off_opt_err == 0
                @test wgt_opt_err ≤ tol
            end
        end
    end

    @testset "Derivatives" begin
        for ker in (BSpline{1}(),
                    BSpline{2}(),
                    BSpline{3}(),
                    BSpline{4}(),
                    CatmullRomSpline(),
                    LanczosKernel{4}(),
                    LanczosKernel{6}())
            der = ker'
            @test derivative_type(ker) === typeof(der)
            @test primitive_type(der) === typeof(ker)
            @test values(der) === values(ker)
            @test der === derivative_type(ker)(values(ker)...)
            @test ker === primitive_type(der)(values(der)...)

            s = 0.123 # step chosen so as to avoid integral values of x
            sup = support(ker)
            x = infimum(sup)+s:s:supremum(sup)-s
            @test der.(x) ≈ simple_derivative(ker, x)

        end
    end
end

end # module
