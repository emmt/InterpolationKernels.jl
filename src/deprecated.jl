@static if VERSION ≥ v"1.9"
    @deprecate(Base.Float16(ker::Kernel), convert(Kernel{Float16}, ker), false)
    @deprecate(Base.Float32(ker::Kernel), convert(Kernel{Float32}, ker), false)
    @deprecate(Base.Float64(ker::Kernel), convert(Kernel{Float64}, ker), false)
    @deprecate(Base.BigFloat(ker::Kernel), convert(Kernel{BigFloat}, ker), false)
else
    @deprecate(Float16(ker::Kernel), convert(Kernel{Float16}, ker), false)
    @deprecate(Float32(ker::Kernel), convert(Kernel{Float32}, ker), false)
    @deprecate(Float64(ker::Kernel), convert(Kernel{Float64}, ker), false)
    @deprecate(BigFloat(ker::Kernel), convert(Kernel{BigFloat}, ker), false)
end

# @deprecate on object calls only work with Julia ≥ 1.9
@static if VERSION ≥ v"1.9"
    @deprecate((ker::BSpline)(A::AbstractArray), map(ker, A), false)
    @deprecate((ker::BSplinePrime)(A::AbstractArray), map(ker, A), false)
    @deprecate((ker::CubicSpline)(A::AbstractArray), map(ker, A), false)
    @deprecate((ker::CubicSplinePrime)(A::AbstractArray), map(ker, A), false)
    @deprecate((ker::CardinalCubicSpline)(A::AbstractArray), map(ker, A), false)
    @deprecate((ker::CardinalCubicSplinePrime)(A::AbstractArray), map(ker, A), false)
    @deprecate((ker::CatmullRomSpline)(A::AbstractArray), map(ker, A), false)
    @deprecate((ker::CatmullRomSplinePrime)(A::AbstractArray), map(ker, A), false)
    @deprecate((ker::LanczosKernel)(A::AbstractArray), map(ker, A), false)
    @deprecate((ker::LanczosKernelPrime)(A::AbstractArray), map(ker, A), false)
else
    (ker::BSpline)(A::AbstractArray) = map(ker, A)
    (ker::BSplinePrime)(A::AbstractArray) = map(ker, A)
    (ker::CubicSpline)(A::AbstractArray) = map(ker, A)
    (ker::CubicSplinePrime)(A::AbstractArray) = map(ker, A)
    (ker::CardinalCubicSpline)(A::AbstractArray) = map(ker, A)
    (ker::CardinalCubicSplinePrime)(A::AbstractArray) = map(ker, A)
    (ker::CatmullRomSpline)(A::AbstractArray) = map(ker, A)
    (ker::CatmullRomSplinePrime)(A::AbstractArray) = map(ker, A)
    (ker::LanczosKernel)(A::AbstractArray) = map(ker, A)
    (ker::LanczosKernelPrime)(A::AbstractArray) = map(ker, A)
end
