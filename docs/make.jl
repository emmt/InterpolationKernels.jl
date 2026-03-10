using Documenter

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "src")))

using InterpolationKernels

DEPLOYDOCS = (get(ENV, "CI", nothing) == "true")

makedocs(
    sitename = "Interpolation kernels for Julia",
    format = Documenter.HTML(
        canonical = "https://github.com/emmt/InterpolationKernels.jl",
        edit_link = "master",
        prettyurls = DEPLOYDOCS,
    ),
    authors = "Éric Thiébaut and contributors",
    pages = [
        "index.md",
        "basics.md",
        "kernels.md",
        "interpolation.md",
        "library.md"]
)

if DEPLOYDOCS
    deploydocs(
        repo = "github.com/emmt/InterpolationKernels.jl.git",
        devbranch = "master",
    )
end
