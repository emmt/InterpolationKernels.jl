using Documenter

push!(LOAD_PATH, "../src/")
using InterpolationKernels

DEPLOYDOCS = (get(ENV, "CI", nothing) == "true")

makedocs(
    sitename = "Interpolation kernels for Julia",
    format = Documenter.HTML(
        prettyurls = DEPLOYDOCS,
    ),
    authors = "Éric Thiébaut and contributors",
    pages = ["index.md", "install.md", "basics.md", "kernels.md", "boundaries.md",
         "interpolation.md", "library.md"]
)

if DEPLOYDOCS
    deploydocs(
        repo = "github.com/emmt/ArrayTools.jl.git",
    )
end
