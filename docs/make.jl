using Documenter
using DeformableBodies
using DeformableBodies: Quaternions, velocity

include(joinpath(@__DIR__, "logo.jl"))

# Make logo and icon
drawlogo()
make_logoicon("favicon.ico", "logo.svg")

# Set the right metadata for doctests
DocMeta.setdocmeta!(DeformableBodies, :DocTestSetup, :(using DeformableBodies); recursive=true)

makedocs(
    sitename = "DeformableBodies.jl",
    authors  = "Iago Leal de Freitas",
    format   = Documenter.HTML(
        prettyurls = !("local" in ARGS),    # Deactivate on local builds
        assets     = ["assets/favicon.ico"],
        canonical  = "https://iagoleal.github.io/DeformableBodies.jl/dev/"
    ),
    modules  = [DeformableBodies, DeformableBodies.Quaternions],
    pages    = [
        "Introduction" => "index.md",
        "Tutorial" => "tutorial-cubecopter.md",
        "Guides" => ["guide-quaternions.md"],
        "Reference" => [
            "Quaternions" => "ref-quaternions.md",
            "DeformableBodies" => "ref-models.md"
        ]
   ]
)

# Deploy site to Github Pages
if !("local" in ARGS)
    deploydocs(
        repo = "github.com/iagoleal/DeformableBodies.jl.git"
    )
end
