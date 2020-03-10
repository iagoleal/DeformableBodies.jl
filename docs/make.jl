using Documenter
using DeformableBodies
using DeformableBodies: Quaternions, velocity

DocMeta.setdocmeta!(DeformableBodies, :DocTestSetup, :(using DeformableBodies); recursive=true)

makedocs(
    sitename = "DeformableBodies.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [DeformableBodies, DeformableBodies.Quaternions],
    pages = [
        "Introduction" => "index.md",
        "Tutorial" => "tutorial-cubecopter.md",
        "theory.md",
        "Reference" => [
            "Quaternions" => "ref-quaternions.md",
            "DeformableBodies" => "ref-models.md"
        ]
   ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
