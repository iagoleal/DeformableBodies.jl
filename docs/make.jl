using Documenter
using DeformableBodies
using DeformableBodies: Quaternions, velocity

assets_dir = joinpath(@__DIR__, "src/assets")

# Use Luxor do draw logo
import Luxor
for data in ("logo" => "black", "logo-dark" => "white")
    # Only create files if they do not already exist
    if ! isfile(joinpath(assets_dir, "$(data.first).svg"))
        Luxor.Drawing(460, 320, joinpath(assets_dir, "$(data.first).svg"))
        Luxor.origin()
        Luxor.translate(-20,25)
        O = Luxor.Point(0,0)
        axis = Luxor.Point(150, -150)
        Luxor.setcolor(data.second)
        Luxor.arrow(-0.5*axis, 1.2*axis, linewidth=7., arrowheadlength=25)
        Luxor.transform([1, 0, tand(45.), 1, 0, 0])
        Luxor.juliacircles()
        axis = Luxor.Point(250, -125)
        Luxor.setcolor(Luxor.julia_blue)
        Luxor.arrow(1.1*axis, 35, 0, 1.7*π, arrowheadangle=π/4, linewidth=4, arrowheadlength=15)
        Luxor.finish()
        println("Sucessfully created asset: $(data.first).")
    end
end
# Turn logo into icon
if !isfile(joinpath(assets_dir, "favicon.ico")) && isfile(joinpath(assets_dir, "logo.svg"))
    run(`convert -background none $(joinpath(assets_dir, "logo.svg")) -define icon:auto-resize $(joinpath(assets_dir, "favicon.ico"))`)
    println("Sucessfully created asset: ico.")
end
DocMeta.setdocmeta!(DeformableBodies, :DocTestSetup, :(using DeformableBodies); recursive=true)

makedocs(
    sitename = "DeformableBodies.jl",
    authors  = "Iago Leal de Freitas",
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets     = ["assets/favicon.ico"]
    ),
    modules  = [DeformableBodies, DeformableBodies.Quaternions],
    pages    = [
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
