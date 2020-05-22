import Luxor

assets_dir = joinpath(@__DIR__, "src/assets")

# Use Luxor do draw logo
function drawlogo()
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
end

# Turn logo into icon
function make_logoicon(icofile, logofile)
    if !isfile(joinpath(assets_dir, icofile)) && isfile(joinpath(assets_dir, logofile))
        run(`convert -background none $(joinpath(assets_dir, logofile)) -define icon:auto-resize $(joinpath(assets_dir, icofile))`)
        println("Sucessfully created asset: ico.")
    end
end
