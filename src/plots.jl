import Plots

# Plot a single body
function plotbody(r; title="")
    local pos = hcat((x.pos for x in r)...)
    x,y,z = pos[1,:], pos[2,:], pos[3,:]
    Plots.scatter(x, y, z, seriestype=:scatter,
                           title=title)
end

"""
    plotmodel(model, SoR; kw...)

Receives a [`Model`](@ref) and returns a [`Animation`](@ref) from its data.
The argument `SoR` means "system of reference" and accepts one of the following symbols: `:bodyframe`, `:inertialframe`, `:both`.

The aditional arguments currently available are
- `fps`: frames per second.
- `duration`: duration of animation in seconds.
- `saveas`: filename to save animation, supported extensions are gif, mp4 or mov. By default, file is not saved.
- `backend`: Which Plots.jl backend to use.
"""
function plotmodel( m::Model
                    , SoR=:bodyframe # System of Reference
                    ; fps = 15
                    , duration=5.0
                    , saveas=nothing
                    , backend=nothing
                    )
    if backend != nothing
        backend()
    end
    local frames = convert(Int, ceil(fps*duration))

    # Build the animation
    anime = Plots.Animation()
    for t in range(m.t_min, m.t_max, length=frames)
        plt =
            if SoR == :bodyframe
                plotbody(m.bodyframe(t), title="Body Frame")
            elseif SoR == :inertialframe
                plotbody(m.inertialframe(t), title="Inertial Frame")
            elseif SoR == :both
                l = @Plots.layout [a b]
                plt1 = plotbody(m.bodyframe(t), title="Body Frame")
                plt2 = plotbody(m.inertialframe(t), title="Inertial Frame")
                Plots.plot(plt1, plt2, layout=l)
            else
                error("Unknown plot type :" * string(SoR) * ".\n Try one of the following: :original, :inertial, :both.")
                return
            end
        # Push plot into animation
        Plots.frame(anime, plt)
    end

    if saveas isa String
        local extension = match(r"\.[0-9a-z]+$", saveas)
        if extension != nothing
            extension = extension.match
        end
        if extension == ".gif"
            Plots.gif(anime, saveas, fps = fps)
        elseif extension == ".mp4"
            Plots.mp4(anime, saveas, fps = fps)
        elseif extension == ".mov"
            Plots.mov(anime, saveas, fps = fps)
        else
            Plots.gif(anime, saveas * ".gif", fps = fps)
        end
    end

    return anime
end
