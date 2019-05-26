import Plots

# Plot a single body
function plotbody!(plt, body)
    local pos = hcat((x.pos for x in body)...)
    Plots.plot!(plt, pos[1,:], pos[2,:], pos[3,:],
                seriestype=:scatter,
               )
    return plt
end

# Plot lines connecting different point masses
function plotbodylines!(plt, body, lines; linecolor=:black)
    for l in lines
        # No point in drawing a line with length zero...
        l[1] == l[2] && continue
        x  = [body[l[1]].pos body[l[2]].pos]
        Plots.plot!(plt, x[1,:], x[2,:], x[3,:],
                    seriestype=:path3d,
                    linecolor=linecolor,
                    leg=nothing
                   )
    end
    return plt
end

# Entire process for each single plot
@inline function plot_process(data, bodylines; title="")
    plt = Plots.plot3d(title=title)
    if bodylines != nothing
        plotbodylines!(plt, data, bodylines,
                       linecolor=:black
                      )
    end
    plotbody!(plt, data)
    return plt
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
- `bodylines`: Array of points pairs. Storages the data about what points should be linked.
"""
function plotmodel( m::Model
                    , SoR=:both # System of Reference
                    ; fps=24
                    , duration=m.t_max
                    , saveas=nothing
                    , backend=nothing
                    , bodylines=nothing
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
                plot_process(m.bodyframe(t), bodylines;
                             title = "Body Frame")
            elseif SoR == :inertialframe
                plot_process(m.inertialframe(t), bodylines;
                             title = "Inertial Frame")
            elseif SoR == :both
                plt_bd = plot_process(m.bodyframe(t), bodylines;
                                      title = "Body Frame")
                plt_it = plot_process(m.inertialframe(t), bodylines;
                                      title = "Inertial Frame")
                Plots.plot(plt_bd, plt_it,
                           layout=(1,2),
                           link=:all
                          )
            else
                error("Unknown plot type :" * string(SoR) * ".\n Try one of the following: :original, :inertial, :both.")
                return
            end
        # Push plot into animation
        Plots.frame(anime, plt)
    end

    saveas isa String && saveanimation(anime, saveas, fps=fps)
    return anime
end

function saveanimation(anime, saveas::String; fps::Int=30)
    print(fps)
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
