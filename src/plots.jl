import Plots
import Plots: RecipesBase

# This macro allows us to directly plot any array of PointMasses
@Plots.recipe function f(body::Vector{PointMass{T}}; bodylines = nothing, markersize_is_mass=true) where T
    grid       --> true
    linecolor  --> :black
    legend     :=  :none
    # Each scatter size depends on the mass of the point
    markersize --> (markersize_is_mass ? 4.0*mass.(body) : 4.0)
    delete!(plotattributes, :markersize_is_mass)

    # Plot each line connecting the Point Massses _before_ the masses themselves
    if bodylines != nothing
        seriestype := :path
        for l in bodylines
            # No point in drawing a line with length zero...
            l[1] == l[2] && continue
            local x = hcat(pos(body[l[1]]), pos(body[l[2]]))
            @Plots.series begin
                (x[1,:], x[2,:], x[3,:])
            end
        end
    end
    delete!(plotattributes, :bodylines)

    # Plot body
    seriestype := :scatter

    local x = hcat(pos.(body)...)
    @Plots.series begin
        (x[1,:], x[2,:], x[3,:])
    end
end

"""
    plotmodel(model, SoR; kw...)

Receive a [`Model`](@ref) and return a `Plots.Animation` from its data.

# Arguments

- `SoR`: means "system of reference" and accepts one of the following symbols: `:bodyframe`, `:inertialframe`, `:both`.
- `fps`: frames per second. Default is `24`.
- `duration`: length of animation in seconds. Default is the total timespan stored in `m`.
- `saveas`: filename to save animation, supported extensions are gif, mp4 or mov. If left blank, file is not saved.
- `bodylines`: Array of points pairs. Stores the data who says which points should be linked. Default is empty.
- `markersize_is_mass`: Says if attribute `markersize` should be weighted by the mass of each particle. Default is `true`.

Additionally, any keyword argument supported by `Plots.plot`
is accepted and will be repassed to the plot.

# Examples

```julia
julia> plotmodel(m, :both, fps=15, saveas="example.gif", color=:green,viewangle=(30,60))
```
"""
function plotmodel( m::Model
                  , SoR # System of Reference
                  ; fps=24
                  , duration=m.timespan[2]-m.timespan[1]
                  , saveas=nothing
                  , bodylines=nothing
                  , markersize_is_mass=true
                  , args... # Additional keyword arguments to be passed to Plots.plot
                  )
    if !(SoR in [:inertialframe, :bodyframe, :both])
        error("Unknown plot type :" * string(SoR) * ".\n Try one of the following: :bodyframe, :inertialframe, :both.")
        return nothing
    end
    local frames = convert(Int, ceil(fps*duration))
    xlim, ylim, zlim = getframelimits(m, SoR, frames)

    # Build the animation
    anime = Plots.Animation()
    for t in range(m.timespan[1], m.timespan[2], length=frames)
        plt =
            if SoR == :inertialframe
                Plots.plot(m.inertialframe(t),
                           bodylines = bodylines,
                           xlim  = xlim,
                           ylim  = ylim,
                           zlim  = zlim,
                           title = "Inertial Frame";
                           args...
                          )
            elseif SoR == :bodyframe
                Plots.plot(m.bodyframe(t),
                           bodylines = bodylines,
                           xlim  = xlim,
                           ylim  = ylim,
                           zlim  = zlim,
                           title = "Body Frame";
                           args...
                          )
            elseif SoR == :both
                plt_it = Plots.plot(m.inertialframe(t),
                                    bodylines = bodylines,
                                    xlim  = xlim,
                                    ylim  = ylim,
                                    zlim  = zlim,
                                    title = "Inertial Frame";
                                    args...
                                   )
                plt_bd = Plots.plot(m.bodyframe(t),
                                    bodylines = bodylines,
                                    xlim  = xlim,
                                    ylim  = ylim,
                                    zlim  = zlim,
                                    title = "Body Frame";
                                    args...
                                   )
                Plots.plot(plt_bd, plt_it, layout=(1,2), link=:all)
            end
        # Push plot into animation
        Plots.frame(anime, plt)
    end

    saveas isa String && saveanimation(anime, saveas, fps=fps)
    return anime
end

# Return limits of coordinates over all instants / frames for each axis
function getframelimits(m::Model, SoR, frames)
    local lim_min = [Inf, Inf, Inf]
    local lim_max = [-Inf, -Inf, -Inf]
    for t in range(m.timespan[1], m.timespan[2], length=frames)
        if SoR == :bodyframe
            x = m.bodyframe(t)
        elseif SoR == :inertialframe
            x = m.inertialframe(t)
        elseif SoR == :both
            x = vcat(m.bodyframe(t), m.inertialframe(t))
        end
        min_t = minimum(hcat(pos.(x)...), dims=2)
        max_t = maximum(hcat(pos.(x)...), dims=2)
        lim_min = min.(lim_min, min_t)
        lim_max = max.(lim_max, max_t)
    end
    return ((lim_min[1], lim_max[1]), (lim_min[2], lim_max[2]), (lim_min[3], lim_max[3]))
end

"""
    saveanimation(anime, saveas; fps=30)

Receive an `Plots.Animation` and save it as a file.

Supported formats are 'gif', 'mp4' and 'mov'.
The extension is automatically detected from `saveas`
and, in case of ambiguity, defaults to '.gif'.
"""
function saveanimation(anime, saveas::String; fps::Int=30)
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
