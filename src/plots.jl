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
            local x = hcat(body[l[1]].pos, body[l[2]].pos)
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

Receives a [`Model`](@ref) and returns a [`Animation`](@ref) from its data.
The argument `SoR` means "system of reference" and accepts one of the following symbols: `:bodyframe`, `:inertialframe`, `:both`.

The aditional arguments currently available are
- `fps`: frames per second.
- `duration`: length of animation in seconds. Default is `model.t_max`.
- `saveas`: filename to save animation, supported extensions are gif, mp4 or mov. By default, file is not saved.
- `bodylines`: Array of points pairs. Stores the data who says which points should be linked.
- `markersize_is_mass`: Says if attribute `markersize` should be weighted by the mass of each particle. Default is `true`.

Additionally, any keyword argument supported by [`Plots.plot`](@ref)
is accepted and will be repassed to the plot.

# Examples
```jldoctest
julia> plotmodel(m, :both, fps=15, saveas="example.gif", color=:green,viewangle=(30,60))
```
"""
function plotmodel( m::Model
                  , SoR=:both # System of Reference
                  ; fps=24
                  , duration=m.t_max
                  , saveas=nothing
                  , bodylines=nothing
                  , args... # Additional keyword arguments to be passed to Plots.plot
                  )
    if !(SoR in [:inertialframe, :bodyframe, :both])
        error("Unknown plot type :" * string(SoR) * ".\n Try one of the following: :bodyframe, :inertialframe, :both.")
        return nothing
    end

    local frames = convert(Int, ceil(fps*duration))

    # Get maximum of coordinates on all iterations
    lim_max = [-Inf, -Inf, -Inf]
    lim_min = [Inf, Inf, Inf]
    for t in range(m.t_min, m.t_max, length=frames)
        if SoR == :bodyframe
            x = m.bodyframe(t)
        elseif SoR == :inertialframe
            x = m.inertialframe(t)
        elseif SoR == :both
            x = vcat(m.bodyframe(t), m.inertialframe(t))
        end
        max_t = maximum(hcat(pos.(x)...), dims=2)
        min_t = minimum(hcat(pos.(x)...), dims=2)
        lim_max = max.(lim_max, max_t)
        lim_min = min.(lim_min, min_t)
    end

    # Build the animation
    anime = Plots.Animation()
    for t in range(m.t_min, m.t_max, length=frames)
        plt =
            if SoR == :bodyframe
                Plots.plot(m.bodyframe(t),
                           bodylines = bodylines,
                           xlim  = (lim_min[1], lim_max[1]),
                           ylim  = (lim_min[2], lim_max[2]),
                           zlim  = (lim_min[3], lim_max[3]),
                           title = "Body Frame";
                           args...
                          )
            elseif SoR == :inertialframe
                Plots.plot(m.inertialframe(t),
                           bodylines = bodylines,
                           xlim  = (lim_min[1], lim_max[1]),
                           ylim  = (lim_min[2], lim_max[2]),
                           zlim  = (lim_min[3], lim_max[3]),
                           title = "Inertial Frame";
                           args...
                          )
            elseif SoR == :both
                plt_bd = Plots.plot(m.bodyframe(t),
                                    bodylines = bodylines,
                                    xlim  = (lim_min[1], lim_max[1]),
                                    ylim  = (lim_min[2], lim_max[2]),
                                    zlim  = (lim_min[3], lim_max[3]),
                                    title = "Body Frame";
                                    args...
                                   )
                plt_it = Plots.plot(m.inertialframe(t),
                                    bodylines = bodylines,
                                    xlim  = (lim_min[1], lim_max[1]),
                                    ylim  = (lim_min[2], lim_max[2]),
                                    zlim  = (lim_min[3], lim_max[3]),
                                    title = "Inertial Frame";
                                    args...
                                   )
                Plots.plot(plt_bd, plt_it,
                           layout=(1,2),
                           link=:all
                          )
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
