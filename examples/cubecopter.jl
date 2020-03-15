################################
# Example of model: Cubecopter #
################################
#=
   8-points standing cube coupled
   with a 2-point moving helix
=#
using DeformableBodies

const r_0 = centralize(
    [ PointMass(1., [ 1.,  1.,  0.])
    , PointMass(1., [ 1., -1.,  0.])
    , PointMass(1., [-1.,  1.,  0.])
    , PointMass(1., [-1., -1.,  0.])
    , PointMass(1., [ 1.,  1., -1.])
    , PointMass(1., [ 1., -1., -1.])
    , PointMass(1., [-1.,  1., -1.])
    , PointMass(1., [-1., -1., -1.])
    , PointMass(.5, [ .0, -1.,  1.])
    , PointMass(.5, [ .0,  1.,  1.])
    ])

# Define vector of trajectories
bodies = Function[]
for x in r_0[1:end-2]
    push!(bodies, let y=x; t -> y;end)
end

# Other parts move according to given rule
const tmax = 5.0
const ω    = 2*π/tmax
const e_1  = [0., 0., 1.]
for x in r_0[end-1:end]
    push!(bodies, let x=x; t -> rotate(x, axis=e_1, angle=ω*t); end)
end

# Create the model
model = Model( bodies
             , 0.
             , tmax
             , one(Quaternion)
             , zeros(3)
             )
println(model)
# And solve it
_, rotations, momentum = solve!(model)

# Lines to connect different parts of the body
# (Only important for plotting)
bodylines = Tuple[]
for i = 1:length(r_0), j = i:length(r_0)
    if count(a -> first(a) == last(a), zip(pos(r_0[i]),pos(r_0[j]))) == 2
        push!(bodylines, (i,j))
    end
end

println("\nTrajectories are stored on the following variables:")
println("Internal frame: model.bodyframe")
println("Inertial frame: model.inertialframe")

println("\nSaving animation to file \"cubecopter.gif\"")
# Plotting time!
anim = plotmodel( model, :both
                , bodylines=bodylines
                , duration=tmax
                , markercolor=:purple
                , size=(800,400)
                , saveas="cubecopter.gif"
                )
