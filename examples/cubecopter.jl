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

model = Model( bodies
             , 0.
             , 7.
             , one(Quaternion)
             , zeros(3)
             )

_, rotations, momentum = solve!(model)

println("Trajectories are stored on the following variables:")
println("Body's internal frame: model.trajectories")
println("Inertial frame       : model.inertial")
