######################################################
# Example of model: falling cat as a deformable body #
######################################################
using DeformableBodies

# Initial values for the body
const r_0 = centralize(
    [ PointMass(3.,  [0.,   0.,   0.])   # torso center
    , PointMass(2.,  [0.,  -1.,   0.])   # torso back
    , PointMass(1.,  [ .5, -1.,   1.])   # left back leg
    , PointMass(1.,  [-.5, -1.,   1.])   # right back leg
    , PointMass(.2,  [0.,  -1.5, -0.5])  # tail
    , PointMass(2.,  [0.,   1.,   0.])   # torso front
    , PointMass(1.,  [ .5,  1.,   1.])   # left front leg
    , PointMass(1.,  [-.5,  1.,   1.])   # right front leg
    , PointMass(1.5, [0.,   1.2, -0.2])  # head
   ])

# Parameters
const tmax    = 15.0                  # Time of complete motion
const tmax_r1 = tmax/15               # Time of bending part
const tmax_r2 = tmax - tmax_r1        # Time of unbending part
const θmax    = -π/6.0                # Bending angle
const ω_1     = θmax/tmax_r1          # (un)Bending frequency
const ω_2     = 6.4*π/tmax              # Rotation frequency
const ω_3     = ω_2/2                 # Rotation frequency

function cat(t)
    body = copy(r_0)
    # Bend the torso
    for i in 2:9
        if i in (2,3,4,5)
            ax = [1.,0.,0.]
        elseif i in (6,7,8,9)
            ax = -[1.,0.,0.]
        end
        body[i] = rotate(body[i], axis=ax, angle=ω_1*min(t,tmax_r1), center=pos(body[1]))
    end
    # Total rotation
    body = map(x -> rotate(x, axis=[0,1,0], angle = ω_3*min(t,tmax), center=center_of_mass(body)), body)
    # Legs rotation
    for i in 2:9
        if i in (2,3,4,5)
            ax = pos(body[2]) - pos(body[1])
        elseif i in (6,7,8,9)
            ax = pos(body[1]) - pos(body[6])
        end
        body[i] = rotate(body[i], axis=ax, angle=ω_2 * min(t,tmax), center=pos(body[1]))
    end
    # Unbend torso
    if t >= tmax_r2
        for i in 2:9
            if i in (2,3,4,5)
                ax = [1.,0.,0.]
            elseif i in (6,7,8,9)
                ax = -[1.,0.,0.]
            end
            body[i] = rotate(body[i], axis=ax, angle=ω_1*(min(t,tmax)-tmax_r2), center=pos(body[1]))
        end
    end
    return body
end

# Create a model
model = Model( cat
             , 0.
             , tmax + 1.5
             , one(Quaternion)
             , zeros(3)
             )
println(model)
# And solve it
rotbodies, R, L = solve!(model)

# Lines to connect different parts of the body
# (Only important for plotting)
bodylines = [(1, 2),
             (2, 3),
             (2, 4),
             (2, 5),
             (1, 6),
             (6, 7),
             (6, 8),
             (6, 9)]

println("\nTrajectories are stored on the following variables:")
println("Cat's internal frame: model.bodyframe")
println("Inertial frame      : model.inertialframe")

# Save the image as a gif
println("\nSaving animation to \"gato.gif\"")
anim = plotmodel(model, :both,
                 bodylines=bodylines,
                 size=(800,400),
                 saveas="gato.gif")
