#################################################################
# Rigid Body: Points on a sphere with non-zero angular momentum #
#################################################################
using DeformableBodies
import LinearAlgebra: norm

# Uniform Distribution on Sphere
function unif_sphere(m)
    v = randn(m)
    return v / norm(v)
end

# Number of points
const npts = 100
# On bodyframe,
# body is made of npts points uniformly distributed on the unit sphere
const r_0 = [PointMass(rand()*10, unif_sphere(3)) for i in 1:npts]

# Create a model
model = Model( t -> r_0 # Constant trajectories
             , 0.
             , 5.
             , axistoquaternion(unif_sphere(3), rand()*2π) # Random initial rotation
             , randn(3) + 3 # Random normally distributed (~N(3,1)) angular momentum
             )
# And solve it
rotbodies, R, L = solve!(model)

println("Trajectories are stored on the following variables:")
println("Cat's internal frame: model.bodyframe")
println("Inertial frame      : model.inertialframe\n")

# Save the image as a mp4
println("\nPlotting result to \"sphere.mp4\"")
anim = plotmodel(model, :both,
                 saveas="sphere.mp4")
