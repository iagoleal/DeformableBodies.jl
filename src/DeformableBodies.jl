module DeformableBodies

# Quaternions definitions
include("Quaternions.jl")
using .Quaternions

export Quaternion
export imagq,
       axis,
       normalize,
       normalize!,
       axis2quaternion,
       rotate

# Physical utilities
include("physics.jl")

export PointMass
export pos,
       mass,
       centralize,
       centralize!,
       center_of_mass,
       inertia_tensor,
       angular_momentum

# Model definitions
include("model.jl")

export Model
export solve!

# Plotting utilities
include("plots.jl")

export plotmodel
end
