module DeformableBodies

# Quaternions definitions
include("Quaternions.jl")
using .Quaternions

export Quaternion
export components,
       imagq,
       axis,
       normalize,
       axistoquaternion,
       quaterniontomatrix,
       matrixtoquaternion,
       rotate

# Physical utilities
include("physics.jl")

export PointMass
export pos,
       mass,
       centralize,
       center_of_mass,
       inertia_tensor,
       angular_momentum

# Model definitions
include("model.jl")

export Model
export solve!

# Plotting utilities
include("plots.jl")

export plotmodel,
       saveanimation
end
