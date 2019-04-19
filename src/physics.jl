# A PointMass is a mass at a given point on space (R^3).
struct PointMass{T} <:Any where T <: Real
    mass::T
    pos::Vector{T}
    function PointMass{T}(m::T, x::Vector) where T <: Real
        if m <= 0
            error("Error: negative mass given.")
        end
        if length(x) != 3
            error("Error: PointMass must three dimensional.")
        end
        new(m, x)
    end
end

# Alternative constructors
PointMass(m::T, x::Vector{T}) where {T <: Real} = PointMass{T}(m,x)

# Cross Product of 3-dimensional vectors
@inline cross(a,b) =
    [ a[2]*b[3] - a[3]*b[2]
    , a[3]*b[1] - a[1]*b[3]
    , a[1]*b[2] - a[2]*b[1]
    ]
@inline a × b = cross(a,b)

function center_of_mass(xs::AbstractArray{PointMass{T}, 1}) where T <: Real
    cm = zeros(3)
    total_mass = 0.
    for x in xs
        total_mass += x.mass
        cm += x.mass * x.pos
    end
    return cm / total_mass
end

function inertia_tensor(xs::AbstractArray{PointMass{T},1}) where T <: Real
    I = zeros(3,3)
    for x in xs
        I[1,1] += x.mass*( x.pos[2]^2 + x.pos[3]^2 )
        I[2,2] += x.mass*( x.pos[1]^2 + x.pos[3]^2 )
        I[3,3] += x.mass*( x.pos[1]^2 + x.pos[2]^2 )

        I[1,2] += -x.mass*x.pos[1]*x.pos[2]
        I[1,3] += -x.mass*x.pos[1]*x.pos[3]
        I[2,3] += -x.mass*x.pos[2]*x.pos[3]
    end
    I[2,1] = I[1,2]
    I[3,1] = I[1,3]
    I[3,2] = I[2,3]
    return I
end

@inline angular_momentum(xs::AbstractArray{PointMass{T},1}, vs::AbstractArray) where T<:Real =
    sum( x.mass * (cross(x.pos, v)) for (x,v) in zip(xs,vs) )

function centralize(xs::AbstractArray{PointMass{T},1}) where T<: Real
    cm = center_of_mass(xs)
    map(r -> PointMass(r.mass, r.pos - cm), xs)
end

# Velocity of a set of trajectories
function velocity(bodies::Union{Function, AbstractVector{Function}}, t; ε=1e-6)
    map( (a,b) -> (a.pos - b.pos)/ ε, (t+ε) .|> bodies , (t-ε) .|> bodies )
end

# Extend Quaternion.rotate to deal with point masses
import .Quaternions: rotate

@inline rotate(x::PointMass; angle=0, axis=[0,0,0]) = let a = angle, b =axis
        PointMass(x.mass, rotate(x.pos, angle=a, axis=b))
end
