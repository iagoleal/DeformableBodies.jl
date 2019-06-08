# A PointMass is a mass at a given point on space (R^3).
"""
    PointMass(m, x)

Wrapper around a mass
and a position on ``R^3``.
"""
struct PointMass{T} <: Any where {T<:Real}
    mass::T
    pos::Vector{T}
    function PointMass{T}(m, x) where {T<:Real}
        m > 0 || error("Error: negative mass given.")
        length(x) == 3 || error("Error: PointMass must three dimensional.")
        return new(m, x)
    end
end

# Alternative constructors
PointMass(m::Real, x::Vector{<:Real}) = PointMass{typeof(promote(m,x[1])[1])}(m,x)

"""
    pos(p::PointMass)

Returns position of a [`PointMass`](@ref).
"""
pos(p::PointMass)  = p.pos

"""
    mass(p::PointMass)

Returns mass of a [`PointMass`](@ref).
"""
mass(p::PointMass) = p.mass

# Cross Product of 3-dimensional vectors
@inline cross(a,b) =
    [ a[2]*b[3] - a[3]*b[2]
    , a[3]*b[1] - a[1]*b[3]
    , a[1]*b[2] - a[2]*b[1]
    ]
@inline a × b = cross(a,b)

# Velocity of a set of trajectories
velocity(xs, t; ε=1e-6) = map((a,b) -> (a.pos - b.pos)/(2*ε), xs(t+ε), xs(t-ε))

"""
    center_of_mass(xs)

Gets a system of [`PointMass`](@ref)es and returns their
center of mass via
```math
cm(x) = \\frac{1}{\\sum m_i}\\sum m_i x_i.
```
"""
function center_of_mass(xs::AbstractVector{PointMass{T}}) where {T<:Real}
    cm = zeros(3)
    total_mass = 0.0
    for x in xs
        total_mass += x.mass
        cm += x.mass * x.pos
    end
    return cm / total_mass
end

"""
    inertia_tensor(xs)

Gets a system of [`PointMass`](@ref)es and returns their
inertia tensor via
```math
I = \\sum m_i \\langle x_i, x_i \\rangle \\text{id} - x_i \\otimes x_i.
```
"""
function inertia_tensor(xs::AbstractVector{PointMass{T}}) where {T<:Real}
    id = one(Array{T}(undef, 3, 3))
    return sum(x.mass * (x.pos'x.pos * id - kron(x.pos',x.pos)) for x in xs)
end
#= function inertia_tensor(xs::AbstractArray{PointMass{T},1}) where T <: Real =#
    #= I = zeros(3,3) =#
    #= for x in xs =#
    #=     I[1,1] += x.mass*( x.pos[2]^2 + x.pos[3]^2 ) =#
    #=     I[2,2] += x.mass*( x.pos[1]^2 + x.pos[3]^2 ) =#
    #=     I[3,3] += x.mass*( x.pos[1]^2 + x.pos[2]^2 ) =#

    #=     I[1,2] += -x.mass*x.pos[1]*x.pos[2] =#
    #=     I[1,3] += -x.mass*x.pos[1]*x.pos[3] =#
    #=     I[2,3] += -x.mass*x.pos[2]*x.pos[3] =#
    #= end =#
    #= I[2,1] = I[1,2] =#
    #= I[3,1] = I[1,3] =#
    #= I[3,2] = I[2,3] =#
    #= return I =#
#= end =#

"""
    angular_momentum(xs, vs)

Receives a system of [`PointMass`](@ref)es and their velocities, and returns their
angular momentum vector via
```math
L = \\sum m_i x_i \times v_i.
```
"""
@inline function angular_momentum(xs, vs)
    return sum(x.mass * (cross(x.pos, v)) for (x,v) in zip(xs,vs))
end

"""
    centralize(xs)

Receives a system of [`PointMass`](@ref)es and returns
the same system translated such that their center of mass
is fixed on the origin.
"""
function centralize(xs::AbstractVector{PointMass{T}}) where {T<:Real}
    cm = center_of_mass(xs)
    return map(r -> PointMass(r.mass, r.pos - cm), xs)
end

function centralize!(xs::AbstractVector{PointMass{T}}) where {T<:Real}
    cm = center_of_mass(xs)
    for i in 1:length(xs)
        xs[i] -= cm
    end
    return xs
end

# Extend Quaternion.rotate to deal with point masses
using .Quaternions: rotate

@inline function Quaternions.rotate(x::PointMass; angle=0, axis=[0,0,0])
    return PointMass(x.mass, rotate(x.pos; angle=angle, axis=axis))
end
