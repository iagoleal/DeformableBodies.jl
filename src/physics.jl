using LinearAlgebra: I, cross, ×

"""
    PointMass(m, x)

Wrapper over a mass and a position on ``R^3``.
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

# Access functions
"""
    pos(p::PointMass)

Return position of a [`PointMass`](@ref).
"""
@inline pos(p::PointMass) = p.pos

"""
    mass(p::PointMass)

Return mass of a [`PointMass`](@ref).
"""
@inline mass(p::PointMass) = p.mass

"""
    velocity(xs, t; ε=1e-6)

Numerically approximate the velocity for a set `xs` of trajectories at time `t`.
The variable `ε` denotes the desired precision.
"""
velocity(xs, t; ε=1e-6) = map((a,b) -> (pos(a) - pos(b))/(2*ε), xs(t+ε), xs(t-ε))

@doc raw"""
    center_of_mass(xs)

Receive a system of [`PointMass`](@ref)es
and return their center of mass through formula
```math
cm(x) = \frac{1}{\sum m_i}\sum m_i x_i.
```
"""
@inline center_of_mass(xs) = sum(mass.(xs) .* pos.(xs)) / sum(mass.(xs))

@doc raw"""
    inertia_tensor(xs)

Receive a system of [`PointMass`](@ref)es
and return their inertia tensor through formula
```math
I(x) = \sum m_i \langle x_i, x_i \rangle \mathrm{id} - x_i \otimes x_i.
```
"""
@inline inertia_tensor(xs) =
    sum(mass(x) * (pos(x)'pos(x) * I - kron(pos(x)',pos(x))) for x in xs)

@doc raw"""
    angular_momentum(xs, vs)

Receive a system of [`PointMass`](@ref)es and their velocities,
and return their angular momentum vector through formula
```math
L(x) = \sum m_i x_i \times v_i.
```
"""
@inline angular_momentum(xs, vs) =
    sum(mass(x) * (cross(pos(x), v)) for (x,v) in zip(xs,vs))

"""
    centralize(xs)

Receive a system of [`PointMass`](@ref)es
and translate it such that the center of mass
is fixed on the origin.
"""
@inline centralize(xs) = map(r -> PointMass(mass(r), pos(r) - center_of_mass(xs)), xs)

# Extend Quaternions.rotate to deal with point masses
using .Quaternions: rotate

@inline Quaternions.rotate(x::PointMass, q::Quaternion) =
    PointMass(mass(x), Quaternions.rotate(pos(x), q))

@inline Quaternions.rotate(x::PointMass; angle, axis) =
     PointMass(mass(x), Quaternions.rotate(pos(x); angle=angle, axis=axis))
