using LinearAlgebra: I, cross, ×

"""
    PointMass(m, x)

Wrapper over a mass and a position on ``R^3``.

This type overloads [`Quaternions.rotate`](@ref)
to rotate only its position.
```jldoctest
julia> a = PointMass(10, [1., 0, 0])
PointMass{Float64}(10.0, [1.0, 0.0, 0.0])

julia> rotate(a; axis=[0., 0., 1.], angle=π/2)
PointMass{Float64}(10.0, [2.220446049250313e-16, 1.0, 0.0])
```
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
    velocity(xs, t; δ=1e-6)

Numerically approximate the velocity for a set `xs` of trajectories at time `t`.
The variable `δ` denotes the step for the finite differences interval.
"""
@inline velocity(xs, t; δ=1e-6) =
    map((a,b) -> (pos(a) - pos(b))/(2*δ), xs(t+δ), xs(t-δ))

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

@inline Quaternions.rotate(x::PointMass, q::Quaternion, center=zeros(3)) =
    PointMass(mass(x), Quaternions.rotate(pos(x), q, center))

@inline Quaternions.rotate(x::PointMass; angle, axis, center=zeros(3)) =
    PointMass(mass(x), Quaternions.rotate(pos(x); axis=axis, angle=angle, center=center))
