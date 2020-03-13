"""
Submodule of `DeformableBodies.jl` implementing quaternion algebra.

Quaternions are a 4-dimensional normed division algebra which extend the complex numbers.
They may be used as a representation of rotations on 3-dimensional space.

# Exports
$(EXPORTS)

"""
module Quaternions

using Base: promote_rule, widen, float, big, show
using Base: ==, +, -, *, /, \
using Base: isreal, isinteger, isfinite, isnan, isinf
using Base: real, imag, conj, conj!, abs2, abs, inv, exp, log, angle

using LinearAlgebra: I
using DocStringExtensions

export Quaternion
export components,
       imagq,
       axis,
       normalize,
       axistoquaternion,
       quaterniontomatrix,
       matrixtoquaternion,
       rotate

# We define quaternions to be a set of 4 real numbers.
# Afterwards, we will define more structure such that
# they form a normed division involution algebra.
"""
    Quaternion{T<:Real} <: Number

Quaternion type with components of type `T`.

This type overloads all the arithmetic operations
as well as the methods defined for Complex numbers
that still make sense for Quaternions.
"""
struct Quaternion{T} <:Number where {T<:Real}
    q::Array{T, 1}
    function Quaternion{T}(x::AbstractVector) where T <: Real
        if length(x) == 4
            return new(x)
        elseif length(x) == 3
            return new([zero(T); x])
        else
            error("Quaternions are four-dimensional!\nDimension given: $(length(x))")
        end
    end
end

Quaternion{T}(t::Real, x::Real = 0, y::Real = 0, z::Real = 0) where {T<:Real} =
    Quaternion{T}([t,x,y,z])

Quaternion(x::AbstractVector{T}) where {T<:Real} = Quaternion{T}(x)
Quaternion(t::Real, x::Real = 0, y::Real = 0, z::Real = 0) = Quaternion([t,x,y,z])
Quaternion(t::Real, v::AbstractVector{<:Real}) = Quaternion([t; v])

Quaternion{T}(q::Quaternion) where {T<:Real} = Quaternion{T}(q.q)
Quaternion(q::Quaternion) = q

Quaternion{T}(z::Complex{T}) where T = Quaternion{T}(real(z), imag(z))
Quaternion(z::Complex) = Quaternion(real(z), imag(z))

Base.promote_rule(::Type{Quaternion{T}}, ::Type{S}) where {T,S<:Real} =
    Quaternion{promote_type(T,S)}
Base.promote_rule(::Type{Quaternion{T}}, ::Type{Quaternion{S}}) where {T,S<:Real} =
    Quaternion{promote_type(T,S)}

Base.widen(::Type{Quaternion{T}}) where {T} = Quaternion{widen(T)}

Base.float(::Type{Quaternion{T}}) where {T<:AbstractFloat} = Quaternion{T}
Base.float(::Type{Quaternion{T}}) where {T} = Quaternion{float(T)}

Base.big(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{big(T)}
Base.big(q::Quaternion{T}) where {T<:Real} = Quaternion{big(T)}(q)

##############
# Components #
##############

# Interface method for list conversion
"""
    components(q)

Return an array with the components of a [`Quaternion`](@ref).

# Example
```jldoctest
julia> components(Quaternion(1.0, 2.0, 3.0, 4.0))
4-element Array{Float64,1}:
 1.0
 2.0
 3.0
 4.0
```
"""
@inline components(q::Quaternion) = q.q

# Analogously to the Complex numbers,
# we define the real and imaginary parts of a Quaternion
@inline Base.real(q::Quaternion) = first(components(q))
@inline Base.imag(q::Quaternion) = components(q)[2:end]

"""
    imagq(q)

Return imaginary part of [`Quaternion`](@ref)
as a [`Quaternion`](@ref) with no real part.

# Examples
```jldoctest
julia> a = Quaternion(1,2,3,4)
1 + 2i + 3j + 4k

julia> imagq(a)
0 + 2i + 3j + 4k
```
"""
@inline imagq(q::Quaternion) = Quaternion(imag(q))

# Or a radius / angle decomposition

# Geometricaly, abs(q) represents the length of a quaternion
# when it is viewed as a vector in R^4
@inline Base.abs2(q::Quaternion) = sum(x ^ 2 for x in components(q))
@inline Base.abs(q::Quaternion)  = (sqrt ∘ abs2)(q)

# Abs of imaginary part of q
@inline _vecnorm(q::Quaternion) = abs(imagq(q))

@inline Base.angle(q::Quaternion) = 2 * atan(_vecnorm(q), real(q))

"""
    axis(q)

Return the unit vector on the direction
of the imaginary part of a [`Quaternion`](@ref).

# Examples
```jldoctest
julia> Quaternion(10,1,1,0.5)
10.0 + 1.0i + 1.0j + 0.5k

julia> axis(Quaternion(10,1,1,0.5))
3-element Array{Float64,1}:
 0.6666666666666666
 0.6666666666666666
 0.3333333333333333
```
"""
axis(q::Quaternion) = imag(q) / _vecnorm(q)

# Quaternions form an involution algebra.
@inline Base.conj(q::Quaternion) = Quaternion(real(q), -imag(q))

###############
# Comparisons #
###############

@inline Base.:(==)(a::Quaternion, b::Quaternion) = (components(a) == components(b))

@inline Base.isreal(q::Quaternion) = iszero(imag(q))
@inline Base.isinteger(q::Quaternion) = isreal(q) & isinteger(real(q))
@inline Base.isfinite(q::Quaternion) = all(isfinite, components(q))
@inline Base.isnan(z::Quaternion) = any(isnan, components(q))
@inline Base.isinf(z::Quaternion) = any(isinf, components(q))

#######################
# Quaternion algebra  #
#######################

# As it is essentially R^4, the Quaternions have a natural vector space structure.
@inline Base.:+(q::Quaternion) = Quaternion(+components(q))
@inline Base.:-(q::Quaternion) = Quaternion(-components(q))

@inline Base.:+(a::Quaternion, b::Quaternion) = Quaternion(components(a) .+ components(b))
@inline Base.:-(a::Quaternion, b::Quaternion) = Quaternion(components(a) .- components(b))

@inline Base.:*(k::Real, q::Quaternion) = Quaternion(k * components(q))
@inline Base.:*(q::Quaternion, k::Real) = k * q
# And also division by a real number by multiplication with inverse
@inline Base.:/(q::Quaternion, k::Real) = q * inv(k)
@inline Base.:\(k::Real, q::Quaternion) = q / k

# Analogously to Complex numbers, Quaternions also form an algebra.
# The difference in this case is that multiplication is not commutative.
@inline function Base.:*(a::Quaternion, b::Quaternion)
    x = components(a)
    y = components(b)
    return Quaternion(
        [ x[1]*y[1] - x[2]*y[2] - x[3]*y[3] - x[4]*y[4]
        , x[1]*y[2] + x[2]*y[1] + x[3]*y[4] - x[4]*y[3]
        , x[1]*y[3] + x[3]*y[1] + x[4]*y[2] - x[2]*y[4]
        , x[1]*y[4] + x[4]*y[1] + x[2]*y[3] - x[3]*y[2]
        ])
end

# Inverse quaternion: inv(q) = q^{-1}
@inline Base.inv(q::Quaternion) = conj(q) / abs2(q)

@inline Base.:(/)(a::Quaternion, b::Quaternion) = a * inv(b)
@inline Base.:(\)(a::Quaternion, b::Quaternion) = inv(a) * b


# When we talk about rotations,
# the norm of a quaternion is generally irrelevant.
# All that we need is it's projection onto the sphere.
"""
    normalize(q)

Return a [`Quaternion`](@ref) with the same direction as `q` but unit norm.

# Examples
```jldoctest
julia> q = Quaternion([1., 1., 1., 1.])
1.0 + 1.0i + 1.0j + 1.0k

julia> a = normalize(q)
0.5 + 0.5i + 0.5j + 0.5k

julia> abs(a)
1.0
```
"""
function normalize(q::Quaternion; tolerance=1e-8)
    local magnitude2 = abs2(q)
    if abs(magnitude2 - 1.0) > tolerance
        return q / sqrt(magnitude2)
    end
    return q
end

#############
# Functions #
#############

function Base.exp(q::Quaternion)
    nv = _vecnorm(q)
    return exp(real(q)) * Quaternion(cos(nv), sin(nv) * axis(q))
end

function Base.log(q::Quaternion)
    nq = abs(q)
    nv = _vecnorm(q)
    return Quaternion(log(nq), acos(real(q) / nq) * axis(q))
end

########
# Show #
########

# For aesthetic purposes, we improve how quaternions are printed:
function Base.show(io::IO, q::Quaternion)
    print(io, real(q))
    iscompact = get(io, :compact, false)
    for (x, e_i) in zip(imag(q), ['i', 'j', 'k'])
        if signbit(x) && !isnan(x)
            print(io, iscompact ? "-" : " - ")
        else
            print(io, iscompact ? "+" : " + ")
        end

        print(io, abs(x))

        if !(isa(x,Integer) && !isa(x,Bool) || isa(x,AbstractFloat) && isfinite(x))
            print(io, "*")
        end
        print(io, e_i)
    end
end

# Puts parentheses around quaternion on expressions
function show_unquoted(io::IO, q::Quaternion, ::Int, prec::Int)
    if operator_precedence(:+) <= prec
        print(io, "(")
        show(io, q)
        print(io, ")")
    else
        show(io, q)
    end
end

#############
# Rotations #
#############

"""
    axistoquaternion(axis, angle)

Receive an axis `v` and angle `θ`
and return the [`Quaternion`](@ref)
corresponding to a rotation of `θ` around `v`.
"""
function axistoquaternion end

function axistoquaternion(ax, angle::Real)
    if length(ax) != 3
        error("Error: Axis must be a 3-dimensional vector.\nDimension given is $(length(ax))")
    end
    return cos(angle/2) + normalize(Quaternion(ax)) * sin(angle/2)
end

@inline axistoquaternion(;axis, angle) = axistoquaternion(axis,angle)

"""
    quaterniontomatrix(q::Quaternion)

Return the rotation matrix associated with a [`Quaternion`](@ref).
"""
function quaterniontomatrix(q::Quaternion)
    q = normalize(q)
    t = real(q)
    v = imag(q)
    v_cross =  _crossmatrix(v)
    return kron(v',v) + t^2*I + 2*t*v_cross + v_cross^2
end

# Skew-symmetric matrix equivalent to cross product from the left
@inline _crossmatrix(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]

"""
    matrixtoquaternion(R)

Given a rotation matrix `R`,
return a quaternion `q` such that `rotate(v,q) = R*v` for all `v`.

The matrix `R` is assumed to be orthogonal
but, for efficiency reasons, no check is made to guarantee that.

Since there are, in general, two quaternions representing the same rotation matrix,
it is not guaranteed that `matrixtoquaternion ∘ quaterniontomatrix` equals the identity.
"""
function matrixtoquaternion(r)
    if     r[2,2] > -r[3,3] && r[1,1] > -r[2,2] && r[1,1] > -r[3,3]
        z = 1 + r[1,1] + r[2,2] + r[3,3]
        q = 0.5 / sqrt(z) * Quaternion(z, r[3,2]-r[2,3], r[1,3]-r[3,1], r[2,1]-r[1,2])
    elseif r[2,2] < -r[3,3] && r[1,1] >  r[2,2] && r[1,1] >  r[3,3]
        z = 1 + r[1,1] - r[2,2] - r[3,3]
        q = 0.5 / sqrt(z) * Quaternion(r[3,2]-r[2,3], z, r[2,1]+r[1,2], r[3,1]+r[1,3])
    elseif r[2,2] >  r[3,3] && r[1,1] <  r[2,2] && r[1,1] < -r[3,3]
        z = 1 - r[1,1] + r[2,2] - r[3,3]
        q = 0.5 / sqrt(z) * Quaternion(r[1,3]-r[3,1], r[2,1]+r[1,2], z, r[3,2]+r[2,3])
    else
        z = 1 - r[1,1] - r[2,2] + r[3,3]
        q = 0.5 / sqrt(z) * Quaternion(r[2,1]-r[1,2], r[3,1]+r[1,3], r[3,2]+r[2,3], z)
    end
    return q
end
# https://arxiv.org/pdf/math/0701759.pdf

"""
    rotate(v::Vector, q::Quaternion)
    rotate(v::Vector; axis, angle)

Rotate a vector `v` by a quaternion `q`.
The quaternion may be given directly or as an axis and an angle.
"""
function rotate end

function rotate(v, q::Quaternion)
    if length(v) != 3
        error("Error: Quaternions only rotate 3-dimensional vectors.\nDimension given is $(length(v))")
    end
    q = normalize(q)
    return imag(q * Quaternion(0, v) * conj(q))
end

@inline function rotate(v; axis, angle)
    return rotate(v, axistoquaternion(axis, angle))
end

end
