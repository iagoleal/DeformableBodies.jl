module Quaternions

import Base: promote_rule, widen, float, big, show
import Base: ==, +, -, *, /, \
import Base: isreal, isinteger, isfinite, isnan, isinf
import Base: real, imag, conj, conj!, abs2, abs, inv, exp, log, angle

export Quaternion
export imagq,
       axis,
       normalize,
       normalize!,
       axis2quaternion,
       rotate

# We define quaternions to be a set of 4 real numbers.
# Afterwards, we will define more structure such that
# they form a normed division involution algebra.
"""
    Quaternion{T<:Real} <: Number

    Quaternion type with components of type `T`.
"""
struct Quaternion{T} <:Number where T <: Real
    q::Array{T, 1}
    function Quaternion{T}(x::AbstractVector) where T <: Real
        if length(x) == 4
            new(x)
        elseif length(x) == 3
            new([zero(T); x])
        else
            error("Quaternions are four-dimensional!\nDimension given: $(length(x))")
        end
    end
end

Quaternion{T}(t::Real, x::Real = 0, y::Real = 0, z::Real = 0) where {T <: Real} =
    Quaternion{T}([t,x,y,z])

Quaternion(x::AbstractVector{T}) where {T <: Real} = Quaternion{T}(x)
Quaternion(t::Real, x::Real = 0, y::Real = 0, z::Real = 0) = Quaternion([t,x,y,z])
Quaternion(t::Real, v::AbstractVector{<:Real}) = Quaternion([t; v])

Quaternion{T}(q::Quaternion) where T <: Real = Quaternion{T}(q.q)
Quaternion(q::Quaternion) = q

Quaternion{T}(z::Complex{T}) where T = Quaternion{T}(real(z), imag(z))
Quaternion(z::Complex) = Quaternion(real(z), imag(z))


Base.promote_rule(::Type{Quaternion{T}}, ::Type{S}) where {T<:Real,S<:Real} =
    Quaternion{promote_type(T,S)}
Base.promote_rule(::Type{Quaternion{T}}, ::Type{Quaternion{S}}) where {T<:Real,S<:Real} =
    Quaternion{promote_type(T,S)}

Base.widen(::Type{Quaternion{T}}) where {T} = Quaternion{widen(T)}

Base.float(::Type{Quaternion{T}}) where {T<:AbstractFloat} = Quaternion{T}
Base.float(::Type{Quaternion{T}}) where {T} = Quaternion{float(T)}

Base.big(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{big(T)}
Base.big(q::Quaternion{T}) where {T<:Real} = Quaternion{big(T)}(q)

##############
# Components #
##############

# Analogous to the Complex numbers,
# we define the real part and imaginary part of q
@inline Base.real(q::Quaternion) = q.q[1]
@inline Base.imag(q::Quaternion) = q.q[2:end]

"""
    imagq(q)

Returns imaginary part of quaternion `q`
as a quaternion.

```jldoctest
julia> a = Quaternion(1,2,3,4)
1 + 2i + 3j + 4k

julia> imagq(a)
0 + 2i + 3j + 4k
```
"""
@inline imagq(q::Quaternion) = Quaternion(q.q[2:end])

# Or a radius / angle decomposition

# Geometricaly, abs(q) represents the length of a quaternion
# when it is viewed as a vector in R^4
@inline Base.abs2(q::Quaternion) = sum( x ^ 2 for x in q.q)
@inline Base.abs(q::Quaternion)  = (sqrt ∘ abs2)(q)

# Abs of imaginary part of q
@inline _vecnorm(q::Quaternion) = (sqrt ∘ sum)( x ^ 2 for x in imag(q))

@inline Base.angle(q::Quaternion) = 2 * atan(_vecnorm(q), real(q))

function axis(q::Quaternion)
    v = imag(q)
    return v / _vecnorm(v)
end

# Quaternions form an involution algebra.
# The conjugation operation is defined as $conj(re, im) = (re, -im)$
@inline Base.conj(q::Quaternion)  = Quaternion(real(q), -imag(q))
@inline function Base.conj!(q::Quaternion)
    q.q[2:end] .*= -1
    return q
end

###############
# Comparisons #
###############

@inline Base.:(==)(a::Quaternion, b::Quaternion) = (a.q == b.q)

@inline Base.isreal(q::Quaternion) = iszero(imag(q))
@inline Base.isinteger(q::Quaternion) = isreal(q) & isinteger(real(q))
@inline Base.isfinite(q::Quaternion) = all(isfinite, q.q)
@inline Base.isnan(z::Quaternion) = any(isnan, q.q)
@inline Base.isinf(z::Quaternion) = any(isinf, q.q)

#######################
# Quaternion algebra  #
#######################


# As it is essentially R^4, the Quaternions have a natural vector space structure.
#
@inline Base.:+(q::Quaternion) = Quaternion(+q.q)
@inline Base.:-(q::Quaternion) = Quaternion(-q.q)

@inline Base.:+(a::Quaternion, b::Quaternion) = Quaternion(a.q .+ b.q)
@inline Base.:-(a::Quaternion, b::Quaternion) = Quaternion(a.q .- b.q)

@inline Base.:*(k::Real, q::Quaternion) = Quaternion(k * q.q)
@inline Base.:*(q::Quaternion, k::Real) = k * q
# And also division by a real number by multiplication with inverse
@inline Base.:/(q::Quaternion, k::Real) = q * inv(k)
@inline Base.:\(k::Real, q::Quaternion) = q / k

# Analogously to Complex numbers,
# Quaternions also form an algebra.
# The difference in this case is
# that multiplication is not commutative.
@inline Base.:*(q::Quaternion, w::Quaternion) =
    Quaternion(
    [ q.q[1]*w.q[1] - q.q[2]*w.q[2] - q.q[3]*w.q[3] - q.q[4]*w.q[4]
    , q.q[1]*w.q[2] + q.q[2]*w.q[1] + q.q[3]*w.q[4] - q.q[4]*w.q[3]
    , q.q[1]*w.q[3] + q.q[3]*w.q[1] + q.q[4]*w.q[2] - q.q[2]*w.q[4]
    , q.q[1]*w.q[4] + q.q[4]*w.q[1] + q.q[2]*w.q[3] - q.q[3]*w.q[2]
    ])

# Inverse quaternion: inv(q) = q^{-1}
@inline Base.inv(q::Quaternion) = conj(q) / abs2(q)

@inline Base.:(/)(q1::Quaternion, q2::Quaternion) = q1 * inv(q2)
@inline Base.:(\)(q1::Quaternion, q2::Quaternion) = inv(q1) * q2


# When we talk about rotations,
# the norm of a quaternion is generally irrelevant.
# All that we need is it's projection onto the sphere.
"""
    normalize(q)

Divides `q` by its norm to return a unit quaternion ``q / |q|``.

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
    magnitude2 = abs2(q)
    if abs( magnitude2 - 1.0 ) < tolerance
        return q
    end
    return q / sqrt(magnitude2)
end
function normalize!(q::Quaternion; tolerance=1e-8)
    magnitude2 = abs2(q)
    if abs( magnitude2 - 1.0 ) > tolerance
        magnitude = sqrt(magnitude2)
        q.q ./= magnitude
    end
    return q
end

#############
# Functions #
#############

function Base.exp(q::Quaternion)
    nv = _vecnorm(q)
    exp(real(q)) * Quaternion(cos(nv), sin(nv) * axis(q))
end

function Base.log(q::Quaternion)
    nq = abs(q)
    nv = _vecnorm(q)
    Quaternion( log(nq), acos(real(q) / nq) * axis(q) )
end

########
# Show #
########

# For aesthetic purposes, we improve how quaternions are printed:
function Base.show(io::IO, q::Quaternion)
    print(io, real(q))
    iscompact = get(io, :compact, false)
    for (x, e_i) in zip(q.q[2:end], ['i', 'j', 'k'])
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

# Any rotatation may be represented by a quaternion.
# If we have a direction on the sphere and an angle θ,
# it gives us an rotation matrix R(v,θ) which
# may be represented according to
"""
    axis2quaternion(axis, angle)

Receives an axis `v` and angle `θ`
and returns the quaternion
who corresponds to a rotation of `θ` around `v`.
"""
function axis2quaternion(ax::AbstractVector,angle::Real)
    if length(ax) != 3
        error("Error: Axis must be a 3-dimensional vector.\nDimension given is $(length(ax))")
    end
    cos(angle/2) + normalize(Quaternion(ax)) * sin(angle/2)
end

"""
    rotate(q::Quaternion, v::Vector)
    rotate(v::Vector; axis, angle)

Rotate a vector `v` by a quaternion `q`.
The quaternion may be given directly
or as an axis and an angle.
"""
function rotate(q::Quaternion, v::AbstractVector)
    if length(v) != 3
        error("Error: Quaternions only rotate 3-dimensional vectors.\nDimension given is $(length(v))")
    end
    q = normalize(q)
    imag( q * Quaternion(0, v) * conj(q) )
end
@inline function rotate( v::AbstractVector
                       ; angle=0
                       , axis=[0, 0, 1]
                       )
    rotate(axis2quaternion(axis, angle), v)
end

end
