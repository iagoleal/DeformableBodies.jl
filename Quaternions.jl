module Quaternions

export Quaternion,
       normalize, normalize!,
       axis2Quaternion,
       rotate

# We define quaternions to be a set of 4 real numbers.
# Aftewards, we will define more structure such that
# they form a normed division involution algebra.
struct Quaternion{T} <:Number where T <: Real
    q::Array{T, 1}
    function Quaternion{T}(x::Array{T, 1}) where T <: Real
        if length(x) != 4
            error("Quaternions are four-dimensional!")
        end
        new(x)
    end
end

# and add some others constructors for easiness.
Quaternion(x::Array{T, 1}) where {T <: Real}   = Quaternion{T}(x)
Quaternion(t::Real, x::Real, y::Real, z::Real) = Quaternion([t,x,y,z])
Quaternion(t::Real, v::Array{<:Real, 1 })      = Quaternion([t, v[1], v[2], v[3]])

#######################
# Quaternion algebra  #
#######################

# As it is essentially R^4,
# the Quaternions have a natural vector space structure.
# That is, we define a sum between Quaternions
# as the sum of their components
@inline Base.:+(a::Quaternion, b::Quaternion) = Quaternion(a.q + b.q)

# and define a neuter element given by the function
@inline Base.zero(a::Type{Quaternion{T}}) where T = Quaternion{T}([zero(T),zero(T),zero(T),zero(T)])
@inline Base.zero(a::Type{Quaternion})            = zero(Quaternion{Int64})
@inline Base.zero(a::Quaternion{T})       where T = zero(Quaternion{T})

# We also give the quaternions a scalar mutiplication structure given by
@inline Base.:*(k::Real, q::Quaternion) = Quaternion(k * q.q)

# To guarantee commutativity, we define that $k * q = q * k$.
@inline Base.:*(q::Quaternion, k::Real) = k * q

# And also division by a real number by multiplication with inverse
@inline Base.:/(q::Quaternion, k::Real) = q * (1/k)

# The other axioms of a vector space follow accordingly.

# Analogously to Complex numbers,
# Quaternions also form an algebra.
# The difference in this case is
# that multiplication is not commutative.

@inline Base.:*(q1::Quaternion, q2::Quaternion) =
    Quaternion(
    [ q1.q[1]*q2.q[1] - q1.q[2]*q2.q[2] - q1.q[3]*q2.q[3] - q1.q[4]*q2.q[4]
    , q1.q[1]*q2.q[2] + q1.q[2]*q2.q[1] + q1.q[3]*q2.q[4] - q1.q[4]*q2.q[3]
    , q1.q[1]*q2.q[3] + q1.q[3]*q2.q[1] + q1.q[4]*q2.q[2] - q1.q[2]*q2.q[4]
    , q1.q[1]*q2.q[4] + q1.q[4]*q2.q[1] + q1.q[2]*q2.q[3] - q1.q[3]*q2.q[2]
    ])

# The multiplication of quaternions has an identity element given by
@inline Base.one(a::Type{Quaternion{T}}) where T = Quaternion{T}([one(T),zero(T),zero(T),zero(T)])
@inline Base.one(a::Type{Quaternion})            = one(Quaternion{Int64})
@inline Base.one(a::Quaternion{T})       where T = one(Quaternion{T})

# Analogous to the Complex numbers,
# we define the real part of q by
@inline Base.real(q::Quaternion) = q.q[1]
# and the imaginary part by
@inline Base.imag(q::Quaternion) = q.q[2:end]

# Quaternions form an involution algebra.
# The conjugation operation is defined as $conj(re + im) = re - im$,
# just as complex numbers
@inline Base.conj(q::Quaternion) = Quaternion(real(q), -imag(q))
@inline function Base.conj!(q::Quaternion)
    for i = 2:length(q.q)
        q.q[i] *= -1
    end
    return q
end

# It is also possible to talk about the norm of a quaternion.
# Geometricaly, abs(q) represents the length of a quaternion
# when it is viewed as a vector in R^4
@inline Base.abs2(q::Quaternion) = sum(q.q .^ 2)
@inline Base.abs(q::Quaternion)  = sqrt(sum(q.q .^ 2))

# Every quaternion has a inverse given by
@inline Base.inv(q::Quaternion) = conj(q) / abs2(q)
# such that they form a normed division algebra.

# When we talk about rotations,
# the norm of a quaternion is generally irrelevant.
# All that we need is it's projection on the sphere.
function normalize(q::Quaternion; tolerance=1e-5)
    magnitude2 = abs2(q)
    if abs( magnitude2 - 1.0 ) < tolerance
        return q
    end
    return q / sqrt(magnitude2)
end

function normalize!(q::Quaternion; tolerance=1e-5)
    magnitude2 = abs2(q)
    if abs( magnitude2 - 1.0 ) > tolerance
        magnitude = sqrt(magnitude2)
        for i in 1:length(q.q)
            q.q[i] = q.q[i]/magnitude
        end
    end
    return q
end

########
# Show #
########

# For aesthetic purposes,
# we improve how quaternions are printed:
function Base.show(io::IO, q::Quaternion)
    output="$(q.q[1])"
    for (x, basis_idx) in zip(q.q[2:end], ['i', 'j', 'k'])
        output *= (x >= 0 ? " + " : " - " ) * "$(abs(x))$(basis_idx)"
    end
    print(io, output)
end

#############
# Rotations #
#############

# Any rotatation may be represented by a quaternion.
# If we have a direction on the sphere and an angle θ,
# it gives us an rotation matrix R(v,θ) which
# may be represented according to
function axis2quaternion(axis::Array{<:Real, 1},angle::Real)
    v = imag(normalize(Quaternion(0,axis)))
    Quaternion(cos(angle/2), v*sin(angle/2))
end

# To apply the rotation represented by a quaternion to a vector,
# we use the following function
function rotate(q::Quaternion, v::Array{<:Real, 1})
    if length(v) != 3
        error("Error: Quaternions only rotate 3-dimensional vectors.\nDimension given is $(length(v))")
    end
    q = normalize(q)
    imag( q * Quaternion(0, v) * conj(q) )
end

# An alias to rotate a vector v
# by an angle around an axis
@inline rotate(v::Array{<:Real,1}; angle=0, axis=[0,0,0]) = rotate(axis2quaternion(axis, angle), v)
end
