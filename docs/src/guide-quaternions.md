# Representing rotations by quaternions

```@setup quat
using DeformableBodies
```

The most usual ways to represent a rotation in three-dimensional space
are [rotations matrices](https://en.wikipedia.org/wiki/Rotation_matrix)
or [Euler angles](https://en.wikipedia.org/wiki/Euler_angles).
In this package, nevertheless,
they are represented using unit quaternions.

Quaternions are a four-dimensional non-commutative algebra
with the property that any rotation may be represented by a quaternion of norm 1.
Some of their advantages consist in that we only need 4 parameters to represent a quaternion
(as opposed to the 9 needed by rotation matrices)
and that they do not suffer from the phenomenon of
[gimbal lock](https://en.wikipedia.org/wiki/Gimbal_lock)
(as opposed to Euler angles).

If a quaternion ``q`` represents a rotation matrix ``R``,
their action on a vector is defined as
```math
R v = q v q^{-1}.
```

!!! warning
    This package assumes that your coordinate system follows the right-hand rule.
    If, for whatever reason, you want to use a left-handed system,
    ``q`` and ``q^{-1}`` must be transposed on the above formula.

## The `rotate` function

All rotations on `DeformableBodies.jl` are done by the function [`rotate`](@ref).
Its signature consists of

```julia
rotate(v::Vector, q::Quaternion, center::Vector)
```

This function rotates the vector `v` by the rotation represented by `q`
in the frame of reference centered on `center`.

If the third argument is left blank, it defaults to the origin.

```@repl quat
v = [3,0,0]
q = Quaternion(1,2,3,4)
rotate(v, q, [0,0,0])
rotate(v, q)
```

## Representing the identity rotation
The identity rotation is the matrix ``I`` satisfying ``I v = v``
for all ``v``.
This is represented by the quaternion ``1``,
which can be gotten using Julia's multiple dispatch via the expression
`one(Quaternion)`.

```@repl quat
q = one(Quaternion)
rotate([1,2,3], q)
```

## Composition of rotations

The composition of rotations translates to the quaternion world
as the product of quaternions.
Thus, it is the same to rotate a vector by ``q_1`` and then by ``q_2``
and to rotate by ``q_2 * q_1``.

```@repl quat
v = [9, 0, 0]
q1 = Quaternion(1,2,3,4)
q2 = Quaternion(4,3,2,1)
rotate(v, q2*q1)
rotate(rotate(v, q1), q2)
```

!!! note
    Rotations are non-commutative.
    Applying ``q_2`` after ``q_1`` is the same as multiplying
    ``q_2`` to the **left** of ``q_1``.

## Axis-angle representation

An important property of three-dimensional space it that every rotation fixes a line.
Therefore, they accept an axis-angle representation.
That is, a rotation ``R`` is defined by a unit vector ``\hat{n}`` and an angle ``\theta``
such that ``R`` rotates a vector counterclockwise by an amount of ``\theta`` around the line defined by ``\hat{n}``.

To get the quaternion representing an rotation of `θ` around `n`,
use the method [`axistoquaternion`](@ref).
You may pass an unormalized vector and the method will normalize it.

```@repl quat
q = axistoquaternion(axis = [0,0,2], angle = π/2)
rotate([1,0,0], q)
```

In fact, this combination is so useful that the function `rotate`
is overloaded to directly convert from an axis-angle pair.

```@repl quat
rotate([1,0,0], axis=[0,0,1], angle=π/2)
```

!!! info
    This version of [`rotate`](@ref) also accepts the optional argument `center`,
    which defaults to the origin.
    ```@repl quat
    rotate([1,0,0], axis=[0,0,1], angle=π/2, center=[0,1,0])
    ```

## Converting between matrices and quaternions

Sometimes a rotation may already come represented as a matrix
or the transition rotations from [`solve!`](@ref) may be needed in matrix form some application.
To deal with these cases,
the package provides two auxiliary functions
[`matrixtoquaternion`](@ref) and [`quaterniontomatrix`](@ref).

To convert a rotation matrix to a unit quaternion, do
```@repl quat
R = [cos(π/4) -sin(π/4) 0;
     sin(π/4)  cos(π/4) 0;
     0         0        1]
R * [sqrt(2), 0, 0]
q = matrixtoquaternion(R)
rotate([sqrt(2), 0, 0], q)
```
Some minor differences may occur due to floating-point rounding errors.

!!! danger
    The function [`matrixtoquaternion`](@ref)
    assumes that the input is a _rotation matrix_
    but, for efficency reasons, no check is done in this regard.
    If you do not make sure beforehand that the matrix is orthogonal,
    bad things may happen.

To convert a quaternion to a matrix simply do
```@repl quat
q = Quaternion(1,2,3,4)
R = quaterniontomatrix(q)
rotate([3,0,0], q)
R * [3,0,0]
```

The function [`quaterniontomatrix`](@ref)
works for every quaternion,
and does not require the input to be a unit quaternion.

!!! note
    There is a **unique** rotation matrix representing a given quaternion
    but there are **two** unit quaternions representing the same matrix.

    This means that `quaterniontomatrix ∘ matrixtoquaternion` equals the identity
    (disconsidering floating-point rouding errors)
    but the opposite is not in general true.
    For a simple example.
    ```@repl quat
    q = Quaternion(-1.0)
    (matrixtoquaternion ∘ quaterniontomatrix)(q)
    ```

    Nevertheless, both these quaternions produce the same rotations.

## Rotating a PointMass

All the previous functionalities only require the submodule [`Quaternions`](@ref)
and work directly with vectors.
Nevertheless, the models on `DeformableBodies.jl` are constructed with respect
to the type [`PointMass`](@ref).
To help with that,
the function rotate is overloaded to directly rotate the position of a `PointMass`
without interfering with its mass.

```julia
rotate(p::PointMass, q::Quaternion, center::Vector)
rotate(p::PointMass; axis::Vector, angle::Real, center::Vector)
```

The usage is identical to the version for vectors including the fact that
the argument `center` must be a Vector and not another PointMass.

An usual application consists in rotating a body around its center of mass.
```@repl quat
body = [ PointMass(rand(), rand(3)) for i in 1:5 ]
center_of_mass(body)
q = Quaternion(1,2,3,4)
rotated = [ rotate(p, q, center_of_mass(body)) for p in body ]
center_of_mass(rotated)
```
