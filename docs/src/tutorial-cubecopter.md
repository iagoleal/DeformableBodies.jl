# Tutorial: The Cubecopter

In this page we model the movement of a simple deformable body
on a reference frame moving with it
and then solve the problem to find its trajectory from the point of view
of an inertial frame.

```@setup cube
using DeformableBodies
```

The example body will be a cube attached to rotating helix,
from now on comradely referred to as the _cubecopter_.
Each particle is represented by a [`PointMass`](@ref),
and the full body is an array of point masses.

```@example cube
r_0 = [
      PointMass(1., [ 1.,  1.,  0.])
    , PointMass(1., [ 1., -1.,  0.])
    , PointMass(1., [-1.,  1.,  0.])
    , PointMass(1., [-1., -1.,  0.])
    , PointMass(1., [ 1.,  1., -1.])
    , PointMass(1., [ 1., -1., -1.])
    , PointMass(1., [-1.,  1., -1.])
    , PointMass(1., [-1., -1., -1.])
    , PointMass(.5, [-1.,  0.,  .5])
    , PointMass(.5, [ 1.,  0.,  .5])
    ]
nothing; #hide
```

For visualization purposes,
the `plot` function from `Plots.jl`
is overloaded to work with arrays of point masses.


```@example cube
using Plots: plot
using Plots: gr; gr(size=(400,400), markercolor=:purple); #hide
plot(r_0)
```

It may be hard to visualize a lot of points scattered on a 3D graph.
Therefore, the `plot` function is also overloaded to accept a parameter `bodylines`
connecting some particles.
It **does not** interferes on the dynamics and serves only for aiding visualization.

Below is a nasty hack for calculating the cubecopter's edges.
You do not have to understand it.
The important part is that the variable `edges` is an array of *indexes*
representing which points should be connected.

```@example cube
edges = Tuple[]
for i = 1:length(r_0), j = i:length(r_0)
    if count(a -> first(a) == last(a), zip(pos(r_0[i]),pos(r_0[j]))) == 2
        push!(edges, (i,j))
    end
end
```

```@example cube
plot(r_0, bodylines=edges)
```

Much better, right?

Now, the array `r_0` represents a stationary body.
What we want is a trajectory over time,
represented by an array of functions.
Our chosen trajectory will be a stationary cube with an helix rotating counter-clockwise over it.

For the cube part, we represent these as constant functions.

```@example cube
trajectory = Function[]
for x in r_0[1:end-2]
    push!(trajectory, let x=x; t -> x;end)
end
nothing; #hide
```

While the helix movement is done using the function [`rotate`](@ref)
over a fixed axis with an angle varying over time.

```@example cube
const ω = 2*π/5.0            # Angular velocity
const z_axis  = [0., 0., 1.] # Axis of rotation is orthogonal to helix

for x in r_0[end-1:end]
    push!(trajectory, let x=x; t -> rotate(x, axis=z_axis, angle=ω*t); end)
end
nothing; #hide
```

With this,
we complete the code for the trajectory on the _body frame_.
Let's animate it so we can see how the trajectory behaves.

```@example cube
using Plots: Animation, frame, gif
anime = Animation()
for t in 0.:0.1:7.
    frame(anime, plot([ x(t) for x in trajectory], bodylines=edges))
end
gif(anime)
```

It is now time to define the problem's model and find the trajectory on the _inertial frame_.
First, we need to define some initial data.

```@example cube
tstart  = 0.                # Starting time
tend    = 5.                # Ending time
Rstart  = one(Quaternion)   # Initial rotation
L_cm    = zeros(3)          # Angular momentum as viewed from center of mass
nothing; #hide
```

The variables `tstart` and `tend` define
respectively the starting and ending time for the dynamics.
The other two variables are the initial data necessary to solve the differential equation.
The term `one(Quaternion)` is our way to represent the identity rotation,
meaning that at the starting time,
the inertial frame coincides with the body frame.
The variable `L_cm` is set to zero,
meaning that the total angular momentum is zero
from the point of view of the center of mass.

With this information, we are ready to define our [`Model`](@ref).

```@example cube
model = Model( trajectory
             , timespan = (tstart, tend)
             , q_0  = Rstart
             , L_cm = L_cm
             )
nothing; #hide
```

Now that everything is setted up,
finding the inertial frame trajectory is as simple as writing

```@example cube
solve!(model)
nothing; #hide
```

To visualize the final result,
we can use the function [`plotmodel`](@ref)
to plot and save an animation of the dynamics of
both frames side by side.


```@setup cube
# Ugly hack to make sure that no output comes from Plots.gif
gr(size=(800,400))
plotmodel( model
         , :both
         , saveas="cubecopter.gif"
         , bodylines=edges          # Only for visualization
         , duration=5.0             # Duration in seconds
         )
```
```julia
plotmodel( model
         , :both
         , saveas="cubecopter.gif"
         , bodylines=edges          # Only for visualization
         , duration=5.0             # Duration in seconds
         )
```

![See how it spins!](cubecopter.gif)

In the inertial frame,
the cube rotates in the opposite direction to the helix
guaranteeing that the total angular momentum is conserved.
