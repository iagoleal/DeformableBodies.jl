# Deformable body dynamics simulator

## Dependencies
- Julia 1.1.0 (or greater)
- DifferentialEquations.jl
- Plots.jl + some plotting backend

## Usage

Define a new _Model_ and store it in a variable.
The constructors receives:
- An array of functions f: R -> PointMass;
- A initial time;
- A final time;
- A Quaternion who gives the initial rotation;
- A vector who gives the initial angular momentum on the inertial frame.

E.g.
```julia
# xs = ... some functions
model = Model(xs, 0., 10., one(Quaternion), zeros(3))
```

Upon solving,
the program will return
the trajectories and angular momenta
on the inertial frame.

```julia
ys, Ls = solve(model)
```

Now compare the trajectories
in `xs` and `ys` to see the physics happening!

The folder `examples/`
contains some examples of how to construct a model.

## Mathematics behind the program
_Geometry goes here._

## TODO
- Write better documentation;
- Plotting utilities;
- Find trajectory who optimizes energy from a parametrized family.

## Disclaimer
No real cats were harmed during the development of this program.
