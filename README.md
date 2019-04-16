# Deformable body dynamics simulator

## Dependencies
- Julia 1.1.0 (or greater)
- DifferentialEquations.jl
- Plots.jl + some plotting backend

## How to use

Define a new _Model_ and run

```julia
prob = constructProblem(model)
solution = ODE.solve(prob)
```

## Mathematics behind the program
_Geometry goes here._

## TODO

- Write better documentation;
- Plotting utilities;
- Organize main file;
- Find trajectory who optimizes energy from a parametrized family.

## Disclaimer
No real cats were harmed during the development of this program.
