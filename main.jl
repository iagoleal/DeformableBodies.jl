#!/usr/bin/env julia
import DifferentialEquations.OrdinaryDiffEq
ODE = OrdinaryDiffEq

# Mathematical and physical definitions and methods
include("Quaternions.jl")
using .Quaternions

include("physics.jl")
# Plotting utilities
#= include("plots.jl") =#

# Struct to store the problem data before solving
struct Model
    trajectories::Array{Function,1} # Particle Trajectories
    t_min::Float64                  # Initial time
    t_max::Float64                  # Ending time
    q_0  ::Quaternion{Float64}      # Initial Rotation
    p_0  ::Array{Float64}           # Initial Angular Momentum
end

function eq_of_motion(du,u,bodies,t)
    ε = 1e-5
    q = Quaternion(u[1:4])
    p = u[5:end]
    particles = t .|> bodies
    v = [(x.pos .- y.pos) ./ ε for (x,y) in zip( (t+ε) .|> bodies, (t-ε) .|> bodies(t-ε))]
    Iinv = inv(inertia_tensor(particles))
    L = angular_momentum(particles, v)
    ω = Iinv*(p - L)

    dq = 0.5 * q * Quaternion(0, ω)
    dp = (Iinv*p) × (p-L)

    du[1:4] .= dq.q
    du[5:end] .= dp
    return du
end

@inline constructProblem(m::Model) =
    ODE.ODEProblem(eq_of_motion
                  , vcat(m.q_0.q, m.p_0)
                  , (m.t_min, m.t_max)
                  , m.trajectories
                  )

# Example model
# include("examples/falling-cat.jl")

# prob = constructProblem(model)
# sol = ODE.solve(prob)
