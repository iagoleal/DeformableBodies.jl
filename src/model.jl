#!/usr/bin/env julia
import DifferentialEquations.OrdinaryDiffEq
ODE = OrdinaryDiffEq

# Struct to store the problem data before solving
struct Model
    trajectories::Array{Function,1}   # Particle Trajectories
    t_min       ::Float64             # Initial time
    t_max       ::Float64             # Ending time
    q_0         ::Quaternion{Float64} # Initial Rotation
    p_0         ::Array{Float64,1}    # Initial Angular Momentum
end

function eq_of_motion!(du,u,bodies,t)
    q = Quaternion(u[1:4])
    p = u[5:end]
    # Change particles to SoR with CM at origin
    r = t .|> bodies |> centralize
    # Velocities must be calculated using a fixed SoR, so we don't centralize here
    v = velocity(bodies, t)
    Iinv = (inv ∘ inertia_tensor)(r)
    L = angular_momentum(r, v)
    ω = Iinv*(p - L)

    dq = 0.5 * q * Quaternion(ω)
    dp = (Iinv*p) × (p-L)

    du[1:4] .= dq.q
    du[5:end] .= dp
    return du
end

@inline construct_problem(m::Model) =
    ODE.ODEProblem(eq_of_motion!
                  , vcat(m.q_0.q, m.p_0)
                  , (m.t_min, m.t_max)
                  , m.trajectories
                  )

function solve(m::Model)
    prob = construct_problem(m)
    solution = ODE.solve( prob
                        #= , alg_hints=[:auto] =#
                        , ODE.Tsit5() # Solver
                        , reltol=1e-8
                        , abstol=1e-8
                        )

    # Evolution of rotations
    R(t) = Quaternion(solution(t)[1:4])
    rotbodies = [t -> PointMass(x(t).mass, rotate(R(t), x(t).pos)) for x in m.trajectories ]
    # Evolution of angular momentum
    momentum(t) = solution(t)[5:end]

    return rotbodies, momentum
end
