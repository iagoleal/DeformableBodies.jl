import OrdinaryDiffEq
const ODE = OrdinaryDiffEq

"""
    Model

Stores the data of a deformable body problem
before and after solving.
"""
mutable struct Model
    bodyframe    ::Function            # Trajectories on body frame
    t_min        ::Float64             # Initial time
    t_max        ::Float64             # Ending time
    q_0          ::Quaternion{Float64} # Initial Rotation
    L_cm         ::Array{Float64,1}    # Center of Mass Angular Momentum
    inertialframe::Function            # Trajectories on inertial frame
    Model(trjs, t0, tf, q_0, L_cm) = new(trjs, t0, tf, normalize(q_0), L_cm)
end

function Model(trajectories::Array{Function,1}, t_min, t_max, q_0, L_cm)
    trajs = t -> [xi(t) for xi in trajectories]
    return Model(trajs, t_min, t_max, normalize(q_0), L_cm)
end

function eq_of_motion!(du,u,trajectories,t)
    q = Quaternion(u[1:4])
    Π = u[5:end]
    # Change particles to SoR with CM at origin
    r = centralize(trajectories(t))
    # Velocities must be calculated using a fixed SoR, so we don't centralize here
    v = velocity(trajectories, t)
    Iinv = (inv ∘ inertia_tensor)(r)
    L = angular_momentum(r, v)
    ω = Iinv * (Π - L)

    dq = 0.5 * q * Quaternion(ω)
    dΠ = Π × (Iinv * (Π - L))

    du[1:4]   .= components(dq)
    du[5:end] .= dΠ
    return du
end

@inline construct_problem(m::Model) =
    ODE.ODEProblem(eq_of_motion!
                   , vcat(components(m.q_0), rotate(m.L_cm, conj(m.q_0)))
                  , (m.t_min, m.t_max)
                  , m.bodyframe
                  )

"""
    solve!(m::Model; reltol, abstol, solver)

Receives a [`Model`](@ref), calculates the
trajectory of the body on a inertial frame
and stores it in the variable `m.inertialframe`.
"""
function solve!(m::Model; reltol=1e-8, abstol=1e-8, solver=ODE.Tsit5())
    prob = construct_problem(m)
    solution = ODE.solve(prob, solver, reltol=reltol, abstol=abstol)
    # Evolution of rotations
    R(t) = Quaternion(solution(t)[1:4])
    # Evolution of angular momentum
    momentum(t) = solution(t)[5:end]
    # Store solution on model
    m.inertialframe = t -> [PointMass(mass(x), rotate(pos(x), R(t))) for x in m.bodyframe(t)]
    return m.inertialframe, R, momentum
end

# Printing Models
using Base: show

function Base.show(io::IO, m::Model)
    print(io, "Model(",m.t_min,", ",m.t_max,", ",m.q_0,", ",m.L_cm,")")
end

function Base.show(io::IO, ::MIME"text/plain", m::Model)
    println(io, "Deformable Body Model")
    println(io, "Initial Time: ", m.t_min)
    println(io, "Final Time  : ", m.t_max)
    println(io, "Initial Data are\n    Rotation: ", m.q_0, "\n    Angular Momentum: ", m.L_cm)
    print(io, "Number of points: ", length(m.bodyframe(m.t_min)))
end
