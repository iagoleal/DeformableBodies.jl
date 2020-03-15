import OrdinaryDiffEq
const ODE = OrdinaryDiffEq

"""
    Model

**Fields**

- `bodyframe` -- Function representing trajectory on reference frame rotating with the body.
- `timespan`  -- Tuple containing starting and ending time of motion.
- `q_0`       -- Quaternion representing initial rotation.
- `L_cm`      -- Vector of angular momentum in relation to center of mass.

Only Initialized after solve:

- `inertialframe` -- Function representing trajectory on inertial reference frame.
- `rotation`      -- Rotations that exchange between `bodyframe` and `inertialframe`.
- `momentum`      -- Internal angular momentum.

Store the data of a deformable body problem before and after solving.
"""
mutable struct Model
    bodyframe    ::Function               # Trajectories on body frame
    timespan     ::Tuple{Float64,Float64} # Initial and Ending time
    q_0          ::Quaternion{Float64}    # Initial Rotation
    L_cm         ::Array{Float64,1}       # Center of Mass Angular Momentum
    inertialframe::Function               # Trajectories on inertial frame
    rotation     ::Function               # Rotations from body to inertial frame
    momentum     ::Function               # Evolution of momentum variables
    Model(trjs, tspan, q_0, L_cm) = new(trjs, tspan, normalize(q_0), L_cm)
end

function Model(trajectories::Array{Function,1}, tspan, q_0, L_cm)
    trajs = t -> [xi(t) for xi in trajectories]
    return Model(trajs, tspan, normalize(q_0), L_cm)
end

function Model(trajectories; timespan, q_0, L_cm)
    return Model(trajectories, timespan, q_0, L_cm)
end

function eq_of_motion!(du,u,trajectories,t)
    q    = Quaternion(u[1:4])        # Rotation variable
    Π    = u[5:end]                  # Momentum variable
    r    = trajectories(t)           # Position
    v    = velocity(trajectories, t) # Velocity
    Iinv = (inv ∘ inertia_tensor)(r) # Inverse of moment of inertia
    L    = angular_momentum(r, v)    # Angular momentum
    ω    = (Iinv * (Π - L))          # Angular velocity

    dq = 0.5 * q * Quaternion(ω)
    dΠ = Π × ω

    du[1:4]   .= components(dq)
    du[5:end] .= dΠ
    return du
end

@inline construct_problem(m::Model) =
    ODE.ODEProblem(eq_of_motion!
                  , vcat(components(m.q_0), rotate(m.L_cm, conj(m.q_0))) # Initial conditions
                  , m.timespan               # Timespan
                  , centralize ∘ m.bodyframe # Movement must be centered on cm
                  )

"""
    solve!(m::Model; reltol=1e-8, abstol=1e-8, solver=Tsit5())

Receive a [`Model`](@ref), calculate the
trajectory of the body on an inertial frame
and store it in the variable `m.inertialframe`.
"""
function solve!(m::Model; reltol=1e-8, abstol=1e-8, solver=ODE.Tsit5())
    prob = construct_problem(m)
    solution = ODE.solve(prob, solver, reltol=reltol, abstol=abstol)
    # Evolution of rotations
    m.rotation = t -> Quaternion(solution(t)[1:4])
    # Evolution of angular momentum
    m.momentum = t -> solution(t)[5:end]
    # Store solution on model
    m.inertialframe = t -> map(x -> rotate(x, m.rotation(t), center_of_mass(m.bodyframe(t))), m.bodyframe(t))
    return m.inertialframe, m.rotation, m.momentum
end

# Printing Models
using Base: show

function Base.show(io::IO, m::Model)
    print(io, "Model(", m.timespan, ", ", ", ", m.q_0, ", ", m.L_cm, ")")
end

function Base.show(io::IO, ::MIME"text/plain", m::Model)
    println(io, "Deformable Body Model")
    println(io, "Initial Time: ", m.timespan[1])
    println(io, "Final Time  : ", m.timespan[2])
    println(io, "Initial Data are\n    Rotation: ", m.q_0, "\n    Angular Momentum: ", m.L_cm)
    print(io, "Number of points: ", length(m.bodyframe(m.timespan[1])))
end
