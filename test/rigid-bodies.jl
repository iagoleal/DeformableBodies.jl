using DeformableBodies
using Test
using Random
import LinearAlgebra: norm

# Uniform Distribution on Sphere
function unif_sphere(m)
    v = randn(m)
    return v / norm(v)
end

approx_zero_array(A) = all(isapprox.(A, 0.0, atol=1e-6))

@testset "Rigid Bodies" for i in 1:5
    @testset "Zero angular momentum" begin
        npts = rand(10:100)
        # A rigid body contained on a ball
        r_0 = [PointMass(rand()*10, unif_sphere(3)) for i in 1:npts]
        bigV_r0 = vcat(pos.(r_0)...)
        model = Model( t -> r_0 # Constant trajectories
                     , (0., 1.)
                     , quaternion(1)
                     , [0.,0.,0.]
                     )
        rotbodies, R, Π = solve!(model)

        for t in range(0., 1., length=30)
            bigV_rt = vcat(pos.(model.inertialframe(t))...)
            @test R(t) ≈ 1 atol=1e-7
            @test approx_zero_array(Π(t))
            @test isapprox(bigV_r0, bigV_rt, atol=1e-6)
        end
    end
    @testset "Gauge Invariance" begin
        npts = rand(10:100)
        # A rigid body contained on a ball
        r_0 = [PointMass(rand()*10, unif_sphere(3)) for i in 1:npts]
        cm = center_of_mass(r_0)

        # Apply random rotations along the trajectory around cm
        freqs = randn(4)
        phases = 2π*rand(4)
        rots(t) = Quaternion(cos.(t*freqs + phases))
        rs = t -> [rotate(x, rots(t), cm) for x in r_0]

        model = Model( rs
                     , (0., 1.)
                     , quaternion(1)
                     , [0.,0.,0.]
                     )
        rotbodies, R, Π = solve!(model)

        bigV_r0 = vcat(pos.(model.inertialframe(0.0))...)
        for t in range(0., 1., length=30)
            bigV_rt = vcat(pos.(model.inertialframe(t))...)
            @test isapprox(Π(t), zeros(3), atol=1e-6)
            @test isapprox(bigV_r0, bigV_rt, atol=1e-6)
        end
    end
end
