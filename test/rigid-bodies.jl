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
end

