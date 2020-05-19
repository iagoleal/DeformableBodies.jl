using DeformableBodies.Quaternions
using Test

@testset "Quaternions" begin
    @testset "Instantiation" begin
        @test quaternion(Int) === Quaternion{Int}
        @test quaternion(Float64) === Quaternion{Float64}
        @test quaternion(Rational) === Quaternion{Rational}
        @test quaternion(Quaternion{Int}) === Quaternion{Int}
        @test quaternion(Quaternion{Float64}) === Quaternion{Float64}

        @test_throws MethodError quaternion(Complex{Int})
    end
    @testset "Hamilton Equations" begin
        i = Quaternion(0,1,0,0)
        j = Quaternion(0,0,1,0)
        k = Quaternion(0,0,0,1)

        @test i^2 == -1
        @test j^2 == -1
        @test k^2 == -1
        @test i*j*k == -1
        @test i*j == k
        @test k*i == j
        @test j*k == i
    end
end
