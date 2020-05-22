using DeformableBodies
using Test
using Random

filepath(x) = joinpath(dirname(@__FILE__), x)

# Fix random number seed, for reproducibility
Random.seed!(12342352154)


@testset "DeformableBodies Package" begin
    @info "Testing Quaternions"
    include(filepath("quaternions.jl"))
    @info "Testing dynamics for Rigid Bodies"
    include(filepath("rigid-bodies.jl"))
end
