using PlatonicSolid
using Test

const radius = 0.5
const N = 128

const fraction = Dict(
    :tetra  =>   sqrt(8)/3 * ( radius /      sqrt(3/2)                          )^3,
    :cubic  =>           8 * ( radius /        sqrt(3)                          )^3,
    :octa   => sqrt(128)/3 * ( radius /        sqrt(2)                          )^3,
    :dodeca =>   61.304952 * ( radius / (      sqrt(3) * (1+sqrt(5))/2)         )^3, 
    :icosa  =>   17.453560 * ( radius / ((1+sqrt(5))/2 * sqrt(3-(1+sqrt(5))/2)) )^3
);

@testset "PlatonicSolid.jl" begin
    for poly in keys(fraction)
        voxel = Array{Float64}(undef, N, N, N)
        generate!(poly, voxel, radius, axis=[rand() rand() rand();], angle=[rand()*π/2])

        @test isapprox(sum(voxel)/length(voxel), fraction[poly], rtol=3e-2)
    end
end
