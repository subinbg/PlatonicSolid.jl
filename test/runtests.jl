push!(LOAD_PATH, "/Users/subeen/Dropbox (Personal)/spack2")
using PlatonicSolid
using Test
using MAT

@testset "PlatonicSolid.jl" begin
    for poly in PlatonicSolid._available_prefixes
        voxel = generate(poly, Float64, 64, 0.5)

        matwrite(string(poly)*".mat", Dict("data"=>voxel))
    end
end
