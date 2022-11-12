module PlatonicSolid

const _available_face_numbers = [4,6,8,12,20]
const _available_prefixes = [:tetra, :cubic, :octa, :dodeca, :icosa]

export generate!

include("utils.jl")
include("structures.jl")


end
