"""
Get coordinates of vertices and face normals for regular polyhedra.      
`radius` is the radius of a smallest sphere that encloses a polyhedron.
- https://en.wikipedia.org/wiki/Regular_dodecahedron
- https://userpages.umbc.edu/~squire/reference/polyhedra.shtml
"""

const _tetra_faces = [
    (0,1,2), (0,2,3), (0,3,1), (1,2,3)
]

"""
Tetrahedron
"""
function tetra(radius)
    phi = -19.471220333 * π / 180
    theta = 0.0

    vertices = Matrix{Float64}(undef, 3, 4)
    vertices[:,1] = [0,0,radius]
    for idx in 2:4
        sph_coord!(view(vertices,:,idx), radius, theta, phi)
        theta += π*2/3
    end
    faces = transform_face_idx(_tetra_faces)
    return vertices, faces
end

const _cubic_faces = [
    (0,1,2,3), (4,5,6,7), (0,1,5,4), (1,2,6,5), (2,3,7,6), (3,0,4,7)
]

"""
Cube
"""
function cubic(radius)
    phi = acos(sqrt(3/4-0.5^2)/sqrt(3)*2)
    theta = 0.0

    vertices = Matrix{Float64}(undef, 3, 8)
    for idx in 1:4
        sph_coord!(view(vertices,:,idx), radius, theta, phi)
        sph_coord!(view(vertices,:,idx+4), radius, theta, -phi)
        theta += π/2
    end
    faces = transform_face_idx(_cubic_faces)
    return vertices, faces
end

const _octa_faces = [
    (0,1,2), (0,2,3), (0,3,4), (0,4,1),
    (5,1,2), (5,2,3), (5,3,4), (5,4,1)
]

"""
Octahedron
"""
function octa(radius)
    theta = 0.0
    vertices = Matrix{Float64}(undef, 3, 6)
    sph_coord!(view(vertices,:,1), radius, 0.0, π/2)
    sph_coord!(view(vertices,:,6), radius, 0.0, -π/2)
    for idx in 2:5
        sph_coord!(view(vertices,:,idx), radius, theta, 0.0)
        theta += π/2
    end
    faces = transform_face_idx(_octa_faces)
    return vertices, faces
end


const _dodeca_faces = [
    (0,1,2,3,4), (0,1,6,10,5), (1,2,7,11,6),
    (2,3,8,12,7), (3,4,9,13,8), (4,0,5,14,9),
    (15,16,11,6,10), (16,17,12,7,11), (17,18,13,8,12),
    (18,19,14,9,13), (19,15,10,5,14), (15,16,17,18,19)
]

"""
Dodecahedron
"""
function dodeca(radius)
    phi_a = 52.62263590 * π / 180 
    phi_b = 10.81231754 * π / 180

    vertices = Matrix{Float64}(undef, 3, 20)
    theta = 0.0
    for idx in 1:5
        sph_coord!(view(vertices,:,idx), radius, theta, phi_a)
        sph_coord!(view(vertices,:,idx+5), radius, theta, phi_b)
        sph_coord!(view(vertices,:,idx+10), radius, theta, -phi_b)
        sph_coord!(view(vertices,:,idx+15), radius, theta, -phi_a)
        theta += 72*π/180
    end
    faces = transform_face_idx(_dodeca_faces)
    return vertices, faces
end


const _icosa_faces = [
    (0,1,2), (0,2,3), (0,3,4), (0,4,5),
    (0,5,1), (11,6,7), (11,7,8), (11,8,9),
    (11,9,10), (11,10,6), (1,2,6), (2,3,7), 
    (3,4,8), (4,5,9), (5,1,10), (6,7,2),
    (7,8,3), (8,9,4), (9,10,5), (10,6,1)
]

"""
Icosahedron
"""
function icosa(radius)
    phi = 26.56505 * π / 180 

    vertices = Matrix{Float64}(undef, 3, 12)
    theta = 0.0
    for idx in 1:6
        sph_coord!(view(vertices,:,idx), radius, theta, phi)
        sph_coord!(view(vertices,:,idx+6), radius, theta+36*π/180, -phi)
        theta += 72*π/180
    end
    faces = transform_face_idx(_icosa_faces)
    return vertices, faces
end

function generate(poly::Union{Symbol,Int}, ::Type{T}, N::Int, radius;
    translation::Vector=[0,0,0], axis::Vector=[1,0,0], angle=0.0) where T

    @assert 0.0 <= radius <= 1.0 "Radius should be in [0,1]"
    radius_N = radius * N
    shape = get_shape(poly)
    @eval func = $shape
    vertices, faces = func(radius_N)
    rmtx = rotation_matrix(axis, angle, type=T)

    generate(T, N, faces, vertices, translation, rmtx)
end

function generate(::Type{T}, N::Int, 
    faces::Matrix, vertices::Matrix, 
    translation::Vector, rmtx::Matrix) where T

    vertices_rot = rmtx*vertices
    generate(T, N, faces, vertices_rot, translation)
end

function generate(::Type{T}, N::Int, 
    faces::Matrix, vertices::Matrix, translation::Vector) where T

    planedist = planeDistance(faces, vertices)
    polyhedron = ones(T, N, N, N)

    @inbounds for x3 in 1:N # @simd 
        x3r = x3 - 1 - N÷2 + translation[3]
        for x2 in 1:N
            x2r = x2 - 1 - N÷2 + translation[2]
            for x1 in 1:N
                x1r = x1 - 1 - N÷2 + translation[1]
                
                for f in 1:size(faces,2)
                    vertices_face = vertices[:,faces[:,f]]
                    if ~check_inside([x1r,x2r,x3r], vertices_face)
                        polyhedron[x1,x2,x3] = zero(T)
                        break
                    end
                end
            end
        end
    end

    return polyhedron
end