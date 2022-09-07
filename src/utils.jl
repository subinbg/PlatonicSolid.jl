function swap(A, i, j)
    x = similar(A)
    copyto!(x, A)
    x[i], x[j] = x[j], x[i]
    return x
end

"""
Counter-clockwise rotation matrix w.r.t. to an axis    
Ref: Eq. (20), http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf
"""
function rotation_matrix(axis::AbstractVector, angle; type::Type{T}=Float64) where T
    @assert length(axis) == 3 "Wrong axis: $axis"

    nv = axis ./ sqrt(sum(axis.^2))
    cos_t = cos(angle)
    sin_t = sin(angle)

    rmtx = Matrix{T}(undef, 3, 3)
    for idx in 1:3
        cps = circshift(nv, 1-idx)
        cms = circshift(swap(nv, 2, 3), idx-1)
        entry = [
            cos_t + nv[idx]^2 * (1 - cos_t),
            cps[1] * cps[2] * (1 - cos_t) + cps[3] * sin_t,
            cms[1] * cms[2] * (1 - cos_t) - cms[3] * sin_t
        ]
        rmtx[:,idx] = circshift(entry, idx-1)
    end

    return rmtx
end

""" Mathematical convention; theta-azimuthal, phi-polar """
function sph_coord!(arr::AbstractArray,r,theta,phi) 
    arr[1] = r*cos(theta)*cos(phi)
    arr[2] = r*sin(theta)*cos(phi)
    arr[3] = r*sin(phi)
end

function transform_face_idx(face_idx::Vector{NTuple{N,Int}}) where N
    n_vertices_per_face = length(face_idx[1])
    n_faces = length(face_idx)

    faces = Matrix{Int}(undef, n_vertices_per_face, n_faces)
    for idx in 1:n_faces
        faces[:,idx] .= face_idx[idx] .+ 1
    end
    return faces
end

function planeDistance(faces::Matrix, vertices::Matrix)
    face = faces[:,1]
    vertices_in_face = vertices[:,face]
    midpoint = sum(vertices_in_face, dims=2) / length(face)
    sqrt(sum(midpoint.^2))
end

function get_shape(poly::Union{Symbol,Int})
    if poly in _available_face_numbers
        idx = findfirst(x->x==poly, _available_face_numbers)
        return _available_prefixes[idx]
    elseif poly in _available_prefixes
        return poly
    else
        ArgumentError("Unknown option: $poly")
    end
end

"""
Line-plane intersection point    
- p_l1, p_l2: points on a line
- p_on_pl: point on a plane
- p_normal: normal to a plane (need not be normalized)
- https://stackoverflow.com/a/18543221
"""
function line_plane_intersection(p_l1::Vector, p_l2::Vector, p_on_pl::Vector, p_normal, epsilon=1e-6)
    u = p_l2 .- p_l1
    inprod = sum(p_normal .* u)
    if abs(inprod) > epsilon
        w = p_l1 .- p_on_pl
        fac = -sum(p_normal .* w) / inprod
        return p_l1 .+ (u .* fac)
    end
    # throw(ArgumentError("Plane and line are parallel"))
    return p_l2
end

function check_inside(point::Vector, vertices::Matrix)
    origin = similar(point)
    fill!(origin, 0)

    p_normal = sum(vertices, dims=2) / size(vertices, 2)
    p_on_pl = vertices[:,1]

    intersection = line_plane_intersection(origin, point, p_on_pl, p_normal)
    sqrt(sum(point.^2)) <= sqrt(sum(intersection.^2))
end