
using FileIO, LinearAlgebra, SparseArrays, Printf, DelimitedFiles, Meshes, GeometryBasics, MeshIO, Plots, LaTeXStrings
import GLMakie as Mke

function matern_sdf(k, nu, l)
    if iszero(nu)
        # squared exponential spectral density function
        return exp(-k*l^2)
    else
        # Matern spectral density function, normalized so g(0) = 1
        return abs(1 + k^2*l)^(-nu/4 - 1/2)
    end
end

function geometrybasics_to_meshes(mesh)
    vertices = Meshes.Point3.(
        getfield.(MeshIO.decompose(MeshIO.Point3{Float64}, mesh), :data)
        )
    connections = Meshes.connect.([
        getfield(GeometryBasics.value.(face), :data) for face in GeometryBasics.faces(mesh)
        ])
    return Meshes.SimpleMesh(vertices, connections)
end

function load_sparse_matrix(data, colind, rowptr)
    colind .+= 1
    rowind   = vcat(
        [fill(i, rowptr[i+1] - rowptr[i]) for i=1:length(rowptr)-1]...
        )
    return sparse(rowind, colind, data)
end

function load_mesh(mesh_file)
    gbmesh   = load(mesh_file)
    vertices = MeshIO.decompose(MeshIO.Point3{Float64}, gbmesh)
    mesh     = geometrybasics_to_meshes(gbmesh)

    return gbmesh, mesh, vertices, length(vertices)
end

function plot_nodal_values(gbmesh, vals; meshname="bunny0", cmap=:turbo, wframe=false)
    fig = Mke.Figure(size=(1000,1000))
    lscene = Mke.LScene(fig[1, 1])
    lscene.show_axis = false

    Mke.mesh!(
        lscene,
        gbmesh, 
        color=vals,
        colormap=cmap,
        showfacets=false
    )
    if wframe
        Mke.wireframe!(
            lscene,
            gbmesh,
            linewidth=0.5,
            alpha=0.5,
            color=:black
        )
    end

    cam = Mke.Camera3D(lscene.scene, center=false)
    if occursin("armadillo", meshname)
        cam.lookat[] = [0, 0, 0]
        cam.eyeposition[] = [0, 0, -100]
        cam.upvector[]    = [0, 1, 0]
        Mke.update_cam!(lscene.scene, cam)
    elseif occursin("bunny", meshname)
        cam.lookat[] = [-0.5, 1, 0]
        cam.eyeposition[] = [-1, 2, 5]
        cam.upvector[]    = [0, 1, 0]
        Mke.update_cam!(lscene.scene, cam)
    end

    return fig
end