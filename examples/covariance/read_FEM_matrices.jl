using Gmsh, FileIO, LinearAlgebra, SparseArrays, Printf, DelimitedFiles, Meshes, GeometryBasics, MeshIO, Plots, LaTeXStrings
import GLMakie as Mke

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

## load mesh

mesh_file = "../../../../butterfly-LBO-models/sphere0.obj"

gbmesh   = load(mesh_file)
vertices = MeshIO.decompose(MeshIO.Point3{Float64}, gbmesh)
n = length(vertices)

mesh = geometrybasics_to_meshes(gbmesh)

## load matrices

M = load_sparse_matrix(
    reinterpret(Float64, read("M_data.bin")), 
    reinterpret(Int64, read("M_colind.bin")), 
    reinterpret(Int64, read("M_rowptr.bin"))
    )
L = load_sparse_matrix(
    reinterpret(Float64, read("L_data.bin")), 
    reinterpret(Int64, read("L_colind.bin")), 
    reinterpret(Int64, read("L_rowptr.bin"))
    )

## compute eigendecomposition

F = eigen(M \ Matrix(L))
Phi = F.vectors
Lam = F.values

# Phi = reshape(reinterpret(Float64, read("Phi.bin")), :, n)'#[:, 2:end]'
# Lam = reinterpret(Float64, read("Lam.bin"))#[2:end]

## define covariance

k  = 1e-5
nu = Inf

if isinf(nu)
    # squared exponential spectral density function
    g(l) = exp(-k*l^2)
else
    # Matern spectral density function, normalized so g(0) = 1
    g(l) = abs(k^(-2) + l)^(-nu/4 - 0.5) * k^(nu/2 + 1)
end

# Plot everything

fig = Mke.Figure(resolution=(1000,1000))
scene = Mke.LScene(fig[1, 1], scenekw=(center=false,))

#

function plot_nodal_values(scene, vals)
    viz!(
        scene, 
        mesh, 
        color=vals,
        showfacets=false
    )
    cam = Mke.Camera3D(scene.scene)
    cam.lookat[] = [-0.5, 1, 0]
    cam.eyeposition[] = [0, 2, 5] * 1.2
    cam.upvector[]    = [0, 1, 0]
    Mke.update_cam!(scene.scene, cam)
    return scene
end

scene = plot_nodal_values(scene, Phi*Diagonal(g.(Lam))*randn(length(Lam)))
save(@sprintf("/Users/beckman/Downloads/sample_kappa%1.0e_nu%1.0e.png", k, nu), scene.scene)

scene = plot_nodal_values(scene, Phi*Diagonal(g.(Lam).^2)*Phi[1029,:])
save(@sprintf("/Users/beckman/Downloads/covariance_kappa%1.0e_nu%1.0e.png", k, nu), scene.scene)

#

sdf = g.(Lam).^2
sdf[sdf .< 1e-20] .= NaN

pl = plot(
    sdf, 
    label="", 
    legend=:right,
    ylims=[1e-16, 1],
    yscale=:log10,
    title=@sprintf(
        "Spectral density\n%s, verts=%i\nκ=%1.1e, ν=%1.1e", 
        split(mesh_file, "/")[end], n, k, nu
        ),
    xlabel=L"\ell"
    )

savefig(pl, @sprintf("/Users/beckman/Downloads/sdf_kappa%1.0e_nu%1.0e.png", k, nu))

#

normM(v, M) = sqrt(dot(v, M*v))

C = Phi * Diagonal(g.(Lam).^2) * Phi'

# rs    = round.(Int64, 10 .^ range(0, stop=log10(n), length=15)[1:end-1])
rs      = round.(Int64, range(2, stop=n, length=20)[1:end-1])
errs    = zeros(2, length(rs))
normC   = norm(g.(Lam).^2)
opnormC = opnorm(C)
L2C     = maximum([sqrt(dot(C[:,j], M*C[:,j])) for j=1:n])
for (i, r) in enumerate(rs)
    @printf("computing covariance error for rank %i (%i of %i)\n", r, i, length(rs))
    errs[1,i] = norm(g.(Lam[1+r:end]).^2) / normC
    # errs[1,i] = opnorm(C - Phi[:,1:r] * Diagonal(g.(Lam[1:r]).^2) * Phi[:,1:r]') / opnormC
    # errs[2,i] = maximum([
    #     normM(
    #         (r > n/2) ? 
    #         Phi[:,r+1:end] * Diagonal(g.(Lam[r+1:end]).^2) * Phi[j,r+1:end] :
    #         C[:,j] - Phi[:,1:r] * Diagonal(g.(Lam[1:r]).^2) * Phi[j,1:r],
    #         M
    #         ) for j=1:round(Int64, n/10):n
    #     ]) / L2C
    nf = min(r, 50)
    mb = dot((r-nf+1:r), Lam[r-nf+1:r]) / norm((r-nf+1:r))^2
    ls = mb[1]*(1:length(Lam)) # .+ mb[2]
    errs[2,i] = norm(g.(ls[1+r:end]).^2) / norm(g.(ls).^2)
end

#

errs[errs .< 1e-20] .= NaN
pl = plot(
    rs, 
    errs', 
    labels=reshape([
        L"\parallel C - C_{\ell} \parallel", 
        # L"\max_i \parallel Ce_i - C_{\ell}e_i \parallel_{L^2(T)}", 
        "error estimator"],:,2), 
    legend=:right,
    marker=3, markerstrokewidth=0,
    linestyle=[:solid :dash], linewidth=2,
    yscale=:log10,
    title=@sprintf(
        "Relative truncation error\n%s, verts=%i\nκ=%1.1e, ν=%1.1e", 
        split(mesh_file, "/")[end], n, k, nu
        ),
    xlabel=L"number of computed eigenmodes $\ell$",
    dpi=200
    )

savefig(pl, @sprintf("/Users/beckman/Downloads/est_truncation_kappa%1.0e_nu%1.0e.png", k, nu))