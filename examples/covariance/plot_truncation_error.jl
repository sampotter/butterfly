using Pkg
Pkg.instantiate()

using FileIO, LinearAlgebra, SparseArrays, Printf, DelimitedFiles, Meshes, GeometryBasics, MeshIO, Plots, LaTeXStrings
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
println("Loading mesh...")

mesh_file = ARGS[1]
meshname  = split(split(mesh_file, "/")[end], ".")[1]

gbmesh   = load(mesh_file)
vertices = MeshIO.decompose(MeshIO.Point3{Float64}, gbmesh)
n = length(vertices)

mesh = geometrybasics_to_meshes(gbmesh)

## load matrices
println("Loading FEM matrices...")

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

## load or compute eigendecomposition
println("Loading eigendecomposition...")

Phi = reshape(reinterpret(Float64, read("Phi.bin")), n, :)
Lam = reinterpret(Float64, read("Lam.bin"))

# F = eigen(M \ Matrix(L))
# Phi = F.vectors
# Lam = F.values

## define covariance

k  = parse(Float64, ARGS[2])
nu = parse(Float64, ARGS[3])

if iszero(nu)
    # squared exponential spectral density function
    g(l) = exp(-k*l^2)
else
    # Matern spectral density function, normalized so g(0) = 1
    g(l) = abs(1 + k^2*l)^(-nu/4 - 1/2)
end

## plot everything
println("Plotting sample and covariance...")

fig = Mke.Figure()
scene = Mke.LScene(fig[1, 1], scenekw=(center=false,))
scene.show_axis = false

function plot_nodal_values!(scene, vals)
    viz!(
        scene, 
        mesh, 
        color=vals,
        showfacets=false
    )
    cam = Mke.Camera3D(scene.scene)
    if occursin("armadillo", meshname)
        cam.lookat[] = [0, 0, 0]
        cam.eyeposition[] = [0, 0, -100]
        cam.upvector[]    = [0, 1, 0]
        Mke.update_cam!(scene.scene, cam)
    elseif occursin("bunny", meshname)
        cam.lookat[] = [-0.5, 1, 0]
        cam.eyeposition[] = [-1, 2, 5]
        cam.upvector[]    = [0, 1, 0]
        Mke.update_cam!(scene.scene, cam)
    end
end

plot_nodal_values!(scene, Phi*Diagonal(g.(Lam))*randn(length(Lam)))
save(@sprintf("output/%s_sample_kappa%1.0e_nu%1.0e.png", meshname, k, nu), scene.scene)

plot_nodal_values!(scene, Phi*Diagonal(g.(Lam).^2)*Phi[450,:])
save(@sprintf("output/%s_covariance_kappa%1.0e_nu%1.0e.png", meshname, k, nu), scene.scene)

##

npl     = 5
rs      = round.(Int64, range(2, stop=length(Lam)-1, length=100))
errs    = zeros(2, length(rs))
normC   = norm(g.(Lam).^2)
pl_lam = plot(
    1:length(Lam), Lam, label="true", line=(:black, 2),
    title=@sprintf(
        "Eigenvalues\n%s, verts=%i\nκ=%1.1e, ν=%1.1e", 
        split(mesh_file, "/")[end], n, k, nu
        ),
    xlabel=L"\ell"
    )
pl_g = plot(
    1:length(Lam), g.(Lam), label="true", line=(:black, 2),
    yscale=:log10,
    legend=:bottomleft,
    ylims=[max(minimum(g.(Lam)), 1e-20), 1.0],
    title=@sprintf(
        "Spectrum of covariance\n%s, verts=%i\nκ=%1.1e, ν=%1.1e", 
        split(mesh_file, "/")[end], n, k, nu
        ),
    xlabel=L"\ell"
    )
for (i, r) in enumerate(rs)
    @printf("computing covariance error for rank %i (%i of %i)\n", r, i, length(rs))
    # compute true truncation error
    errs[1,i] = norm(g.(Lam[1+r:end]).^2) / normC
    # estimate truncation error using least squares
    nf   = min(r, 100)
    xhat = sum((r-nf+1:r)) / nf
    yhat = sum(Lam[r-nf+1:r]) / nf
    m = dot((r-nf+1:r) .- xhat, Lam[r-nf+1:r] .- yhat) / norm((r-nf+1:r) .- xhat)^2
    b = yhat - m*xhat
    ls = collect(m*(1:length(Lam)) .+ b)
    ls[1:r] .= Lam[1:r] # don't use linear fit for already computed eigenvalues
    errs[2,i] = norm(g.(ls[1+r:end]).^2) / norm(g.(ls).^2)
    if i % div(length(rs), npl) == 0
        scatter!(pl_lam, r-nf+1:r, Lam[r-nf+1:r], c=i, label="", markerstrokewidth=0, markersize=2)
        plot!(pl_lam, 
            r-nf+1:length(Lam), ls[r-nf+1:end], label="rank $r",
            line=(palette(:default)[div(i-1,div(length(rs), npl))+1], 1, :dash)
            )
        plot!(pl_g,   
            r-nf+1:length(Lam), g.(ls[r-nf+1:end]), label="rank $r",
            line=(palette(:default)[div(i-1,div(length(rs), npl))+1], 1, :dash)
            )
    end
end
savefig(pl_lam, @sprintf("output/%s_est_eigenvalues.png", meshname))
savefig(pl_g,   @sprintf("output/%s_est_covspectrum_kappa%1.0e_nu%1.0e.png", meshname, k, nu))

##

errs[errs .< 1e-20] .= NaN
pl = plot(
    rs, 
    errs', 
    labels=reshape([
        L"\parallel C - C_{\ell} \parallel", 
        "error estimator"],:,2), 
    legend=:bottomleft,
    linestyle=[:solid :dash], linewidth=2,
    yscale=:log10,
    title=@sprintf(
        "Relative truncation error\n%s, verts=%i\nκ=%1.1e, ν=%1.1e", 
        split(mesh_file, "/")[end], n, k, nu
        ),
    xlabel=L"number of computed eigenmodes $\ell$",
    dpi=200
    )

savefig(pl, @sprintf("output/%s_est_truncation_kappa%1.0e_nu%1.0e.png", meshname, k, nu))