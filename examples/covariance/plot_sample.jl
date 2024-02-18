using FileIO, LinearAlgebra, SparseArrays, Printf, DelimitedFiles, Meshes, GeometryBasics, MeshIO, Plots, LaTeXStrings
import GLMakie as Mke

include("utils.jl")

## load mesh
println("Loading mesh...")

mesh_file = ARGS[1]
meshname  = split(split(mesh_file, "/")[end], ".")[1]
gbmesh, _, _, n = load_mesh(mesh_file)

## read other parameters

k    = parse(Float64, ARGS[2])
nu   = parse(Float64, ARGS[3])
g(l) = matern_sdf(k, nu, l)

method = ARGS[4]

if method == "lbo"
    tol = parse(Float64, ARGS[5])
    filestring = @sprintf("lbo_tol%.0e_kappa%.1e_nu%.1e", tol, k, nu)
elseif method == "cheb"
    p   = parse(Int64, ARGS[5]) 
    filestring = @sprintf("cheb_p%i_kappa%.1e_nu%.1e", p, k, nu)
elseif method == "exact"
    filestring = @sprintf("exact_kappa%.1e_nu%.1e", k, nu)

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

    ## load eigendecomposition
    println("Loading eigendecomposition...")

    Phi = reshape(reinterpret(Float64, read("Phi.bin")), n, :)
    Lam = reinterpret(Float64, read("Lam.bin"))
end

## plot sample and marginal covariance
println("Plotting sample and marginal covariance...")

if method == "exact"
    sample   = Phi*Diagonal(g.(Lam))*randn(length(Lam))
    marg_cov = Phi*Diagonal(g.(Lam).^2)*Phi[1,:]
    wframe   = true
else
    sample   = reinterpret(Float64, read("z_$filestring.bin"))
    marg_cov = reinterpret(Float64, read("c_$filestring.bin"))
    wframe   = false
end

scene = plot_nodal_values(
    gbmesh, sample, 
    meshname=meshname, cmap=:turbo, wframe=wframe
    )
save(@sprintf("output/%s_sample_%s.png", meshname, filestring), scene)

scene = plot_nodal_values(
    gbmesh, marg_cov, 
    meshname=meshname, cmap=:inferno, wframe=wframe
    )
save(@sprintf("output/%s_marginal_%s.png", meshname, filestring), scene)