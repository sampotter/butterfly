using Gmsh, FileIO, Plots, LaTeXStrings, Printf, DelimitedFiles
import MeshIO: decompose, Point3

include("sphere_cov.jl")

mesh_file = ARGS[1]
k  = parse(Float64, ARGS[2])
nu = parse(Float64, ARGS[3])
P_file   = ARGS[4]
tol_file = ARGS[5]

Ps   = Int64.(readdlm(P_file)[1:end-1])
tols = readdlm(tol_file)[1:end-1]

mesh = load(mesh_file)

# from the mesh, compute geodesic distances at which to evaluate the covariance 
vertices = decompose(Point3{Float64}, mesh)
z        = getindex.(vertices, 3)
ts       = acos.(z)
perm     = sortperm(ts)

# Matern power spectral density function
g(l) = abs(k^2 + l)^(-nu/4 - 0.5)

# compute true covariance
c_true = sphere_cov.(g, ts, rtol=1e-15, atol=0, maxterms=500000, blksz=1000, verb=1)

# read chebyshev and butterfly covariances from file
c_chebs = zeros(length(ts), length(Ps))
for (i, P) in enumerate(Ps)
    cov_file = "output/c_cheb_p$P.bin"
    c_chebs[:, i] .= reinterpret(Float64, read(cov_file))
end
c_lbos = zeros(length(ts), length(tols))
freqs_vec = []
for (i, tol) in enumerate(tols)
    cov_file = @sprintf("output/c_lbo_tol%.0e.bin", tol)
    c_lbos[:, i] .= reinterpret(Float64, read(cov_file))
    freqs_file = @sprintf("output/freqs_tol%.0e.bin", tol)
    push!(freqs_vec, reinterpret(Float64, read(freqs_file)))
end
freqs = hcat(freqs_vec...)

L = size(freqs, 1)
true_eigs = zeros(2*L)
i = 1; l = 1
while i < L
    true_eigs[i:i+2l+1] .= l*(l+1)
    global i += 2l + 1
    global l += 1
end
true_eigs = true_eigs[1:L]

## Plot freqs

pl = plot(
    freqs[2:end,:].^2,
    title="$(split(mesh_file, '/')[end]), " * L"\kappa" * " = $k, " * L"\nu" * " = $nu",
    xlabel=L"\ell",
    yscale=:log10, dpi=300, legend=:bottomright, 
    ylabel=L"\tilde{\omega}_\ell^2", 
    labels=reshape([@sprintf("tol = %.0e", tol) for tol in tols], 1, length(tols)),
    # colormap=:rainbow
    )
display(pl)

##

pl = plot(
    abs.((freqs[2:end,:].^2 .- true_eigs[1:end-1]) ./ true_eigs[1:end-1]),
    title="$(split(mesh_file, '/')[end]), " * L"\kappa" * " = $k, " * L"\nu" * " = $nu",
    xlabel=L"\ell",
    yscale=:log10, dpi=300, legend=:bottomright, 
    ylabel=L"|\tilde{\omega}_\ell^2 - \lambda_\ell| / \lambda_\ell", 
    labels=reshape([@sprintf("tol = %.0e", tol) for tol in tols], 1, length(tols)),
    # colormap=:rainbow
    )
display(pl)
# savefig(pl, "output/eigenvalue_comparison.png")

## Plot covariances

pl = plot(
    ts[perm], c_true[perm], 
    c=:black, 
    la=0.0, marker=(:circle,1),
    label="true", title="$(split(mesh_file, '/')[end]), " * L"\kappa" * " = $k, " * L"\nu" * " = $nu",
    xlabel=L"\theta", ylabel=L"c(\theta)", 
    yscale=:log10, 
    ylims=(1e-16, 1e0),
    dpi=300
    )

plot!(
    pl, ts[perm], c_chebs[perm, :], 
    linestyle=:dash, 
    la=0.0, marker=(:circle,1), markerstrokewidth=0,
    labels=reshape(["P = $P" for P in Ps], 1, length(Ps)),
    palette=cgrad(:darkrainbow, length(Ps), categorical = true)[1:end]
    )

plot!(
    pl, ts[perm], c_lbos[perm, :], 
    la=0.0, marker=(:circle,1), markerstrokewidth=0,
    labels=reshape([@sprintf("tol = %.0e", tol) for tol in tols], 1, length(tols)),
    palette=cgrad(:lightrainbow, length(tols), categorical = true)[1:end]
    )
# display(pl)
savefig(pl, "output/cov_comparison.png")

## Plot errors

pl = plot(
    title="$(split(mesh_file, '/')[end]), " * L"\kappa" * " = $k, " * L"\nu" * " = $nu",
    xlabel=L"\theta",
    yscale=:log10, dpi=300, legend=:bottomright, 
    ylabel=L"|c_{true}(\theta) - c_{cheb}(\theta)|", 
    colormap=:inferno
    )

plot!(
    pl, ts[perm], 
    abs.(c_true[perm] .- c_chebs[perm,:]), 
    la=0.0, marker=(:circle,1), markerstrokewidth=0,
    linestyle=:dash, labels=reshape(["P = $P" for P in Ps], 1, length(Ps)),
    palette=cgrad(:darkrainbow, length(Ps), categorical = true)[1:end]
    )

plot!(
    pl, ts[perm], 
    abs.(c_true[perm] .- c_lbos[perm,:]),
    la=0.0, marker=(:circle,1), markerstrokewidth=0,
    labels=reshape([@sprintf("tol = %.0e", tol) for tol in tols], 1, length(tols)),
    palette=cgrad(:lightrainbow, length(tols), categorical = true)[1:end]
    )
# display(pl)
savefig(pl, "output/absolute_error_comparison.png")

pl = plot(
    title="$(split(mesh_file, '/')[end]), " * L"\kappa" * " = $k, " * L"\nu" * " = $nu",
    xlabel=L"\theta",
    yscale=:log10, dpi=300, legend=:bottomright, 
    ylabel=L"|c_{true}(\theta) - c_{cheb}(\theta)| \ / \ |c_{true}(\theta)|", 
    colormap=:inferno
    )

plot!(
    pl, ts[perm], 
    abs.((c_true[perm] .- c_chebs[perm,:]) ./ c_true[perm]), 
    la=0.0, marker=(:circle,1), markerstrokewidth=0,
    linestyle=:dash, labels=reshape(["P = $P" for P in Ps], 1, length(Ps)),
    palette=cgrad(:darkrainbow, length(Ps), categorical = true)[1:end]
    )

plot!(
    pl, ts[perm], 
    abs.((c_true[perm] .- c_lbos[perm,:]) ./ c_true[perm]),
    la=0.0, marker=(:circle,1), markerstrokewidth=0,
    labels=reshape([@sprintf("tol = %.0e", tol) for tol in tols], 1, length(tols)),
    palette=cgrad(:lightrainbow, length(tols), categorical = true)[1:end]
    )
# display(pl)
savefig(pl, "output/relative_error_comparison.png")
