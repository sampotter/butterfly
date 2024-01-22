using FileIO, LinearAlgebra, SparseArrays, Printf, DelimitedFiles, Meshes, GeometryBasics, MeshIO, Plots, LaTeXStrings
import GLMakie as Mke

include("utils.jl")

## load mesh
println("Loading mesh...")

mesh_file = ARGS[1]
meshname  = split(split(mesh_file, "/")[end], ".")[1]
_, mesh, vertices, n = load_mesh(mesh_file)

## read other parameters

k  = parse(Float64, ARGS[2])
nu = parse(Float64, ARGS[3])

ps   = Int64.(readdlm("ps.txt")[1:end-1])
tols = Float64.(readdlm("tols.txt")[1:end-1])

num_ps   = length(ps)
num_tols = length(tols)-1 # last is reference

cheb_matvecs = Vector{Matrix{Float64}}(undef, num_ps)
lbo_matvecs  = Vector{Matrix{Float64}}(undef, num_tols)

## read covariance matvecs

for (i, p) in enumerate(ps)
    filestring = @sprintf("cheb_p%i_kappa%.1e_nu%.1e", p, k, nu)
    cheb_matvecs[i] = reshape(reinterpret(Float64, read("matvecs_$filestring.bin")), :, n)
end

for (i, tol) in enumerate(tols)
    filestring = @sprintf("lbo_tol%.0e_kappa%.1e_nu%.1e", tol, k, nu)
    if i==length(tols)
        global ref_matvecs = reshape(reinterpret(Float64, read("matvecs_$filestring.bin")), :, n)
    else
        lbo_matvecs[i] = reshape(reinterpret(Float64, read("matvecs_$filestring.bin")), :, n)
    end
end

## read timings

perf = readdlm(@sprintf("performance_kappa%.1e_nu%.1e.txt", k, nu))

lbo_timings  = Float64.(perf[2:2+num_tols-1, 3])
cheb_timings = Float64.(perf[num_tols+4:end, 2])

## estimate errors using Hutchinson

function hutchinson_frobenius(matvecs)
    m  = size(matvecs, 2) 
    zs = norm.(eachcol(matvecs)).^2
    mean_z = sum(zs) / m
    var_z  = sum((zs .- mean_z).^2) / (m-1)

    tr     = mean_z
    tr_var = var_z / m

    return tr, tr_var
end

nmv = size(lbo_matvecs[1], 2)

cheb_frobsq = Vector{Tuple{Float64, Float64}}(undef, num_ps)
lbo_frobsq  = Vector{Tuple{Float64, Float64}}(undef, num_tols)

ref_frobsq = hutchinson_frobenius(ref_matvecs)

for i in eachindex(cheb_matvecs)
    cheb_frobsq[i] = hutchinson_frobenius(cheb_matvecs[i] - ref_matvecs[:, 1:nmv])
end

for i in eachindex(lbo_matvecs)
    lbo_frobsq[i]  = hutchinson_frobenius(lbo_matvecs[i] -  ref_matvecs[:, 1:nmv])
end

cheb_opnorms = sqrt.(getindex.(cheb_frobsq, 1) / ref_frobsq[1])
lbo_opnorms  = sqrt.(getindex.(lbo_frobsq, 1) / ref_frobsq[1])

@show ref_frobsq
@show cheb_opnorms
@show lbo_opnorms

##

pl = plot(
    title=@sprintf(
        "%s, verts=%i\nκ=%1.1e, ν=%1.1e", 
        split(mesh_file, "/")[end], n, k, nu
        ),
    yscale=:log10,
    xlabel="time per sample (s)",
    ylabel=L"||C - \tilde{C}||"
    )
# plot!(pl,
#     cheb_timings, cheb_opnorms, 
#     label="Chebyshev", line=(:red, 2, :dash), 
#     marker=4, markerstrokewidth=0, markercolor=:red
#     )
plot!(pl,
    lbo_timings, lbo_opnorms, 
    label="Butterfly", line=(:blue, 2, :dash), 
    marker=4, markerstrokewidth=0, markercolor=:blue
    )
pl
# savefig(pl, @sprintf("output/%s_est_error_kappa%.1e_nu%.1e.png", meshname, k, nu))

# npl     = 5
# rs      = round.(Int64, range(2, stop=length(Lam)-1, length=100))
# errs    = zeros(2, length(rs))
# normC   = norm(g.(Lam).^2)
# pl_lam = plot(
#     1:length(Lam), Lam, label="true", line=(:black, 2),
#     title=@sprintf(
#         "Eigenvalues\n%s, verts=%i\nκ=%1.1e, ν=%1.1e", 
#         split(mesh_file, "/")[end], n, k, nu
#         ),
#     xlabel=L"\ell"
#     )
# pl_g = plot(
#     1:length(Lam), g.(Lam), label="true", line=(:black, 2),
#     yscale=:log10,
#     legend=:bottomleft,
#     ylims=[max(minimum(g.(Lam)), 1e-20), 1.0],
#     title=@sprintf(
#         "Spectrum of covariance\n%s, verts=%i\nκ=%1.1e, ν=%1.1e", 
#         split(mesh_file, "/")[end], n, k, nu
#         ),
#     xlabel=L"\ell"
#     )
# for (i, r) in enumerate(rs)
#     @printf("computing covariance error for rank %i (%i of %i)\n", r, i, length(rs))
#     # compute true truncation error
#     errs[1,i] = norm(g.(Lam[1+r:end]).^2) / normC
#     # estimate truncation error using least squares
#     nf   = min(r, 100)
#     xhat = sum((r-nf+1:r)) / nf
#     yhat = sum(Lam[r-nf+1:r]) / nf
#     m = dot((r-nf+1:r) .- xhat, Lam[r-nf+1:r] .- yhat) / norm((r-nf+1:r) .- xhat)^2
#     b = yhat - m*xhat
#     ls = collect(m*(1:length(Lam)) .+ b)
#     ls[1:r] .= Lam[1:r] # don't use linear fit for already computed eigenvalues
#     errs[2,i] = norm(g.(ls[1+r:end]).^2) / norm(g.(ls).^2)
#     if i % div(length(rs), npl) == 0
#         scatter!(pl_lam, r-nf+1:r, Lam[r-nf+1:r], c=i, label="", markerstrokewidth=0, markersize=2)
#         plot!(pl_lam, 
#             r-nf+1:length(Lam), ls[r-nf+1:end], label="rank $r",
#             line=(palette(:default)[div(i-1,div(length(rs), npl))+1], 1, :dash)
#             )
#         plot!(pl_g,   
#             r-nf+1:length(Lam), g.(ls[r-nf+1:end]), label="rank $r",
#             line=(palette(:default)[div(i-1,div(length(rs), npl))+1], 1, :dash)
#             )
#     end
# end
# savefig(pl_lam, @sprintf("output/%s_est_eigenvalues.png", meshname))
# savefig(pl_g,   @sprintf("output/%s_est_covspectrum_kappa%.1e_nu%.1e.png", meshname, k, nu))

# ##

# errs[errs .< 1e-20] .= NaN
# pl = plot(
#     rs, 
#     errs', 
#     labels=reshape([
#         L"\parallel C - C_{\ell} \parallel", 
#         "error estimator"],:,2), 
#     legend=:bottomleft,
#     linestyle=[:solid :dash], linewidth=2,
#     yscale=:log10,
#     title=@sprintf(
#         "Relative truncation error\n%s, verts=%i\nκ=%1.1e, ν=%1.1e", 
#         split(mesh_file, "/")[end], n, k, nu
#         ),
#     xlabel=L"number of computed eigenmodes $\ell$",
#     dpi=200
#     )

# savefig(pl, @sprintf("output/%s_est_truncation_kappa%.1e_nu%.1e.png", meshname, k, nu))