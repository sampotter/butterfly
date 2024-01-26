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

function hutchinson_frobenius(matvecs; corr_matvecs=nothing)
    n  = size(matvecs, 2) 
    zs = norm.(eachcol(matvecs)).^2
    mean_z = sum(zs) / n
    var_z  = sum((zs .- mean_z).^2) / (n-1)

    tr     = mean_z
    tr_var = var_z / n

    if !isnothing(corr_matvecs)
        m  = size(corr_matvecs, 2) 
        ws = norm.(eachcol(corr_matvecs)).^2
        mean_w = sum(ws) / m
        var_w  = sum((ws .- mean_w).^2) / (m-1)

        rho = sum((zs .- mean_z).*(ws .- mean_w)) / sqrt(var_z*var_w)

        return tr, tr_var, rho
    else
        return tr, tr_var
    end
end

nmv = size(lbo_matvecs[1], 1)

cheb_frobsq = Vector{Tuple{Float64, Float64, Float64}}(undef, num_ps)
lbo_frobsq  = Vector{Tuple{Float64, Float64, Float64}}(undef, num_tols)

ref_frobsq = hutchinson_frobenius(ref_matvecs)

for i in eachindex(cheb_matvecs)
    cheb_frobsq[i] = hutchinson_frobenius(
        cheb_matvecs[i] - ref_matvecs[:, 1:nmv],
        corr_matvecs=ref_matvecs
        )
end

for i in eachindex(lbo_matvecs)
    lbo_frobsq[i]  = hutchinson_frobenius(
        lbo_matvecs[i] -  ref_matvecs[:, 1:nmv],
        corr_matvecs=ref_matvecs
        )
end

# cheb_opnorms = sqrt.(getindex.(cheb_frobsq, 1) / ref_frobsq[1])
lbo_opnorms = sqrt.(getindex.(lbo_frobsq, 1)  / ref_frobsq[1])

# lbo_logvars = 0.25*(
#     getindex.(lbo_frobsq, 2) ./ getindex.(lbo_frobsq, 1).^2 .+
#     ref_frobsq[2] / ref_frobsq[1]^2
#     )
# lbo_2std = (
#     abs.(lbo_opnorms .* exp.(-2sqrt.(lbo_logvars)) .- lbo_opnorms),
#          lbo_opnorms .* exp.( 2sqrt.(lbo_logvars)) .- lbo_opnorms
# )

lbo_sqvars = (
    getindex.(lbo_frobsq, 2).^2 ./ ref_frobsq[1].^2 .+
    getindex.(lbo_frobsq, 1).^2 .* ref_frobsq[2].^2 ./ ref_frobsq[1].^4 .-
    2*getindex.(lbo_frobsq, 1) .* getindex.(lbo_frobsq, 3) ./ ref_frobsq[1].^3
    ) ./ nmv
lbo_2std = (
    abs.(lbo_opnorms .* exp.(-2sqrt.(lbo_logvars)) .- lbo_opnorms),
         lbo_opnorms .* exp.( 2sqrt.(lbo_logvars)) .- lbo_opnorms
)

@show ref_frobsq
# @show cheb_opnorms
@show lbo_opnorms

# @show cheb_1std
@show lbo_2std

##

pl = plot(
    title=@sprintf(
        "%s, verts=%i\nκ=%1.1e, ν=%1.1e", 
        split(mesh_file, "/")[end], n, k, nu
        ),
    # yscale=:log10,
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
    yerror=lbo_2std,
    label="Butterfly", line=(:blue, 1, :dash, 0.5), 
    marker=4, msw=1, markercolor=:blue
    )
savefig(pl, @sprintf("output/%s_est_error_kappa%.1e_nu%.1e.png", meshname, k, nu))