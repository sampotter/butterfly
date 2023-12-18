using Gmsh, FileIO, Plots, Plots.PlotMeasures
import MeshIO: decompose, Point3

mesh_file = ARGS[1]
Ps        = ARGS[2:end]

mesh = load(mesh_file)

vertices = decompose(Point3{Float64}, mesh)
x = getindex.(vertices, 1)
y = getindex.(vertices, 2)
z = getindex.(vertices, 3)
theta = acos.(z)
rho   = replace(
    v -> isnan(v) ? 0 : v, 
    sign.(y) .* acos.(x ./ sqrt.(x.^2 .+ y.^2))
    )

for P in Ps
    cov_file    = "output/c_cheb_p$P.bin"
    sample_file = "output/z_cheb_p$P.bin"
    fig_file    = "output/fig_p$P.png"

    c_cheb = reinterpret(Float64, read(cov_file))
    z_cheb = reinterpret(Float64, read(sample_file))

    gr(size=(1000,700))
    p1 = scatter(
        theta, rho, marker_z=c_cheb,
        label="", xlabel="theta", ylabel="rho", title="Covariance",
        markerstrokewidth=0, markersize=2,
        c=:reds
        )
    p2 = scatter(
        theta, rho, marker_z=z_cheb,
        label="", xlabel="theta", ylabel="rho", title="Sample",
        markerstrokewidth=0, markersize=2,
        c=:lightrainbow
        )
    p = plot(p1, p2, layout=@layout [a b])
    savefig(p, fig_file)
end

# plot(
#     scatter(
#     x=getindex.(vertices,1),
#     y=getindex.(vertices,2),
#     z=getindex.(vertices,3),
#     mode="markers",
#     marker=attr(
#         size=4,
#         color=z_cheb,
#         colorscale="Viridis",
#         opacity=0.8
#     ),
#     type="scatter3d"
#     )
# )

##

