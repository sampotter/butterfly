
using Gmsh, FileIO

## Generate mesh on the unit sphere

gmsh.initialize()

# add the unit sphere to the Julia model
gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)

# synchronize Julia model with the Gmsh environment
gmsh.model.occ.synchronize()

# set length factor to change mesh coarseness
L = parse(Float64, ARGS[1])
gmsh.option.set_number("Mesh.CharacteristicLengthFactor", L)

# generate the mesh
gmsh.model.mesh.generate(2)

# save mesh
i = ARGS[2]
filename = "./sphere$i"
gmsh.write(filename * ".msh")

# visualize
# gmsh.fltk.run()

gmsh.finalize()

## Convert to obj file

save(filename * ".obj", load(filename * ".msh"))

rm(filename * ".msh")
