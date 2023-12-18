# Generate meshes ./sphere$i.obj where $i is the refinement number
# I have to open these meshes and re-export them in Meshlab for ./bf_lbo to run
# properly... not sure what's up with that

Ls=(0.36 0.25 0.18 0.127 0.0889 0.0628 0.044 0.031  0.022  0.0155 0.011)
Vs=(1059 2111 3998 7918  16117  31875  64626 129216 256603 514632 1018122)

declare -i N=${#Ls[@]}-1

for i in {0..$N}
do
    echo ""
    julia mesh_sphere.jl ${Ls[i+1]} $i
done
