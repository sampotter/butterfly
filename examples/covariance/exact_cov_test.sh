#!/bin/zsh
set -o nullglob

MESH=$1
MESHBASE=${MESH##*/}
MESHNAME=${MESHBASE%.*}

# KAPPA, NU pairs to test
PARAMS=(1e-6,0.0 1e-1,4.0 3e-1,0.5)

# set up julia
julia --project=. -e "using Pkg; Pkg.instantiate()"

# clean this directory and make output directory
rm -f -- *.bin
mkdir -p output

# compute and store L, M, and Phi if they haven't already been computed
if test -f output/"${MESHNAME}_L_rowptr.bin"; 
then
    for filename in "L_rowptr" "L_colind" "L_data" "M_rowptr" "M_colind" "M_data" "Phi" "Lam"
    do
        mv output/"${MESHNAME}_$filename.bin" "$filename.bin"
    done
elif ! test -f "Phi.bin";
then
    ./save_FEM_matrices $MESH
fi

# estimate and plot error
for params in $PARAMS
do 
    IFS="," read KAPPA NU <<< "$params"
    julia --project=. plot_sample.jl         $MESH $KAPPA $NU exact
    julia --project=. plot_exact_cov_test.jl $MESH $KAPPA $NU
done

# move generated files to output folder
for file in *.bin
do
    mv $file output/"${MESHNAME}_$file"
done