#!/bin/zsh
set -o nullglob

MESH=$1
MESHBASE=${MESH##*/}
MESHNAME=${MESHBASE%.*}

KAPPA=1e-5
NU=0
NUM_SAMPLES=1000
FRACTION=1.0

REFTOL=1e-8
TOLS=(1e-2 1e-4 1e-6)
LOGPS=(4 8)

# clean this directory and make output directory
rm -f -- *.bin *.csv *.txt
:> tols.csv
:> Ps.csv
mkdir -p output

# compute reference solution using butterfly
cmd="./lbo_cov $MESH $KAPPA $NU $NUM_SAMPLES $REFTOL $FRACTION" > "lbo_ref_tol${REFTOL}_log.txt"
echo "\n$cmd"; eval $cmd

# compute covariances for various Chebyshev orders
for logP in ${LOGPS[@]}
do
    P="$((2**$logP))"
    cmd="./cheb_cov $MESH $P $KAPPA $NU $NUM_SAMPLES" > "cheb_p${P}_log.txt"
    echo "\n$cmd"; eval $cmd
    echo -n "$P " >> Ps.csv
done

# compute covariance for various butterfly tolerances
for TOL in ${TOLS[@]}
do
    cmd="./lbo_cov $MESH $KAPPA $NU $NUM_SAMPLES $TOL $FRACTION" > "lbo_tol${TOL}_log.txt"
    echo "\n$cmd"; eval $cmd
    echo -n "$TOL " >> tols.csv
done

# move generated files to output folder
for file in *.bin *.txt
do
    mv $file output/"${MESHNAME}_$file"
done
mv Ps.csv output
mv tols.csv output