#!/bin/zsh
set -o nullglob

MESH=$1
MESHBASE=${MESH##*/}
MESHNAME=${MESHBASE%.*}

KAPPA=1e-4
NU=0
NUM_SAMPLES=1000

REFTOL=1e-6
TOLS=(1e-1 3e-2 1e-2 3e-3 1e-3 3e-4 1e-4 $REFTOL)
LOGPS=(4 6)

# set up julia
julia --project=. -e "using Pkg; Pkg.instantiate()"

# clean this directory and make output directory
rm -f -- *.bin *.txt
:> tols.txt
:> ps.txt
mkdir -p output

echo "TOL\tPRE(s)\t\tSAMP(s)\t\tCOMP(MB)\tUNCOMP(MB)\tUNTRUNC(MB)" >> $(printf "performance_kappa%.1e_nu%.1e.txt" $KAPPA $NU)

# compute covariance for various butterfly tolerances
for TOL in ${TOLS[@]}
do
    # compute relevant quantities with butterfly covariance
    cmd="./lbo_cov $MESH $KAPPA $NU $NUM_SAMPLES $TOL" > $(printf "log_lbo_tol%.1e_kappa%.1e_nu%.1e.txt" $TOL $KAPPA $NU) 
    echo "\n$cmd"; eval $cmd
    echo -n "$TOL " >> tols.txt

    # plot sample
    julia --project=. plot_sample.jl $MESH $KAPPA $NU lbo $TOL
done

echo "\nP\tSAMP(s)" >> $(printf "performance_kappa%.1e_nu%.1e.txt" $KAPPA $NU)

# compute covariances for various Chebyshev orders
for logP in ${LOGPS[@]}
do
    P="$((2**$logP))"

    # compute relevant quantities with Chebyshev covariance
    cmd="./cheb_cov $MESH $KAPPA $NU $NUM_SAMPLES $P" > $(printf "log_cheb_p%i_kappa%.1e_nu%.1e.txt" $P $KAPPA $NU)
    echo "\n$cmd"; eval $cmd
    echo -n "$P " >> ps.txt

    # plot sample
    julia --project=. plot_sample.jl $MESH $KAPPA $NU cheb $P
done

# plot error
julia --project=. plot_fast_cov_test.jl $MESH $KAPPA $NU

# move generated files to output folder
for file in *.bin *.txt
do
    mv $file output/"${MESHNAME}_$file"
done
