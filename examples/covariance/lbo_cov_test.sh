
MESH=$1
KAPPA=$2
NU=$3
NUM_SAMPLES=$4
FRACTION=1.0
LOGPS=(4 6 8 10)
TOLS=(1e-2 1e-4 1e-6)

mkdir -p output
:> ps.txt
:> tols.txt

# compute covariances for various Chebyshev orders
for logP in ${LOGPS[@]}
do
    let P=2**$logP
    cmd="./cheb_cov $MESH $P $KAPPA $NU $NUM_SAMPLES"
    printf "\n" 
    echo $cmd
    # eval $cmd
    echo -n "$P " >> ps.txt
done

# compute covariance for various butterfly tolerances
for TOL in ${TOLS[@]}
do
    cmd="./lbo_cov $MESH $KAPPA $NU $NUM_SAMPLES $TOL $FRACTION"
    printf "\n" 
    echo $cmd
    eval $cmd
    echo -n "$TOL " >> tols.txt
done

# move generated files to output folder
for file in *.bin *.txt
do
    mv $file output
done

# # make covariance and sample plots for each P
# cmd="julia plot_test_output.jl $MESH ps.txt"
# printf "\n"
# echo $cmd
# eval $cmd

# # make comparison plots against analytic covariance 
# cmd="julia covariance_sphere.jl $MESH $KAPPA $NU output/ps.txt output/tols.txt"
# printf "\n"
# echo $cmd
# eval $cmd