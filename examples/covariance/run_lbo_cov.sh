#!/usr/bin/env bash

OBJ_PATH=../../../../butterfly-LBO-models/armadillo3.obj
KAPPA=0.1
NU=0 # Gaussian spectral decay
NUM_SAMPLES=100
TOL=1e-2
FRACTION=1.0
ROW_TREE_OFFSET=0
FREQ_TREE_OFFSET=8

./lbo_cov $OBJ_PATH $KAPPA $NU $NUM_SAMPLES $TOL $FRACTION