#!/usr/bin/env bash

OBJ_PATH=dragon0.obj
TOL=1e-3
P=0.0625
SPACE_TREE_OFFSET=0
FREQ_TREE_OFFSET=5

echo "running bf_lbo example:"
echo "- obj file: $OBJ_PATH"
echo "- tol: $TOL"
echo "- space tree offset: $SPACE_TREE_OFFSET"
echo "- freq tree offset: $FREQ_TREE_OFFSET"

set -x

./bf_lbo $OBJ_PATH $TOL $P $SPACE_TREE_OFFSET $FREQ_TREE_OFFSET

set +x
