#!/usr/bin/env bash

SHAPE=$1
NUM_POINTS=$2
K=$3
LAYER_POT=$4
SRC_DEPTH=$5
SRC_INDEX=$6
TGT_DEPTH=$7
TGT_INDEX=$8

rm -f *.bin
rm -f *.txt
rm -rf *factor*
rm -rf plots

./make_${SHAPE}_test_data.py ${NUM_POINTS}
./bf_one_block ${K} \
			   ${SRC_DEPTH} ${SRC_INDEX} \
			   ${TGT_DEPTH} ${TGT_INDEX} \
			   ${SHAPE}${NUM_POINTS}_points.bin \
			   ${LAYER_POT} \
			   ${SHAPE}${NUM_POINTS}_normals.bin
./plot_bf_one_block.py
