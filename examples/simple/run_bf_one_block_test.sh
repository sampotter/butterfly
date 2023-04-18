#!/usr/bin/env bash

usage() {
	>&2 echo "usage: $0 <shape> <n> <k> <layerPot> <srcDepth> <srcIndex> <tgtDepth> <tgtIndex>"
}

SHAPE="ellipse"
NUM_POINTS=4096
K=100
LAYER_POT="Sp"
SRC_DEPTH=2
SRC_INDEX=1
TGT_DEPTH=2
TGT_INDEX=10

if (( $# >= 1 )); then
	SHAPE=$1
fi

if (( $# >= 2 )); then
	NUM_POINTS=$2
fi

if (( $# >= 3 )); then
	K=$3
fi

if (( $# >= 4 )); then
	LAYER_POT=$4
fi

if (( $# >= 5 )); then
	SRC_DEPTH=$5
fi

if (( $# >= 6 )); then
	SRC_INDEX=$6
fi

if (( $# >= 7 )); then
	TGT_DEPTH=$7
fi

if (( $# >= 8 )); then
	TGT_INDEX=$8
fi

echo "running bf_one_block example:"
echo "- shape: $SHAPE"
echo "- number of discretization points: n = $NUM_POINTS"
echo "- wavenumber: k = $K"
echo "- layer potential: $LAYER_POT"
echo "- source node: depth = $SRC_DEPTH, index = $SRC_INDEX"
echo "- target node: depth = $TGT_DEPTH, index = $TGT_INDEX"

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
