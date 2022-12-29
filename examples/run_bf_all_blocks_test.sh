#!/usr/bin/env bash

usage() {
	>&2 echo "usage: $0 <shape> <n> <k> <layerPot> <blocks.txt>"
}

SHAPE="ellipse"
NUM_POINTS=4096
K=100
LAYER_POT="Sp"
BLOCKS_PATH="./blocks.txt"

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
	BLOCKS_PATH=$5
fi

echo "running bf_all_blocks example:"
echo "- shape: $SHAPE"
echo "- number of discretization points: n = $NUM_POINTS"
echo "- wavenumber: k = $K"
echo "- layer potential: $LAYER_POT"
echo "- saving blocks to $BLOCKS_PATH"

set -x

./make_${SHAPE}_test_data.py ${NUM_POINTS}

./bf_all_blocks ${K} \
				${SHAPE}${NUM_POINTS}_points.bin \
				${BLOCKS_PATH} \
				${LAYER_POT} \
				${SHAPE}${NUM_POINTS}_normals.bin

./plot_blocks.py ${BLOCKS_PATH}

./plot_bf_all_blocks.py

set +x
