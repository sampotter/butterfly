#!/usr/bin/env bash

usage() {
	>&2 echo "usage: $0 <shape> <n> <k> <blocks.txt>"
}

SHAPE="ellipse"
NUM_POINTS=4096
K=100
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
	BLOCKS_PATH=$4
fi

echo "running helm2_bie example:"
echo "- shape: $SHAPE"
echo "- number of discretization points: n = $NUM_POINTS"
echo "- wavenumber: k = $K"
echo "- saving blocks to $BLOCKS_PATH"

./make_${SHAPE}_test_data.py ${NUM_POINTS}

./helm2_bie ${K} \
			${SHAPE}${NUM_POINTS}_points.bin \
			${SHAPE}${NUM_POINTS}_normals.bin \
			${SHAPE}${NUM_POINTS}_weights.bin \
			${SHAPE}${NUM_POINTS}_sources.bin \
			${SHAPE}${NUM_POINTS}_targets.bin \
			${BLOCKS_PATH}
