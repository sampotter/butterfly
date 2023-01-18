#!/usr/bin/env bash

OBJ_PATH=dragon0.obj
TOL=1e-3
INIT_DEPTH=0

echo "running bf_lbo example:"
echo "- obj file: $OBJ_PATH"
echo "- tol: $TOL"
echo "- initial depth: $INIT_DEPTH"

set -x

./bf_lbo $OBJ_PATH $TOL $INIT_DEPTH

set +x
