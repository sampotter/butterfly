# Example: `simple`

This directory contains several simple tests with shell scripts to run them.

## `run_bf_one_block_test.sh`

This script will set up a problem and butterfly /one/ of its off-diagonal blocks. It will fail if the block doesn't correspond to a pair of well-separated nodes. Afterwards, it will make some plots which are helpful for debugging.

## `run_bf_all_blocks_test.sh`

Sets up a problem and builds the multilevel butterfly factorization of the kernel matrix for a given layer potential. Runs some basic tests to check for error in the matrix-vector product. Afterwards, makes debugging plots, including diagrammatic plots of the multilevel butterfly factorization and the pointwise error in the approximate kernel matrix.

## `run_helm2_bie.sh`

Solves a simple exterior Helmholtz problem which has been discretized using Kapur-Rokhlin quadrature. The problem is solved three ways:

- using the dense system matrix and Gaussian elimination
- again using the dense system matrix, but using GMRES instead of GE
- using the multilevel butterfly factorization of the system matrix and GMRES

Timings and errors are reported, both for the matrix-vector products of the approximated system matrix and the errors in the computed field.
