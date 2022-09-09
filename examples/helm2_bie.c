#include "bf/geom.h"
#include <assert.h>
#include <stdio.h>

#include <openblas/cblas.h>

#include <bf/error_macros.h>
#include <bf/fac.h>
#include <bf/helm2.h>
#include <bf/layer_pot.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_block_dense.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_dense_complex.h>
#include <bf/mat_diag_real.h>
#include <bf/mat_util.h>
#include <bf/points.h>
#include <bf/quadrature.h>
#include <bf/quadtree.h>
#include <bf/rand.h>
#include <bf/util.h>
#include <bf/vectors.h>

int main(int argc, char const *argv[]) {
  if (argc != 8) {
    printf("usage: %s <K> <points.bin> <normals.bin> <weights.bin> "
           "<sources.bin> <targets.bin> <blocks.txt>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  BEGIN_ERROR_HANDLING();

  BfReal K = atoi(argv[1]);

  BfMat *X = bfMatFromFile(argv[2], -1, 2, BF_DTYPE_REAL);
  HANDLE_ERROR();

  BfMat *N = bfMatFromFile(argv[3], -1, 2, BF_DTYPE_REAL);
  HANDLE_ERROR();

  BfMat *w = bfMatFromFile(argv[4], -1, 1, BF_DTYPE_REAL);
  HANDLE_ERROR();

  BfMat *Xsrc = bfMatFromFile(argv[5], -1, 2, BF_DTYPE_REAL);
  HANDLE_ERROR();

  BfMat *Xtgt = bfMatFromFile(argv[6], -1, 2, BF_DTYPE_REAL);
  HANDLE_ERROR();

  /** Make sure everything is compatibly sized */

  BfSize n = bfMatGetNumRows(X);

  /** Convert matrices to points */

  BfPoints2 const *X_points = bfPoints2ConstViewFromMat(X);
  BfPoints2 const *X_source_points = bfPoints2ConstViewFromMat(Xsrc);
  BfPoints2 const *X_target_points = bfPoints2ConstViewFromMat(Xtgt);
  BfVectors2 const *N_vectors = bfVectors2ConstViewFromMat(N);

  /* Build a quadtree on the points sampling the boundary */
  bfToc();
  BfQuadtree tree;
  bfInitQuadtreeFromMat(&tree, X);
  HANDLE_ERROR();
  printf("built quadtree [%0.2fs]\n", bfToc());

  /** Discretize Helmholtz BIE with Neumann BCs using Kapur-Rokhlin
   *  quadrature. */

  /* Set up the LHS of the problem */
  BfMat *phi_in = bf_hh2_get_dGdN(X_source_points, X_points, K, N_vectors);
  HANDLE_ERROR();

  /* Evaluate the normal derivative of the SLP on the boundary (i.e.,
   * compute S_k') */
  BfMat *A_dense = bf_hh2_get_dGdN(X_points, X_points, K, N_vectors);
  HANDLE_ERROR();

  /* Scale the columns by the discretization weights */
  bfMatScaleCols(A_dense, w);

  /* Apply the KR quadrature correction to the system matrix */
  bf_apply_KR_correction(A_dense, 6);

  /* Perturb by one-half the identity to get a second-kind IE */
  BfMat *oneHalfEye;
  {
    BfMatDiagReal *_ = bfMatDiagRealNew();
    bfMatDiagRealInit(_, n, n);
    bfMatDiagRealSetConstant(_, 1./2);
    oneHalfEye = bfMatDiagRealToMat(_);
  }
  bfMatAddInplace(A_dense, oneHalfEye);

  /* Solve the discretized BIE */
  BfMat *sigma = bfMatSolve(A_dense, phi_in);

  BfMat *G_eval = bfMatDenseComplexToMat(bfGetHelm2KernelMatrix(X_points, X_target_points, K));

  bfMatScaleCols(G_eval, w);

  BfMat *phi_dense = bfMatMul(G_eval, sigma);

  BfMat *phi_exact = bfMatDenseComplexToMat(bfGetHelm2KernelMatrix(X_source_points, X_target_points, K));

  BfMat *error_l2_dense = bfMatColDists(phi_dense, phi_exact);

  bfMatPrint(stdout, error_l2_dense);

  END_ERROR_HANDLING() {}
}
