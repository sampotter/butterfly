#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <bf/blas.h>
#include <bf/error_macros.h>
#include <bf/fac.h>
#include "bf/geom.h"
#include <bf/helm2.h>
#include <bf/layer_pot.h>
#include <bf/mat_block_coo.h>
#include <bf/mat_block_dense.h>
#include <bf/mat_block_diag.h>
#include <bf/mat_coo_complex.h>
#include <bf/mat_dense_complex.h>
#include <bf/mat_diag_real.h>
#include <bf/mat_solve.h>
#include <bf/mat_zero.h>
#include <bf/points.h>
#include <bf/quadrature.h>
#include <bf/quadtree.h>
#include <bf/rand.h>
#include <bf/util.h>
#include <bf/vectors.h>

static int const KR_order = 6;

struct K_helm2_wkspc {
  BfPoints2 const *points;
  BfVectors2 const *normals;
  BfReal K;
};

BfComplex K_helm2(BfSize i, BfSize j, void *aux) {
  struct K_helm2_wkspc *wkspc = aux;
  BfReal const *xsrc = &wkspc->points->data[i][0];
  BfReal const *xtgt = &wkspc->points->data[j][0];
  BfReal const *ntgt = &wkspc->normals->data[j][0];
  return bfHelm2GetKernelValue(
    xsrc, xtgt, ntgt, wkspc->K, BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE);
}

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

  BfVec *w = bfVecFromFile(argv[4], -1, BF_DTYPE_REAL);
  HANDLE_ERROR();

  BfMat *Xsrc = bfMatFromFile(argv[5], -1, 2, BF_DTYPE_REAL);
  HANDLE_ERROR();

  BfMat *Xtgt = bfMatFromFile(argv[6], -1, 2, BF_DTYPE_REAL);
  HANDLE_ERROR();

  /* Number of iterations */
  BfSize const numIter = 256;

  /* GMRES tolerance */
  BfReal const tol = 1e-10;

  /** Make sure everything is compatibly sized */

  BfSize n = bfMatGetNumRows(X);

  /** Convert matrices to points */

  BfPoints2 const *X_points = bfPoints2ConstViewFromMat(X);
  BfPoints2 const *X_source_points = bfPoints2ConstViewFromMat(Xsrc);
  BfPoints2 const *X_target_points = bfPoints2ConstViewFromMat(Xtgt);
  BfVectors2 const *N_vectors = bfVectors2ConstViewFromMat(N);

  /** Setup... miscellaneous things */

  /* Build a quadtree on the points sampling the boundary */
  bfToc();
  BfQuadtree quadtree;
  bfQuadtreeInit(&quadtree, X_points, N_vectors);
  HANDLE_ERROR();
  printf("built quadtree [%0.2fs]\n", bfToc());

  BfTree *tree = bfQuadtreeToTree(&quadtree);

  BfPerm const *perm = bfTreeGetPermConst(tree);

  /* Compute the reverse of the quadtree permutation */
  BfPerm revPerm = bfPermGetReversePerm(perm);

  /* Set up the LHS of the problem */
  // BfMat *phi_in = bf_hh2_get_dGdN(X_source_points, X_points, K, N_vectors);
  BfMat *phi_in = bfGetHelm2KernelMatrix(
    X_source_points, X_points, N_vectors, K, BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE);
  HANDLE_ERROR();

  /* Permute the LHS for the butterfly factorized version of A*/
  BfMat *phi_in_perm = bfMatCopy(phi_in);
  bfMatPermuteRows(phi_in_perm, &revPerm);

  /* One-half times the identity matrix---used to set up BIEs below */
  BfMat *oneHalfEye;
  { BfMatDiagReal *_ = bfMatDiagRealNew();
    bfMatDiagRealInit(_, n, n);
    bfMatDiagRealSetConstant(_, 1./2);
    oneHalfEye = bfMatDiagRealToMat(_); }

  /* Workspace for evaluating Helmholtz kernel when applying KR
   * quadrature corrections */
  struct K_helm2_wkspc K_wkspc = {
    .points = X_points,
    .normals = N_vectors,
    .K = K
  };

  /** Set up the dense system matrix */

  bfToc();

  /* Compute S_k' (normal derivative of single-layer potential) */
  BfMat *A_dense = bfGetHelm2KernelMatrix(
    X_points, X_points, N_vectors, K, BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE);
  HANDLE_ERROR();

  /* Perturb by the KR correction */
  bf_apply_KR_correction(A_dense, KR_order, K_helm2, (void *)&K_wkspc);

  /* Scale the columns by the trapezoid rule weights */
  bfMatScaleCols(A_dense, w);

  /* Perturb by one-half the identity to get a second-kind IE */
  bfMatAddInplace(A_dense, oneHalfEye);

  printf("assembled dense system matrix [%0.2fs]\n", bfToc());

  /** Set up the butterfly-factorized system matrix */

  bfToc();

  BfMat *A_BF = bfFacHelm2MakeMultilevel(
    &quadtree, K, BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE);
  HANDLE_ERROR();

  /* Perturb by the KR correction */
  bf_apply_KR_correction_quadtree(A_BF, KR_order, tree, K_helm2, (void *)&K_wkspc);

  /* Scale the columns by the trapezoid rule weights */
  BfVec *w_perm = bfVecCopy(w);
  bfVecPermute(w_perm, &revPerm);
  bfMatScaleCols(A_BF, w_perm);

  /* Perturb by one-half the identity to get a second-kind IE */
  bfMatAddInplace(A_BF, oneHalfEye);

  printf("assembled BF'd system matrix [%0.2fs]\n", bfToc());

  /* Write blocks to disk: */
  FILE *fp = fopen(argv[7], "w");
  bfPrintBlocks(A_BF, 2, fp);
  fclose(fp);
  printf("wrote blocks to %s\n", argv[7]);

  /** Copy the BF'd system matrix and do an MVP to verify that it was
   * copied correctly (recursive copying logic is a bit
   * involved---worth validating this here) */

  /** Verify and time the matrix multiplications */

  bfToc();
  BfMat *y_test_dense = bfMatMul(A_dense, phi_in);
  printf("did test dense MVP [%0.2fs]\n", bfToc());

  bfToc();
  BfMat *y_test_BF = bfMatMul(A_BF, phi_in_perm);
  bfMatPermuteRows(y_test_BF, perm);
  printf("did test butterfly MVP [%0.2fs]\n", bfToc());

  BfVec *err_MVP = bfMatColDists(y_test_dense, y_test_BF);
  BfVec *y_col_norms = bfMatColNorms(y_test_dense);
  BfReal max_rel_err_MVP = bfVecNormMax(err_MVP)/bfVecNormMax(y_col_norms);
  bfVecDelete(&err_MVP);
  bfVecDelete(&y_col_norms);
  printf("MVP rel. l2 errors:\n");
  printf("- BF: %g\n", max_rel_err_MVP);

  /** Solve the system using different methods */

  /* Solve dense system using LU decomposition */
  BfMat *A_dense_LU_copy = bfMatCopy(A_dense);
  bfToc();
  BfMat *sigma_dense_LU = bfMatSolveLU(A_dense_LU_copy, phi_in);
  printf("solved using Gaussian elimination [%0.2fs]\n", bfToc());
  bfMatDelete(&A_dense_LU_copy);

  /* Solve dense system using GMRES */
  bfToc();
  BfSize num_iter_dense_GMRES;
  BfMat *sigma_dense_GMRES = bfMatSolveGMRES(
    A_dense, phi_in, NULL, tol, numIter, &num_iter_dense_GMRES);
  printf("solved using GMRES (dense): %lu iter. [%0.2fs]\n",
         num_iter_dense_GMRES, bfToc());

  /* Solve butterfly-factorized system using GMRES */
  bfToc();
  BfSize num_iter_BF_GMRES;
  BfMat *sigma_BF_GMRES = bfMatSolveGMRES(
    A_BF, phi_in_perm, NULL, tol, numIter, &num_iter_BF_GMRES);
  bfMatPermuteRows(sigma_BF_GMRES, perm);
  printf("solved using GMRES (BF): %lu iter. [%0.2fs]\n",
         num_iter_BF_GMRES, bfToc());

  /** Evaluate solution and check errors */

  /* Set up evaluation matrix */
  BfMat *G_eval = bfGetHelm2KernelMatrix(
    X_points, X_target_points, NULL, K, BF_LAYER_POTENTIAL_SINGLE);
  bfMatScaleCols(G_eval, w);

  BfMat *phi_exact = bfGetHelm2KernelMatrix(
    X_source_points, X_target_points, NULL, K, BF_LAYER_POTENTIAL_SINGLE);

  BfMat *phi_dense_LU = bfMatMul(G_eval, sigma_dense_LU);
  BfMat *phi_dense_GMRES = bfMatMul(G_eval, sigma_dense_GMRES);
  BfMat *phi_BF_GMRES = bfMatMul(G_eval, sigma_BF_GMRES);

  /* Compute errors */

  BfVec *error_l2_dense_LU_vec = bfMatColDists(phi_dense_LU, phi_exact);
  BfReal error_l2_dense_LU = *(BfReal *)bfVecGetEltPtr(error_l2_dense_LU_vec, 0);

  BfVec *error_l2_dense_GMRES_vec = bfMatColDists(phi_dense_GMRES, phi_exact);
  BfReal error_l2_dense_GMRES = *(BfReal *)bfVecGetEltPtr(error_l2_dense_GMRES_vec, 0);

  BfVec *error_l2_BF_GMRES_vec = bfMatColDists(phi_BF_GMRES, phi_exact);
  BfReal error_l2_BF_GMRES = *(BfReal *)bfVecGetEltPtr(error_l2_BF_GMRES_vec, 0);

  BfVec *phi_l2_norm_vec = bfMatColNorms(phi_exact);
  BfReal phi_l2_norm = *(BfReal *)bfVecGetEltPtr(phi_l2_norm_vec, 0);

  printf("rel l2 errors:\n");
  printf("- dense LU: %g\n", error_l2_dense_LU/phi_l2_norm);
  printf("- dense GMRES: %g\n", error_l2_dense_GMRES/phi_l2_norm);
  printf("- BF GMRES: %g\n", error_l2_BF_GMRES/phi_l2_norm);

  END_ERROR_HANDLING() {}
}
