#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/fac_helm2.h>
#include <bf/helm2.h>
#include <bf/linalg.h>
#include <bf/mat_diag_real.h>
#include <bf/points.h>
#include <bf/util.h>
#include <bf/vectors.h>

#include <stdlib.h>

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
    xsrc, xtgt, NULL, ntgt, wkspc->K, BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE);
}

int main(int argc, char const *argv[]) {
  if (argc != 8) {
    printf("usage: %s <K> <points.bin> <normals.bin> <weights.bin> "
           "<sources.bin> <targets.bin> <blocks.txt>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  BF_ERROR_BEGIN();

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
  BfPerm *revPerm = bfPermGetReversePerm(perm);

  /* Set up the LHS of the problem */
  BfMat *phi_in = bfGetHelm2KernelMatrix(
    X_source_points, X_points, NULL, N_vectors, K,
    BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE, NULL, NULL);
  HANDLE_ERROR();

  /* Permute the LHS for the butterfly factorized version of A*/
  BfMat *phi_in_perm = bfMatCopy(phi_in);
  bfMatPermuteRows(phi_in_perm, revPerm);

  /* One-half times the identity matrix---used to set up BIEs below */
  BfMat *oneHalfEye = bfMatDiagRealToMat(bfMatDiagRealNewConstant(n, n, 1./2));

  /* Workspace for evaluating Helmholtz kernel when applying KR
   * quadrature corrections */
  struct K_helm2_wkspc K_wkspc = {
    .points = X_points,
    .normals = N_vectors,
    .K = K
  };

  /** Set up the butterfly-factorized system matrix */

  bfToc();

  BfMat *A_BF = bfFacHelm2MakeMultilevel(
    &quadtree, &quadtree, K, BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE, NULL, NULL);
  HANDLE_ERROR();

  /* Perturb by the KR correction */
  bf_apply_KR_correction_quadtree(A_BF, KR_order, tree, K_helm2, (void *)&K_wkspc);

  /* Scale the columns by the trapezoid rule weights */
  BfVec *w_perm = bfVecCopy(w);
  bfVecPermute(w_perm, revPerm);
  bfMatScaleCols(A_BF, w_perm);

  /* Perturb by one-half the identity to get a second-kind IE */
  bfMatAddInplace(A_BF, oneHalfEye);

  printf("assembled BF'd system matrix [%0.2fs]\n", bfToc());

  /* Write blocks to disk: */
  FILE *fp = fopen(argv[7], "w");
  bfPrintBlocks(A_BF, 2, fp);
  fclose(fp);
  printf("wrote blocks to %s\n", argv[7]);

  /** Solve the system using different methods */

  /* Solve butterfly-factorized system using GMRES */
  bfToc();
  BfSize num_iter_BF_GMRES;
  BfMat *sigma_BF_GMRES = bfSolveGMRES(
    A_BF, phi_in_perm, NULL, tol, numIter, &num_iter_BF_GMRES, NULL);
  bfMatPermuteRows(sigma_BF_GMRES, perm);
  printf("solved using GMRES (BF): %lu iter. [%0.2fs]\n",
         num_iter_BF_GMRES, bfToc());

  /** Evaluate solution and check errors */

  /* Set up evaluation matrix */
  BfMat *G_eval = bfGetHelm2KernelMatrix(
    X_points, X_target_points, NULL, NULL, K, BF_LAYER_POTENTIAL_SINGLE, NULL, NULL);
  bfMatScaleCols(G_eval, w);

  BfMat *phi_exact = bfGetHelm2KernelMatrix(
    X_source_points, X_target_points, NULL, NULL, K, BF_LAYER_POTENTIAL_SINGLE, NULL, NULL);

  BfMat *phi_BF_GMRES = bfMatMul(G_eval, sigma_BF_GMRES);

  /* Compute errors */

  BfVec *error_l2_BF_GMRES_vec = bfMatColDists(phi_BF_GMRES, phi_exact);
  BfReal error_l2_BF_GMRES = *(BfReal *)bfVecGetEltPtr(error_l2_BF_GMRES_vec, 0);

  BfVec *phi_l2_norm_vec = bfMatColNorms(phi_exact);
  BfReal phi_l2_norm = *(BfReal *)bfVecGetEltPtr(phi_l2_norm_vec, 0);

  printf("rel l2 errors:\n");
  printf("- BF GMRES: %g\n", error_l2_BF_GMRES/phi_l2_norm);

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMatDelete(&X);
  bfMatDelete(&N);
  bfVecDelete(&w);
  bfMatDelete(&Xsrc);
  bfMatDelete(&Xtgt);

  bfPoints2DeinitAndDealloc((BfPoints2 **)&X_points);
  bfPoints2DeinitAndDealloc((BfPoints2 **)&X_source_points);
  bfPoints2DeinitAndDealloc((BfPoints2 **)&X_target_points);
  bfVectors2DeinitAndDealloc((BfVectors2 **)&N_vectors);

  bfQuadtreeDeinit(&quadtree);
  bfPermDeinitAndDealloc(&revPerm);

  bfMatDelete(&phi_in);
  bfMatDelete(&phi_in_perm);
  bfMatDelete(&oneHalfEye);

  bfMatDelete(&A_BF);
  bfVecDelete(&w_perm);

  bfMatDelete(&sigma_BF_GMRES);

  bfMatDelete(&G_eval);
  bfMatDelete(&phi_exact);
  bfMatDelete(&phi_BF_GMRES);

  bfVecDelete(&error_l2_BF_GMRES_vec);
  bfVecDelete(&phi_l2_norm_vec);
}
