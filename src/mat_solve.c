#include <bf/mat_solve.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

BfMat *bfMatSolveGMRES(BfMat const *A, BfMat const *B, BfMat *X0,
                       BfReal tol, BfSize numIter) {
  BEGIN_ERROR_HANDLING();

  /* Solution of the system */
  BfMat *X = NULL;

  /* Intermediate matrices */
  BfMat *R = NULL;
  BfMat **V = NULL;
  BfMat **H = NULL;
  BfMat *S = NULL;

  /* Givens rotations for solving the least squares problem */
  BfMat **J = NULL;

  BfSize const m = numIter;

  bool converged = false;

  /* Make sure B is a dense matrix */
  if (!bfMatInstanceOf(B, BF_TYPE_MAT_DENSE_REAL) &&
      !bfMatInstanceOf(B, BF_TYPE_MAT_DENSE_COMPLEX))
    RAISE_ERROR(BF_ERROR_TYPE_ERROR);

  /* Make sure that either X0 wasn't passed or is a dense matrix */
  if (X0 != NULL)
    if (!bfMatInstanceOf(X0, BF_TYPE_MAT_DENSE_REAL) &&
        !bfMatInstanceOf(X0, BF_TYPE_MAT_DENSE_COMPLEX))
      RAISE_ERROR(BF_ERROR_TYPE_ERROR);

  /* Make sure m is positive */
  if (numIter == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Get order of system (n) and check compatibility */
  BfSize n = bfMatGetNumRows(A);
  if (n != bfMatGetNumCols(A))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
  if (n != bfMatGetNumRows(B))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
  if (X0 != NULL && n != bfMatGetNumRows(X0))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Get number of RHSs (p) and check compatibility */
  BfSize p = bfMatGetNumCols(B);
  if (X0 != NULL && p != bfMatGetNumCols(X0))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  V = malloc((m + 1)*sizeof(BfMatDenseComplex *));
  if (V == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  H = malloc(m*sizeof(BfMatDenseComplex *));
  if (H == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  J = malloc(m*p*sizeof(BfMat *));
  if (J == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* If an initial iterate wasn't passed, set it to zero */
  if (X0 == NULL) {
    X0 = bfMatZerosLike(B, n, p);
    HANDLE_ERROR();
  }

  /* Compute residual for first iteration */
  BfMat *Y = bfMatMul(A, X0);
  HANDLE_ERROR();
  R = bfMatSub(B, Y);
  HANDLE_ERROR();
  bfMatDelete(&Y);

  /* Compute the norm of the residual for each righthand side */
  BfVec *R_col_norms = bfMatColNorms(R);
  BfReal beta = bfVecNormMax(R_col_norms);

  /* Set the first basis vector to the normalized residual */
  V[0] = bfMatCopy(R);
  bfVecRecipInplace(R_col_norms);
  bfMatScaleCols(V[0], R_col_norms);

  /* Initialize S to ||r||_2*e1---need to take the reciprocal of
   * `R_col_norms` again first */
  S = bfMatZerosLike(B, m + 1, p);
  bfVecRecipInplace(R_col_norms);
  bfMatSetRow(S, 0, R_col_norms);

  BfSize i = 0;
  for (i = 0; i < m; ++i) {
    BfMat *W = bfMatMul(A, V[i]);

    H[i] = bfMatEmptyLike(W, i + 2, p);

    for (BfSize k = 0; k <= i; ++k) {
      BfVec *Hk = bfMatColDots(W, V[k]);
      bfMatSetRow(H[i], k, Hk);

      BfMat *Vk_Hk = bfMatCopy(V[k]);
      bfMatScaleCols(Vk_Hk, Hk);
      bfMatSubInplace(W, Vk_Hk);

      bfMatDelete(&Vk_Hk);
      bfVecDelete(&Hk);
    }

    BfVec *W_col_norms = bfMatColNorms(W);
    bfMatSetRow(H[i], i + 1, W_col_norms);

    V[i + 1] = bfMatCopy(W);
    bfVecRecipInplace(W_col_norms);
    bfMatScaleCols(V[i + 1], W_col_norms);
    bfMatDelete(&W);
    bfVecDelete(&W_col_norms);

    for (BfSize k = 0; k < i; ++k) {
      for (BfSize j = 0; j < p; ++j) {
        BfVec *h = bfMatGetColRangeView(H[i], 0, k + 2, j);
        bfVecSolveInplace(h, J[p*k + j]);
        bfVecDelete(&h);
      }
    }

    for (BfSize j = 0; j < p; ++j) {
      BfVec *h = bfMatGetColView(H[i], j);
      J[p*i + j] = bfVecGetGivensRotation(h, i, i + 1);
      bfVecSolveInplace(h, J[p*i + j]);
      bfVecDelete(&h);
    }

    /* Apply most recent set of Givens rotations to S */
    for (BfSize j = 0; j < p; ++j) {
      BfVec *s = bfMatGetColView(S, j);
      BfVec *sSub = bfVecGetSubvecView(s, 0, i + 2);
      bfVecSolveInplace(sSub, J[p*i + j]);
      bfVecDelete(&sSub);
      bfVecDelete(&s);
    }

    BfVec *sLastRow = bfMatGetRowView(S, i + 1);
    if (bfVecNormMax(sLastRow) < tol*beta)
      converged = true;
    bfVecDelete(&sLastRow);

    if (converged)
      break;
  }

  /* Construct the solution */
  X = bfMatEmptyLike(B, n, p);
  for (BfSize j = 0; j < p; ++j) {
    /* Extract the version of H for this RHS */
    BfMat *Hj = bfMatZerosLike(H[0], i + 1, i + 1);
    for (BfSize k = 0; k < i + 1; ++k) {
      BfVec *h = bfMatGetColView(H[k], j);
      bfMatSetColRange(Hj, k, 0, h->size - 1, h);
    }
    assert(bfMatIsUpperTri(Hj));

    /* Extract the version of V for this RHS */
    BfMat *Vj = bfMatEmptyLike(V[0], n, i + 1);
    for (BfSize k = 0; k < i + 1; ++k) {
      BfVec *v = bfMatGetColView(V[k], j);
      bfMatSetCol(Vj, k, v);
    }

    /* Solve for the coefficients representing x - x0 in the V-basis */
    BfVec *s = bfMatGetColRangeView(S, 0, i + 1, j);
    BfVec *y = bfMatBackwardSolveVec(Hj, s);
    BfVec *x0 = bfMatGetColView(X0, j);
    BfVec *x = bfMatMulVec(Vj, y);
    bfVecAddInplace(x, x0);
    bfMatSetCol(X, j, x);
  }

  END_ERROR_HANDLING() {}

  return X;
}
