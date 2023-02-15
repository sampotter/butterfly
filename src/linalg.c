#include <bf/linalg.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/lu.h>
#include <bf/mat_dense_real.h>
#include <bf/util.h>
#include <bf/vec_real.h>

#include <arpack.h>

BfMat *bfMatSolveGMRES(BfMat const *A, BfMat const *B, BfMat *X0,
                       BfReal tol, BfSize maxNumIter, BfSize *numIter) {
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

  BfSize const m = maxNumIter;

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
  if (maxNumIter == 0)
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

  if (numIter != NULL)
    *numIter = i;

  return X;
}

static BfSize estimateNcv(BfSize nev, BfSize N) {
  BfSize ncv = 2*nev + 1;
  if (ncv < 20)
    ncv = 20;
  if (ncv > N)
    ncv = N;
  return ncv;
}

BfReal bfGetMaxEigenvalue(BfMat const *L, BfMat const *M) {
  /* TODO: this is a work in progress! This does NOT work for any type
   * of BfMat yet. Just real ones... */

  BEGIN_ERROR_HANDLING();

  a_int const N = bfMatGetNumRows(L);
  char const which[] = "LM";
  a_int const nev = 1;
  a_int const ncv = estimateNcv(nev, N);
  char bmat = 'G'; /* Solve (G)eneralized eigenvalue problem */
  BfReal const tol = 0;
  a_int const ldv = N;
  a_int const lworkl = 3*ncv*ncv + 6*ncv;
  a_int const rvec = 0; /* only computing eigenvalues */
  char const howmny[] = "A";

  BfLu *M_lu = bfLuNew();
  HANDLE_ERROR();

  bfLuInit(M_lu, M);
  HANDLE_ERROR();

  double *resid = malloc(N*sizeof(BfReal));
  if (resid == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  double *V = malloc(N*ncv*sizeof(BfReal));
  if (V == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  a_int *select = calloc(ncv, sizeof(a_int));

  a_int iparam[11] = {
    [0] = 1, /* compute exact shifts */
    [2] = 10*N, /* max number of iterations */
    // [3] = 1, /* only value allowed */
    [6] = 2, /* mode: A*x = lam*M*x */
  };

  a_int ipntr[11];

  double *workd = malloc(3*N*sizeof(BfReal));
  if (workd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (a_int i = 0; i < 3*N; ++i)
    workd[i] = 0;

  double *workl = malloc(lworkl*sizeof(BfReal));
  if (workl == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *workev = malloc(3*ncv*sizeof(BfReal));
  if (workev == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  a_int ido = 0;
  a_int info = 0;

dnaupd:
  dnaupd_c(&ido, &bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
           workd, workl, lworkl, &info);
  if (ido == 1 || ido == -1) {
    assert(ipntr[0] > 0);
    assert(ipntr[1] > 0);

    BfVecReal x;
    bfVecRealInitView(&x, N, BF_DEFAULT_STRIDE, &workd[ipntr[0] - 1]);
    BfVec *tmp = bfMatMulVec(L, bfVecRealToVec(&x));
    BfVecReal *y = bfVecToVecReal(bfLuSolve(M_lu, tmp));

    memcpy(&workd[ipntr[1] - 1], y->data, N*sizeof(BfReal));

    bfVecDelete(&tmp);
    bfVecRealDeinitAndDealloc(&y);

    goto dnaupd;
  } else if (ido == 2) {
    assert(ipntr[0] > 0);
    assert(ipntr[1] > 0);

    BfVecReal x;
    bfVecRealInitView(&x, N, BF_DEFAULT_STRIDE, &workd[ipntr[0] - 1]);

    BfVecReal *y = bfVecToVecReal(bfMatMulVec(M, bfVecRealToVec(&x)));

    memcpy(&workd[ipntr[1] - 1], y->data, N*sizeof(BfReal));

    bfVecRealDeinitAndDealloc(&y);

    goto dnaupd;
  }

  if (info < 0 || iparam[4] < nev)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfReal dr[2], di[2];

  dneupd_c(
    rvec, /* == 0 -> not computing Ritz vectors */
    howmny, /* == "A" -> compute all requested eigenvalues */
    select, /* used as internal workspace since howmny == "A" */
    dr, /* will contain real part of first nev + 1 eigenvalues */
    di,       /* ... imaginary part ... */
    NULL, /* not referenced since rvec == 0 */
    0,              /* ditto */
    0.0, /* not referenced since mode == 2 */
    0.0,           /* ditto */
    workev, /* internal workspace */

    /* dnaupd parameters: don't modify before calling dseupd */
    &bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
    workd, workl, lworkl, &info);

  if (info != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfReal eigmax = dr[0];

  if (fabs(di[0]) > 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {
    eigmax = NAN;
  }

  bfLuDeinitAndDealloc(&M_lu);

  free(resid);
  free(V);
  free(workd);
  free(workl);
  free(workev);

  return eigmax;
}

void bfGetShiftedEigs(BfMat const *A, BfMat const *M, BfReal sigma, BfSize k,
                      BfMat **PhiPtr, BfVecReal **LambdaPtr) {
  /* TODO: this is a work in progress! This does NOT work for any type
   * of BfMat yet. Just real ones... */

  BEGIN_ERROR_HANDLING();

  a_int const N = bfMatGetNumRows(A);
  char const which[] = "LM";
  a_int const nev = k;
  a_int const ncv = estimateNcv(nev, N);
  char bmat = 'G'; /* Solve (G)eneralized eigenvalue problem */
  BfReal const tol = 0;
  a_int const ldv = N;
  a_int const lworkl = 3*ncv*ncv + 6*ncv;
  a_int const rvec = 1; /* computing eigenvalues and eigenvectors */
  char const howmny[] = "A";

  BfMat *A_minus_sigma_M = bfMatCopy(M);
  HANDLE_ERROR();

  bfMatScale(A_minus_sigma_M, -sigma);
  HANDLE_ERROR();

  bfMatAddInplace(A_minus_sigma_M, A);
  HANDLE_ERROR();

  BfLu *A_minus_sigma_M_lu = bfLuNew();
  HANDLE_ERROR();

  bfLuInit(A_minus_sigma_M_lu, A_minus_sigma_M);
  HANDLE_ERROR();

  double *resid = malloc(N*sizeof(BfReal));
  if (resid == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  double *V = malloc(N*ncv*sizeof(BfReal));
  if (V == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  a_int *select = calloc(ncv, sizeof(a_int));

  a_int iparam[11] = {
    [0] = 1, /* compute exact shifts */
    [2] = 10*N, /* max number of iterations */
    [6] = 3, /* shift-invert mode */
  };

  a_int ipntr[11];

  double *workd = malloc(3*N*sizeof(BfReal));
  if (workd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (a_int i = 0; i < 3*N; ++i)
    workd[i] = 0;

  double *workl = malloc(lworkl*sizeof(BfReal));
  if (workl == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *workev = malloc(3*ncv*sizeof(BfReal));
  if (workev == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  a_int ido = 0;
  a_int info = 0;

dnaupd:
  dnaupd_c(&ido, &bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
           workd, workl, lworkl, &info);
  if (ido == 1 || ido == -1) {
    assert(ipntr[0] > 0);
    assert(ipntr[1] > 0);

    BfVecReal x;
    bfVecRealInitView(&x, N, BF_DEFAULT_STRIDE, &workd[ipntr[0] - 1]);
    BfVec *tmp = bfMatMulVec(M, bfVecRealToVec(&x));
    BfVecReal *y = bfVecToVecReal(bfLuSolve(A_minus_sigma_M_lu, tmp));

    memcpy(&workd[ipntr[1] - 1], y->data, N*sizeof(BfReal));

    bfVecDelete(&tmp);
    bfVecRealDeinitAndDealloc(&y);

    goto dnaupd;
  } else if (ido == 2) {
    assert(ipntr[0] > 0);
    assert(ipntr[1] > 0);

    BfVecReal x;
    bfVecRealInitView(&x, N, BF_DEFAULT_STRIDE, &workd[ipntr[0] - 1]);

    BfVecReal *y = bfVecToVecReal(bfMatMulVec(M, bfVecRealToVec(&x)));

    memcpy(&workd[ipntr[1] - 1], y->data, N*sizeof(BfReal));

    bfVecRealDeinitAndDealloc(&y);

    goto dnaupd;
  }

  if (info < 0 || iparam[4] < nev)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfReal *dr = malloc((nev + 1)*sizeof(BfReal));
  if (dr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *di = malloc((nev + 1)*sizeof(BfReal));
  if (di == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *z = malloc(N*(nev + 1)*sizeof(BfReal));
  if (z == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfSize ldz = N;

  BfReal sigmar = sigma;
  BfReal sigmai = 0.0;

  dneupd_c(
    rvec, /* == 0 -> not computing Ritz vectors */
    howmny, /* == "A" -> compute all requested eigenvalues */
    select, /* used as internal workspace since howmny == "A" */
    dr, /* will contain real part of first nev + 1 eigenvalues */
    di,       /* ... imaginary part ... */
    z,
    ldz,
    sigmar,
    sigmai,
    workev, /* internal workspace */

    /* dnaupd parameters: don't modify before calling dseupd */
    &bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
    workd, workl, lworkl, &info);

  if (info != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* TODO: evals could obviously be complex but right now we're only
   * interested in real evals and haven't introduced any other
   * safeguards to deal with complex case yet. So just treat complex
   * eigenvalues as an error for now. */
  for (BfSize i = 0; i < k; ++i)
    if (fabs(di[i]) != 0)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfSize *J = malloc(k*sizeof(BfSize));
  if (J == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  bfRealArgsort(dr, k, J);

  if (LambdaPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (*LambdaPtr != NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfVecReal *Lambda = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(Lambda, k);
  HANDLE_ERROR();

  for (BfSize j = 0; j < k; ++j)
    *(Lambda->data + j*Lambda->stride) = dr[J[j]];

  if (PhiPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (*PhiPtr != NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatDenseReal *Phi = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInit(Phi, N, k);
  HANDLE_ERROR();

  for (BfSize j = 0; j < k; ++j) {
    BfVecReal *col = bfVecRealNew();
    HANDLE_ERROR();

    bfVecRealInitView(col, N, BF_DEFAULT_STRIDE, z + N*J[j]);
    HANDLE_ERROR();

    bfMatDenseRealSetCol(bfMatDenseRealToMat(Phi), j, bfVecRealToVec(col));

    bfVecRealDeinitAndDealloc(&col);
  }

  *LambdaPtr = Lambda;
  *PhiPtr = bfMatDenseRealToMat(Phi);

  END_ERROR_HANDLING() {
    bfVecRealDeinitAndDealloc(&Lambda);
    bfMatDenseRealDeinitAndDealloc(&Phi);

    *LambdaPtr = NULL;
    *PhiPtr = NULL;
  }

  free(J);
  free(z);
  free(dr);
  free(di);

  bfLuDeinitAndDealloc(&A_minus_sigma_M_lu);
  bfMatDelete(&A_minus_sigma_M);

  free(resid);
  free(V);
  free(workd);
  free(workl);
  free(workev);
}

void bfGetEigenband(BfMat const *A, BfMat const *M, BfReal lam0, BfReal lam1,
                    BfReal sigma, BfMat **PhiPtr, BfVecReal **LambdaPtr) {
  BEGIN_ERROR_HANDLING();

  BfMat *Phi = NULL;
  BfVecReal *Lambda = NULL;

  if (PhiPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (*PhiPtr != NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (LambdaPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (*LambdaPtr != NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize k = 8;

  bool useLeft = isfinite(lam0);
  bool useRight = isfinite(lam1);

get_shifted_eigs:
  bfGetShiftedEigs(A, M, sigma, k, &Phi, &Lambda);
  HANDLE_ERROR();

  BfReal const *lam = Lambda->data;
  BfReal lamMin = lam[0];
  BfReal lamMax = lam[k - 1];
  if ((useLeft && lam0 < lamMin) || (useRight && lamMax < lam1)) {
    k *= 2;
    bfMatDelete(&Phi);
    bfVecRealDeinitAndDealloc(&Lambda);
    goto get_shifted_eigs;
  }

  /* Find the first eigenpair in the band */
  BfSize j0 = 0;
  if (useLeft) {
    while (j0 < k && lam[j0] < lam0)
      ++j0;
    assert(lam0 <= lam[j0]);
  }
  assert(j0 < k);

  /* Find the last eigenpair in the band */
  BfSize j1 = 0;
  if (useRight) {
    while (j1 < k && lam[j1] <= lam1)
      ++j1;
    assert(lam[j1 - 1] < lam1);
    assert(j1 == k || lam1 <= lam[j1]);
  }
  assert(j1 <= k);

  if (0 < j0 || j1 < k) {
    /* Prune unnecessary eigenvectors */
    BfMat *oldPhi = Phi;
    Phi = bfMatGetColRangeCopy(oldPhi, j0, j1);
    HANDLE_ERROR();

    /* Prune unnecessary eigenvalues */
    BfVecReal *oldLambda = Lambda;
    Lambda = bfVecToVecReal(bfVecGetSubvecCopy(bfVecRealToVec(oldLambda), j0, j1));
    HANDLE_ERROR();

    /* Free old eigenvectors and values */
    bfMatDelete(&oldPhi);
    bfVecRealDeinitAndDealloc(&oldLambda);
  }

  *PhiPtr = Phi;
  *LambdaPtr = Lambda;

  END_ERROR_HANDLING() {
    bfMatDelete(&Phi);
    bfVecRealDeinitAndDealloc(&Lambda);

    *PhiPtr = NULL;
    *LambdaPtr = NULL;
  }
}

bool bfGetTruncatedSvd(BfMat const *mat, BfMat **U, BfMatDiagReal **S, BfMat **VTrans,
                       BfTruncSpec truncSpec, BfBackend backend) {
  BEGIN_ERROR_HANDLING();

  bool success = true;

  if (backend != BF_BACKEND_LAPACK)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  (void)mat;
  (void)U;
  (void)S;
  (void)VTrans;
  (void)truncSpec;
  (void)backend;

  assert(false);

  END_ERROR_HANDLING() {}

  return success;
}
