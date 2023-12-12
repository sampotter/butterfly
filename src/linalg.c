#include <bf/linalg.h>

#include <math.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/disjoint_interval_list.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/logging.h>
#include <bf/lu_csr_real.h>
#include <bf/mat_dense_real.h>
#include <bf/mem.h>
#include <bf/ptr_array.h>
#include <bf/real_array.h>
#include <bf/timer.h>
#include <bf/util.h>
#include <bf/vec_real.h>

#if BF_DEBUG
#include <bf/vec_complex.h>
#endif

#include <arpack.h>

BfSize bfTruncSpecGetNumTerms(BfTruncSpec const *truncSpec, BfMatDiagReal const *S) {
  BfSize k = 0;
  if (truncSpec->usingTol) {
    while (k < S->numElts && S->data[k] >= truncSpec->tol*S->data[0])
      ++k;
  } else {
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
  }
  return k;
}

/* Solve the linear system with matrix `A`, RHS `B`, and initial guess
 * `X0`. Tolerance and maximum number of iterations specified by `tol`
 * and `maxNumIter`, respectively. The final iteration count will be
 * passed back in `numIter` if `numIter != NULL`.
 *
 * For left preconditioned GMRES, set the preconditioner in `M`.  Note
 * that the residual will be determined from the *preconditioned*
 * residual vectors. See Saad for more details.
 *
 * TODO: right and split preconditioned GMRES. */
BfMat *bfSolveGMRES(BfMat const *A, BfMat const *B, BfMat *X0,
                    BfReal tol, BfSize maxNumIter, BfSize *numIter,
                    BfMat const *M) {
  BF_ERROR_BEGIN();

  /* Solution of the system */
  BfMat *X = NULL;

  /* Intermediate matrices */
  BfMat *R = NULL;
  BfMat **V = NULL;
  BfMat **H = NULL;
  BfMat *S = NULL;

  /* Givens rotations for solving the least squares problem */
  BfMat **J = NULL;

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

  bool leftPrecond = M != NULL;

  /* If we're using a left preconditioner, make sure it has a
   * compatible shape. */
  if (leftPrecond) {
    if (n != bfMatGetNumRows(M))
      RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
    if (n != bfMatGetNumCols(M))
      RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
  }

  /* Get number of RHSs (P) and check compatibility */
  BfSize numRhs = bfMatGetNumCols(B);
  if (X0 != NULL && numRhs != bfMatGetNumCols(X0))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  V = bfMemAllocAndZero(maxNumIter + 1, sizeof(BfMatDenseComplex *));
  if (V == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  H = bfMemAllocAndZero(maxNumIter, sizeof(BfMatDenseComplex *));
  if (H == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  J = bfMemAllocAndZero(numRhs*maxNumIter, sizeof(BfMat *));
  if (J == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* If an initial iterate wasn't passed, set it to zero */
  if (X0 == NULL) {
    X0 = bfMatZerosLike(B, n, numRhs);
    HANDLE_ERROR();
  }

  /* Compute residual for first iteration */
  BfMat *Y = bfMatMul(A, X0);
  HANDLE_ERROR();
  R = bfMatSub(B, Y);
  HANDLE_ERROR();
  bfMatDelete(&Y);
  if (leftPrecond) {
    BfMat *_ = bfMatSolve(M, R);
    bfMatDelete(&R);
    R = _;
  }

  /* Compute the norm of the residual for each righthand side */
  BfVec *RNorm = bfMatColNorms(R);

  // TODO: beta is different for each RHS
  BfReal beta = bfVecNormMax(RNorm);

  /* Set the first basis vector to the normalized residual */
  V[0] = bfMatCopy(R);
  bfMatDivideCols(V[0], RNorm);

  /* Initialize S to ||r||_2*e1---need to take the reciprocal of
   * `R_col_norms` again first */
  S = bfMatZerosLike(B, maxNumIter + 1, numRhs);
  bfMatSetRow(S, 0, RNorm);

  /* The iteration count/the current column: */
  BfSize j = 0;

  for (j = 0; j < maxNumIter; ++j) {
    BfMat *W = bfMatMul(A, V[j]);

    /* Apply left preconditioner if we're using one: */
    if (leftPrecond) {
      BfMat *_ = bfMatSolve(M, W);
      bfMatDelete(&W);
      W = _;
    }

    /** Run modified Gram-Schmidt: */

    /* Compute W column norms before applying Gram-Schmidt: */
    BfVec *WNormBefore = bfMatColNorms(W); // 433

    /* Allocate space for next column of H: */
    H[j] = bfMatEmptyLike(W, j + 2, numRhs);

    for (BfSize i = 0; i <= j; ++i) {
      BfVec *Hij = bfMatColDots(V[i], W); // 435
      bfMatSetRow(H[j], i, Hij);          // "

      BfMat *Hij_Vi = bfMatCopy(V[i]);
      bfMatScaleCols(Hij_Vi, Hij); // 436
      bfMatSubInplace(W, Hij_Vi);  // "

      bfMatDelete(&Hij_Vi);
      bfVecDelete(&Hij);
    }

    BfVec *WNorm = bfMatColNorms(W); // 438
    bfMatSetRow(H[j], j + 1, WNorm); // 439

#if BF_DEBUG // TODO: deal with breakdown
    {
      BF_ASSERT(numRhs == 1); // TODO: handle numRhs > 1
      BfReal _ = bfVecNormMax(WNorm)/bfVecNormMax(WNormBefore);
      BF_ASSERT(_ >= tol);
    }
#endif

    V[j + 1] = bfMatCopy(W); // 440
    bfMatDivideCols(V[j + 1], WNorm); // 447 & 448

    bfMatDelete(&W);
    bfVecDelete(&WNorm);
    bfVecDelete(&WNormBefore);

    /** Use Givens rotations to reduce H to upper triangular form: */

    for (BfSize i = 0; i < j; ++i) {
      for (BfSize p = 0; p < numRhs; ++p) {
        BfVec *h = bfMatGetColRangeView(H[j], 0, i + 2, p);
        bfVecMulInplace(h, J[numRhs*i + p]);
        bfVecDelete(&h);
      }
    }

    for (BfSize p = 0; p < numRhs; ++p) {
      BfVec *h = bfMatGetColView(H[j], p);
      J[numRhs*j + p] = bfVecGetGivensRotation(h, j, j + 1);
      bfVecMulInplace(h, J[numRhs*j + p]);
      bfVecDelete(&h);
    }

    /* Apply most recent set of Givens rotations to S */
    for (BfSize p = 0; p < numRhs; ++p) {
      BfVec *s = bfMatGetColView(S, p);
      BfVec *sSub = bfVecGetSubvecView(s, 0, j + 2);
      bfVecMulInplace(sSub, J[numRhs*j + p]);
      bfVecDelete(&sSub);
      bfVecDelete(&s);
    }

    BfVec *sLastRow = bfMatGetRowView(S, j + 1);
    BfReal residual = bfVecNormMax(sLastRow)/beta;

    /** Deal with convergence: */

    if (residual < tol)
      converged = true;

    bfVecDelete(&sLastRow);

    if (converged)
      break;
  }

  /* Construct the solution */
  X = bfMatEmptyLike(B, n, numRhs);
  for (BfSize p = 0; p < numRhs; ++p) {
    /* Extract the triangularized version of the upper Hessenberg
     * matrix H for the current RHS: */
    BfMat *Hp = bfMatZerosLike(H[0], j, j);
    for (BfSize i = 0; i < j; ++i) {
      BfVec *h = bfMatGetColView(H[i], p);
#if BF_DEBUG
      if (bfVecGetType(h) == BF_TYPE_VEC_COMPLEX) {
        BfVecComplex *_ = bfVecToVecComplex(h);
        BfComplex z = *(_->data + (i + 1)*_->stride);
        BF_ASSERT(cabs(z) <= 1e-15);
      }
#endif
      bfMatSetColRange(Hp, i, 0, i + 1, h);
      bfVecDelete(&h);
    }

    /* Extract the version of V for this RHS */
    BfMat *Vp = bfMatEmptyLike(V[0], n, j);
    for (BfSize i = 0; i < j; ++i) {
      BfVec *v = bfMatGetColView(V[i], p);
      bfMatSetCol(Vp, i, v);
      bfVecDelete(&v);
    }

    /* Solve for the coefficients representing x - x0 in the V-basis */
    BfVec *s = bfMatGetColRangeView(S, 0, j, p);
    BfVec *y = bfMatBackwardSolveVec(Hp, s);
    BfVec *x0 = bfMatGetColView(X0, p);
    BfVec *x = bfMatMulVec(Vp, y);
    bfVecAddInplace(x, x0);
    bfMatSetCol(X, p, x);

    bfMatDelete(&Hp);
    bfMatDelete(&Vp);
    bfVecDelete(&s);
    bfVecDelete(&y);
    bfVecDelete(&x0);
    bfVecDelete(&x);
  }

  if (numIter != NULL)
    *numIter = j;

  BF_ERROR_END() {
    BF_DIE();
  }

  for (BfSize i = 0; i < maxNumIter + 1; ++i)
    if (V[i] != NULL)
      bfMatDelete(&V[i]);

  for (BfSize i = 0; i < maxNumIter; ++i)
    if (H[i] != NULL)
      bfMatDelete(&H[i]);

  for (BfSize i = 0; i < numRhs*maxNumIter; ++i)
    if (J[i] != NULL)
      bfMatDelete(&J[i]);

  bfVecDelete(&RNorm);
  bfMatDelete(&R);

  bfMemFree(V);
  bfMemFree(H);
  bfMemFree(J);

  bfMatDelete(&S);

  if (X0 != NULL)
    bfMatDelete(&X0);

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

  BF_ERROR_BEGIN();

  BfTimer timerFunc;
  bfTimerReset(&timerFunc);

  double *resid = NULL;
  double *V = NULL;
  a_int *select = NULL;
  double *workd = NULL;
  double *workl = NULL;
  BfReal *workev = NULL;

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

  BfLuCsrReal *M_lu = bfLuCsrRealNew();
  HANDLE_ERROR();

  bfLuCsrRealInit(M_lu, M);
  HANDLE_ERROR();

  resid = bfMemAlloc(N, sizeof(BfReal));
  if (resid == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  V = bfMemAlloc(N, ncv*sizeof(BfReal));
  if (V == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  select = bfMemAllocAndZero(ncv, sizeof(a_int));
  HANDLE_ERROR();

  a_int iparam[11] = {
    [0] = 1, /* compute exact shifts */
    [2] = 10*N, /* max number of iterations */
    // [3] = 1, /* only value allowed */
    [6] = 2, /* mode: A*x = lam*M*x */
  };

  a_int ipntr[11];

  workd = bfMemAllocAndZero(3*N, sizeof(BfReal));
  HANDLE_ERROR();

  workl = bfMemAlloc(lworkl, sizeof(BfReal));
  HANDLE_ERROR();

  workev = bfMemAlloc(3*ncv, sizeof(BfReal));
  HANDLE_ERROR();

  a_int ido = 0;
  a_int info = 0;

dnaupd:
  dnaupd_c(&ido, &bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
           workd, workl, lworkl, &info);
  if (ido == 1 || ido == -1) {
    BF_ASSERT(ipntr[0] > 0);
    BF_ASSERT(ipntr[1] > 0);

    BfVecReal x;
    bfVecRealInitView(&x, N, BF_DEFAULT_STRIDE, &workd[ipntr[0] - 1]);
    BfVec *tmp = bfMatMulVec(L, bfVecRealToVec(&x));
    BfVecReal *y = bfVecToVecReal(bfLuCsrRealSolveVec(M_lu, tmp));

    bfMemCopy(y->data, N, sizeof(BfReal), &workd[ipntr[1] - 1]);

    bfVecDelete(&tmp);
    bfVecRealDeinitAndDealloc(&y);

    goto dnaupd;
  } else if (ido == 2) {
    BF_ASSERT(ipntr[0] > 0);
    BF_ASSERT(ipntr[1] > 0);

    BfVecReal x;
    bfVecRealInitView(&x, N, BF_DEFAULT_STRIDE, &workd[ipntr[0] - 1]);

    BfVecReal *y = bfVecToVecReal(bfMatMulVec(M, bfVecRealToVec(&x)));

    bfMemCopy(y->data, N, sizeof(BfReal), &workd[ipntr[1] - 1]);

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

  BF_ERROR_END() {
    eigmax = NAN;
  }

  bfLuCsrRealDeinitAndDealloc(&M_lu);

  bfMemFree(resid);
  bfMemFree(V);
  bfMemFree(select);
  bfMemFree(workd);
  bfMemFree(workl);
  bfMemFree(workev);

  return eigmax;
}

void bfGetShiftedEigs(BfMat const *A, BfMat const *M, BfReal sigma, BfSize k,
                      BfMat **PhiTransposePtr, BfVecReal **LambdaPtr) {
  /* TODO: this is a work in progress! This does NOT work for any type
   * of BfMat yet. Just real ones... */

  BF_ERROR_BEGIN();

  if (!isfinite(sigma))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfTimer timerFunc;
  bfTimerReset(&timerFunc);

  double *resid = NULL;
  double *V = NULL;
  a_int *select = NULL;
  double *workd = NULL;
  double *workl = NULL;
  BfReal *workev = NULL;
  BfReal *dr = NULL;
  BfReal *di = NULL;
  BfReal *z = NULL;

  BfSize *J = NULL;
  BfVecReal *Lambda = NULL;

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

  BfTimer timer;

  bfTimerReset(&timer);

  BfMat *A_minus_sigma_M = bfMatCopy(M);
  HANDLE_ERROR();

  bfMatScale(A_minus_sigma_M, -sigma);
  HANDLE_ERROR();

  bfMatAddInplace(A_minus_sigma_M, A);
  HANDLE_ERROR();

  BfLuCsrReal *A_minus_sigma_M_lu = bfLuCsrRealNew();
  HANDLE_ERROR();

  bfLuCsrRealInit(A_minus_sigma_M_lu, A_minus_sigma_M);
  HANDLE_ERROR();

  bfLogInfo("bfGetShiftedEigs: created and factorized A - sigma*M [%.1fs]\n",
            bfTimerGetElapsedTimeInSeconds(&timer));

  resid = bfMemAlloc(N, sizeof(BfReal));
  if (resid == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  V = bfMemAlloc(N*ncv, sizeof(BfReal));
  if (V == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  select = bfMemAllocAndZero(ncv, sizeof(a_int));
  HANDLE_ERROR();

  a_int iparam[11] = {
    [0] = 1, /* compute exact shifts */
    [2] = 10*N, /* max number of iterations */
    [6] = 3, /* shift-invert mode */
  };

  a_int ipntr[11];

  workd = bfMemAlloc(3*N, sizeof(BfReal));
  if (workd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (a_int i = 0; i < 3*N; ++i)
    workd[i] = 0;

  workl = bfMemAlloc(lworkl, sizeof(BfReal));
  if (workl == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  workev = bfMemAlloc(3*ncv, sizeof(BfReal));
  if (workev == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  a_int ido = 0;
  a_int info = 0;

  BfSize numSolve = 0;
  BfSize numMvp = 0;

  BfReal totalSolveTime = 0;
  BfReal totalMvpTime = 0;

dnaupd:
  dnaupd_c(&ido, &bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
           workd, workl, lworkl, &info);
  if (ido == 1 || ido == -1) {
    BF_ASSERT(ipntr[0] > 0);
    BF_ASSERT(ipntr[1] > 0);

    bfTimerReset(&timer);

    BfVecReal x;
    bfVecRealInitView(&x, N, BF_DEFAULT_STRIDE, &workd[ipntr[0] - 1]);
    BfVec *tmp = bfMatMulVec(M, bfVecRealToVec(&x));
    BfVecReal *y = bfVecToVecReal(bfLuCsrRealSolveVec(A_minus_sigma_M_lu, tmp));

    totalSolveTime += bfTimerGetElapsedTimeInSeconds(&timer);

    ++numSolve;

    bfMemCopy(y->data, N, sizeof(BfReal), &workd[ipntr[1] - 1]);

    bfVecDelete(&tmp);
    bfVecRealDeinitAndDealloc(&y);

    goto dnaupd;
  } else if (ido == 2) {
    BF_ASSERT(ipntr[0] > 0);
    BF_ASSERT(ipntr[1] > 0);

    BfVecReal x;
    bfVecRealInitView(&x, N, BF_DEFAULT_STRIDE, &workd[ipntr[0] - 1]);

    bfTimerReset(&timer);

    BfVecReal *y = bfVecToVecReal(bfMatMulVec(M, bfVecRealToVec(&x)));

    totalMvpTime += bfTimerGetElapsedTimeInSeconds(&timer);

    ++numMvp;

    bfMemCopy(y->data, N, sizeof(BfReal), &workd[ipntr[1] - 1]);

    bfVecRealDeinitAndDealloc(&y);

    goto dnaupd;
  }

  bfLogInfo("bfGetShiftedEigs: finished iterating [%.1fs]\n",
            bfTimerGetElapsedTimeInSeconds(&timer));
  bfLogInfo("bfGetShiftedEigs: did %lu solves [%.1fs]\n", numSolve, totalSolveTime);
  bfLogInfo("bfGetShiftedEigs: did %lu MVPs [%.1fs]\n", numMvp, totalMvpTime);

  if (info < 0 || iparam[4] < nev)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  dr = bfMemAlloc(nev + 1, sizeof(BfReal));
  if (dr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  di = bfMemAlloc(nev + 1, sizeof(BfReal));
  if (di == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  z = bfMemAlloc(N*(nev + 1), sizeof(BfReal));
  if (z == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfSize ldz = N;

  BfReal sigmar = sigma;
  BfReal sigmai = 0.0;

  bfTimerReset(&timer);

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

  bfLogInfo("bfGetShiftedEigs: extracted eigs [%.1fs]\n",
            bfTimerGetElapsedTimeInSeconds(&timer));

  if (info != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* TODO: evals could obviously be complex but right now we're only
   * interested in real evals and haven't introduced any other
   * safeguards to deal with complex case yet. So just treat complex
   * eigenvalues as an error for now. */
  for (BfSize i = 0; i < k; ++i)
    if (fabs(di[i]) != 0)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  J = bfMemAlloc(k, sizeof(BfSize));
  if (J == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  bfRealArgsort(dr, k, J);

  if (LambdaPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (*LambdaPtr != NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  Lambda = bfVecRealNew();
  HANDLE_ERROR();

  bfVecRealInit(Lambda, k);
  HANDLE_ERROR();

  for (BfSize j = 0; j < k; ++j)
    *(Lambda->data + J[j]*Lambda->stride) = dr[j];

  if (PhiTransposePtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (*PhiTransposePtr != NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMatDenseReal *PhiTranspose = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInit(PhiTranspose, k, N);
  HANDLE_ERROR();

  for (BfSize j = 0; j < k; ++j) {
    BfVecReal *phi = bfVecRealNew();
    HANDLE_ERROR();

    bfVecRealInitView(phi, N, BF_DEFAULT_STRIDE, z + N*j);
    HANDLE_ERROR();

    bfMatDenseRealSetRow(bfMatDenseRealToMat(PhiTranspose), J[j], bfVecRealToVec(phi));

    bfVecRealDeinitAndDealloc(&phi);
  }

  *LambdaPtr = Lambda;
  *PhiTransposePtr = bfMatDenseRealToMat(PhiTranspose);

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMemFree(J);
  bfMemFree(z);
  bfMemFree(dr);
  bfMemFree(di);

  bfLuCsrRealDeinitAndDealloc(&A_minus_sigma_M_lu);
  bfMatDelete(&A_minus_sigma_M);

  bfMemFree(resid);
  bfMemFree(V);
  bfMemFree(select);
  bfMemFree(workd);
  bfMemFree(workl);
  bfMemFree(workev);

  bfLogInfo("bfGetShiftedEigs: done [%.1fs]\n",
            bfTimerGetElapsedTimeInSeconds(&timerFunc));
}

static void getEigenband_doubling(BfMat const *A, BfMat const *M, BfInterval const *interval,
                                  BfMat **PhiTransposePtr, BfVecReal **LambdaPtr) {
  BF_ERROR_BEGIN();

  BfMat *PhiTranspose = NULL;
  BfVecReal *Lambda = NULL;

  BfSize k = 8;

  BfReal const sigma = bfIntervalIsFinite(interval) ?
    bfIntervalGetMidpoint(interval) :
    bfIntervalGetFiniteEndpoint(interval);

get_shifted_eigs:
  bfLogInfo("bfGetEigenband: k = %lu\n", k);

  bfGetShiftedEigs(A, M, sigma, k, &PhiTranspose, &Lambda);
  HANDLE_ERROR();

  BfReal const *lam = Lambda->data;

  /* Keep doubling the number of eigenvalues until we cover `interval`
   * with `currentInterval`: */
  BfInterval currentInterval = {
    .endpoint = {lam[0], lam[k - 1]},
    .closed = {true, true}
  };
  if (!bfIntervalContainsInterval(&currentInterval, interval)) {
    k *= 2;
    bfMatDelete(&PhiTranspose);
    bfVecRealDeinitAndDealloc(&Lambda);
    goto get_shifted_eigs;
  }

  /* Find the first eigenpair in the band */
  BfSize j0 = 0;
  while (j0 < k && !bfIntervalContainsPoint(interval, lam[j0])) ++j0;

  /* Find the last eigenpair in the band */
  BfSize j1 = k;
  while (j1 > j0 && !bfIntervalContainsPoint(interval, lam[j1 - 1])) --j1;

  if (0 < j0 || j1 < k) {
    /* Prune unnecessary eigenvectors */
    BfMat *oldPhiTranspose = PhiTranspose;
    PhiTranspose = bfMatGetRowRangeCopy(oldPhiTranspose, j0, j1);
    HANDLE_ERROR();

    /* Prune unnecessary eigenvalues */
    BfVecReal *oldLambda = Lambda;
    Lambda = bfVecToVecReal(bfVecGetSubvecCopy(bfVecRealToVec(oldLambda), j0, j1));
    HANDLE_ERROR();

    /* Free old eigenvectors and values */
    bfMatDelete(&oldPhiTranspose);
    bfVecRealDeinitAndDealloc(&oldLambda);
  }

  *PhiTransposePtr = PhiTranspose;
  *LambdaPtr = Lambda;

  BF_ERROR_END() {
    bfMatDelete(&PhiTranspose);
    bfVecRealDeinitAndDealloc(&Lambda);

    *PhiTransposePtr = NULL;
    *LambdaPtr = NULL;
  }
}

static BfInterval getPairsCoveringInterval(BfMat const *A, BfMat const *M, BfInterval const *interval, BfRealArray *LamData, BfRealArray *PhiTransposeData) {
  BF_ERROR_BEGIN();

  BfMat *coverPhiTranspose = NULL;
  BfVecReal *coverLam = NULL;
  BfReal const sigma = bfIntervalIsFinite(interval) ?
    bfIntervalGetMidpoint(interval) : bfIntervalGetFiniteEndpoint(interval);
  BfSize const k = 8;
  bfGetShiftedEigs(A, M, sigma, k + 2, &coverPhiTranspose, &coverLam);
  HANDLE_ERROR();

  BfReal const *lam = coverLam->data;

  BfSize i0 = 0;
  while (i0 < k + 2 && !bfIntervalContainsPoint(interval, lam[i0])) ++i0;

  BfSize i1 = k + 2;
  while (i1 > 0 && !bfIntervalContainsPoint(interval, lam[i1 - 1])) --i1;

  BfInterval cover = {.closed = {false, false}};

  if (i0 == 0 && i1 == k + 2) {
    cover.endpoint[0] = (lam[0] + lam[1])/2;
    cover.endpoint[1] = (lam[k] + lam[k + 1])/2;
    i0 = 1;
    i1 = k + 1;
  } else if (i0 == 0 && i1 == 1) {
    cover = *interval;
  } else if (i0 == k + 1 && i1 == k + 2) {
    cover = *interval;
  } else if (i0 == 0 && i1 < k + 2) {
    cover.endpoint[0] = (lam[0] + lam[1])/2;
    cover.endpoint[1] = interval->endpoint[1];
    cover.closed[1] = true;
  } else if (0 < i0 && i1 == k + 2) {
    cover.endpoint[0] = interval->endpoint[0];
    cover.endpoint[1] = (lam[k] + lam[k + 1])/2;
    cover.closed[0] = true;
    i1 = k + 1;
  } else if (i1 < i0) {
    cover.endpoint[0] = BF_INFINITY;
    cover.endpoint[1] = -BF_INFINITY;
  } else if (0 < i0 && i1 < k + 2) {
    cover.endpoint[0] = interval->endpoint[0];
    cover.endpoint[1] = interval->endpoint[1];
    cover.closed[0] = true;
    cover.closed[1] = true;
  } else {
    BF_DIE();
  }

  for (BfSize i = i0; i < i1; ++i) {
    BfReal lam = bfVecRealGetElt(coverLam, i);
    BF_ASSERT(bfIntervalContainsPoint(interval, lam));

    bfRealArrayAppend(LamData, lam);
    HANDLE_ERROR();

    BfVecReal *phi = bfVecToVecReal(bfMatGetRowView(coverPhiTranspose, i));
    HANDLE_ERROR();

    BfRealArray *phiArray = bfVecRealGetArrayView(phi);
    HANDLE_ERROR();

    bfRealArrayExtend(PhiTransposeData, phiArray);
    HANDLE_ERROR();

    bfVecRealDeinitAndDealloc(&phi);
    bfRealArrayDeinitAndDealloc(&phiArray);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMatDelete(&coverPhiTranspose);
  bfVecRealDeinitAndDealloc(&coverLam);

  return cover;
}

static void getEigenband_covering(BfMat const *A, BfMat const *M, BfInterval const *interval,
                                  BfMat **PhiPtr, BfVecReal **LambdaPtr) {
  BF_ERROR_BEGIN();

  BfRealArray *evals = bfRealArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  BfRealArray *evecs = bfRealArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  BfDisjointIntervalList *intervals = bfDisjointIntervalListNewEmpty();
  HANDLE_ERROR();

  bfDisjointIntervalListAdd(intervals, interval);
  HANDLE_ERROR();

  while (!bfDisjointIntervalListIsEmpty(intervals)) {
    BfInterval const *nextInterval = bfDisjointIntervalListGetFirstPtrConst(intervals);

    BfInterval cover = getPairsCoveringInterval(A, M, nextInterval, evals, evecs);
    HANDLE_ERROR();

    bfDisjointIntervalListRemove(intervals, bfIntervalIsEmpty(&cover) ? nextInterval : &cover);
  }

  BfSize numRows = bfMatGetNumRows(A);

  BfSize numCols = bfRealArrayGetSize(evecs);
  BF_ASSERT(numCols % numRows == 0);
  numCols /= numRows;

  BF_ASSERT(numCols == bfRealArrayGetSize(evals));
  for (BfSize j = 0; j < numCols; ++j) {
    BfReal eval = bfRealArrayGetValue(evals, j);
    BF_ASSERT(bfIntervalContainsPoint(interval, eval));
  }

  BfMatDenseReal *Phi = bfMatDenseRealNew();
  HANDLE_ERROR();

  bfMatDenseRealInit(Phi, numRows, numCols);
  HANDLE_ERROR();

  for (BfSize j = 0; j < numCols; ++j) {
    BfVec *col = bfRealArrayGetSubvecView(evecs, j*numRows, (j + 1)*numRows);
    bfMatDenseRealSetCol(Phi, j, col);
    bfVecDelete(&col);
  }

  BfPerm *perm = bfRealArrayArgsort(evals);
  HANDLE_ERROR();

  BfVecReal *Lambda = bfVecRealNewFromRealArray(evals, BF_POLICY_STEAL);
  HANDLE_ERROR();

  bfRealArrayPermute(evals, perm);
  bfMatDenseRealPermuteCols(Phi, perm);

  *PhiPtr = bfMatDenseRealToMat(Phi);
  *LambdaPtr = Lambda;

  BF_ERROR_END() {
    BF_DIE();
  }

  bfDisjointIntervalListDeinitAndDealloc(&intervals);
}

void bfGetEigenband(BfMat const *A, BfMat const *M, BfInterval const *interval,
                    BfEigenbandMethod method, BfMat **PhiTransposePtr, BfVecReal **LambdaPtr) {
  BF_ERROR_BEGIN();

  if (PhiTransposePtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (*PhiTransposePtr != NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (LambdaPtr == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (*LambdaPtr != NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (method == BF_EIGENBAND_METHOD_DOUBLING) {
    getEigenband_doubling(A, M, interval, PhiTransposePtr, LambdaPtr);
    HANDLE_ERROR();
  }

  else if (method == BF_EIGENBAND_METHOD_COVERING) {
    getEigenband_covering(A, M, interval, PhiTransposePtr, LambdaPtr);
    HANDLE_ERROR();
  }

  else RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BF_ERROR_END() {
    BF_DIE();
  }
}

bool bfGetTruncatedSvd(BfMat const *mat, BfMat **UPtr, BfMatDiagReal **SPtr, BfMat **VTPtr,
                       BfTruncSpec const *truncSpec, BfBackend backend) {
  BF_ERROR_BEGIN();

  bool truncated = true;

  if (backend != BF_BACKEND_LAPACK)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  /* Cast to `MatDenseReal` if we can, otherwise we'll need to convert
   * and allocate later. */
  BfMatDenseReal const *matDenseReal = NULL;
  bool shouldDeleteMatDenseReal = false;
  if (bfMatGetType(mat) == BF_TYPE_MAT_DENSE_REAL) {
    matDenseReal = bfMatConstToMatDenseRealConst(mat);
  } else {
    matDenseReal = bfMatDenseRealNewFromMatrix(mat);
    HANDLE_ERROR();

    shouldDeleteMatDenseReal = true;
  }

  BfMatDenseReal *U = NULL;
  BfMatDiagReal *S = NULL;
  BfMatDenseReal *VT = NULL;
  bfMatDenseRealSvd(matDenseReal, &U, &S, &VT);
  HANDLE_ERROR();

  /* Find the number of terms in the truncated SVD: */
  BfSize k = bfTruncSpecGetNumTerms(truncSpec, S);

  /* Should actually truncate? If so, do it: */
  truncated = k < S->numElts;
  if (truncated) {
    /** Truncate U: */

    BfMatDenseReal *Uk = bfMatToMatDenseReal(bfMatDenseRealGetColRangeCopy(U, 0, k));
    HANDLE_ERROR();

    bfMatDenseRealDeinitAndDealloc(&U);

    U = Uk;

    /** Truncate S: */

    BfMatDiagReal *Sk = bfMatDiagRealNew();
    HANDLE_ERROR();

    bfMatDiagRealInit(Sk, k, k);
    HANDLE_ERROR();

    for (BfSize i = 0; i < k; ++i)
      Sk->data[i] = S->data[i];

    bfMatDiagRealDeinitAndDealloc(&S);

    S = Sk;

    /** Truncate VT: */

    BfMatDenseReal *VkT = bfMatToMatDenseReal(bfMatDenseRealGetRowRangeCopy(VT, 0, k));
    HANDLE_ERROR();

    bfMatDenseRealDeinitAndDealloc(&VT);

    VT = VkT;
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  if (shouldDeleteMatDenseReal)
    bfMatDenseRealDeinitAndDealloc((BfMatDenseReal **)&matDenseReal);

  *UPtr = bfMatDenseRealToMat(U);
  *SPtr = S;
  *VTPtr = bfMatDenseRealToMat(VT);

  return truncated;
}
