#include <bf/linalg.h>

#include <math.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/lu_csr_real.h>
#include <bf/mat_dense_real.h>
#include <bf/mem.h>
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

  V = bfMemAlloc((maxNumIter + 1), sizeof(BfMatDenseComplex *));
  if (V == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  H = bfMemAlloc(maxNumIter, sizeof(BfMatDenseComplex *));
  if (H == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  J = bfMemAlloc(maxNumIter, numRhs*sizeof(BfMat *));
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
    }

    /* Extract the version of V for this RHS */
    BfMat *Vp = bfMatEmptyLike(V[0], n, j);
    for (BfSize i = 0; i < j; ++i) {
      BfVec *v = bfMatGetColView(V[i], p);
      bfMatSetCol(Vp, i, v);
    }

    /* Solve for the coefficients representing x - x0 in the V-basis */
    BfVec *s = bfMatGetColRangeView(S, 0, j, p);
    BfVec *y = bfMatBackwardSolveVec(Hp, s);
    BfVec *x0 = bfMatGetColView(X0, p);
    BfVec *x = bfMatMulVec(Vp, y);
    bfVecAddInplace(x, x0);
    bfMatSetCol(X, p, x);
  }

  END_ERROR_HANDLING() {}

  if (numIter != NULL)
    *numIter = j;

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

  BfLuCsrReal *M_lu = bfLuCsrRealNew();
  HANDLE_ERROR();

  bfLuCsrRealInit(M_lu, M);
  HANDLE_ERROR();

  double *resid = bfMemAlloc(N, sizeof(BfReal));
  if (resid == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  double *V = bfMemAlloc(N, ncv*sizeof(BfReal));
  if (V == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  a_int *select = bfMemAllocAndZero(ncv, sizeof(a_int));

  a_int iparam[11] = {
    [0] = 1, /* compute exact shifts */
    [2] = 10*N, /* max number of iterations */
    // [3] = 1, /* only value allowed */
    [6] = 2, /* mode: A*x = lam*M*x */
  };

  a_int ipntr[11];

  double *workd = bfMemAlloc(3*N, sizeof(BfReal));
  if (workd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (a_int i = 0; i < 3*N; ++i)
    workd[i] = 0;

  double *workl = bfMemAlloc(lworkl, sizeof(BfReal));
  if (workl == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *workev = bfMemAlloc(3*ncv, sizeof(BfReal));
  if (workev == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

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

  END_ERROR_HANDLING() {
    eigmax = NAN;
  }

  bfLuCsrRealDeinitAndDealloc(&M_lu);

  bfMemFree(resid);
  bfMemFree(V);
  bfMemFree(workd);
  bfMemFree(workl);
  bfMemFree(workev);

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

  BfLuCsrReal *A_minus_sigma_M_lu = bfLuCsrRealNew();
  HANDLE_ERROR();

  bfLuCsrRealInit(A_minus_sigma_M_lu, A_minus_sigma_M);
  HANDLE_ERROR();

  double *resid = bfMemAlloc(N, sizeof(BfReal));
  if (resid == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  double *V = bfMemAlloc(N*ncv, sizeof(BfReal));
  if (V == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  a_int *select = bfMemAllocAndZero(ncv, sizeof(a_int));

  a_int iparam[11] = {
    [0] = 1, /* compute exact shifts */
    [2] = 10*N, /* max number of iterations */
    [6] = 3, /* shift-invert mode */
  };

  a_int ipntr[11];

  double *workd = bfMemAlloc(3*N, sizeof(BfReal));
  if (workd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (a_int i = 0; i < 3*N; ++i)
    workd[i] = 0;

  double *workl = bfMemAlloc(lworkl, sizeof(BfReal));
  if (workl == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *workev = bfMemAlloc(3*ncv, sizeof(BfReal));
  if (workev == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

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
    BfVec *tmp = bfMatMulVec(M, bfVecRealToVec(&x));
    BfVecReal *y = bfVecToVecReal(bfLuCsrRealSolveVec(A_minus_sigma_M_lu, tmp));

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

  BfReal *dr = bfMemAlloc(nev + 1, sizeof(BfReal));
  if (dr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *di = bfMemAlloc(nev + 1, sizeof(BfReal));
  if (di == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *z = bfMemAlloc(N*(nev + 1), sizeof(BfReal));
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

  BfSize *J = bfMemAlloc(k, sizeof(BfSize));
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

  bfMemFree(J);
  bfMemFree(z);
  bfMemFree(dr);
  bfMemFree(di);

  bfLuCsrRealDeinitAndDealloc(&A_minus_sigma_M_lu);
  bfMatDelete(&A_minus_sigma_M);

  bfMemFree(resid);
  bfMemFree(V);
  bfMemFree(workd);
  bfMemFree(workl);
  bfMemFree(workev);
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
    BF_ASSERT(lam0 <= lam[j0]);
  }
  BF_ASSERT(j0 < k);

  /* Find the last eigenpair in the band */
  BfSize j1 = k;
  if (useRight) {
    while (j1 > j0 && lam1 <= lam[j1 - 1])
      --j1;
    BF_ASSERT(lam[j1 - 1] < lam1);
    BF_ASSERT(j1 == k || lam1 <= lam[j1]);
  }
  BF_ASSERT(j0 <= j1 && j1 <= k);

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

bool bfGetTruncatedSvd(BfMat const *mat, BfMat **UPtr, BfMatDiagReal **SPtr, BfMat **VTPtr,
                       BfTruncSpec const *truncSpec, BfBackend backend) {
  BEGIN_ERROR_HANDLING();

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

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  if (shouldDeleteMatDenseReal)
    bfMatDenseRealDeinitAndDealloc((BfMatDenseReal **)&matDenseReal);

  *UPtr = bfMatDenseRealToMat(U);
  *SPtr = S;
  *VTPtr = bfMatDenseRealToMat(VT);

  return truncated;
}
