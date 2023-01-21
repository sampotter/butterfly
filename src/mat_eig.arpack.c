#include <bf/mat.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <arpack.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/lu.h>
#include <bf/vec_real.h>

BfReal bfMatGetEigMaxGen(BfMat const *L, BfMat const *M) {
  /* TODO: this is a work in progress! This does NOT work for any type
   * of BfMat yet. Just real ones... */

  BEGIN_ERROR_HANDLING();

  a_int const N = bfMatGetNumRows(L);
  char const which[] = "LM";
  a_int const nev = 1;
  a_int const ncv = N < 20 ? N : 20; /* TODO: how many is best? */
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

  a_int *select = malloc(ncv*sizeof(a_int));

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
    NULL, /* not reversed since rvec == 0 */
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

void bfMatGetEigBandGen(BfMat const *A, BfMat const *M, BfReal lam0, BfReal lam1,
                        BfMat **Phi, BfMat **Lam) {
  (void)A;
  (void)M;
  (void)lam0;
  (void)lam1;
  (void)Phi;
  (void)Lam;

  assert(false);
}
