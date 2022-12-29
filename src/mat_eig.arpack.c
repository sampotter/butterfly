#include <bf/mat.h>

#include <stdlib.h>
#include <string.h>

#include <arpack/arpack.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/vec_real.h>

BfReal bfMatGetEigMaxGenSym(BfMat const *A, BfMat const *M, BfMat **Lhandle) {
  /* TODO: this is a work in progress! This does NOT work for any type
   * of BfMat yet. Just real ones... */

  BEGIN_ERROR_HANDLING();

  a_int const N = bfMatGetNumRows(A);
  char const which[] = "LM";
  a_int const nev = 1;
  a_int const ncv = 20; /* TODO: how many is best? */
  char const bmat[] = "G";
  BfReal const tol = 0;
  a_int const ldv = N;
  a_int const lworkl = ncv*(ncv + 8);
  a_int const rvec = 0; /* only computing eigenvalues */
  char const howmny[] = "A";

  BfMat *L = NULL;
  if (Lhandle == NULL) {
    L = bfMatCholesky(M);
  } else {
    L = *Lhandle;
  }

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
    [4] = 1, /* only value allowed */
    [6] = 1, /* mode: A*x = lam*M*x */
  };

  a_int ipntr[11];

  double *workd = malloc(3*N*sizeof(BfReal));
  if (workd == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  double *workl = malloc(lworkl*sizeof(BfReal));
  if (workl == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  a_int ido = 0, info = 0;
  do {
    dsaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, &info);

    BfVecReal x;
    bfVecRealInitView(&x, N, 0, &workd[ipntr[0] - 1]);

    BfVec *tmp1 = bfMatForwardSolveVec(bfMatTrans(L), bfVecRealToVec(&x));
    BfVec *tmp2 = bfMatMulVec(A, tmp1);
    BfVecReal *y = bfVecToVecReal(bfMatBackwardSolveVec(L, tmp2));

    bfVecDelete(&tmp1);
    bfVecDelete(&tmp2);

    memcpy(&workd[ipntr[1] - 1], y->data, N*sizeof(BfReal));
  } while ((ido == 1) || (ido == -1));

  if (info < 0 || iparam[4] < nev)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfReal eigmax;

  dseupd_c(rvec, howmny, select, &eigmax,
           NULL, 0, 0, /* <-- these three parameters aren't referenced */
           /* dsaupd parameters: don't modify before calling dseupd */
           bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
           workd, workl, lworkl, &info);

  if (Lhandle == NULL)
    bfMatDelete(&L);

  if (info < 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  END_ERROR_HANDLING() {}

  free(resid);
  free(V);
  free(workd);
  free(workl);

  return eigmax;
}

void bfMatGetEigBandGenSym(BfMat const *A, BfMat const *M, BfReal lam0, BfReal lam1,
                           BfMat **Phi, BfMat **Lam) {
  (void)A;
  (void)M;
  (void)lam0;
  (void)lam1;
  (void)Phi;
  (void)Lam;

  assert(false);
}
