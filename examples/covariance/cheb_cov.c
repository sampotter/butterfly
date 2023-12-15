#include <bf/assert.h>
#include <bf/cheb.h>
#include <bf/const.h>
#include <bf/linalg.h>
#include <bf/mat_csr_real.h>
#include <bf/mat_product.h>
#include <bf/rand.h>
#include <bf/trimesh.h>
#include <bf/util.h>
#include <bf/vec_real.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static BfReal kappa = BF_NAN;
static BfReal nu = BF_NAN;

static BfReal gamma_(BfReal lambda) {
  return pow(fabs(kappa*kappa + lambda), -nu/4 - 1./2);
}

static BfVec *chebmul(BfCheb const *cheb, BfMat const *S, BfVec *w) {
  BfReal const *c = cheb->c;

  BF_ASSERT(cheb->a == 0);
  BfReal lamMax = cheb->b;

  BfSize n = bfMatGetNumRows(S);
  BfVec *x = bfVecRealToVec(bfVecRealNewWithValue(n, 0));

  BfVec *y2 = bfVecCopy(w);

  bfVecDaxpy(x, c[0], y2);

  BfVec *y1 = bfMatMulVec(S, w);
  bfVecDscal(y1, 2/lamMax);
  bfVecDaxpy(y1, -1, w);

  bfVecDaxpy(x, c[1], y1);

  for (BfSize k = 2; k < cheb->order; ++k) {
    BfVec *y = bfMatMulVec(S, y1);
    bfVecDscal(y, 4/lamMax);
    bfVecDaxpy(y, -2, y1);
    bfVecDaxpy(y, -1, y2);

    bfVecDaxpy(x, c[k], y);

    bfVecDelete(&y2);
    y2 = y1;
    y1 = y;
  }

  if (y1 != NULL) bfVecDelete(&y1);
  if (y2 != NULL) bfVecDelete(&y2);

  return x;
}

static BfVec *sample_z(BfCheb const *cheb, BfMat const *S, BfMat const *MLumpSqrtInv) {
  BfSize n = bfMatGetNumRows(S);
  BfVec *w = bfVecRealToVec(bfVecRealNewRandn(n));
  BfVec *x = chebmul(cheb, S, w);
  BfVec *z = bfMatMulVec(MLumpSqrtInv, x);
  bfVecDelete(&w);
  bfVecDelete(&x);
  return z;
}

static BfVec *get_c(BfCheb const *cheb, BfMat const *S, BfMat const *MLumpSqrtInv) {
  BfSize n = bfMatGetNumRows(S);

  BfVecReal *e = bfVecRealNewWithValue(n, 0);
  e->data[0] = 1;

  BfVec *tmp = bfMatMulVec(MLumpSqrtInv, bfVecRealToVec(e));
  BfVec *tmp1 = chebmul(cheb, S, tmp);
  bfVecDelete(&tmp);
  tmp = chebmul(cheb, S, tmp1);
  bfVecDelete(&tmp1);
  tmp1 = bfMatMulVec(MLumpSqrtInv, tmp);
  bfVecDelete(&tmp);

  bfVecRealDelete(&e);

  return tmp1;
}

static void getS(BfMat *L, BfMat *M, BfMat **SHandle, BfMat **MLumpSqrtInvHandle) {
  BfSize n = bfMatGetNumRows(L);

  BfMatCsrReal const *MCsr = bfMatConstToMatCsrRealConst(M);

  BfMatDiagReal *MLumpSqrtInv = bfMatDiagRealNew();
  bfMatDiagRealInit(MLumpSqrtInv, n, n);

  for (BfSize i = 0; i < n; ++i) {
    MLumpSqrtInv->data[i] = 0;
    for (BfSize j = MCsr->rowptr[i]; j < MCsr->rowptr[i + 1]; ++j)
      MLumpSqrtInv->data[i] += MCsr->data[j];
    MLumpSqrtInv->data[i] = 1/sqrt(MLumpSqrtInv->data[i]);
  }

  BfMatProduct *S = bfMatProductNew();
  bfMatProductInit(S);

  // TODO: should add a policy argument to PostMultiply...
  bfMatProductPostMultiply(S, bfMatGetView(bfMatDiagRealToMat(MLumpSqrtInv)));
  bfMatProductPostMultiply(S, bfMatGetView(L));
  bfMatProductPostMultiply(S, bfMatGetView(bfMatDiagRealToMat(MLumpSqrtInv)));

  *SHandle = bfMatProductToMat(S);
  *MLumpSqrtInvHandle = bfMatDiagRealToMat(MLumpSqrtInv);
}

int main(int argc, char const *argv[]) {
  if (argc != 6) {
    printf("usage: %s mesh.obj p kappa nu num_samples\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  BfSize p = atoi(argv[2]);
  kappa = atof(argv[3]);
  nu = atof(argv[4]);
  BfSize numSamples = atoi(argv[5]);

  bfSeed(0); // must seed before using PRNG

  BfTrimesh *trimesh = bfTrimeshNewFromObjFile(argv[1]);
  printf("triangle mesh with %lu verts\n", bfTrimeshGetNumVerts(trimesh));

  bfToc();

  BfMat *L = NULL, *M = NULL;
  bfTrimeshGetLboFemDiscretization(trimesh, &L, &M);
  printf("set up FEM discretization [%0.1fs]\n", bfToc());

  BfReal lamMax = bfGetMaxEigenvalue(L, M);
  printf("maximum eigenvalue: lambda = %g [%0.1fs]\n", lamMax, bfToc());

  BfCheb gammaCheb;
  bfChebInitWithDegree(&gammaCheb, p, 0, lamMax);
  bfChebInterp(&gammaCheb, gamma_, NULL);

  /** Evaluate gamma and our Chebyshev approximation of gamma on a
   ** grid, and write the grid, the true gamma values, and the
   ** polynomial values to disk for plotting. */

  BfSize N = 1000;
  FILE *fp = NULL;

  char filename[50];
  sprintf(filename, "lambda_p%i.bin", (int)p);
  fp = fopen(filename, "w");
  for (BfSize i = 0; i <= N; ++i) {
    BfReal lam = (i*lamMax)/N;
    fwrite(&lam, sizeof(BfReal), 1, fp);
  }
  fclose(fp);

  sprintf(filename, "gamma_p%i.bin", (int)p);
  fp = fopen(filename, "w");
  for (BfSize i = 0; i <= N; ++i) {
    BfReal y = gamma_((i*lamMax)/N);
    fwrite(&y, sizeof(BfReal), 1, fp);
  }
  fclose(fp);

  sprintf(filename, "gamma_cheb_p%i.bin", (int)p);
  fp = fopen(filename, "w");
  for (BfSize i = 0; i <= N; ++i) {
    BfReal y = bfChebEval(&gammaCheb, (i*lamMax)/N);
    fwrite(&y, sizeof(BfReal), 1, fp);
  }
  fclose(fp);

  /** Get S = sqrt(M)^-T*L*sqrt(M)^-1 (although the "T" does nothing
   ** here since M is lumped hence diagonal). */

  BfMat *S = NULL, *MLumpSqrtInv = NULL;
  getS(L, M, &S, &MLumpSqrtInv);

  /** Sample z once and write it out to disk for plotting. */

  BfVec *z = sample_z(&gammaCheb, S, MLumpSqrtInv);
  sprintf(filename, "z_cheb_p%i.bin", (int)p);
  bfVecSave(z, filename);
  bfVecDelete(&z);

  /** Time how long it takes to sample z numSamples times. */

  bfToc();
  for (BfSize _ = 0; _ < numSamples; ++_) {
    z = sample_z(&gammaCheb, S, MLumpSqrtInv);
    bfVecDelete(&z);
  }
  printf("drew %lu samples [%0.1fs]\n", numSamples, bfToc());

  /** Evaluate the covariance function with respect to a fixed point
   ** on the mesh. */

  BfVec *c = get_c(&gammaCheb, S, MLumpSqrtInv);
  sprintf(filename, "c_cheb_p%i.bin", (int)p);
  bfVecSave(c, filename);

  /** Clean up: */

  bfTrimeshDeinitAndDealloc(&trimesh);
  bfMatDelete(&L);
  bfMatDelete(&M);
  bfChebDeinit(&gammaCheb);
  bfVecDelete(&c);
  bfMatDelete(&S);
  bfMatDelete(&MLumpSqrtInv);
}
