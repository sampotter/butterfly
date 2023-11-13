#include "cmocka_setup.h"

#include <stdio.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/interval.h>
#include <bf/lbo.h>
#include <bf/linalg.h>
#include <bf/mat_dense_real.h>
#include <bf/trimesh.h>
#include <bf/vec_real.h>

static BfEigenbandMethod const eigenbandMethods[2] = {
  BF_EIGENBAND_METHOD_DOUBLING,
  BF_EIGENBAND_METHOD_COVERING
};

void test_bfGetEigenband(void **state) {
  (void)state;

  BfMat *PhiTrue = bfMatDenseRealToMat(bfMatDenseRealNewFromCsv("./sphere_Phi.txt"));
  BfVecReal *LamTrue = bfVecRealNewFromCsv("./sphere_Lam.txt");

  BfTrimesh *trimesh = bfTrimeshNewFromObjFile("./sphere.obj");

  BfMat *L = NULL;
  BfMat *M = NULL;
  bfLboGetFemDiscretization(trimesh, &L, &M);

  BfReal lamMin = 50;
  BfReal lamMax = 100;
  BfInterval interval = {.endpoint = {lamMin, lamMax}, .closed = {true, true}};

  BfMat *PhiTranspose = NULL;
  BfVecReal *Lam = NULL;

  for (BfSize i = 0; i < 2; ++i) {
    bfGetEigenband(L, M, &interval, eigenbandMethods[i], &PhiTranspose, &Lam);

    BfSize numEigs = bfVecRealToVec(Lam)->size;
    BF_ASSERT(bfMatGetNumRows(PhiTranspose) == numEigs);

    BF_ASSERT(bfVecRealDistMax(LamTrue, bfVecRealToVec(Lam))/bfVecRealNormMax(LamTrue) <= 1e-15);

    for (BfSize j = 0; j < numEigs; ++j) {
      BfVec *phi = bfMatGetRowView(PhiTranspose, j);
      BfVec *phiTrue = bfMatGetColView(PhiTrue, j);

      /* When we solve the GEP using the two different methods, we can
       * have at most a sign difference between corresponding
       * eigenvectors since the GEP is symmetric. To compute the "max
       * norm" between the two eigenvectors matrices, we compute the max
       * norm distance twice: once with a sign flip, and once
       * without. We use the minimum of the two distances for the
       * comparison. */
      BfReal dists[2];
      dists[0] = bfVecDistMax(phi, phiTrue);
      bfVecDscal(phi, -1);
      dists[1] = bfVecDistMax(phi, phiTrue);
      BfReal dist = fmin(dists[0], dists[1]);
      BF_ASSERT(dist <= 1.12e-11);

      bfVecDelete(&phi);
      bfVecDelete(&phiTrue);
    }

    bfMatDelete(&PhiTranspose);
    bfVecRealDeinitAndDealloc(&Lam);
  }

  bfMatDelete(&M);
  bfMatDelete(&L);
  bfTrimeshDeinitAndDealloc(&trimesh);
  bfVecRealDeinitAndDealloc(&LamTrue);
  bfMatDelete(&PhiTrue);
}

int main(int argc, char const *argv[]) {
  printf("%d\n", argc);
  for (int i = 0; i < argc; ++i)
    printf("%s\n", argv[i]);

  struct CMUnitTest const tests[] = {
    cmocka_unit_test(test_bfGetEigenband)
  };

  return cmocka_run_group_tests(tests, NULL, NULL);
}
