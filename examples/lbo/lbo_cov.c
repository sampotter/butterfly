#include <bf/assert.h>
#include <bf/const.h>
#include <bf/fac_span.h>
#include <bf/fac_streamer.h>
#include <bf/interval_tree.h>
#include <bf/lbo.h>
#include <bf/linalg.h>
#include <bf/logging.h>
#include <bf/mat_dense_real.h>
#include <bf/octree.h>
#include <bf/rand.h>
#include <bf/util.h>
#include <bf/vec_real.h>

#include <math.h>
#include <stdlib.h>

#include "lbo.h"

static BfReal kappa = BF_NAN;
static BfReal nu = BF_NAN;

static BfReal gamma_(BfReal lambda) {
  return pow(fabs(kappa*kappa + lambda), -nu/4 - 1./2);
}

static BfReal gammaFromFreq(BfReal omega) {
  return gamma_(pow(omega, 2));
}

static BfVec *sample_z(BfMat const *Phi, BfMat const *GammaLam, BfPerm const *rowPerm) {
  BfSize n = bfMatGetNumCols(GammaLam);
  BfVec *w = bfVecRealToVec(bfVecRealNewRandn(n));
  BfVec *x = bfMatMulVec(GammaLam, w);
  BfVec *z = bfMatMulVec(Phi, x);
  bfVecPermute(z, rowPerm);
  bfVecDelete(&w);
  bfVecDelete(&x);
  return z;
}

static BfVec *get_c(BfMat const *Phi, BfMat const *GammaLam, BfMat const *M, BfPerm const *rowPerm, BfPerm const *revRowPerm) {
  BfSize n = bfMatGetNumCols(M);
  BfVec *e = bfVecRealToVec(bfVecRealNewStdBasis(n, 0));
  BfVec *tmp1 = bfMatMulVec(M, e);
  bfVecPermute(tmp1, revRowPerm);
  BfVec *tmp2 = bfMatRmulVec(Phi, tmp1);
  tmp1 = bfMatMulVec(GammaLam, tmp2);
  bfVecDelete(&tmp2);
  tmp2 = bfMatMulVec(GammaLam, tmp1);
  bfVecDelete(&tmp1);
  BfVec *z = bfMatMulVec(Phi, tmp2);
  bfVecPermute(z, rowPerm);
  bfVecDelete(&tmp2);
  return z;
}

int main(int argc, char const *argv[]) {
  if (argc < 5) {
    printf("usage: %s mesh.obj kappa nu num_samples [tol] [p] [rowTreeOffset] [freqTreeDepth]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  bfSeed(0);
  bfSetLogLevel(BF_LOG_LEVEL_INFO);

  char const *objPath = argv[1];
  kappa = atof(argv[2]);
  nu = atof(argv[3]);
  BfSize numSamples = atoi(argv[4]);

  BfReal tol = argc > 5 ? strtod(argv[5], NULL) : 1e-3;
  BfReal p = argc > 6 ? strtod(argv[6], NULL) : 0.0625;
  BfSize rowTreeOffset = argc > 7 ? strtoull(argv[7], NULL, 10) : 0;
  BfSize freqTreeDepth = argc > 8 ? strtoull(argv[8], NULL, 10) : BF_SIZE_BAD_VALUE;

  BfTrimesh trimesh;
  bfTrimeshInitFromObjFile(&trimesh, objPath);

  BfSize numVerts = bfTrimeshGetNumVerts(&trimesh);
  BfSize numEigs = (BfSize)(p*numVerts);

  printf("triangle mesh with %lu verts\n", numVerts);
  printf("streaming %lu eigenpairs\n", numEigs);

  BfOctree octree;
  bfOctreeInit(&octree, trimesh.verts, NULL, /* maxLeafSize: */ 1);

  BfTree *rowTree = bfOctreeToTree(&octree);
  BfSize rowTreeMaxDepth = bfTreeGetMaxDepth(rowTree);
  printf("row tree with depth %lu\n", rowTreeMaxDepth);

  BfPerm const *rowPerm = bfTreeGetPermConst(rowTree);
  BfPerm *revRowPerm = bfPermGetReversePerm(rowPerm);

  if (freqTreeDepth == BF_SIZE_BAD_VALUE)
    freqTreeDepth = rowTreeMaxDepth - 3;

  BfMat *L, *M;
  bfLboGetFemDiscretization(&trimesh, &L, &M);
  printf("set up FEM discretization [%0.1fs]\n", bfToc());

  bfMatCsrRealDump(bfMatToMatCsrReal(L), "L_rowptr.bin", "L_colind.bin", "L_data.bin");
  bfMatCsrRealDump(bfMatToMatCsrReal(M), "M_rowptr.bin", "M_colind.bin", "M_data.bin");

  BfReal lamMax = bfGetMaxEigenvalue(L, M);
  printf("maximum eigenvalue: lambda = %g [%0.1fs]\n", lamMax, bfToc());

  BfIntervalTree *freqTree = bfIntervalTreeNew();
  bfIntervalTreeInitEmpty(freqTree, 0, sqrt(lamMax), 2, freqTreeDepth);

  BfPoints1 *freqs = bfPoints1New();
  bfPoints1InitEmpty(freqs, BF_ARRAY_DEFAULT_CAPACITY);

  BfFacSpec spec = {
    .rowTree = rowTree,
    .colTree = bfIntervalTreeToTree(freqTree),
    .rowTreeInitDepth = rowTreeOffset,
    .colTreeInitDepth = freqTreeDepth, // TODO: this is unused!
    .tol = tol,
    .minNumRows = 20,
    .minNumCols = 20,
  };

  BfFacStreamer *facStreamer = bfFacStreamerNew();
  bfFacStreamerInit(facStreamer, &spec);

  while (!bfFacStreamerIsDone(facStreamer)) {
    feedFacStreamerNextEigenband(facStreamer, freqs, L, M);
    if (freqs->size >= numEigs) break;
  }

  printf("finished streaming BF (actually factorized %lu eigenpairs)\n", freqs->size);

  bfPoints1Save(freqs, "freqs.bin");

  BfPoints1 *gammaLam = bfPoints1Copy(freqs);
  bfPoints1Map(gammaLam, gammaFromFreq);

  BfFacSpan *facSpan = bfFacStreamerGetFacSpan(facStreamer);
  BfMat *Phi = bfFacSpanGetMat(facSpan, BF_POLICY_VIEW);
  BfMat *GammaLam = bfMatDiagRealToMat(
    bfMatDiagRealNewFromData(gammaLam->size, gammaLam->size, gammaLam->data));

  /** Sample z once and write it out to disk for plotting. */

  BfVec *z = sample_z(Phi, GammaLam, rowPerm);
  bfVecSave(z, "z_lbo.bin");
  bfVecDelete(&z);

  /** Time how long it takes to sample z numSamples times. */

  bfToc();
  for (BfSize _ = 0; _ < numSamples; ++_) {
    z = sample_z(Phi, GammaLam, rowPerm);
    bfVecDelete(&z);
  }
  printf("drew %lu samples [%0.1fs]\n", numSamples, bfToc());

  /** Evaluate the covariance function with respect to a fixed point
   ** on the mesh. */

  BfVec *c = get_c(Phi, GammaLam, M, rowPerm, revRowPerm);
  bfVecSave(c, "c_lbo.bin");

  /* Extract dense Phi: */

  BfSize m = bfMatGetNumRows(Phi);
  BfSize n = bfMatGetNumCols(Phi);
  BfMat *PhiDense = bfMatDenseRealToMat(bfMatDenseRealNewWithValue(m, n, BF_NAN));
  for (BfSize j = 0; j < n; ++j) {
    /* Get jth standard basis vector: */
    BfVec *e = bfVecRealToVec(bfVecRealNewStdBasis(n, j));
    BfVec *phi = bfMatMulVec(Phi, e);
    bfVecPermute(phi, rowPerm);
    bfMatSetCol(PhiDense, j, phi);
    bfVecDelete(&phi);
    bfVecDelete(&e);
  }
  bfMatSave(PhiDense, "PhiDense.bin");

  /* Clean up */
  bfTreeDelete(&rowTree);
  bfOctreeDeinit(&octree);
  bfMatDelete(&M);
  bfMatDelete(&L);
  bfTrimeshDeinit(&trimesh);
}
