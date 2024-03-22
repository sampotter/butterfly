#include <bf/assert.h>
#include <bf/const.h>
#include <bf/fac_span.h>
#include <bf/fac_streamer.h>
#include <bf/interval_tree.h>
#include <bf/lbo.h>
#include <bf/linalg.h>
#include <bf/logging.h>
#include <bf/mat_csr_real.h>
#include <bf/mat_dense_real.h>
#include <bf/octree.h>
#include <bf/rand.h>
#include <bf/trimesh.h>
#include <bf/util.h>
#include <bf/vec_real.h>

#include <math.h>
#include <stdlib.h>

static BfReal kappa = BF_NAN;
static BfReal nu = BF_NAN;

static BfReal gamma_(BfReal lambda) {
  if (nu == 0) {
    // squared exponential spectral density function
    return exp(-kappa*lambda*lambda);
  } else {
    // Matern spectral density function, normalized so g(0) = 1
    return pow(fabs(1 + kappa*kappa*lambda), -nu/4 - 1./2);
  }
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

static BfVec *cov_matvec(BfVec *v, BfMat const *Phi, BfMat const *GammaLam, BfPerm const *rowPerm, BfPerm const *revRowPerm) {
  BfVec *tmp1 = v;
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
    printf("usage: %s L_rowptr.bin L_colind.bin L_data.bin M_rowptr.bin M_colind.bin M_data.bin points.bin kappa nu num_samples [tol] [fraction] [rowTreeOffset] [freqTreeDepth]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  bfSeed(0);
  bfSetLogLevel(BF_LOG_LEVEL_WARN);

  char const *LRowptrPath = argv[1];
  char const *LColindPath = argv[2];
  char const *LDataPath = argv[3];
  char const *MRowptrPath = argv[4];
  char const *MColindPath = argv[5];
  char const *MDataPath = argv[6];
  char const *PointsPath = argv[7];
  kappa = atof(argv[8]);
  nu = atof(argv[9]);
  BfSize numSamples = atoi(argv[10]);

  BfReal tol = argc > 11 ? strtod(argv[11], NULL) : 1e-3;
  BfReal fraction = argc > 12 ? strtod(argv[12], NULL) : 1.0;
  BfSize rowTreeOffset = argc > 13 ? strtoull(argv[13], NULL, 10) : 0;
  BfSize freqTreeDepth = argc > 14 ? strtoull(argv[14], NULL, 10) : BF_SIZE_BAD_VALUE;

  BfPoints3 *points = bfPoints3NewFromBinaryFile(PointsPath);

  BfSize numPoints = bfPoints3GetSize(points);
  BfSize numEigs = fraction*numPoints;

  BfOctree octree;
  bfOctreeInit(&octree, points, NULL, /* maxLeafSize: */ 1);

  BfTree *rowTree = bfOctreeToTree(&octree);
  BfSize rowTreeMaxDepth = bfTreeGetMaxDepth(rowTree);
  printf("row tree with depth %lu\n", rowTreeMaxDepth);

  BfPerm const *rowPerm = bfTreeGetPermConst(rowTree);
  BfPerm *revRowPerm = bfPermGetReversePerm(rowPerm);

  if (freqTreeDepth == BF_SIZE_BAD_VALUE)
    freqTreeDepth = rowTreeMaxDepth - 3;

  BfMat *L = bfMatCsrRealToMat(
    bfMatCsrRealNewFromBinaryFiles(LRowptrPath, LColindPath, LDataPath));
  BfMat *M = bfMatCsrRealToMat(
    bfMatCsrRealNewFromBinaryFiles(MRowptrPath, MColindPath, MDataPath));
  printf("loaded stiffness and mass matrices [%0.1fs]\n", bfToc());

  BfReal lamMax = bfGetMaxEigenvalue(L, M);
  printf("maximum eigenvalue: lambda = %g [%0.1fs]\n", lamMax, bfToc());

  BfIntervalTree *freqTree = bfIntervalTreeNew();
  bfIntervalTreeInitEmpty(freqTree, 0, sqrt(lamMax), 2, freqTreeDepth);
  BfSize freqTreeMaxDepth = bfTreeGetMaxDepth(bfIntervalTreeToTree(freqTree));
  printf("freq tree with depth %lu\n", freqTreeMaxDepth);

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

  BfSize nfit = 100; // number of eigenvalues to fit for extrapolation
  BfReal err_est = 1.0;
  while (!bfFacStreamerIsDone(facStreamer) && err_est > tol) {
    bfLboFeedFacStreamerNextEigenband(facStreamer, freqs, L, M);
    if (freqs->size >= numEigs) break;

    // Don't try to extrapolate if we don't have enough frequencies:
    if (freqs->size <= nfit) continue;

    BfReal numer = 0;
    BfReal denom = 0;
    BfReal xhat  = 0;
    BfReal yhat  = 0;
    for (BfSize i = freqs->size - nfit; i < freqs->size; ++i) {
      BfReal lam = pow(freqs->data[i], 2);
      xhat  += i;
      yhat  += lam;
    }
    xhat /= nfit;
    yhat /= nfit;
    for (BfSize i = freqs->size - nfit; i < freqs->size; ++i) {
      BfReal lam = pow(freqs->data[i], 2);
      numer += (i - xhat)*(lam - yhat);
      denom += pow(i - xhat, 2);
    }
    BfReal m = numer/denom;
    BfReal b = yhat - m*xhat;

    numer = 0;
    for (BfSize i = freqs->size; i < numPoints; ++i) {
      numer += pow(gamma_(m*i + b), 4);
    }
    denom = numer;
    for (BfSize i = 0; i < freqs->size; ++i) {
      BfReal lam = pow(freqs->data[i], 2);
      denom += pow(gamma_(lam), 4);
    }
    err_est = sqrt(numer)/sqrt(denom);
    printf("truncation error estimate after %lu eigenpairs is %.2e\n", freqs->size, err_est);
  }
  double precomp_time = bfToc();

  printf("finished streaming BF (actually factorized %lu eigenpairs) [%0.1fs]\n", freqs->size, precomp_time);

  BfFacSpan *facSpan = bfFacStreamerGetFacSpan(facStreamer);
  BfMat *Phi = bfFacSpanGetMat(facSpan, BF_POLICY_VIEW);

  BfReal numBytesCompressed = bfMatNumBytes(Phi);
  BfReal numBytesUncompressed = sizeof(BfReal)*numPoints*freqs->size;
  BfReal numBytesUntruncated  = sizeof(BfReal)*numPoints*numPoints;

  printf("  compressed   size: %.1f MB\n", numBytesCompressed/pow(1024, 2));
  printf("  uncompressed size: %.1f MB\n", numBytesUncompressed/pow(1024, 2));
  printf("  untruncated  size: %.1f MB\n", numBytesUntruncated/pow(1024, 2));
  printf("  compression  rate: %.1f\n", numBytesUncompressed/numBytesCompressed);

  char filename[50];
  sprintf(filename, "freqs_tol%.1e.bin", tol);
  bfPoints1Save(freqs, filename);

  BfPoints1 *gammaLam = bfPoints1Copy(freqs);
  bfPoints1Map(gammaLam, gammaFromFreq);

  BfMat *GammaLam = bfMatDiagRealToMat(
    bfMatDiagRealNewFromData(gammaLam->size, gammaLam->size, gammaLam->data));

  /** Sample z once and write it out to disk for plotting. */

  BfVec *z = sample_z(Phi, GammaLam, rowPerm);
  sprintf(filename, "z_lbo_tol%.0e_kappa%.1e_nu%.1e.bin", tol, kappa, nu);
  bfVecSave(z, filename);
  bfVecDelete(&z);

  /** Time how long it takes to sample z numSamples times. */

  bfToc();
  for (BfSize _ = 0; _ < numSamples; ++_) {
    z = sample_z(Phi, GammaLam, rowPerm);
    bfVecDelete(&z);
  }
  double sampling_time = bfToc();
  printf("drew %lu samples [%0.1fs]\n", numSamples, sampling_time);

  // save factorization time and memory sizes to file
  char line[100];
  FILE *fptr;
  sprintf(filename, "performance_kappa%.1e_nu%.1e.txt", kappa, nu);
  sprintf(
    line,
    "%.1e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\n",
    tol, precomp_time, sampling_time/numSamples,
    numBytesCompressed/pow(1024, 2),
    numBytesUncompressed/pow(1024, 2),
    numBytesUntruncated/pow(1024, 2)
    );
  fptr = fopen(filename, "a");
  fprintf(fptr, line);
  fclose(fptr);

  /** Evaluate the covariance function with respect to a fixed point
   ** on the mesh. */

  BfVec *e = bfVecRealToVec(bfVecRealNewStdBasis(numPoints, 0));
  BfVec *c = cov_matvec(e, Phi, GammaLam, rowPerm, revRowPerm);
  sprintf(filename, "c_lbo_tol%.0e_kappa%.1e_nu%.1e.bin", tol, kappa, nu);
  bfVecSave(c, filename);

  // /* Extract dense Phi: */

  // BfSize m = bfMatGetNumRows(Phi);
  // BfSize n = bfMatGetNumCols(Phi);
  // BfMat *PhiDense = bfMatDenseRealToMat(bfMatDenseRealNewWithValue(m, n, BF_NAN));
  // for (BfSize j = 0; j < n; ++j) {
  //   /* Get jth standard basis vector: */
  //   BfVec *e = bfVecRealToVec(bfVecRealNewStdBasis(n, j));
  //   BfVec *phi = bfMatMulVec(Phi, e);
  //   bfVecPermute(phi, rowPerm);
  //   bfMatSetCol(PhiDense, j, phi);
  //   bfVecDelete(&phi);
  //   bfVecDelete(&e);
  // }
  // bfMatSave(PhiDense, "PhiDense.bin");

  /* Compute and store covariance matrix vector products */

  bfSeed(0);
  BfSize s = numSamples;
  BfMatDenseReal *matvecs = bfMatDenseRealNewZeros(numPoints, s);

  printf("computing %lu matvecs with covariance\n", s);
  for (BfSize j = 0; j < s; ++j) {
      BfVec *x = bfVecRealToVec(bfVecRealNewRandn(numPoints));

      // apply covariance matrix
      BfVec *tmp1 = cov_matvec(x, Phi, GammaLam, rowPerm, revRowPerm);

      // Set column of results:
      bfMatDenseRealSetCol(matvecs, j, tmp1);

      // Clean up:
      bfVecDelete(&x);
      bfVecDelete(&tmp1);
  }

  sprintf(filename, "matvecs_lbo_tol%.0e_kappa%.1e_nu%.1e.bin", tol, kappa, nu);
  bfMatDenseRealSave(matvecs, filename);

  /* Clean up */
  // bfTreeDelete(&rowTree); // This segfaults...
  // bfOctreeDeinit(&octree);
  bfMatDelete(&M);
  bfMatDelete(&L);
  bfPoints3DeinitAndDealloc(&points);
}
