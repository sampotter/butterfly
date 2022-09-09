#include <bf/helm2.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <bf/bessel.h>
#include <bf/circle.h>
#include <bf/const.h>
#include <bf/error_macros.h>
#include <bf/points.h>
#include <bf/vectors.h>

BfSize bfHelm2RankEstForTwoCircles(BfCircle2 const *circ1,
                                   BfCircle2 const *circ2,
                                   BfReal k, BfReal C, BfReal eps)
{
  assert(k > 0);
  assert(C > 0);
  assert(eps > 0);

  BfReal r1 = circ1->r;
  BfReal r2 = circ2->r;
  BfReal R = bfPoint2Dist(circ1->center, circ2->center);
  assert(r1 > 0 && r2 > 0 && R > 0);

  BfReal d = R - r1 - r2;
  assert(d > 0);

  /* Reading Michielssen & Boag, seems like there should be a factor
   * of two pi in front of this, but that leads to a rank estimate
   * which is way higher than necessary. */
  BfReal p = k*r1*r2/d - C*log10(eps);
  assert(p > 0);

  BfSize rank = (BfSize)ceil(p);
  assert(rank > 0);
  return rank;
}

BfComplex
bfHelm2GetKernelValue(BfReal r, BfReal K) {
  return r == 0 ? NAN : I*bf_H0(K*r)/4;
}

BfMatDenseComplex *
bfGetHelm2KernelMatrix(BfPoints2 const *srcPts, BfPoints2 const *tgtPts, BfReal K)
{
  BEGIN_ERROR_HANDLING();

  BfSize m = tgtPts->size; /* number of rows */
  BfSize n = srcPts->size; /* number of columns */
  BfReal *r = NULL;

  /* length m*n array of pairwise dists in row major order */
  r = bfPoints2PairwiseDists(tgtPts, srcPts);
  HANDLE_ERROR();

  BfMatDenseComplex *kernelMat = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(kernelMat, m, n);
  HANDLE_ERROR();

  BfSize k = 0;
  for (BfSize i = 0; i < m; ++i) {
    for (BfSize j = 0; j < n; ++j) {
      kernelMat->data[k] = r[k] == 0 ? 0 : I*bf_H0(K*r[k])/4;
      ++k;
    }
  }

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&kernelMat);

  free(r);

  return kernelMat;
}

BfMat *bf_hh2_get_dGdN(BfPoints2 const *Xsrc,
                       BfPoints2 const *Xtgt,
                       BfReal K,
                       BfVectors2 const *Ntgt) {
  BEGIN_ERROR_HANDLING();

  BfSize m = Xtgt->size; /* number of rows */
  BfSize n = Xsrc->size; /* number of columns */
  BfReal *r = NULL;

  /* length m*n array of pairwise dists in row major order */
  r = bfPoints2PairwiseDists(Xtgt, Xsrc);
  HANDLE_ERROR();

  BfMatDenseComplex *kernelMat = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(kernelMat, m, n);
  HANDLE_ERROR();

  BfSize k = 0;
  for (BfSize i = 0; i < m; ++i) {
    BfReal const *xsrc = Xsrc->data[i];
    for (BfSize j = 0; j < n; ++j) {
      if (r[k] == 0) {
        kernelMat->data[k] = 0;
      } else {
        BfReal const *xtgt = Xtgt->data[j];
        BfReal const *ntgt = Ntgt->data[j];
        BfReal dot = ntgt[0]*(xtgt[0] - xsrc[0]) + ntgt[1]*(xtgt[1] - xsrc[1]);
        BfReal scale = (I/4)*K*bf_H1(K*r[k])/r[k];
        kernelMat->data[k] = scale*dot;
      }
      ++k;
    }
  }

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&kernelMat);

  free(r);

  return bfMatDenseComplexToMat(kernelMat);
}

BfMatDenseComplex *
bfHelm2GetReexpansionMatrix(BfPoints2 const *srcPtsOrig,
                            BfPoints2 const *srcPtsEquiv,
                            BfPoints2 const *tgtPts, BfReal K)
{
  BEGIN_ERROR_HANDLING();

  /* compute the kernel matrix mapping charges on the original sources
   * points to potentials on the original target points */
  BfMatDenseComplex *Z_orig = bfGetHelm2KernelMatrix(srcPtsOrig, tgtPts, K);
  HANDLE_ERROR();

  /* compute the kernel matrix mapping charges on the source
   * circle to potentials on the target circle */
  BfMatDenseComplex *Z_equiv = bfGetHelm2KernelMatrix(srcPtsEquiv, tgtPts, K);
  HANDLE_ERROR();

  /* set the "shift matrix" to Z_equiv\Z_orig */
  BfMatDenseComplex *Z_shift = bfMatDenseComplexDenseComplexLstSq(Z_equiv, Z_orig);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfMatDenseComplexDeinit(Z_shift);
    bfMatDenseComplexDealloc(&Z_shift);
  }

  bfMatDenseComplexDeinitAndDealloc(&Z_orig);
  bfMatDenseComplexDeinitAndDealloc(&Z_equiv);

  return Z_shift;
}
