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

BfSize bfHelm2RankEstForTwoCircles(BfCircle const *circ1,
                                   BfCircle const *circ2,
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

BfComplex get_S_value(BfPoint2 const xsrc, BfPoint2 const xtgt, BfReal K) {
  BfReal r = bfPoint2Dist(xsrc, xtgt);
  return r == 0 ? NAN : I*bf_H0(K*r)/4;
}

BfComplex get_Sp_value(BfPoint2 const xsrc, BfPoint2 const xtgt, BfVector2 const ntgt, BfReal K) {
  BfReal r = bfPoint2Dist(xsrc, xtgt);
  if (r < 1e-15)
    return 0;
  BfReal dot = ntgt[0]*(xtgt[0] - xsrc[0]) + ntgt[1]*(xtgt[1] - xsrc[1]);
  BfComplex scale = (I/4)*K*bf_H1(K*r)/r;
  return scale*dot;
}

BfComplex get_D_value(BfPoint2 const xsrc, BfPoint2 const xtgt, BfVector2 const nsrc, BfReal K) {
  BfReal r = bfPoint2Dist(xsrc, xtgt);
  if (r < 1e-15)
    return 0;
  BfReal dot = nsrc[0]*(xtgt[0] - xsrc[0]) + nsrc[1]*(xtgt[1] - xsrc[1]);
  BfComplex scale = (I/4)*K*bf_H1(K*r)/r;
  return scale*dot;
}

BfComplex bfHelm2GetKernelValue(BfPoint2 const xsrc, BfPoint2 const xtgt,
                                BfVector2 const nsrc, BfVector2 const ntgt,
                                BfReal K, BfLayerPotential layerPot) {
  BEGIN_ERROR_HANDLING();

  BfComplex z;

  switch (layerPot) {
  case BF_LAYER_POTENTIAL_SINGLE:
    z = get_S_value(xsrc, xtgt, K);
    break;
  case BF_LAYER_POTENTIAL_PV_DOUBLE:
    z = get_D_value(xsrc, xtgt, nsrc, K);
    break;
  case BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE:
    z = get_Sp_value(xsrc, xtgt, ntgt, K);
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {
    z = NAN;
  }

  return z;
}

static BfMat *
get_S_kernel_matrix(BfPoints2 const *Xsrc, BfPoints2 const *Xtgt, BfReal K) {
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
    for (BfSize j = 0; j < n; ++j) {
      kernelMat->data[k] = r[k] == 0 ? 0 : I*bf_H0(K*r[k])/4;
      ++k;
    }
  }

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&kernelMat);

  free(r);

  return bfMatDenseComplexToMat(kernelMat);
}

static BfMat *
get_Sp_kernel_matrix(BfPoints2 const *Xsrc, BfPoints2 const *Xtgt,
                       BfVectors2 const *Ntgt, BfReal K) {
  BEGIN_ERROR_HANDLING();

  BfSize m = Xtgt->size; /* number of rows */
  BfSize n = Xsrc->size; /* number of columns */
  BfReal *r = NULL;

  BfMatDenseComplex *kernelMat = NULL;

  if (K <= 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* length m*n array of pairwise dists in row major order */
  r = bfPoints2PairwiseDists(Xtgt, Xsrc);
  HANDLE_ERROR();

  kernelMat = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(kernelMat, m, n);
  HANDLE_ERROR();

  BfSize k = 0;
  for (BfSize i = 0; i < m; ++i) {
    BfReal const *xtgt = Xtgt->data[i];
    BfReal const *ntgt = Ntgt->data[i];
    for (BfSize j = 0; j < n; ++j) {
      if (r[k] == 0) {
        kernelMat->data[k] = 0;
      } else {
        BfReal const *xsrc = Xsrc->data[j];
        BfReal dot = ntgt[0]*(xtgt[0] - xsrc[0]) + ntgt[1]*(xtgt[1] - xsrc[1]);
        BfComplex scale = (I/4)*K*bf_H1(K*r[k])/r[k];
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

static BfMat *
get_D_kernel_matrix(BfPoints2 const *Xsrc, BfPoints2 const *Xtgt,
                    BfVectors2 const *Nsrc, BfReal K) {
  BEGIN_ERROR_HANDLING();

  BfSize m = Xtgt->size; /* number of rows */
  BfSize n = Xsrc->size; /* number of columns */
  BfReal *r = NULL;

  BfMatDenseComplex *kernelMat = NULL;

  if (K <= 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* length m*n array of pairwise dists in row major order */
  r = bfPoints2PairwiseDists(Xtgt, Xsrc);
  HANDLE_ERROR();

  kernelMat = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(kernelMat, m, n);
  HANDLE_ERROR();

  BfSize k = 0;
  for (BfSize i = 0; i < m; ++i) {
    BfReal const *xtgt = Xtgt->data[i];
    for (BfSize j = 0; j < n; ++j) {
      if (r[k] == 0) {
        kernelMat->data[k] = 0;
      } else {
        BfReal const *xsrc = Xsrc->data[j];
        BfReal const *nsrc = Nsrc->data[j];
        BfReal dot = -(nsrc[0]*(xtgt[0] - xsrc[0]) + nsrc[1]*(xtgt[1] - xsrc[1]));
        BfComplex scale = (I/4)*K*bf_H1(K*r[k])/r[k];
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

static BfMat *
get_S_plus_D_kernel_matrix(BfPoints2 const *Xsrc, BfPoints2 const *Xtgt,
                           BfVectors2 const *Nsrc, BfReal K,
                           BfComplex alpha, BfComplex beta) {
  BEGIN_ERROR_HANDLING();

  BfSize m = Xtgt->size; /* number of rows */
  BfSize n = Xsrc->size; /* number of columns */
  BfReal *r = NULL;

  BfMatDenseComplex *kernelMat = NULL;

  if (K <= 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* length m*n array of pairwise dists in row major order */
  r = bfPoints2PairwiseDists(Xtgt, Xsrc);
  HANDLE_ERROR();

  kernelMat = bfMatDenseComplexNew();
  HANDLE_ERROR();

  bfMatDenseComplexInit(kernelMat, m, n);
  HANDLE_ERROR();

  BfSize k = 0;
  for (BfSize i = 0; i < m; ++i) {
    BfReal const *xtgt = Xtgt->data[i];
    for (BfSize j = 0; j < n; ++j) {
      if (r[k] == 0) {
        kernelMat->data[k] = 0;
      } else {
        BfReal const *xsrc = Xsrc->data[j];
        BfReal const *nsrc = Nsrc->data[j];

        /* Compute single-layer potential kernel value: */
        BfComplex S = (I/4)*bf_H0(K*r[k]);

        /* Compute double-layer potential kernel value: */
        BfReal dot = -(nsrc[0]*(xtgt[0] - xsrc[0]) + nsrc[1]*(xtgt[1] - xsrc[1]));
        BfComplex scale = (I/4)*K*bf_H1(K*r[k])/r[k];
        BfComplex D = scale*dot;

        kernelMat->data[k] = alpha*S + beta*D;
      }
      ++k;
    }
  }

  END_ERROR_HANDLING()
    bfMatDenseComplexDeinitAndDealloc(&kernelMat);

  free(r);

  return bfMatDenseComplexToMat(kernelMat);
}

BfMat *
bfGetHelm2KernelMatrix(BfPoints2 const *Xsrc, BfPoints2 const *Xtgt,
                       BfVectors2 const *Nsrc, BfVectors2 const *Ntgt,
                       BfReal K, BfLayerPotential layerPot,
                       BfComplex const *alpha, BfComplex const *beta)
{
  BEGIN_ERROR_HANDLING();

  if (Nsrc != NULL && bfPoints2GetSize(Xsrc) != bfVectors2GetSize(Nsrc))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (Ntgt != NULL && bfPoints2GetSize(Xtgt) != bfVectors2GetSize(Ntgt))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMat *kernelMat = NULL;

  switch (layerPot) {
  case BF_LAYER_POTENTIAL_SINGLE:
    kernelMat = get_S_kernel_matrix(Xsrc, Xtgt, K);
    HANDLE_ERROR();
    break;
  case BF_LAYER_POTENTIAL_PV_DOUBLE:
    kernelMat = get_D_kernel_matrix(Xsrc, Xtgt, Nsrc, K);
    break;
  case BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE:
    kernelMat = get_Sp_kernel_matrix(Xsrc, Xtgt, Ntgt, K);
    HANDLE_ERROR();
    break;
  case BF_LAYER_POTENTIAL_COMBINED_FIELD:
    kernelMat = get_S_plus_D_kernel_matrix(Xsrc, Xtgt, Nsrc, K, *alpha, *beta);
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  END_ERROR_HANDLING() {}

  return kernelMat;
}

BfMat *
bfHelm2GetReexpansionMatrix(BfPoints2 const *srcPtsOrig,
                            BfPoints2 const *srcPtsEquiv,
                            BfVectors2 const *srcNormalsOrig,
                            BfVectors2 const *srcNormalsEquiv,
                            BfPoints2 const *tgtPts,
                            BfReal K, BfLayerPotential layerPot,
                            BfComplex const *alpha, BfComplex const *beta)
{
  BEGIN_ERROR_HANDLING();

  if (srcNormalsOrig != NULL
      && bfPoints2GetSize(srcPtsOrig) != bfVectors2GetSize(srcNormalsOrig))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (srcNormalsEquiv != NULL
      && bfPoints2GetSize(srcPtsEquiv) != bfVectors2GetSize(srcNormalsEquiv))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Computing a reexpansion matrix doesn't make any sense for a layer
   * potential which depends on unit normals at target points */
  if (BF_LAYER_POT_USES_TGT_NORMALS[layerPot])
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* compute the kernel matrix mapping charges on the original sources
   * points to potentials on the original target points */
  BfMat *Z_orig = bfGetHelm2KernelMatrix(srcPtsOrig, tgtPts, srcNormalsOrig, NULL, K, layerPot, alpha, beta);
  HANDLE_ERROR();

  /* compute the kernel matrix mapping charges on the source
   * circle to potentials on the target circle */
  BfMat *Z_equiv = bfGetHelm2KernelMatrix(srcPtsEquiv, tgtPts, srcNormalsEquiv, NULL, K, layerPot, alpha, beta);
  HANDLE_ERROR();

  /* set the "shift matrix" to Z_equiv\Z_orig */
  BfMat *Z_shift = bfMatLstSq(Z_equiv, Z_orig);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatDelete(&Z_shift);

  bfMatDelete(&Z_orig);
  bfMatDelete(&Z_equiv);

  return Z_shift;
}
