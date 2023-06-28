#include "multiple_scattering_context.h"

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/ellipse.h>
#include <bf/fac_helm2.h>
#include <bf/helm2.h>
#include <bf/linalg.h>
#include <bf/lu.h>
#include <bf/mat_block_diag.h>
// #include <bf/mat_diag_real.h>
#include <bf/mat_func.h>
#include <bf/mat_perm.h>
#include <bf/mem.h>
#include <bf/poisson_disk_sampling.h>
#include <bf/rand.h>
// #include <bf/real_array.h>
#include <bf/size_array.h>
#include <bf/util.h>
#include <bf/vec_real.h>
#include <bf/vectors.h>

#include <math.h>

#include <fmm2d/helmholtz.h>

BfComplex K_helm2(BfSize i, BfSize j, void *aux) {
  KHelm2Wkspc *wkspc = aux;

  BfReal const *xsrc = &wkspc->points->data[i][0];
  BfReal const *xtgt = &wkspc->points->data[j][0];
  BfReal const *nsrc = &wkspc->normals->data[i][0];
  BfReal k = wkspc->k;

  BfComplex G = bfHelm2GetKernelValue(
    xsrc, xtgt, NULL, NULL, k, BF_LAYER_POTENTIAL_SINGLE);
  BfComplex dGdN = bfHelm2GetKernelValue(
    xsrc, xtgt, nsrc, NULL, k, BF_LAYER_POTENTIAL_PV_DOUBLE);
  BfComplex K = wkspc->alpha*G + wkspc->beta*dGdN;

  return K;
}

static BfMat *mulFmm(BfMat const *sigma, void const *aux) {
  MultipleScatteringContext const *context = (MultipleScatteringContext const *)aux;

//   /* Assemble discretized CFIE kernel matrix: */
//   BfMat *KDenseTest = bfGetHelm2KernelMatrix(
//     context->X,
//     context->X,
//     context->N,
//     NULL,
//     context->k,
//     BF_LAYER_POTENTIAL_COMBINED_FIELD,
//     &context->alpha,
//     &context->beta);
//   /* Apply KR correction blockwise: */
//   bf_apply_block_KR_correction(
//     KDenseTest,
//     context->ellipseOffsets,
//     context->orderKR,
//     K_helm2,
//     (void *)&context->KWkspc);
//   BfVec *WVec = bfRealArrayGetVecView(context->W);
//   bfMatScaleCols(KDenseTest, WVec);
//   bfVecDelete(&WVec);
//   BfMat *oneHalfEye = bfMatDiagRealToMat(
//     bfMatDiagRealNewConstant(context->n, context->n, 0.5));
//   bfMatAddInplace(KDenseTest, oneHalfEye);
//   bfMatDelete(&oneHalfEye);

  /** Allocate space and do some basic setup: */

  BfSize n = context->X->size;

  BfSize numRhs = 1;

  BF_ASSERT(context->N->size == n);
  BF_ASSERT(bfMatGetNumRows(sigma) == n);
  BF_ASSERT(bfMatGetNumCols(sigma) == numRhs); // TODO: handle multiple righthand sides
  BF_ASSERT(bfMatConstToMatDenseComplexConst(sigma)->rowStride == 1); // XXX

  BfMatDenseComplex *u = bfMatDenseComplexNew();
  bfMatDenseComplexInit(u, n, numRhs);

  BfComplex *uPtr = u->data;

  /** Copy over sigma, scale it by the trapezoid rule weights, and set
   ** up the `charge` and `dipstr` arguments for `hmm2d_s_cd_p`: */

  BfComplex *sigmaPtr = bfMatConstToMatDenseComplexConst(sigma)->data;

  BfComplex chargeScale = context->alpha;
  BfComplex dipstrScale = context->beta;

  BfComplex *charge = bfMemAlloc(n, sizeof(BfComplex));
  BF_ASSERT(charge != NULL);
  for (BfSize i = 0; i < n; ++i)
    charge[i] = chargeScale*sigmaPtr[i];

  BfComplex *dipstr = bfMemAlloc(n, sizeof(BfComplex));
  BF_ASSERT(dipstr != NULL);
  for (BfSize i = 0; i < n; ++i)
    dipstr[i] = dipstrScale*sigmaPtr[i];

  BfReal const *w = context->W->data;
  for (BfSize i = 0; i < n; ++i) charge[i] *= w[i];
  for (BfSize i = 0; i < n; ++i) dipstr[i] *= w[i];

  /** Everything has been set up---call `hfmm2d_s_cd_p` now: */

  int64_t ier = 0;
  hfmm2d_s_cd_p(
    /* eps: */ context->tolFmm,
    /* zk: */ context->k,
    /* ns: */ n,
    /* sources: */ (BfReal const *)&context->X->data[0],
    /* charge: */ charge,
    /* dipstr: */ dipstr,
    /* dipvec: */ (BfReal const *)&context->N->data[0],
    /* pot: */ uPtr,
    /* ier: */ &ier);

  // TODO: error handling
  //
  //   ier == 1, 2, 3, or 4 -> "allocate" failed
  //   ier == 13 -> "computational box too big"
  //
  BF_ASSERT(ier == 0);

  bfMemFree(charge);
  bfMemFree(dipstr);

  /** Apply the KR correction: */

  for (BfSize i = 0; i < context->numEllipses; ++i) {
    BfSize i0 = bfSizeArrayGet(context->ellipseOffsets, i);
    BfSize i1 = bfSizeArrayGet(context->ellipseOffsets, i + 1);

    BfPoints2 *XBlock = bfPoints2GetRangeView(context->X, i0, i1);
    BfVectors2 *NBlock = bfVectors2GetRangeView(context->N, i0, i1);

    KHelm2Wkspc KWkspcBlock = context->KWkspc;
    KWkspcBlock.points = XBlock;
    KWkspcBlock.normals = NBlock;

    BfSize n = i1 - i0;

    bf_accum_with_KR_correction(
      /* order: */ context->orderKR,
      /* K: */ K_helm2,
      /* aux: */ (void *)&KWkspcBlock,
      /* n: */ n,
      /* x: */ &sigmaPtr[i0],
      /* h: */ &w[i0],
      /* y: */ &uPtr[i0]);

    bfMemFree(XBlock);
    bfMemFree(NBlock);
  }

  /** Multiply by I/2: */

  for (BfSize i = 0; i < n; ++i)
    uPtr[i] += sigmaPtr[i]/2;

//   bfMatSave(bfMatDenseComplexToMat(u), "testFmm.bin");
//   bfMatSave(bfMatMul(KDenseTest, sigma), "testDense.bin");

  return bfMatDenseComplexToMat(u);
}

void init(MultipleScatteringContext *context, Opts const *opts) {
  /** Copy over parameters from `opts`: */

  context->k = opts->wavenumber;
  context->minDist = opts->minDist;
  context->axisLow = opts->axisLow;
  context->axisHigh = opts->axisHigh;
  context->alpha = opts->alpha;
  context->beta = opts->beta;
  context->d[0] = opts->d[0];
  context->d[1] = opts->d[1];

  context->h = opts->h;
  context->tol = opts->tol;
  context->tolFmm = BF_NAN;
  context->orderKR = opts->orderKR;

  context->bboxEval = opts->bboxEval;
  context->nxEval = opts->nxEval;
  context->nyEval = opts->nyEval;

  /** Default-initialize "derived" parameters: */

  context->R = BF_NAN;
  context->ellipseBbox = bfGetEmptyBbox2();
  context->ellipseCenters = NULL;
  context->numEllipses = BF_SIZE_BAD_VALUE;
  context->ellipse = NULL;
  context->ellipseOffsets = NULL;

  context->X = NULL;
  context->N = NULL;
  context->W = NULL;
  context->n = BF_SIZE_BAD_VALUE;

  /* Set by `buildQuadtrees`: */
  context->quadtree = NULL;
  context->perm = NULL;
  context->revPerm = NULL;

  /* Set by `setUpKrWorkspace`: */
  context->KWkspc = (KHelm2Wkspc) {
    .points = NULL,
    .normals = NULL,
    .k = BF_NAN,
    .alpha = context->alpha,
    .beta = context->beta,
  };

  context->K = NULL;
  context->KButterfly = NULL;
  context->KButterflyDense = NULL;
  context->KLu = NULL;

  context->KBlockLus = NULL;

  context->M = NULL;
  context->MPerm = NULL;

  context->rhs = NULL;
  context->rhsPerm = NULL;
  context->yTestDense = NULL;
  context->yTestButterfly = NULL;

  context->numIterDense = BF_SIZE_BAD_VALUE;
  context->numIterDensePrecondLeft = BF_SIZE_BAD_VALUE;
  context->numIterButterfly = BF_SIZE_BAD_VALUE;
  context->numIterButterflyPrecondLeft = BF_SIZE_BAD_VALUE;

  context->sigmaLu = NULL;
  context->sigmaDense = NULL;
  context->sigmaDensePrecondLeft = NULL;
  context->sigmaFmm = NULL;
  context->sigmaFmmPrecondLeft = NULL;
  context->sigmaButterfly = NULL;
  context->sigmaButterflyPrecondLeft = NULL;

  context->XEval = NULL;

  context->quadtreeEval = NULL;
  context->permEval = NULL;

  context->KEvalButterfly = NULL;

  context->uIn = NULL;
}


void printInfo(MultipleScatteringContext const *context) {
  puts("Solving multiple scattering problem:");
  printf("- uIn = exp(i*k*d*r), where:\n");
  printf("  * k = %g\n", context->k);
  printf("  * d = (%g, %g)\n", context->d[0], context->d[1]);
  printf("- average points per wavelength: %0.2f\n", (2*BF_PI/context->k)/context->h);
  printf("- using CFIE: alpha*G + beta*dG/dn, where:\n");
  printf("  * alpha = %g + %gi\n", creal(context->alpha), cimag(context->alpha));
  printf("  * beta = %g + %gi\n", creal(context->beta), cimag(context->beta));
  printf("- number of wavelengths across domain: %g\n", 2.0*sqrt(2.0)*context->k/(2*BF_PI));
}

/* Set up problem geometry: randomly sample lots of little
 * well-separated ellipses with different semimajor/minor axes and
 * orientations: */
void setUpGeometry(MultipleScatteringContext *context) {
  printf("Setting up problem geometry...");
  fflush(stdout);

  bfToc();

  /* Set up bounding box of sampling domain---give a little margin to
   * keep ellipses fully inside */
  context->R = 1 - context->minDist;
  context->ellipseBbox = (BfBbox2) {
    .min = {-context->R, -context->R},
    .max = { context->R,  context->R}
  };

  /* Sample ellipse centers */
  context->ellipseCenters = bfPoints2SamplePoissonDisk(
    &context->ellipseBbox,
    context->minDist,
    30);
  bfSavePoints2(context->ellipseCenters, "ellipseCenters.bin");

  context->numEllipses = bfPoints2GetSize(context->ellipseCenters);

  /* Randomly sample ellipses */
  context->ellipse = bfMemAlloc(context->numEllipses, sizeof(BfEllipse));
  for (BfSize i = 0; i < context->numEllipses; ++i) {
    BfEllipse *ellipse = &context->ellipse[i];

    /* Sample ellipse semi-major and semi-minor axes: */
    BfReal a = (context->axisHigh - context->axisLow)*bfRealUniform1() + context->axisLow;
    BfReal b = (context->axisHigh - context->axisLow)*bfRealUniform1() + context->axisLow;
    ellipse->semiMajorAxis = fmax(a, b);
    ellipse->semiMinorAxis = fmin(a, b);

    /* Give ellipse random orientation */
    ellipse->theta = BF_TWO_PI*bfRealUniform1();

    /* Set ellipse center */
    bfMemCopy(context->ellipseCenters->data[i], 1, sizeof(BfPoint2), ellipse->center);
  }

  printf(" sampled %lu ellipses [%0.2fs]\n", context->numEllipses, bfToc());

  FILE *fp = fopen("ellipseData.bin", "w");
  for (BfSize i = 0; i < context->numEllipses; ++i) {
    BfEllipse const *ellipse = &context->ellipse[i];
    fwrite(&ellipse->center[0], sizeof(BfReal), 1, fp);
    fwrite(&ellipse->center[1], sizeof(BfReal), 1, fp);
    fwrite(&ellipse->semiMajorAxis, sizeof(BfReal), 1, fp);
    fwrite(&ellipse->semiMinorAxis, sizeof(BfReal), 1, fp);
    fwrite(&ellipse->theta, sizeof(BfReal), 1, fp);
  }
  fclose(fp);
}

/* Sample discretization points from the ellipses using inverse
 * curvature weighting: */
void setUpDiscretization(MultipleScatteringContext *context) {
  printf("Setting up discretization...");
  fflush(stdout);

  bfToc();

  /* Keep track of the offset to each ellipse---we'll use these
   * indices to figure out where to apply the KR correction and do
   * preconditioning *blockwise*. */
  context->ellipseOffsets = bfSizeArrayNewWithDefaultCapacity();
  bfSizeArrayAppend(context->ellipseOffsets, 0);

  context->X = bfPoints2NewEmpty();
  context->N = bfVectors2NewEmpty();
  context->W = bfRealArrayNewWithDefaultCapacity();
  for (BfSize i = 0; i < context->numEllipses; ++i) {
    BfEllipse const *ellipse = &context->ellipse[i];
    BfReal p = bfEllipseGetPerimeter(ellipse);
    BfSize n = floor(p/context->h) + 1;
    bfEllipseSampleLinspaced(ellipse, n, context->X, NULL, context->N, context->W);
    // bfEllipseSampleEquispaced(ellipse, n, context->X, NULL, context->N);
    // bfEllipseSampleWithInverseCurvatureSpacing(ellipse, n, context->X, NULL, context->N);
    bfSizeArrayAppend(context->ellipseOffsets, bfPoints2GetSize(context->X));
  }
  bfSavePoints2(context->X, "X.bin");
  bfSaveVectors2(context->N, "N.bin");

  context->n = bfPoints2GetSize(context->X);

  printf(" sampled %lu points [%0.2fs]\n", context->n, bfToc());

  FILE *fp = fopen("ellipseOffsets.bin", "w");
  for (BfSize i = 0; i <= context->numEllipses; ++i) {
    BfSize i0 = bfSizeArrayGet(context->ellipseOffsets, i);
    fwrite(&i0, sizeof(BfSize), 1, fp);
  }
  fclose(fp);
}

/* Build quadtree on discretization points and get the forward and
 * reverse permutations. */
void buildQuadtrees(MultipleScatteringContext *context) {
  bfToc();
  context->quadtree = bfQuadtreeNew();
  bfQuadtreeInit(context->quadtree, context->X, context->N);

  BfTree *tree = bfQuadtreeToTree(context->quadtree);
  context->perm = bfTreeGetPermConst(tree);

  context->revPerm = bfMemAlloc(1, sizeof(BfPerm));
  *context->revPerm = bfPermGetReversePerm(context->perm);

  printf("Built quadtrees [%0.2fs]\n", bfToc());

  FILE *fp = fopen("perm.bin", "w");
  fwrite(context->perm->index, sizeof(BfSize), context->perm->size, fp);
  fclose(fp);
}

/* Set up workspace for applying KR corrections: */
void setUpKrWorkspace(MultipleScatteringContext *context) {
  context->KWkspc.points = context->X;
  context->KWkspc.normals = context->N;
  context->KWkspc.k = context->k;
}

/* Set up the RHS for the scattering problem: */
void setUpRhs(MultipleScatteringContext *context) {
  bfToc();

  context->rhs = NULL;
  { BfMatDenseComplex *_ = bfMatDenseComplexNew();
    bfMatDenseComplexInit(_, context->n, 1);
    for (BfSize i = 0; i < context->n; ++i) {
      BfPoint2 x;
      bfPoints2Get(context->X, i, x);
      *(_->data + i*_->rowStride) = -cexp(I*context->k*(context->d[0]*x[0] + context->d[1]*x[1]));
    }
    context->rhs = bfMatDenseComplexToMat(_); }

  context->rhsPerm = bfMatCopy(context->rhs);
  bfMatPermuteRows(context->rhsPerm, context->revPerm);

  printf("Set up RHS for scattering problem [%0.2fs]\n", bfToc());

  bfMatSave(context->rhs, "rhs.bin");
}

/* Set up the dense kernel matrix using the combined field integral
 * equation for the exterior Dirichlet problem: */
void assembleDenseK(MultipleScatteringContext *context) {
  printf("Assembling dense kernel matrix...");
  fflush(stdout);

  bfToc();

  /* Assemble discretized CFIE kernel matrix: */
  context->K = bfGetHelm2KernelMatrix(
    context->X,
    context->X,
    context->N,
    NULL,
    context->k,
    BF_LAYER_POTENTIAL_COMBINED_FIELD,
    &context->alpha,
    &context->beta);

  /* Apply KR correction blockwise: */
  bf_apply_block_KR_correction(
    context->K,
    context->ellipseOffsets,
    context->orderKR,
    K_helm2,
    (void *)&context->KWkspc);

  /* Scale columns by trapezoid rule weights: */
  BfVec *WVec = bfRealArrayGetVecView(context->W);
  bfMatScaleCols(context->K, WVec);
  bfVecDelete(&WVec);

  /* Perturb by I/2: */
  BfMat *oneHalfEye = bfMatDiagRealToMat(
    bfMatDiagRealNewConstant(context->n, context->n, 0.5));
  bfMatAddInplace(context->K, oneHalfEye);
  bfMatDelete(&oneHalfEye);

  printf(" done [%0.2fs]\n", bfToc());

  bfMatSave(context->K, "K.bin");
}

/* Assemble the multilevel BF approximation of the system matrix: */
void assembleButterfliedK(MultipleScatteringContext *context) {
  printf("Assembling butterflied kernel matrix...");
  fflush(stdout);

  bfToc();

  /* Assemble and butterfly compress discretized CFIE kernel matrix: */
  context->KButterfly = bfFacHelm2MakeMultilevel(
    context->quadtree,
    context->quadtree,
    context->k,
    BF_LAYER_POTENTIAL_COMBINED_FIELD,
    &context->alpha,
    &context->beta);

  /* Apply KR correction blockwise: */
  bf_apply_block_KR_correction_quadtree(
    context->KButterfly,
    context->ellipseOffsets,
    context->orderKR,
    bfQuadtreeConstToTreeConst(context->quadtree),
    K_helm2,
    (void *)&context->KWkspc);

  /* Scale columns by trapezoid rule weights: */
  BfVec *WVec = bfRealArrayGetVecView(context->W);
  BfVec *WVecPerm = bfVecCopy(WVec);
  bfVecPermute(WVecPerm, context->revPerm);
  bfMatScaleCols(context->KButterfly, WVecPerm);
  bfVecDelete(&WVecPerm);
  bfVecDelete(&WVec);

  /* Perturb by I/2: */
  BfMat *oneHalfEye = bfMatDiagRealToMat(
    bfMatDiagRealNewConstant(context->n, context->n, 0.5));
  bfMatAddInplace(context->KButterfly, oneHalfEye);
  bfMatDelete(&oneHalfEye);

  printf(" done [%0.2fs]\n", bfToc());

  /* Write BF blocks to disk: */
  FILE *fp = fopen("blocks.txt", "w");
  bfPrintBlocks(context->KButterfly, 2, fp);
  fclose(fp);

  printf("- saved blocks to blocks.txt\n");

  BfSize numBytes = bfMatNumBytes(context->KButterfly);
  printf("- size of factorization: ");
  if (numBytes < 1024) printf("%lu B\n", numBytes);
  else if (numBytes < 1024*1024) printf("%0.1f KB\n", numBytes/1024.0);
  else if (numBytes < 1024*1024*1024) printf("%0.1f MB\n", numBytes/pow(1024, 2));
  else printf("%0.1f GB\n", numBytes/pow(1024, 3));

  BfReal numBytesDense = 8*pow(context->n, 2);
  printf("- compression rate: %0.1f\n", numBytesDense/numBytes);
}

/* Extract components of BF matrix (NOTE: very expensive!!!): */
void extractDenseButterfliedK(MultipleScatteringContext *context) {
  printf("Extracting dense version of butterflied kernel matrix...");
  fflush(stdout);

  bfToc();

  context->KButterflyDense = bfMatEmptyLike(context->K, context->n, context->n);
  for (BfSize j = 0; j < context->n; ++j) {
    /* Get jth standard basis vector: */
    BfMat *ej; {
      BfMatDenseComplex *_ = bfMatDenseComplexNew();
      bfMatDenseComplexInit(_, context->n, 1);
      for (BfSize i = 0; i < context->n; ++i)
        *(_->data + i*_->rowStride) = i == j ? 1 : 0;
      ej = bfMatDenseComplexToMat(_); }
    BfMat *kj = bfMatMul(context->KButterfly, ej);
    BfVec *kjVecView = bfMatGetColView(kj, 0);
    bfMatSetCol(context->KButterflyDense, j, kjVecView);
    bfVecDelete(&kjVecView);
    bfMatDelete(&kj);
    bfMatDelete(&ej);
  }

  printf(" done [%0.2fs]\n", bfToc());

  bfMatSave(context->KButterflyDense, "KBFDense.bin");
}

void assembleFmmK(MultipleScatteringContext *context) {
  BfMatFunc *KFmm = bfMatFuncNew();
  bfMatFuncInit(KFmm, context->n, context->n, (MatMulFunc)mulFmm, context);
  context->KFmm = bfMatFuncToMat(KFmm);
}

/* Set up block Jacobi preconditioner. */
void assemblePreconditioner(MultipleScatteringContext *context) {
  printf("Assembling block Jacobi preconditioner... ");
  fflush(stdout);

  bfToc();

  context->KBlockLus = bfPtrArrayNewWithDefaultCapacity();

  /* Build block Jacobi preconditioner: */
  BfPtrArray *MBlocks = bfPtrArrayNewWithDefaultCapacity();
  for (BfSize i = 0; i < context->numEllipses; ++i) {
    BfSize i0 = bfSizeArrayGet(context->ellipseOffsets, i);
    BfSize i1 = bfSizeArrayGet(context->ellipseOffsets, i + 1);

    /* Get the current of K: if we've already computed the dense
     * kernel matrix, just grab a view of it; otherwise we need to
     * compute the dense block now. */
    BfMat *KBlock = NULL;
    if (context->K) {
      KBlock = bfMatGetBlockView(context->K, i0, i1, i0, i1);
    } else {
      BfPoints2 *XBlock = bfPoints2GetRangeView(context->X, i0, i1);
      BfVectors2 *NBlock = bfVectors2GetRangeView(context->N, i0, i1);

      /* Get dense kernel matrix block for CFIE: */
      KBlock = bfGetHelm2KernelMatrix(
        XBlock, XBlock, NBlock, NULL, context->k,
        BF_LAYER_POTENTIAL_COMBINED_FIELD, &context->alpha, &context->beta);

      /* Set up K workspace for this block: */
      KHelm2Wkspc KWkspcBlock = context->KWkspc;
      KWkspcBlock.points = XBlock;
      KWkspcBlock.normals = NBlock;

      /* Apply KR correction: */
      bf_apply_KR_correction(
        KBlock, context->orderKR, K_helm2, (void *)&KWkspcBlock);

      /* Scale columns by trapezoid rule weights: */
      BfVec *WSubvec = bfRealArrayGetSubvecView(context->W, i0, i1);
      bfMatScaleCols(KBlock, WSubvec);
      bfVecDelete(&WSubvec);

      /* Perturb by I/2: */
      BfSize nBlock = i1 - i0;
      BfMat *oneHalfEye = bfMatDiagRealToMat(bfMatDiagRealNewConstant(nBlock, nBlock, 0.5));
      bfMatAddInplace(KBlock, oneHalfEye);
      bfMatDelete(&oneHalfEye);

      // TODO: need to implement Dealloc for Points2 and Vectors2...
      bfMemFree(XBlock);
      bfMemFree(NBlock);
    }

    BfLu *KBlockLu = bfMatGetLu(KBlock);
    bfPtrArrayAppend(context->KBlockLus, KBlockLu);

    BfMat *MBlock = bfLuGetMatView(KBlockLu);
    bfPtrArrayAppend(MBlocks, MBlock);

    bfMatDelete(&KBlock);
  }
  context->M = bfMatBlockDiagToMat(
    bfMatBlockDiagNewFromBlocks(MBlocks, BF_POLICY_STEAL));

  /* Get the permuted version of M for use with the
   * butterfly-accelerated version: */

  context->MPerm = NULL;
  { BfMat *matPerm = bfMatPermToMat(bfMatPermNewFromPerm(context->perm));
    BfMat *matPermInverse = bfMatGetInverse(matPerm);
    BfMatProduct *matProduct = bfMatProductNew();
    bfMatProductInit(matProduct);
    bfMatProductPostMultiply(matProduct, matPerm);
    bfMatProductPostMultiply(matProduct, bfMatGet(context->M, BF_POLICY_VIEW));
    bfMatProductPostMultiply(matProduct, matPermInverse);
    context->MPerm = bfMatProductToMat(matProduct); }

  printf(" done [%0.2fs]\n", bfToc());

  BfSize numBytes = bfMatNumBytes(context->M);
  printf("- size of block Jacobi preconditioner: ");
  if (numBytes < 1024) printf("%lu B\n", numBytes);
  else if (numBytes < 1024*1024) printf("%0.1f KB\n", numBytes/1024.0);
  else if (numBytes < 1024*1024*1024) printf("%0.1f MB\n", numBytes/pow(1024, 2));
  else printf("%0.1f GB\n", numBytes/pow(1024, 3));

  for (BfSize i = 0; i < bfPtrArraySize(MBlocks); ++i) {
    BfMat *MBlock = bfPtrArrayGet(MBlocks, i);
    bfMatDelete(&MBlock);
  }
  bfPtrArrayDelete(&MBlocks);
}

static BfReal getMaxRelErrForFmmTol(MultipleScatteringContext *context, BfReal tolFmm) {
  BfSize const numTrials = 5;

  BfReal oldTolFmm = context->tolFmm;

  context->tolFmm = tolFmm;

  /* Estimate the relative error between the FMM and BF by doing
   * `numTrials` test multiplications and taking the largest
   * relative maximum error. */
  BfReal maxRelError = -BF_INFINITY;
  for (BfSize _ = 0; _ < numTrials; ++_) {
    /* Sample random test vector: */
    BfMat *x; {
      BfMatDenseComplex *_ = bfMatDenseComplexNew();
      bfMatDenseComplexInit(_, context->n, 1);
      bfComplexRandn(context->n, _->data);
      x = bfMatDenseComplexToMat(_);
    }

    /* Permute test vector for BF multiply: */
    BfMat *xPerm = bfMatCopy(x);
    bfMatPermuteRows(xPerm, context->revPerm);

    /* Do BF multiply: */
    BfMat *yButterfly = bfMatMul(context->KButterfly, xPerm);
    bfMatPermuteRows(yButterfly, context->perm);

    /* Do FMM multiply: */
    BfMat *yFmm = mulFmm(x, context);

    /* Compute relative l2 error: */
    BfVec *tmp1 = bfMatColDists(yButterfly, yFmm);
    BfVec *tmp2 = bfMatColNorms(yButterfly);
    BfVec *tmp3 = bfMatColNorms(yFmm);
    BfReal denom = fmax(bfVecNormMax(tmp2), bfVecNormMax(tmp3));
    BfReal relError = bfVecNormMax(tmp1)/denom;

    maxRelError = fmax(maxRelError, relError);

    bfMatDelete(&x);
    bfMatDelete(&xPerm);
    bfMatDelete(&yButterfly);
    bfMatDelete(&yFmm);
    bfVecDelete(&tmp1);
    bfVecDelete(&tmp2);
    bfVecDelete(&tmp3);
  }

  context->tolFmm = oldTolFmm;

  return maxRelError;
}

void estimateFmmTol(MultipleScatteringContext *context) {
  printf("Estimating FMM tolerance... ");
  fflush(stdout);

  BfReal tolFmm = BF_NAN;

  BfReal tolFmmInit = getMaxRelErrForFmmTol(context, 1e-15);

  BfReal p = ceil(-log10(tolFmmInit)) - 0.5;
  while (p > 0) {
    tolFmm = pow(10, -p);
    BfReal maxRelErr = getMaxRelErrForFmmTol(context, tolFmm);
    if (tolFmm > maxRelErr)
      break;
    p -= 0.5;
  }

  BF_ASSERT(p > 0);

  context->tolFmm = tolFmm;

  printf("done [%0.2fs]\n", bfToc());

  printf("- FMM tolerance: %g\n", context->tolFmm);
}

void computeLu(MultipleScatteringContext *context) {
  printf("Computing LU decomposition of dense kernel matrix...");
  fflush(stdout);

  bfToc();

  context->KLu = bfMatGetLu(context->K);

  printf("done [%0.2fs]\n", bfToc());
}

/* Solve using dense elimination. */
void solveLu(MultipleScatteringContext *context) {
  printf("Solving dense system using LU decomposition...");
  fflush(stdout);

  context->sigmaLu = bfLuSolve(context->KLu, context->rhs);

  printf("done [%0.2fs]\n", bfToc());
}

/* Solve the dense system using GMRES. */
void solveDenseGmres(MultipleScatteringContext *context) {
  printf("Solving dense system using GMRES... ");
  fflush(stdout);

  bfToc();

  context->sigmaDense = bfSolveGMRES(
    context->K,
    context->rhs,
    NULL,
    context->tol,
    context->n,
    &context->numIterDense,
    NULL);

  printf("finished in %lu iterations [%0.2fs]\n", context->numIterDense, bfToc());
}

/* Solve dense system using block LU preconditioned GMRES: */
void solveDensePreconditionedGmres(MultipleScatteringContext *context) {
  printf("Solving dense system using preconditioned GMRES... ");
  fflush(stdout);

  bfToc();

  context->sigmaDensePrecondLeft = bfSolveGMRES(
    context->K,
    context->rhs,
    NULL,
    context->tol,
    context->n,
    &context->numIterDensePrecondLeft,
    context->M);

  printf("finished in %lu iterations [%0.2fs]\n", context->numIterDensePrecondLeft, bfToc());
}

/* Solve butterfly-factorized system using GMRES: */
void solveButterflyGmres(MultipleScatteringContext *context) {
  printf("Solving butterflied system using GMRES... ");
  fflush(stdout);

  bfToc();

  context->sigmaButterfly = bfSolveGMRES(
    context->KButterfly,
    context->rhsPerm,
    NULL,
    context->tol,
    context->n,
    &context->numIterButterfly,
    NULL);

  bfMatPermuteRows(context->sigmaButterfly, context->perm);

  printf("finished in %lu iterations [%0.2fs]\n", context->numIterButterfly, bfToc());
}

/* Solve butterfly-factorized system using GMRES: */
void solveButterflyPreconditionedGmres(MultipleScatteringContext *context) {
  printf("Solving butterflied system using preconditioned GMRES... ");
  fflush(stdout);

  bfToc();

  context->sigmaButterflyPrecondLeft = bfSolveGMRES(
    context->KButterfly,
    context->rhsPerm,
    NULL,
    context->tol,
    context->n,
    &context->numIterButterflyPrecondLeft,
    context->MPerm);

  bfMatPermuteRows(context->sigmaButterflyPrecondLeft, context->perm);

  printf(" finished in %lu iterations [%0.2fs]\n", context->numIterButterflyPrecondLeft, bfToc());
}

/* Solve the FMM system using GMRES. */
void solveFmmGmres(MultipleScatteringContext *context) {
  printf("Solving FMM system using GMRES... ");
  fflush(stdout);

  bfToc();

  context->sigmaFmm = bfSolveGMRES(
    context->KFmm,
    context->rhs,
    NULL,
    context->tol,
    context->n,
    &context->numIterFmm,
    NULL);

  printf("finished in %lu iterations [%0.2fs]\n", context->numIterFmm, bfToc());
}

/* Solve FMM system using block LU preconditioned GMRES: */
void solveFmmPreconditionedGmres(MultipleScatteringContext *context) {
  printf("Solving FMM system using preconditioned GMRES... ");
  fflush(stdout);

  bfToc();

  context->sigmaFmmPrecondLeft = bfSolveGMRES(
    context->KFmm,
    context->rhs,
    NULL,
    context->tol,
    context->n,
    &context->numIterFmmPrecondLeft,
    context->M);

  printf("finished in %lu iterations [%0.2fs]\n", context->numIterFmmPrecondLeft, bfToc());
}

void collectAndPrintStats(MultipleScatteringContext *context) {
  puts("Error and timing info:");

  /** Compare dense and butterfly MVPs: */

  BfMat *yTestDense = NULL;
  BfMat *yTestButterfly = NULL;
  BfMat *yTestFmm = NULL;

  if (context->K) {
    bfToc();
    yTestDense = bfMatMul(context->K, context->rhs);
    printf("- did test dense MVP [%0.2fs]\n", bfToc());
    bfMatDelete(&yTestDense);
  }

  if (context->KButterfly) {
    bfToc();
    yTestButterfly = bfMatMul(context->KButterfly, context->rhsPerm);
    bfMatPermuteRows(yTestButterfly, context->perm);
    bfMatDelete(&yTestButterfly);
    printf("- did test butterfly MVP [%0.2fs]\n", bfToc());
  }

  if (context->KFmm) {
    bfToc();
    yTestFmm = bfMatMul(context->KFmm, context->rhs);
    bfMatDelete(&yTestFmm);
    printf("- did test FMM MVP [%0.2fs]\n", bfToc());
  }

  if (yTestDense && yTestButterfly) {
    BfVec *mvpError = bfMatColDists(yTestDense, yTestButterfly);
    BfVec *yColNorms = bfMatColNorms(yTestDense);
    BfReal mvpRelErrorMax = bfVecNormMax(mvpError)/bfVecNormMax(yColNorms);
    bfVecDelete(&mvpError);
    bfVecDelete(&yColNorms);
    printf("- relative error between dense and butterfly MVPs: %g\n", mvpRelErrorMax);
  }

  if (yTestDense && yTestFmm) {
    BfVec *mvpError = bfMatColDists(yTestDense, yTestFmm);
    BfVec *yColNorms = bfMatColNorms(yTestDense);
    BfReal mvpRelErrorMax = bfVecNormMax(mvpError)/bfVecNormMax(yColNorms);
    bfVecDelete(&mvpError);
    bfVecDelete(&yColNorms);
    printf("- relative error between dense and FMM MVPs: %g\n", mvpRelErrorMax);
  }

  if (yTestButterfly && yTestFmm) {
    BfVec *mvpError = bfMatColDists(yTestButterfly, yTestFmm);
    BfVec *yColNorms = bfMatColNorms(yTestButterfly);
    BfReal mvpRelErrorMax = bfVecNormMax(mvpError)/bfVecNormMax(yColNorms);
    bfVecDelete(&mvpError);
    bfVecDelete(&yColNorms);
    printf("- relative error between butterfly and FMM MVPs: %g\n", mvpRelErrorMax);
  }

  /* All the different ways we've computed sigma: */
  BfMat const *sigma[] = {
    context->sigmaLu,
    context->sigmaDense,
    context->sigmaDensePrecondLeft,
    context->sigmaFmm,
    context->sigmaFmmPrecondLeft,
    context->sigmaButterfly,
    context->sigmaButterflyPrecondLeft,
  };

  /* Names for each method: */
  char const *methodName[] = {
    "LU",
    "dense (GMRES)",
    "dense (preconditioned GMRES)",
    "FMM (GMRES)",
    "FMM (preconditioned GMRES)",
    "butterfly (GMRES)",
    "butterfly (preconditioned GMRES)"
  };

  BfSize numMethods = sizeof(sigma)/sizeof(sigma[0]);
  BF_ASSERT(numMethods == sizeof(methodName)/sizeof(methodName[0]));

  for (BfSize i = 0; i < numMethods; ++i) {
    if (sigma[i] == NULL) continue;
    for (BfSize j = i + 1; j < numMethods; ++j) {
      if (sigma[j] == NULL) continue;
      BfVecReal *l2Dist = bfVecToVecReal(bfMatColDists(sigma[i], sigma[j]));
      BfVecReal *l2Norm_i = bfVecToVecReal(bfMatColNorms(sigma[i]));
      BfVecReal *l2Norm_j = bfVecToVecReal(bfMatColNorms(sigma[j]));
      BfReal l2ErrorRel = l2Dist->data[0]/fmax(l2Norm_i->data[0], l2Norm_j->data[0]);
      printf("- rel l2 error in sigma (%s vs %s): %g\n", methodName[i], methodName[j], l2ErrorRel);
      bfVecRealDeinitAndDealloc(&l2Dist);
      bfVecRealDeinitAndDealloc(&l2Norm_i);
      bfVecRealDeinitAndDealloc(&l2Norm_j);
    }
  }
}

void doPostprocessing(MultipleScatteringContext *context) {
  puts("Doing postprocessing:");

  /** Set up evaluation grid: */

  bfToc();

  context->XEval = bfPoints2NewGrid(&context->bboxEval, context->nxEval, context->nyEval);

  printf("- set up %lu x %lu evaluation grid [%0.2fs]\n", context->nxEval, context->nyEval, bfToc());

  bfSavePoints2(context->XEval, "XEval.bin");

  /** Build quadtree on evaluation grid: */

  context->quadtreeEval = bfQuadtreeNew();
  bfQuadtreeInit(context->quadtreeEval, context->XEval, NULL);

  BfTree const *treeEval = bfQuadtreeConstToTreeConst(context->quadtreeEval);
  context->permEval = bfTreeGetPermConst(treeEval);

  /** Compute uIn on evaluation grid: */

  bfToc();

  context->uIn = NULL;
  { BfSize nEval = context->nxEval*context->nyEval;
    BfMatDenseComplex *_ = bfMatDenseComplexNew();
    bfMatDenseComplexInit(_, nEval, 1);
    for (BfSize i = 0; i < nEval; ++i) {
      BfPoint2 x;
      bfPoints2Get(context->XEval, i, x);
      *(_->data + i*_->rowStride) = cexp(I*context->k*(context->d[0]*x[0] + context->d[1]*x[1]));
    }
    context->uIn = bfMatDenseComplexToMat(_); }

  printf("- computed uIn on evaluation grid [%0.2fs]\n", bfToc());

  bfMatSave(context->uIn, "uIn.bin");

  /** Sample total field for different solution methods: */

  bfToc();

  context->KEvalButterfly = bfFacHelm2MakeMultilevel(
    context->quadtree,
    context->quadtreeEval,
    context->k,
    BF_LAYER_POTENTIAL_COMBINED_FIELD,
    &context->alpha,
    &context->beta);

  BfVec *WVec = bfRealArrayGetVecView(context->W);
  BfVec *WVecPerm = bfVecCopy(WVec);
  bfVecPermute(WVecPerm, context->revPerm);
  bfMatScaleCols(context->KEvalButterfly, WVecPerm);
  bfVecDelete(&WVecPerm);
  bfVecDelete(&WVec);

  printf("- set up butterfly factorized evaluation matrix [%0.2fs]\n", bfToc());

  /* All the different ways we've computed sigma: */
  BfMat const *sigma[] = {
    context->sigmaLu,
    context->sigmaDense,
    context->sigmaDensePrecondLeft,
    context->sigmaFmm,
    context->sigmaFmmPrecondLeft,
    context->sigmaButterfly,
    context->sigmaButterflyPrecondLeft,
  };

  /* Names for each method: */
  char const *methodName[] = {
    "LU",
    "dense (GMRES)",
    "dense (preconditioned GMRES)",
    "FMM (GMRES)",
    "FMM (preconditioned GMRES)",
    "butterfly (GMRES)",
    "butterfly (preconditioned GMRES)"
  };

  /* Filenames for each method: */
  char const *filename[] = {
    "uScatLu.bin",
    "uScatDense.bin",
    "uScatDensePrecond.bin",
    "uScatFmm.bin",
    "uScatFmmPrecond.bin",
    "uScatButterfly.bin",
    "uScatButterflyPrecond.bin",
  };

  BfSize numMethods = sizeof(sigma)/sizeof(sigma[0]);
  BF_ASSERT(numMethods == sizeof(methodName)/sizeof(methodName[0]));

  for (BfSize i = 0; i < numMethods; ++i) {
    if (sigma[i] == NULL) continue;

    bfToc();

    BfMat *sigmaPerm = bfMatCopy(sigma[i]);
    bfMatPermuteRows(sigmaPerm, context->revPerm);

    BfMat *uScat = bfMatMul(context->KEvalButterfly, sigmaPerm);
    bfMatPermuteRows(uScat, context->permEval);

    printf("- computed %s scattered field [%0.2fs]\n", methodName[i], bfToc());

    bfMatSave(uScat, filename[i]);

    bfMatDelete(&sigmaPerm);
    bfMatDelete(&uScat);
  }
}

void deinit(MultipleScatteringContext *context) {
  if (context->ellipseCenters != NULL) {
    bfFreePoints2(context->ellipseCenters);
    bfMemFree(context->ellipseCenters);
  }

  if (context->ellipse != NULL)
    bfMemFree(context->ellipse);

  if (context->ellipseOffsets != NULL)
    bfSizeArrayDeinitAndDealloc(&context->ellipseOffsets);

  if (context->X) {
    bfFreePoints2(context->X);
    bfMemFree(context->X);
  }

  if (context->N) {
    bfVectors2DeinitAndDealloc(&context->N);
    bfMemFree(context->N);
  }

  if (context->W)
    bfRealArrayDeinitAndDealloc(&context->W);

  if (context->quadtree) {
    bfQuadtreeDeinit(context->quadtree);
    bfMemFree(context->quadtree);
  }

  if (context->revPerm)
    bfPermDelete(&context->revPerm);

  if (context->K != NULL)
    bfMatDelete(&context->K);

  if (context->KButterfly != NULL)
    bfMatDelete(&context->KButterfly);

  if (context->KButterflyDense != NULL)
    bfMatDelete(&context->KButterflyDense);

  if (context->KFmm != NULL)
    bfMatDelete(&context->KFmm);

  if (context->KLu != NULL)
    bfLuDelete(&context->KLu);

  if (context->KBlockLus != NULL) {
    for (BfSize i = 0; i < bfPtrArraySize(context->KBlockLus); ++i) {
      BfLu *KBlockLu = bfPtrArrayGet(context->KBlockLus, i);
      bfLuDelete(&KBlockLu);
    }
    bfPtrArrayDelete(&context->KBlockLus);
  }

  if (context->M != NULL)
    bfMatDelete(&context->M);

  if (context->MPerm != NULL)
    bfMatDelete(&context->MPerm);

  if (context->rhs != NULL)
    bfMatDelete(&context->rhs);

  if (context->rhsPerm != NULL)
    bfMatDelete(&context->rhsPerm);

  if (context->yTestDense != NULL)
    bfMatDelete(&context->yTestDense);

  if (context->yTestButterfly != NULL)
    bfMatDelete(&context->yTestButterfly);

  if (context->sigmaLu != NULL)
    bfMatDelete(&context->sigmaLu);

  if (context->sigmaDense != NULL)
    bfMatDelete(&context->sigmaDense);

  if (context->sigmaDensePrecondLeft != NULL)
    bfMatDelete(&context->sigmaDensePrecondLeft);

  if (context->sigmaFmm != NULL)
    bfMatDelete(&context->sigmaFmm);

  if (context->sigmaFmmPrecondLeft != NULL)
    bfMatDelete(&context->sigmaFmmPrecondLeft);

  if (context->sigmaButterfly != NULL)
    bfMatDelete(&context->sigmaButterfly);

  if (context->sigmaButterflyPrecondLeft != NULL)
    bfMatDelete(&context->sigmaButterflyPrecondLeft);

  if (context->XEval != NULL) {
    bfFreePoints2(context->XEval);
    bfMemFree(context->XEval);
  }

  if (context->quadtreeEval != NULL) {
    bfQuadtreeDeinit(context->quadtreeEval);
    bfMemFree(context->quadtreeEval);
  }

  if (context->KEvalButterfly != NULL)
    bfMatDelete(&context->KEvalButterfly);

  if (context->uIn != NULL)
    bfMatDelete(&context->uIn);
}
