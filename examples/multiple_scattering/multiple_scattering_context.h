#pragma once

#include "multiple_scattering_opts.h"

#include <bf/quadtree.h>
// #include <bf/size_array.h>

typedef struct KHelm2Wkspc KHelm2Wkspc;
typedef struct MultipleScatteringContext MultipleScatteringContext;

struct KHelm2Wkspc {
  BfPoints2 const *points;
  BfVectors2 const *normals;
  BfReal k;
  BfComplex alpha;
  BfComplex beta;
};

struct MultipleScatteringContext {
  /** These parameters are determined by the user: */

  /* Problem parameters: */
  BfReal k;
  BfReal minDist;
  BfReal axisLow;
  BfReal axisHigh;
  BfComplex alpha;
  BfComplex beta;
  BfVector2 d;

  /* Discretization parameters: */
  BfReal h;
  BfReal tol;
  BfReal tolFmm;
  BfSize orderKR;

  /* Postprocessing parameters: */
  BfBbox2 bboxEval;
  BfSize nxEval;
  BfSize nyEval;

  /** Everything else is derived from the preceding parameters: */

  /* Set by `setUpGeometry`: */
  BfReal R;
  BfBbox2 ellipseBbox;
  BfPoints2 *ellipseCenters;
  BfSize numEllipses;
  BfEllipse *ellipse;
  BfSizeArray *ellipseOffsets;

  /* Set by `setUpDiscretization`: */
  BfPoints2 *X;
  BfVectors2 *N;
  BfRealArray *W;
  BfSize n;

  /* Set by `buildQuadtrees`: */
  BfQuadtree *quadtree;
  BfPerm const *perm;
  BfPerm *revPerm;

  /* Set by `setUpKrWorkspace`: */
  KHelm2Wkspc KWkspc;

  BfMat *K;
  BfMat *KButterfly;
  BfMat *KButterflyDense;
  BfMat *KFmm;
  BfLu *KLu;

  BfMat *M;
  BfMat *MPerm;

  BfMat *rhs;
  BfMat *rhsPerm;
  BfMat *yTestDense;
  BfMat *yTestButterfly;

  BfSize numIterDense;
  BfSize numIterDensePrecondLeft;
  BfSize numIterFmm;
  BfSize numIterFmmPrecondLeft;
  BfSize numIterButterfly;
  BfSize numIterButterflyPrecondLeft;

  BfMat *sigmaLu;
  BfMat *sigmaDense;
  BfMat *sigmaDensePrecondLeft;
  BfMat *sigmaFmm;
  BfMat *sigmaFmmPrecondLeft;
  BfMat *sigmaButterfly;
  BfMat *sigmaButterflyPrecondLeft;

  BfPoints2 *XEval;

  BfQuadtree *quadtreeEval;
  BfPerm const *permEval;

  BfMat *KEvalButterfly;

  BfMat *uIn;
};

void init(MultipleScatteringContext *context, Opts const *opts);
void deinit(MultipleScatteringContext *context, Opts *opts);

void printInfo(MultipleScatteringContext const *context);
void setUpGeometry(MultipleScatteringContext *context);
void setUpDiscretization(MultipleScatteringContext *context);
void buildQuadtrees(MultipleScatteringContext *context);
void setUpKrWorkspace(MultipleScatteringContext *context);
void setUpRhs(MultipleScatteringContext *context);
void assembleDenseK(MultipleScatteringContext *context);
void assembleButterfliedK(MultipleScatteringContext *context);
void extractDenseButterfliedK(MultipleScatteringContext *context);
void assembleFmmK(MultipleScatteringContext *context);
void estimateFmmTol(MultipleScatteringContext *context);
void assemblePreconditioner(MultipleScatteringContext *context);
void computeLu(MultipleScatteringContext *context);
void solveLu(MultipleScatteringContext *context);
void solveDenseGmres(MultipleScatteringContext *context);
void solveDensePreconditionedGmres(MultipleScatteringContext *context);
void solveButterflyGmres(MultipleScatteringContext *context);
void solveButterflyPreconditionedGmres(MultipleScatteringContext *context);
void solveFmmGmres(MultipleScatteringContext *context);
void solveFmmPreconditionedGmres(MultipleScatteringContext *context);
void collectAndPrintStats(MultipleScatteringContext *context);
void doPostprocessing(MultipleScatteringContext *context);
