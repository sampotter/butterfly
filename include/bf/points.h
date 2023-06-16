#pragma once

#include "geom.h"
#include "mat.h"

BfReal bfPoint2Dist(BfPoint2 const p, BfPoint2 const q);
BfReal bfPoint2Magnitude(BfPoint2 const p);
void bfPoint2SampleUniformlyFromBoundingBox(BfBbox2 const *bbox, BfPoint2 p);
void bfPoint2RotateAboutOrigin(BfPoint2 p, BfReal theta);
void bfPoint2Translate(BfPoint2 p, BfVector2 const u);

BfReal bfPoint3Dist(BfPoint3 const p, BfPoint3 const q);
void bfPoint3Sub(BfPoint3 const v, BfPoint3 const u, BfVector3 uv);
void bfPoint3GetPointOnRay(BfPoint3 const r0, BfVector3 const dr, BfReal t, BfPoint3 rt);
void bfPoint3Copy(BfPoint3 x, BfPoint3 const y);
bool bfPoint3Equal(BfPoint3 const x, BfPoint3 const y);
bool bfPoint3EqualApprox(BfPoint3 const x, BfPoint3 const y, BfReal tol);

struct BfPoints1 {
  BfPoint1 *data;
  BfSize size;
  BfSize capacity;
  bool isView;
};

BfPoints1 *bfPoints1New();
BfPoints1 *bfPoints1Copy();
void bfPoints1InitEmpty(BfPoints1 *points, BfSize capacity);
void bfPoints1InitViewFromVecReal(BfPoints1 *points, BfVecReal const *vecReal);
void bfPoints1Deinit(BfPoints1 *points);
void bfPoints1Dealloc(BfPoints1 **points);
void bfPoints1Delete(BfPoints1 **points);
bool bfPoints1IsSorted(BfPoints1 const *points);
void bfPoints1Append(BfPoints1 *points, BfPoint1 point);
void bfPoints1InsertPointsSorted(BfPoints1 *points, BfPoints1 const *newPoints);
void bfPoints1Map(BfPoints1 *points, BfReal (*func)(BfReal));
void bfPoints1Save(BfPoints1 const *points, char const *path);

struct BfPoints2 {
  BfPoint2 *data;
  BfSize size;
  BfSize capacity;
  bool isView;
};

BfPoints2 *bfPoints2NewEmpty();
BfPoints2 *bfPoints2NewGrid(BfBbox2 const *bbox, BfSize nx, BfSize ny);
BfPoints2 const *bfPoints2ConstViewFromMat(BfMat const *mat);
BfPoints2 const *bfPoints2ConstViewFromMatDenseReal(BfMatDenseReal const *matDenseReal);
BfPoints2 bfGetUninitializedPoints2();
void bfInitEmptyPoints2(BfPoints2 *points, BfSize numPoints);
void bfReadPoints2FromFile(char const *path, BfPoints2 *points);
void bfFreePoints2(BfPoints2 *points);
void bfPoints2Dealloc(BfPoints2 **points);
void bfPoints2DeinitAndDealloc(BfPoints2 **points);
bool bfPoints2Initialized(BfPoints2 const *points);
BfBbox2 bfPoints2GetBoundingBox(BfPoints2 const *points);
void bfGetPointsByIndex(BfPoints2 const *points, BfSize numInds, BfSize const *inds, BfPoints2 *indexedPoints);
void bfPrintPoints2(BfPoints2 const *points);
void bfSavePoints2(BfPoints2 const *points, char const *path);
BfReal *bfPoints2PairwiseDists(BfPoints2 const *X, BfPoints2 const *Y);
void bfPoints2Append(BfPoints2 *points, BfPoint2 const p);
void bfPoints2Extend(BfPoints2 *points, BfPoints2 const *newPoints);
void bfPoints2Get(BfPoints2 const *points, BfSize i, BfPoint2 p);
BfSize bfPoints2GetSize(BfPoints2 const *points);
BfPoints2 *bfPoints2GetRangeView(BfPoints2 *points, BfSize i0, BfSize i1);

struct BfPoints3 {
  BfPoint3 *data;
  BfSize size;
  BfSize capacity;
  bool isView;
};

BfPoints3 *bfPoints3NewFromBinaryFile(char const *path);
BfPoints3 *bfPoints3NewWithDefaultCapacity();
BfPoints3 *bfPoints3Copy(BfPoints3 const *points);
void bfPoints3InitWithDefaultCapacity(BfPoints3 *points);
void bfPoints3Deinit(BfPoints3 *points);
void bfPoints3Dealloc(BfPoints3 **points);
void bfPoints3DeinitAndDealloc(BfPoints3 **points);
BfBoundingBox3 bfPoints3GetBoundingBox(BfPoints3 const *points);
void bfPoints3GetByIndex(BfPoints3 const *points, BfSize numInds, BfSize const *inds, BfPoints3 *indexedPoints);
BfReal const *bfPoints3GetPtrConst(BfPoints3 const *points, BfSize i);
void bfPoints3Append(BfPoints3 *points, BfPoint3 const point);
void bfPoints3Extend(BfPoints3 *points, BfPoints3 const *newPoints);
bool bfPoints3Contains(BfPoints3 const *points, BfPoint3 const point);
BfSize bfPoints3GetSize(BfPoints3 const *points);
BfSize bfPoints3Find(BfPoints3 const *points, BfPoint3 const point);
void bfPoints3Set(BfPoints3 *points, BfSize i, BfPoint3 const point);
void bfPoints3Delete(BfPoints3 *points, BfSize i);
bool bfPoints3AllUnique(BfPoints3 const *points);
void bfPoints3Save(BfPoints3 const *points, char const *path);
void bfPoints3ContainsApprox(BfPoint3 const *points, BfPoint3 const point, BfReal tol);
BfSize bfPoints3FindApprox(BfPoints3 const *points, BfPoint3 const point, BfReal tol);
