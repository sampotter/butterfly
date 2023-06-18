#pragma once

#include <bf/geom.h>
#include <bf/points.h>

typedef struct BfTrimesh {
  BfPoints3 *verts;

  BfSize3 *faces;
  BfSize numFaces;

  BfSize2 *edges;
  BfSize numEdges;

  BfSize *vfOffset, *vf;
  BfSize *vvOffset, *vv;

  bool *isBoundaryEdge;
  bool *isBoundaryVert;

  BfSize2 *boundaryEdges;
  BfSize numBoundaryEdges;
} BfTrimesh;

BfTrimesh *bfTrimeshCopy(BfTrimesh const *trimesh);
BfTrimesh *bfTrimeshNewFromObjFile(char const *objPath);
BfTrimesh *bfTrimeshNewFromVertsAndFaces(BfPoints3 const *verts, BfSize numFaces, BfSize3 const *faces);
void bfTrimeshInitFromBinaryFiles(BfTrimesh *trimesh, char const *vertsPath, char const *facesPath);
void bfTrimeshInitFromObjFile(BfTrimesh *trimesh, char const *objPath);
void bfTrimeshInitFromVertsAndFaces(BfTrimesh *trimesh, BfPoints3 const *verts, BfSize numFaces, BfSize3 const *faces);
void bfTrimeshDeinit(BfTrimesh *trimesh);
void bfTrimeshDealloc(BfTrimesh **trimesh);
void bfTrimeshDeinitAndDealloc(BfTrimesh **trimesh);
BfSize bfTrimeshGetNumVerts(BfTrimesh const *trimesh);
BfSize bfTrimeshGetNumFaces(BfTrimesh const *trimesh);
void bfTrimeshGetVertex(BfTrimesh const *trimesh, BfSize i, BfPoint3 x);
BfReal const *bfTrimeshGetVertPtrConst(BfTrimesh const *trimesh, BfSize i);
void bfTrimeshGetOpFaceVerts(BfTrimesh const *trimesh, BfSize faceIndex, BfSize i, BfSize *i0, BfSize *i1);
BfTrimesh *bfTrimeshGetLevelSetSubmesh(BfTrimesh const *trimesh, BfVecReal const *phi, bool const *permMask, BfSizeArray **permPtr);
BfSizeArray *bfTrimeshGetInteriorInds(BfTrimesh const *trimesh);
BfSize bfTrimeshGetBoundaryFaceIndexByBoundaryEdge(BfTrimesh const *trimesh, BfSize2 const boundaryEdge);
void bfTrimeshSplitEdge(BfTrimesh *trimesh, BfSize i, BfReal lam);
BfSize bfTrimeshGetEdgeIndex(BfTrimesh const *trimesh, BfSize2 const edge);
BfSize bfTrimeshGetFacesIncOnEdge(BfTrimesh const *trimesh, BfSize edgeInd, BfSize2 incFaceInds);
bool bfTrimeshIsOK(BfTrimesh const *trimesh);
void bfTrimeshDumpVerts(BfTrimesh const *trimesh, char const *path);
void bfTrimeshDumpFaces(BfTrimesh const *trimesh, char const *path);
BfRealArray *bfTrimeshGetNormalDeriv(BfTrimesh const *trimesh, BfRealArray const *values, bool computeAtAllVerts);
BfRealArray *bfTrimeshGetFiedler(BfTrimesh const *trimesh);
