#pragma once

#include "geom.h"
#include "types.h"

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
BfPoints3 const *bfTrimeshGetVertsConst(BfTrimesh const *trimesh);
void bfTrimeshGetVertex(BfTrimesh const *trimesh, BfSize i, BfPoint3 x);
BfReal const *bfTrimeshGetVertPtrConst(BfTrimesh const *trimesh, BfSize i);
void bfTrimeshGetOpFaceVerts(BfTrimesh const *trimesh, BfSize faceIndex, BfSize i, BfSize *i0, BfSize *i1);
BfTrimesh *bfTrimeshGetLevelSetSubmesh(BfTrimesh const *trimesh, BfRealArray const *phi, BfReal tol, bool const *permMask, BfSizeArray **permPtr);
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
BfSize bfTrimeshGetNumEdges(BfTrimesh const *trimesh);
BfSizeArray *bfTrimeshGetVertNbs(BfTrimesh const *trimesh, BfSize i);
bool bfTrimeshHasDuplicateFaces(BfTrimesh const *trimesh);
BfSize bfTrimeshGetNumVertexNeighbors(BfTrimesh const *trimesh, BfSize i);
BfSize bfTrimeshGetVertexNeighbor(BfTrimesh const *trimesh, BfSize i, BfSize j);
bool bfTrimeshIsBoundaryVertex(BfTrimesh const *trimesh, BfSize i);
BfSize bfTrimeshGetNumBoundaryEdges(BfTrimesh const *trimesh);
void bfTrimeshGetBoundaryEdge(BfTrimesh const *trimesh, BfSize i, BfSize2 boundaryEdge);
BfSize const *bfTrimeshGetBoundaryEdgeConstPtr(BfTrimesh const *trimesh, BfSize i);
BfSize const *bfTrimeshGetFaceConstPtr(BfTrimesh const *trimesh, BfSize i);
void bfTrimeshGetLboFemDiscretization(BfTrimesh const *trimesh, BfMat **L, BfMat **M);
#ifdef BF_EMBREE
BfSizeArray *bfTrimeshGetVisibility(BfTrimesh const *trimesh, BfSize srcInd, BfSizeArray const *tgtInds);
#endif
BfReal const *bfTrimeshGetFaceCentroidConstPtr(BfTrimesh const *trimesh, BfSize i);
BfReal const *bfTrimeshGetFaceUnitNormalConstPtr(BfTrimesh const *trimesh, BfSize i);
BfReal bfTrimeshGetFaceArea(BfTrimesh const *trimesh, BfSize i);
bool bfTrimeshHasFaceNormals(BfTrimesh const *trimesh);
bool bfTrimeshHasVertexNormals(BfTrimesh const *trimesh);
BfVectors3 *bfTrimeshGetVertexNormalsPtr(BfTrimesh *trimesh);
BfVectors3 *bfTrimeshGetFaceNormalsPtr(BfTrimesh *trimesh);
void bfTrimeshComputeFaceNormalsMatchingVertexNormals(BfTrimesh *trimesh);
