#pragma once

#include <bf/geom.h>
#include <bf/points.h>

typedef struct BfTrimesh {
  BfPoints3 verts;

  BfSize3 *faces;
  BfSize numFaces;

  BfSize *vfOffset, *vf;
  BfSize *vvOffset, *vv;
} BfTrimesh;

void bfTrimeshInitFromBinaryFiles(BfTrimesh *trimesh, char const *vertsPath,
                                  char const *facesPath);
void bfTrimeshInitFromObjFile(BfTrimesh *trimesh, char const *objPath);
void bfTrimeshDeinit(BfTrimesh *trimesh);
BfSize bfTrimeshGetNumVerts(BfTrimesh const *trimesh);
BfSize bfTrimeshGetNumFaces(BfTrimesh const *trimesh);
void bfTrimeshGetVertex(BfTrimesh const *trimesh, BfSize i, BfPoint3 x);
void bfTrimeshGetOpFaceVerts(BfTrimesh const *trimesh, BfSize faceIndex, BfSize i, BfSize *i0, BfSize *i1);
