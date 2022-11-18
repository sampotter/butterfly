#pragma once

#include <bf/geom.h>

typedef struct BfTrimesh {
  BfPoint3 *verts;
  BfSize numVerts;
  BfSize3 *faces;
  BfSize numFaces;
} BfTrimesh;

BfTrimesh *bfTrimeshNew();
void bfTrimeshInitFromBinaryFiles(BfTrimesh *trimesh, char const *vertsPath,
                                  char const *facesPath);
