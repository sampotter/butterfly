#include <bf/mat_csr_real.h>
#include <bf/size_array.h>
#include <bf/trimesh.h>
#include <bf/util.h>

#include <stdlib.h>

int main(void) {
  bfToc();

  BfTrimesh *trimesh = bfTrimeshNewFromObjFile("67p.obj");
  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);
  BfSize numFaces = bfTrimeshGetNumFaces(trimesh);

  if (bfTrimeshHasVertexNormals(trimesh) && !bfTrimeshHasFaceNormals(trimesh))
    bfTrimeshComputeFaceNormalsMatchingVertexNormals(trimesh);

  printf("loaded mesh with %lu verts and %lu faces [%0.2fs]\n", numVerts, numFaces, bfToc());

  BfSizeArray *rowInds = bfSizeArrayNewIota(numFaces);
  BfSizeArray *colInds = bfSizeArrayNewIota(numFaces);
  BfMatCsrReal *mat = bfMatCsrRealNewViewFactorMatrixFromTrimesh(trimesh, rowInds, colInds);

  printf("computed view factor matrix [%0.2fs]\n", bfToc());

  bfMatCsrRealDeinitAndDealloc(&mat);

  return EXIT_SUCCESS;
}
