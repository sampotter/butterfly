#include <bf/trimesh.h>

#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/util.h>

BfTrimesh *bfTrimeshNew() {
  BEGIN_ERROR_HANDLING();

  BfTrimesh *trimesh = malloc(sizeof(BfTrimesh));
  if (trimesh == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  trimesh->verts = NULL;
  trimesh->numVerts = BF_SIZE_BAD_VALUE;

  trimesh->faces = NULL;
  trimesh->numFaces = BF_SIZE_BAD_VALUE;

  END_ERROR_HANDLING() {
    free(trimesh);
    trimesh = NULL;
  }

  return trimesh;
}

void bfTrimeshInitFromBinaryFiles(BfTrimesh *trimesh,
                                  char const *vertsPath,
                                  char const *facesPath) {
  BEGIN_ERROR_HANDLING();

  BfSize numBytes = BF_SIZE_BAD_VALUE;

  if (trimesh->verts != NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (trimesh->faces != NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Read in the number of vertices based on the size of the binary
   * file containing the vertex data */

  numBytes = bfGetFileSizeInBytes(vertsPath);
  HANDLE_ERROR();

  if (numBytes % sizeof(BfPoint3) != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
  trimesh->numVerts = numBytes/sizeof(BfPoint3);

  /* Allocate memory for the vertices and read them in */

  trimesh->verts = malloc(trimesh->numVerts*sizeof(BfPoint3));
  if (trimesh->verts == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  bfReadFileToMemory(vertsPath, numBytes, (BfByte *)trimesh->verts);
  HANDLE_ERROR();

  /* Get the face count */

  numBytes = bfGetFileSizeInBytes(facesPath);
  HANDLE_ERROR();

  if (numBytes % sizeof(BfSize3) != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
  trimesh->numFaces = numBytes/sizeof(BfSize3);

  /* Allocate space and load the faces */

  trimesh->faces = malloc(trimesh->numFaces*sizeof(BfSize3));
  if (trimesh->faces == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  bfReadFileToMemory(facesPath, numBytes, (BfByte *)trimesh->faces);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}
}
