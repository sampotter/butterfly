#include <bf/trimesh.h>

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/points.h>
#include <bf/util.h>

static void initVf(BfTrimesh *trimesh) {
  BEGIN_ERROR_HANDLING();

  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);
  BfSize numFaces = bfTrimeshGetNumFaces(trimesh);

  trimesh->vfOffset = malloc((numVerts + 1)*sizeof(BfSize));
  if (trimesh->vfOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  trimesh->vfOffset[0] = 0;

  for (BfSize i = 0; i < numFaces; ++i)
    for (BfSize j = 0; j < 3; ++j)
      ++trimesh->vfOffset[trimesh->faces[i][j] + 1];

  for (BfSize i = 0; i < numVerts; ++i)
    trimesh->vfOffset[i + 1] += trimesh->vfOffset[i];

  trimesh->vf = malloc(trimesh->vfOffset[numVerts]*sizeof(BfSize));
  if (trimesh->vf == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfSize *offset = calloc(numVerts, sizeof(BfSize));
  if (offset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < numFaces; ++i) {
    for (BfSize j = 0; j < 3; ++j) {
      BfSize k = trimesh->faces[i][j];
      trimesh->vf[trimesh->vfOffset[k] + offset[k]++] = i;
    }
  }

  END_ERROR_HANDLING() {
    free(trimesh->vfOffset);
    trimesh->vfOffset = NULL;

    free(trimesh->vf);
    trimesh->vf = NULL;
  }

  free(offset);
}

/* Insert a new entry into `trimesh->vv`. This assumes that the `vv`
 * blocks are being filled in increasing vertex order. */
static void vvInsertNewSorted(BfTrimesh *trimesh, BfSize i, BfSize i0, BfSize *vvCapacity) {
  BEGIN_ERROR_HANDLING();

  /* Find position of `i0` in `vv` */
  BfSize j = trimesh->vvOffset[i];
  for (; j < trimesh->vvOffset[i + 1]; ++j)
    if (trimesh->vv[j] >= i0)
      break;
  assert (j <= *vvCapacity);

  /* If `vv` already contains `i0`, return */
  if (trimesh->vv[j] == i0)
    return;

  /* If we're at capacity, expand the array */
  if (j == *vvCapacity) {
    *vvCapacity *= 2;
    BfSize *vvNew = realloc(trimesh->vv, *vvCapacity*sizeof(BfSize));
    if (vvNew == NULL)
      RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
    trimesh->vv = vvNew;
  }

  /* Move the rest of the current `vv` block over one element to make
   * room for `i0` (no-op if `j` is at the end of the block) */
  memmove(&trimesh->vv[j + 1], &trimesh->vv[j],
          (trimesh->vvOffset[i + 1] - j)*sizeof(BfSize));

  /* Insert `i0` and update the offset to the next block */
  trimesh->vv[j] = i0;
  ++trimesh->vvOffset[i + 1];

  END_ERROR_HANDLING() {}
}

static void initVv(BfTrimesh *trimesh) {
  BEGIN_ERROR_HANDLING();

  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);

  trimesh->vvOffset = malloc((numVerts + 1)*sizeof(BfSize));
  if (trimesh->vvOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  trimesh->vvOffset[0] = 0;

  BfSize vvCapacity = numVerts;

  trimesh->vv = malloc(vvCapacity*sizeof(BfSize));
  if (trimesh->vv == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < numVerts; ++i) {
    for (BfSize j = trimesh->vfOffset[i]; j < trimesh->vfOffset[i + 1]; ++j) {
      BfSize i0, i1;
      bfTrimeshGetOpFaceVerts(trimesh, trimesh->vf[j], i, &i0, &i1);
      vvInsertNewSorted(trimesh, i, i0, &vvCapacity);
      vvInsertNewSorted(trimesh, i, i1, &vvCapacity);
    }
  }

  // resize to fit

  END_ERROR_HANDLING() {
  }
}

void bfTrimeshInitFromBinaryFiles(BfTrimesh *trimesh,
                                  char const *vertsPath,
                                  char const *facesPath) {
  BEGIN_ERROR_HANDLING();

  bfPoints3InitFromBinaryFile(&trimesh->verts, vertsPath);

  /* Get the face count */

  BfSize numBytes = bfGetFileSizeInBytes(facesPath);
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

  initVf(trimesh);
  HANDLE_ERROR();

  initVv(trimesh);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}
}

void bfTrimeshInitFromObjFile(BfTrimesh *trimesh, char const *objPath) {
  BEGIN_ERROR_HANDLING();

  FILE *fp;
  char *lineptr;
  ssize_t nread;
  size_t n;
  char *saveptr;
  char *tok;

  size_t num_verts = 0;
  size_t num_vert_normals = 0;
  size_t num_face_vertex_indices = 0;
  size_t num_face_vertex_normal_indices = 0;

  fp = fopen(objPath, "r");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* Make first pass over OBJ file to count number of elements of each
   * type and make sure the file isn't malformed */

  do {
    lineptr = NULL;
    nread = getline(&lineptr, &n, fp);

    saveptr = NULL;
    tok = strtok_r(lineptr, " ", &saveptr);

    if (!strcmp(tok, "v")) ++num_verts;
    if (!strcmp(tok, "vn")) ++num_vert_normals;

    if (!strcmp(tok, "f")) {
      ++num_face_vertex_indices;
      ++num_face_vertex_normal_indices;
    }

    free(lineptr);
  } while (nread >= 0);

  /* Allocate space for the vertices and faces */

  bfPoints3InitEmpty(&trimesh->verts, num_verts);
  HANDLE_ERROR();

  trimesh->numFaces = num_face_vertex_indices;

  trimesh->faces = malloc(trimesh->numFaces*sizeof(BfSize3));
  if (trimesh->faces == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* Make second pass over file and parse data */

  rewind(fp);

  size_t v_index = 0;
//   size_t vn_index = 0;
  size_t fv_index = 0;
//   size_t fvn_index = 0;

  do {
    lineptr = NULL;
    nread = getline(&lineptr, &n, fp);

    saveptr = NULL;
    tok = strtok_r(lineptr, " ", &saveptr);

    if (!strcmp(tok, "v")) {
      for (size_t i = 0; i < 3; ++i) {
        tok = strtok_r(NULL, " ", &saveptr);
        trimesh->verts.data[v_index][i] = strtod(tok, NULL);
      }
      ++v_index;
    }

//     if (!strcmp(tok, "vn")) {
//       for (size_t i = 0; i < 3; ++i) {
//         tok = strtok_r(NULL, " ", &saveptr);
//         vert_normals[vn_index][i] = strtod(tok, NULL);
//       }
//       ++vn_index;
//     }

    if (!strcmp(tok, "f")) {
      if (num_vert_normals > 0) {
        char *toks[3], *saveptr_, *tok;
        for (size_t i = 0; i < 3; ++i) {
          toks[i] = strtok_r(NULL, " ", &saveptr);
        }
        for (size_t i = 0; i < 3; ++i) {
          saveptr_ = NULL;
          tok = strtok_r(toks[i], "//", &saveptr_);
          trimesh->faces[fv_index][i] = strtoull(tok, NULL, 10) - 1;
//           tok = strtok_r(NULL, "//", &saveptr_);
//           face_vertex_normal_indices[fvn_index][i] = strtoull(tok, NULL, 10) - 1;
        }
        ++fv_index;
//         ++fvn_index;
      } else {
        assert(false);
      }
    }

    free(lineptr);
  } while (nread >= 0);

  initVf(trimesh);
  HANDLE_ERROR();

  initVv(trimesh);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfTrimeshDeinit(trimesh);

  fclose(fp);
}

void bfTrimeshDeinit(BfTrimesh *trimesh) {
  (void)trimesh;
  assert(false);
}

BfSize bfTrimeshGetNumVerts(BfTrimesh const *trimesh) {
  return trimesh->verts.size;
}

BfSize bfTrimeshGetNumFaces(BfTrimesh const *trimesh) {
  return trimesh->numFaces;
}

void bfTrimeshGetVertex(BfTrimesh const *trimesh, BfSize i, BfPoint3 x) {
  memcpy(x, trimesh->verts.data[i], sizeof(BfPoint3));
}

void bfTrimeshGetOpFaceVerts(BfTrimesh const *trimesh, BfSize faceIndex,
                             BfSize i, BfSize *i0, BfSize *i1) {
  BfSize *F = trimesh->faces[faceIndex];
  assert(F[0] == i || F[1] == i || F[2] == i);
  BfSize k = 0;
  while (F[k] != i) ++k;
  *i0 = F[k];
  while (F[k] != i) ++k;
  *i1 = F[k];
}
