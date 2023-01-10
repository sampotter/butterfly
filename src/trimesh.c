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

  trimesh->vfOffset = calloc(numVerts + 1, sizeof(BfSize));
  if (trimesh->vfOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* Count the number of faces incident on each vertex by making a
   * pass over the vertex array. Accumulate the count into the
   * *following* entry of vfOffset, since we will do a cumulative sum
   * below to get the offsets. */
  for (BfSize i = 0; i < numFaces; ++i)
    for (BfSize j = 0; j < 3; ++j)
      ++trimesh->vfOffset[trimesh->faces[i][j] + 1];

  /* Raise an error if there are any isolated vertices. */
  for (BfSize i = 0; i < numVerts; ++i)
    if (trimesh->vfOffset[i + 1] == 0)
      RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* Cumsum vfOffset to convert counts to offsets. */
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

static void fillWithAdjacentVertsUsingVf(BfTrimesh const *trimesh, BfSize i,
                                         BfSize *nvvPtr, BfSize *capacityPtr,
                                         BfSize **vvPtr) {
  BEGIN_ERROR_HANDLING();

  BfSize nvv = *nvvPtr;
  BfSize capacity = *capacityPtr;
  BfSize *vv = *vvPtr;

  for (BfSize j = trimesh->vfOffset[i]; j < trimesh->vfOffset[i + 1]; ++j) {
    BfSize f = trimesh->vf[j];

    BfSize inds[2];
    bfTrimeshGetOpFaceVerts(trimesh, f, i, &inds[0], &inds[1]);

    for (BfSize k = 0; k < 2; ++k) {
      /* Scan through vv looking for inds[k] */
      BfSize l = 0;
      while (l < nvv && vv[l] < inds[k])
        ++l;

      /* If we found inds[k], continue---otherwise, need to insert
       * inds[k] into vv (in sorted order) */
      if (l < nvv && vv[l] == inds[k])
        continue;

      /* If we're at capacity, need to allocate */
      if (nvv == capacity) {
        BfSize capacityNew = 2*capacity;
        BfSize *vvNew = realloc(vv, capacityNew*sizeof(BfSize));
        if (vvNew == NULL)
          RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
        capacity = capacityNew;
        vv = vvNew;
      }

      /* Shift the latter part of vv back */
      if (l < nvv)
        memmove(&vv[l + 1], &vv[l], (nvv - l)*sizeof(BfSize));

      /* Insert inds[k] */
      vv[l] = inds[k];

      /* Grow the size of the adjacency list */
      ++nvv;
    }
  }

  *nvvPtr = nvv;
  *capacityPtr = capacity;
  *vvPtr = vv;

  END_ERROR_HANDLING() {}
}

static size_t countAdjacentVertsUsingVf(BfTrimesh const *trimesh, BfSize i) {
  BEGIN_ERROR_HANDLING();

  BfSize nvv = 0, capacity = 16;
  BfSize *vv = malloc(capacity*sizeof(BfSize));
  fillWithAdjacentVertsUsingVf(trimesh, i, &nvv, &capacity, &vv);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    nvv = BF_SIZE_BAD_VALUE;
  }

  free(vv);

  return nvv;
}

static void initVv(BfTrimesh *trimesh) {
  BEGIN_ERROR_HANDLING();

  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);

  /* Allocate space for vvOffset */
  trimesh->vvOffset = malloc((numVerts + 1)*sizeof(BfSize));
  if (trimesh->vvOffset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* First entry of vvOffset is just 0 */
  trimesh->vvOffset[0] = 0;

  /* Count the number of vertices neighboring each vertex using vf,
   * which should already be initialized */
  for (BfSize i = 0; i < numVerts; ++i)
    trimesh->vvOffset[i + 1] = countAdjacentVertsUsingVf(trimesh, i);

  /* Transform vvOffset using a cumulative sum over each entry */
  for (BfSize i = 0; i < numVerts; ++i)
    trimesh->vvOffset[i + 1] += trimesh->vvOffset[i];

  /* Allocate space for vv */
  trimesh->vv = malloc(trimesh->vvOffset[numVerts]*sizeof(BfSize));
  if (trimesh->vv == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* Fill vv */
  for (BfSize i = 0; i < numVerts; ++i) {
    BfSize nvv = 0;
    BfSize capacity = trimesh->vvOffset[i + 1] - trimesh->vvOffset[i];
    BfSize *vv = &trimesh->vv[trimesh->vvOffset[i]];
    fillWithAdjacentVertsUsingVf(trimesh, i, &nvv, &capacity, &vv);
  }

  END_ERROR_HANDLING() {}
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
  char *lineptr = NULL;
  ssize_t nread;
  size_t n = 1024;
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
    if (nread < 0)
      break;

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
    if (nread < 0)
      break;

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
  while (F[k] == i) ++k;
  *i0 = F[k];
  assert(*i0 != i);
  while (F[k] == *i0 || F[k] == i) ++k;
  *i1 = F[k];
  assert(*i1 != *i0 && *i1 != i);
}
