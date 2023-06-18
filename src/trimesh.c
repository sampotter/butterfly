#include <bf/trimesh.h>

#include <stdlib.h>
#include <string.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/real_array.h>
#include <bf/util.h>
#include <bf/vectors.h>

#include "macros.h"

BfTrimesh *bfTrimeshCopy(BfTrimesh const *trimesh) {
  BF_ERROR_BEGIN();

  BfTrimesh *trimeshCopy = bfMemAlloc(1, sizeof(BfTrimesh));
  HANDLE_ERROR();

  trimeshCopy->verts = bfPoints3Copy(trimesh->verts);
  HANDLE_ERROR();

  trimeshCopy->faces = bfMemAllocCopy(trimesh->faces, trimesh->numFaces, sizeof(BfSize3));
  HANDLE_ERROR();

  trimeshCopy->numFaces = trimesh->numFaces;

  trimeshCopy->edges = bfMemAllocCopy(trimesh->edges, trimesh->numEdges, sizeof(BfSize2));
  HANDLE_ERROR();

  trimeshCopy->numEdges = trimesh->numEdges;

  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);

  trimeshCopy->vfOffset = bfMemAllocCopy(trimesh->vfOffset, numVerts + 1, sizeof(BfSize));
  HANDLE_ERROR();

  trimeshCopy->vf = bfMemAllocCopy(trimesh->vf, trimesh->vfOffset[numVerts], sizeof(BfSize));
  HANDLE_ERROR();

  trimeshCopy->vvOffset = bfMemAllocCopy(trimesh->vvOffset, numVerts + 1, sizeof(BfSize));
  HANDLE_ERROR();

  trimeshCopy->vv = bfMemAllocCopy(trimesh->vv, trimesh->vvOffset[numVerts], sizeof(BfSize));
  HANDLE_ERROR();

  trimeshCopy->isBoundaryEdge = bfMemAllocCopy(trimesh->isBoundaryEdge, trimesh->numEdges, sizeof(bool));
  HANDLE_ERROR();

  trimeshCopy->isBoundaryVert = bfMemAllocCopy(trimesh->isBoundaryVert, numVerts, sizeof(bool));
  HANDLE_ERROR();

  trimeshCopy->boundaryEdges = bfMemAllocCopy(trimesh->boundaryEdges, trimesh->numBoundaryEdges, sizeof(BfSize2));
  HANDLE_ERROR();

  trimeshCopy->numBoundaryEdges = trimesh->numBoundaryEdges;

  BF_ERROR_END() {
    BF_DIE();
  }

  return trimeshCopy;
}

BfTrimesh *bfTrimeshNewFromObjFile(char const *objPath) {
  BF_ERROR_BEGIN();

  BfTrimesh *trimesh = bfMemAlloc(1, sizeof(BfTrimesh));
  HANDLE_ERROR();

  bfTrimeshInitFromObjFile(trimesh, objPath);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return trimesh;
}

BfTrimesh *bfTrimeshNewFromVertsAndFaces(BfPoints3 const *verts, BfSize numFaces, BfSize3 const *faces) {
  BF_ERROR_BEGIN();

  BfTrimesh *trimesh = bfMemAlloc(1, sizeof(BfTrimesh));
  HANDLE_ERROR();

  bfTrimeshInitFromVertsAndFaces(trimesh, verts, numFaces, faces);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return trimesh;
}

static void initVf(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);
  BfSize numFaces = bfTrimeshGetNumFaces(trimesh);

  trimesh->vfOffset = bfMemAllocAndZero(numVerts + 1, sizeof(BfSize));
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

  trimesh->vf = bfMemAlloc(trimesh->vfOffset[numVerts], sizeof(BfSize));
  if (trimesh->vf == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfSize *offset = bfMemAllocAndZero(numVerts, sizeof(BfSize));
  if (offset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  for (BfSize i = 0; i < numFaces; ++i) {
    for (BfSize j = 0; j < 3; ++j) {
      BfSize k = trimesh->faces[i][j];
      trimesh->vf[trimesh->vfOffset[k] + offset[k]++] = i;
    }
  }

  BF_ERROR_END() {
    bfMemFree(trimesh->vfOffset);
    trimesh->vfOffset = NULL;

    bfMemFree(trimesh->vf);
    trimesh->vf = NULL;
  }

  bfMemFree(offset);
}

static void fillWithAdjacentVertsUsingVf(BfTrimesh const *trimesh, BfSize i,
                                         BfSize *nvvPtr, BfSize *capacityPtr,
                                         BfSize **vvPtr) {
  BF_ERROR_BEGIN();

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
        BfSize *vvNew = bfMemRealloc(vv, capacityNew, sizeof(BfSize));
        if (vvNew == NULL)
          RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
        capacity = capacityNew;
        vv = vvNew;
      }

      /* Shift the latter part of vv back */
      if (l < nvv)
        bfMemMove(&vv[l], nvv - l, sizeof(BfSize), &vv[l + 1]);

      /* Insert inds[k] */
      vv[l] = inds[k];

      /* Grow the size of the adjacency list */
      ++nvv;
    }
  }

  *nvvPtr = nvv;
  *capacityPtr = capacity;
  *vvPtr = vv;

  BF_ERROR_END() {}
}

static size_t countAdjacentVertsUsingVf(BfTrimesh const *trimesh, BfSize i) {
  BF_ERROR_BEGIN();

  BfSize nvv = 0, capacity = 16;
  BfSize *vv = bfMemAlloc(capacity, sizeof(BfSize));
  fillWithAdjacentVertsUsingVf(trimesh, i, &nvv, &capacity, &vv);
  HANDLE_ERROR();

  BF_ERROR_END() {
    nvv = BF_SIZE_BAD_VALUE;
  }

  bfMemFree(vv);

  return nvv;
}

static void initVv(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);

  /* Allocate space for vvOffset */
  trimesh->vvOffset = bfMemAlloc(numVerts + 1, sizeof(BfSize));
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
  trimesh->vv = bfMemAlloc(trimesh->vvOffset[numVerts], sizeof(BfSize));
  if (trimesh->vv == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* Fill vv */
  for (BfSize i = 0; i < numVerts; ++i) {
    BfSize nvv = 0;
    BfSize capacity = trimesh->vvOffset[i + 1] - trimesh->vvOffset[i];
    BfSize *vv = &trimesh->vv[trimesh->vvOffset[i]];
    fillWithAdjacentVertsUsingVf(trimesh, i, &nvv, &capacity, &vv);
  }

  BF_ERROR_END() {}
}

int compar(BfSize2 const elt1, BfSize2 const elt2) {
  return elt1[0] == elt2[0] ? elt1[1] - elt2[1] : elt1[0] - elt2[0];
}

static void initEdges(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  BfSize numEdges = 0;
  BfSize capacity = BF_ARRAY_DEFAULT_CAPACITY;
  BfSize2 *edges = bfMemAlloc(capacity, sizeof(BfSize2));
  HANDLE_ERROR();

  /* Iterate over each pair of adjacent vertices: */
  for (BfSize i0 = 0; i0 < bfTrimeshGetNumVerts(trimesh); ++i0) {
    for (BfSize j = trimesh->vvOffset[i0]; j < trimesh->vvOffset[i0 + 1]; ++j) {
      BfSize i1 = trimesh->vv[j];

      /* Set up the edge (represented as a sorted index pair for
       * uniqueness): */
      BfSize2 edge = {i0, i1};
      SORT2(edge[0], edge[1]);

      /* Check if we already added this edge: */
      bool found = false;
      for (BfSize k = 0; k < numEdges; ++k) {
        if (!compar(edge, edges[k])) {
          found = true;
          break;
        }
      }
      if (found) continue;

      /* Grow the edge array now if we need to: */
      if (numEdges == capacity) {
        capacity *= 2;
        BfSize2 *newEdges = bfMemRealloc(edges, capacity, sizeof(BfSize2));
        HANDLE_ERROR();

        edges = newEdges;
      }

      bfMemCopy(edge, 1, sizeof(BfSize2), edges[numEdges++]);
    }
  }

  trimesh->numEdges = numEdges;

  trimesh->edges = bfMemRealloc(edges, numEdges, sizeof(BfSize2));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void initBoundaryInfo(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  trimesh->isBoundaryEdge = bfMemAllocAndZero(trimesh->numEdges, sizeof(bool));
  HANDLE_ERROR();

  for (BfSize i = 0; i < trimesh->numFaces; ++i) {
    for (BfSize j = 0; j < 3; ++j) {
      BfSize2 edge = {trimesh->faces[i][j], trimesh->faces[i][(j + 1) % 3]};
      SORT2(edge[0], edge[1]);

      BfSize k = 0;
      for (; k < trimesh->numEdges; ++k)
        if (trimesh->edges[k][0] == edge[0] && trimesh->edges[k][1] == edge[1])
          break;

      trimesh->isBoundaryEdge[k] = !trimesh->isBoundaryEdge[k];
    }
  }

  trimesh->isBoundaryVert = bfMemAllocAndZero(bfTrimeshGetNumVerts(trimesh), sizeof(bool));
  HANDLE_ERROR();

  for (BfSize i = 0; i < trimesh->numEdges; ++i)
    if (trimesh->isBoundaryEdge[i])
      for (BfSize j = 0; j < 2; ++j)
        trimesh->isBoundaryVert[trimesh->edges[i][j]] = true;

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void initBoundaryEdges(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  trimesh->numBoundaryEdges = 0;
  for (BfSize i = 0; i < trimesh->numEdges; ++i)
    if (trimesh->isBoundaryEdge[i])
      ++trimesh->numBoundaryEdges;

  trimesh->boundaryEdges = bfMemAlloc(trimesh->numBoundaryEdges, sizeof(BfSize2));
  HANDLE_ERROR();

  BfSize j = 0;
  for (BfSize i = 0; i < trimesh->numEdges; ++i)
    if (trimesh->isBoundaryEdge[i])
      bfMemCopy(trimesh->edges[i], 1, sizeof(BfSize2), trimesh->boundaryEdges[j++]);
  BF_ASSERT(j == trimesh->numBoundaryEdges);

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void initCommon(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  initVf(trimesh);
  HANDLE_ERROR();

  initVv(trimesh);
  HANDLE_ERROR();

  initEdges(trimesh);
  HANDLE_ERROR();

  initBoundaryInfo(trimesh);
  HANDLE_ERROR();

  initBoundaryEdges(trimesh);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void rebuildMesh(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  bfMemFree(trimesh->edges);
  bfMemFree(trimesh->vfOffset);
  bfMemFree(trimesh->vf);
  bfMemFree(trimesh->vvOffset);
  bfMemFree(trimesh->vv);
  bfMemFree(trimesh->isBoundaryEdge);
  bfMemFree(trimesh->isBoundaryVert);
  bfMemFree(trimesh->boundaryEdges);

  initCommon(trimesh);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfTrimeshInitFromBinaryFiles(BfTrimesh *trimesh,
                                  char const *vertsPath,
                                  char const *facesPath) {
  BF_ERROR_BEGIN();

  trimesh->verts = bfPoints3NewFromBinaryFile(vertsPath);
  HANDLE_ERROR();

  /* Get the face count */

  BfSize numBytes = bfGetFileSizeInBytes(facesPath);
  HANDLE_ERROR();

  if (numBytes % sizeof(BfSize3) != 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
  trimesh->numFaces = numBytes/sizeof(BfSize3);

  /* Allocate space and load the faces */

  trimesh->faces = bfMemAlloc(trimesh->numFaces, sizeof(BfSize3));
  if (trimesh->faces == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  bfReadFileToMemory(facesPath, numBytes, (BfByte *)trimesh->faces);
  HANDLE_ERROR();

  initCommon(trimesh);
  HANDLE_ERROR();

  BF_ERROR_END() {}
}

void bfTrimeshInitFromObjFile(BfTrimesh *trimesh, char const *objPath) {
  BF_ERROR_BEGIN();

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

    bfMemFree(lineptr);
  } while (nread >= 0);

  /* Allocate space for the vertices and faces */

  trimesh->verts = bfPoints3NewWithDefaultCapacity();
  HANDLE_ERROR();

  trimesh->numFaces = num_face_vertex_indices;

  trimesh->faces = bfMemAlloc(trimesh->numFaces, sizeof(BfSize3));
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
      BfPoint3 vert;
      for (size_t i = 0; i < 3; ++i) {
        tok = strtok_r(NULL, " ", &saveptr);
        vert[i] = strtod(tok, NULL);
      }
      bfPoints3Append(trimesh->verts, vert);
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
        BF_DIE();
      }
    }

    bfMemFree(lineptr);
  } while (nread >= 0);

  BF_ASSERT(bfTrimeshGetNumVerts(trimesh) == num_verts);

  initCommon(trimesh);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  fclose(fp);
}

void bfTrimeshInitFromVertsAndFaces(BfTrimesh *trimesh, BfPoints3 const *verts, BfSize numFaces, BfSize3 const *faces) {
  BF_ERROR_BEGIN();

  /* Make sure none of the face indices are out of bounds: */
  bool anyFaceIndexOutOfBounds = false;
  for (BfSize i = 0; i < numFaces; ++i)
    for (BfSize j = 0; j < 3; ++j)
      if (faces[i][j] >= verts->size)
        anyFaceIndexOutOfBounds = true;
  if (anyFaceIndexOutOfBounds)
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  /* Make sure all vertices are used in the mesh: */
  bool *vertInMesh = bfMemAllocAndZero(verts->size, sizeof(bool));
  for (BfSize i = 0; i < numFaces; ++i)
    for (BfSize j = 0; j < 3; ++j)
      vertInMesh[faces[i][j]] = true;
  bool allVertsUsed = true;
  for (BfSize i = 0; i < verts->size; ++i)
    if (!vertInMesh[i])
      allVertsUsed = false;
  if (!allVertsUsed)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  trimesh->verts = bfPoints3Copy(verts);
  HANDLE_ERROR();

  trimesh->faces = bfMemAlloc(numFaces, sizeof(BfSize3));
  HANDLE_ERROR();

  bfMemCopy(faces, numFaces, sizeof(BfSize3), trimesh->faces);

  trimesh->numFaces = numFaces;

  initCommon(trimesh);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMemFree(vertInMesh);
}

void bfTrimeshDeinit(BfTrimesh *trimesh) {
  bfPoints3DeinitAndDealloc(&trimesh->verts);

  bfMemFree(trimesh->faces);
  trimesh->faces = NULL;

  trimesh->numFaces = BF_SIZE_BAD_VALUE;

  bfMemFree(trimesh->vfOffset);
  trimesh->vfOffset = NULL;

  bfMemFree(trimesh->vf);
  trimesh->vf = NULL;

  bfMemFree(trimesh->vvOffset);
  trimesh->vvOffset = NULL;

  bfMemFree(trimesh->vv);
  trimesh->vv = NULL;
}

void bfTrimeshDealloc(BfTrimesh **trimesh) {
  (void)trimesh;
  BF_DIE();
}

void bfTrimeshDeinitAndDealloc(BfTrimesh **trimesh) {
  (void)trimesh;
  BF_DIE();
}

BfSize bfTrimeshGetNumVerts(BfTrimesh const *trimesh) {
  return trimesh->verts->size;
}

BfSize bfTrimeshGetNumFaces(BfTrimesh const *trimesh) {
  return trimesh->numFaces;
}

void bfTrimeshGetVertex(BfTrimesh const *trimesh, BfSize i, BfPoint3 x) {
  bfMemCopy(trimesh->verts->data[i], 1, sizeof(BfPoint3), x);
}

BfReal const *bfTrimeshGetVertPtrConst(BfTrimesh const *trimesh, BfSize i) {
  BF_ERROR_BEGIN();

  BfReal const *vert = NULL;

  if (i >= bfTrimeshGetNumVerts(trimesh))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  vert = trimesh->verts->data[i];

  BF_ERROR_END() {
    BF_DIE();
  }

  return vert;
}

void bfTrimeshGetOpFaceVerts(BfTrimesh const *trimesh, BfSize faceIndex,
                             BfSize i, BfSize *i0, BfSize *i1) {
  BfSize *F = trimesh->faces[faceIndex];
  BF_ASSERT(F[0] == i || F[1] == i || F[2] == i);
  BfSize k = 0;
  while (F[k] == i) ++k;
  *i0 = F[k];
  BF_ASSERT(*i0 != i);
  while (F[k] == *i0 || F[k] == i) ++k;
  *i1 = F[k];
  BF_ASSERT(*i1 != *i0 && *i1 != i);
}

BfSize bfTrimeshGetBoundaryFaceIndexByBoundaryEdge(BfTrimesh const *trimesh, BfSize2 const boundaryEdge) {
  BF_ERROR_BEGIN();

  BfSize faceIndex = BF_SIZE_BAD_VALUE;

  if (boundaryEdge[0] >= boundaryEdge[1])
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (boundaryEdge[0] == BF_SIZE_BAD_VALUE || boundaryEdge[1] == BF_SIZE_BAD_VALUE)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize edgeIndex = bfTrimeshGetEdgeIndex(trimesh, boundaryEdge);
  if (edgeIndex == BF_SIZE_BAD_VALUE)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (!trimesh->isBoundaryEdge[edgeIndex])
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  for (BfSize i = 0; i < trimesh->numFaces; ++i) {
    for (BfSize j = 0; j < 3; ++j) {
      BfSize2 faceEdge = {
        trimesh->faces[i][j],
        trimesh->faces[i][(j + 1) % 3]
      };
      SORT2(faceEdge[0], faceEdge[1]);

      if (faceEdge[0] == boundaryEdge[0] && faceEdge[1] == boundaryEdge[1]) {
        faceIndex = i;
        break;
      }
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return faceIndex;
}

void bfTrimeshAddFace(BfTrimesh *trimesh, BfSize3 const face) {
  BF_ERROR_BEGIN();

  BfSize3 *newFaces = bfMemRealloc(trimesh->faces, trimesh->numFaces + 1, sizeof(BfSize3));
  HANDLE_ERROR();

  trimesh->faces = newFaces;

  bfMemCopy(face, 1, sizeof(BfSize3), trimesh->faces[trimesh->numFaces++]);

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfTrimeshDeleteFace(BfTrimesh *trimesh, BfSize faceInd) {
  BF_ERROR_BEGIN();

  if (faceInd >= trimesh->numFaces)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfMemMove(trimesh->faces + faceInd + 1,
            trimesh->numFaces - faceInd - 1,
            sizeof(BfSize3),
            trimesh->faces + faceInd);

  --trimesh->numFaces;

  BfSize3 *newFaces = bfMemRealloc(trimesh->faces, trimesh->numFaces, sizeof(BfSize3));
  HANDLE_ERROR();

  trimesh->faces = newFaces;

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfTrimeshSplitEdge(BfTrimesh *trimesh, BfSize edgeInd, BfReal lam) {
  BF_ERROR_BEGIN();

  if (edgeInd >= trimesh->numEdges)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (lam <= 0 || 1 <= lam)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize const *edge = trimesh->edges[edgeInd];

  BfReal const *x0 = bfTrimeshGetVertPtrConst(trimesh, edge[0]);
  BfReal const *x1 = bfTrimeshGetVertPtrConst(trimesh, edge[1]);

  BfPoint3 dx; bfPoint3Sub(x1, x0, dx);

  BfPoint3 xlam; bfPoint3GetPointOnRay(x0, dx, lam, xlam);

  BfSize newVertInd = bfTrimeshGetNumVerts(trimesh);

  bfPoints3Append(trimesh->verts, xlam);

  BfSize incFaceInds[2] = {BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE};
  BfSize numIncFaces = bfTrimeshGetFacesIncOnEdge(trimesh, edgeInd, incFaceInds);

  for (BfSize i = 0; i < numIncFaces; ++i) {
    BfSize3 face;
    for (BfSize j = 0; j < 3; ++j)
      face[j] = trimesh->faces[incFaceInds[i]][j];

    bfTrimeshDeleteFace(trimesh, incFaceInds[i]);

    for (BfSize j = 0; j < 2; ++j) {
      for (BfSize k = 0; k < 3; ++k) {
        if (face[k] == edge[j]) {
          BfSize old = face[k];
          face[k] = newVertInd;
          bfTrimeshAddFace(trimesh, face);
          face[k] = old;
        }
      }
    }
  }

  rebuildMesh(trimesh);

  BF_ERROR_END() {
    BF_DIE();
  }
}

BfSize bfTrimeshGetEdgeIndex(BfTrimesh const *trimesh, BfSize2 const edge) {
  BF_ERROR_BEGIN();

  BfSize edgeIndex = BF_SIZE_BAD_VALUE;

  if (edge[0] >= edge[1])
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (edge[0] == BF_SIZE_BAD_VALUE || edge[1] == BF_SIZE_BAD_VALUE)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  for (BfSize i = 0; i < trimesh->numEdges; ++i) {
    if (trimesh->edges[i][0] == edge[0] && trimesh->edges[i][1] == edge[1]) {
      edgeIndex = i;
      break;
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return edgeIndex;
}

static bool faceContainsEdge(BfSize const *face, BfSize const *edge) {
  BF_ERROR_BEGIN();

  bool containsEdge = false;

  if (edge[0] >= edge[1])
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (edge[0] == BF_SIZE_BAD_VALUE || edge[1] == BF_SIZE_BAD_VALUE)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  for (BfSize i = 0; i < 3; ++i) {
    BfSize2 faceEdge = {face[i], face[(i + 1) % 3]};
    SORT2(faceEdge[0], faceEdge[1]);

    if ((containsEdge = faceEdge[0] == edge[0] && faceEdge[1] == edge[1]))
      break;
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return containsEdge;
}

BfSize bfTrimeshGetFacesIncOnEdge(BfTrimesh const *trimesh, BfSize edgeInd, BfSize2 incFaceInds) {
  BF_ERROR_BEGIN();

  BfSize numIncFaces = 0;
  for (BfSize i = 0; i < trimesh->numFaces; ++i) {
    if (faceContainsEdge(trimesh->faces[i], trimesh->edges[edgeInd])) {
      /* TODO: assuming our mesh is manifold... */
      if (numIncFaces == 2) RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);
      incFaceInds[numIncFaces++] = i;
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return numIncFaces;
}

void bfTrimeshDumpVerts(BfTrimesh const *trimesh, char const *path) {
  BF_ERROR_BEGIN();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fwrite(trimesh->verts->data, sizeof(BfPoint3), trimesh->verts->size, fp);

  BF_ERROR_END() {
    BF_DIE();
  }

  fclose(fp);
}

void bfTrimeshDumpFaces(BfTrimesh const *trimesh, char const *path) {
  BF_ERROR_BEGIN();

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  fwrite(trimesh->faces, sizeof(BfSize3), trimesh->numFaces, fp);

  BF_ERROR_END() {
    BF_DIE();
  }

  fclose(fp);
}

static BfReal getNormalDerivEdge(BfTrimesh const *trimesh, BfRealArray const *values, BfSize const *boundaryEdge) {
  BF_ERROR_BEGIN();

  BfReal dudn = BF_NAN;

  BfSize faceInd = bfTrimeshGetBoundaryFaceIndexByBoundaryEdge(trimesh, boundaryEdge);
  HANDLE_ERROR();

  BF_ASSERT(faceInd != BF_SIZE_BAD_VALUE);

  BfSize const *face = trimesh->faces[faceInd];

  /* Fetch the index of the interior vertex from `face`. */

  BfSize iInt = BF_SIZE_BAD_VALUE;
  for (BfSize j = 0; j < 3; ++j) {
    bool found = false;
    for (BfSize k = 0; k < 2; ++k)
      if (face[j] == boundaryEdge[k])
        found = true;
    if (!found) {
      iInt = face[j];
      break;
    }
  }
  BF_ASSERT(iInt != BF_SIZE_BAD_VALUE);

  /* Next, we'll orthogonally project the interior vertex onto the
   * boundary edge. The projected point is `xProj` and the convex
   * coefficient giving it as a combination of the two boundary
   * vertices is `lamProj`. */

  BfReal const *xInt = bfTrimeshGetVertPtrConst(trimesh, iInt);
  BfReal const *x0 = bfTrimeshGetVertPtrConst(trimesh, boundaryEdge[0]);
  BfReal const *x1 = bfTrimeshGetVertPtrConst(trimesh, boundaryEdge[1]);

  BfVector3 dx; bfPoint3Sub(xInt, x0, dx);
  BfVector3 dx1; bfPoint3Sub(x1, x0, dx1);

  BfReal dx1Norm = bfVector3Norm(dx1);

  BfVector3 t;
  for (BfSize i = 0; i < 3; ++i) t[i] = dx1[i]/dx1Norm;

  BfReal t_dot_dx = bfVector3Dot(t, dx);

  BfVector3 xProj; bfPoint3GetPointOnRay(x0, t, t_dot_dx, xProj);

  BfReal lamProj = (xProj[0] - x0[0])/dx1[0];

  /* Finally, compute the normal derivative and return it. */

  BfReal uInt = bfRealArrayGetValue(values, iInt);
  BfReal u0 = bfRealArrayGetValue(values, boundaryEdge[0]);
  BfReal u1 = bfRealArrayGetValue(values, boundaryEdge[1]);

  BfReal du = uInt - (1 - lamProj)*u0 - lamProj*u1;
  BfReal dn = bfPoint3Dist(xInt, xProj);

  BF_ASSERT(dn > 0);
  dudn = du/dn;

  BF_ERROR_END() {
    BF_DIE();
  }

  return dudn;
}

/* This function computes the normal derivative of the piecewise
 * linear function defined by the nodal `values` at each of the
 * boundary vertices. This is done by first computing the normal
 * derivative on each edge and then using weighted interpolation to
 * approximate the normal derivative at each boundary vertex, which
 * should approximate the L2 projection.
 *
 * If `computeAtAllVerts` is true, then this will return an array of
 * the same size as `values` which is equal to `NAN` at each interior
 * vertex. */
BfRealArray *bfTrimeshGetNormalDeriv(BfTrimesh const *trimesh, BfRealArray const *values,
                                     bool computeAtAllVerts) {
  BF_ERROR_BEGIN();

  if (!computeAtAllVerts)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  /* Start by computing the boundary edge normal derivatives: */

  BfReal *normalDerivEdge = bfMemAlloc(trimesh->numBoundaryEdges, sizeof(BfReal));
  HANDLE_ERROR();

  for (BfSize i = 0; i < trimesh->numBoundaryEdges; ++i)
    normalDerivEdge[i] = getNormalDerivEdge(trimesh, values, trimesh->boundaryEdges[i]);

  /* Do a weighted interpolation to approximate the normal derivative
   * at each vertex.
   *
   * The idea here is that we want to do something simple which will
   * approximate the L2 projection of the piecewise constant normal
   * derivative defined on the boundary onto the space of piecewise
   * linear functions defined on the boundary. */

  BfReal *h = bfMemAlloc(trimesh->numBoundaryEdges, sizeof(BfReal));
  HANDLE_ERROR();

  for (BfSize i = 0; i < trimesh->numBoundaryEdges; ++i) {
    BfReal const *x0 = trimesh->verts->data[trimesh->boundaryEdges[i][0]];
    BfReal const *x1 = trimesh->verts->data[trimesh->boundaryEdges[i][1]];
    h[i] = bfPoint3Dist(x0, x1);
  }

  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);

  BfReal *denom = bfMemAllocAndZero(numVerts, sizeof(BfReal));
  HANDLE_ERROR();

  for (BfSize i = 0; i < trimesh->numBoundaryEdges; ++i)
    for (BfSize j = 0; j < 2; ++j)
      denom[trimesh->boundaryEdges[i][j]] += h[i];

  BfRealArray *normalDeriv = bfRealArrayNewWithValue(numVerts, BF_NAN);
  HANDLE_ERROR();

  for (BfSize i = 0; i < trimesh->numBoundaryEdges; ++i) {
    for (BfSize j = 0; j < 2; ++j) {
      BfSize k = trimesh->boundaryEdges[i][j];
      if (isnan(normalDeriv->data[k]))
        normalDeriv->data[k] = 0;
      normalDeriv->data[k] += normalDerivEdge[i]*h[i]/denom[k];
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMemFree(normalDerivEdge);
  bfMemFree(h);
  bfMemFree(denom);

  return normalDeriv;
}
