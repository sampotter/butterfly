#include <bf/trimesh.h>

#include <stdlib.h>
#include <string.h>

#include <bf/array.h>
#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/linalg.h>
#include <bf/mat_csr_real.h>
#include <bf/mem.h>
#include <bf/points.h>
#include <bf/real_array.h>
#include <bf/size_array.h>
#include <bf/util.h>
#include <bf/vec_real.h>
#include <bf/vectors.h>

#include "macros.h"

#ifdef BF_EMBREE
#  include <embree4/rtcore.h>
#endif

struct BfTrimesh {
  BfPoints3 *verts;

  BfSize numFaces;
  BfSize3 *faces;

  /*! Array containing `BfSize2`s representing the edge indices. This
   *  array works analogously to `faces`, with the same sorting
   *  scheme. */
  BfArray *edges;

  BfSize *vfOffset, *vf;
  BfSize *vvOffset, *vv;

  /*! An array containing `BfSize2`s representing a mapping from edges
   *  to incident faces. The two components of the `i`th entry are the
   *  indices of faces incident on the `i`th edge. We assume the mesh
   *  is manifold, so there can be at most two incident faces. If the
   *  edge is a boundary edge, then one of the components of the `i`th
   *  entry equals `BF_SIZE_BAD_VALUE`. */
  BfArray *ef;

  bool *isBoundaryEdge;
  bool *isBoundaryVert;

  BfSize2 *boundaryEdges;
  BfSize numBoundaryEdges;

#ifdef BF_EMBREE
  RTCDevice device;
  RTCGeometry geometry;
  RTCScene scene;
  float (*vertexBuffer)[3];
  uint32_t (*indexBuffer)[3];
#endif

  BfPoint3 *faceCentroids;
  BfVector3 *faceUnitNormals;
  BfReal *faceAreas;
};

int comparFace(BfSize const *face1, BfSize const *face2, void *arg) {
  (void)arg;
  if (face1[0] == face2[0] && face1[1] == face2[1])
    return face1[2] - face2[2];
  else if (face1[0] == face2[0])
    return face1[1] - face2[1];
  else
    return face1[0] - face2[0];
}

BfTrimesh *bfTrimeshCopy(BfTrimesh const *trimesh) {
  BF_ERROR_BEGIN();

  BfTrimesh *trimeshCopy = bfMemAlloc(1, sizeof(BfTrimesh));
  HANDLE_ERROR();

  trimeshCopy->verts = bfPoints3Copy(trimesh->verts);
  HANDLE_ERROR();

  trimeshCopy->faces = bfMemAllocCopy(trimesh->faces, trimesh->numFaces, sizeof(BfSize3));
  HANDLE_ERROR();

  trimeshCopy->numFaces = trimesh->numFaces;

  trimeshCopy->edges = bfArrayCopy(trimesh->edges);
  HANDLE_ERROR();

  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);

  trimeshCopy->vfOffset = bfMemAllocCopy(trimesh->vfOffset, numVerts + 1, sizeof(BfSize));
  HANDLE_ERROR();

  trimeshCopy->vf = bfMemAllocCopy(trimesh->vf, trimesh->vfOffset[numVerts], sizeof(BfSize));
  HANDLE_ERROR();

  trimeshCopy->vvOffset = bfMemAllocCopy(trimesh->vvOffset, numVerts + 1, sizeof(BfSize));
  HANDLE_ERROR();

  trimeshCopy->vv = bfMemAllocCopy(trimesh->vv, trimesh->vvOffset[numVerts], sizeof(BfSize));
  HANDLE_ERROR();

  trimeshCopy->ef = bfArrayCopy(trimesh->ef);
  HANDLE_ERROR();

  trimeshCopy->isBoundaryEdge = bfMemAllocCopy(trimesh->isBoundaryEdge, bfTrimeshGetNumEdges(trimesh), sizeof(bool));
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

int comparEdge(BfSize const *edge1, BfSize const *edge2, void *arg) {
  (void)arg;
  return edge1[0] == edge2[0] ? edge1[1] - edge2[1] : edge1[0] - edge2[0];
}

static void initEdges(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  trimesh->edges = bfArrayNewEmpty(sizeof(BfSize2));
  HANDLE_ERROR();

  /* Iterate over each pair of adjacent vertices: */
  for (BfSize i0 = 0; i0 < bfTrimeshGetNumVerts(trimesh); ++i0) {
    BfSizeArray *vertNbs = bfTrimeshGetVertNbs(trimesh, i0);
    HANDLE_ERROR();

    for (BfSize j = 0; j < bfSizeArrayGetSize(vertNbs); ++j) {
      BfSize i1 = bfSizeArrayGet(vertNbs, j);

      BF_DEFINE_EDGE(edge, i0, i1);

      BfSize pos = bfArrayFindSorted(trimesh->edges, edge, (BfCompar)comparEdge);
      if (!BF_SIZE_OK(pos))
        continue;

      /* If `pos` is a valid index, check if the edge at the `pos`th
       * position is the same as `edge`. If it is, we avoid adding it
       * again by continuing. */
      if (pos < bfArrayGetSize(trimesh->edges)) {
        BfSize2 otherEdge;
        bfArrayGet(trimesh->edges, pos, otherEdge);
        if (!comparEdge(edge, otherEdge, NULL))
          continue;
      }

      bfArrayInsert(trimesh->edges, pos, edge);
      HANDLE_ERROR();
    }

    bfSizeArrayDeinitAndDealloc(&vertNbs);
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void initEf(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  BF_DEFINE_EDGE(emptyEdge, BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE);

  trimesh->ef = bfArrayNewWithValue(sizeof(BfSize2), bfTrimeshGetNumEdges(trimesh), emptyEdge);
  HANDLE_ERROR();

  for (BfSize i = 0; i < bfTrimeshGetNumFaces(trimesh); ++i) {
    for (BfSize j = 0; j < 3; ++j) {
      BF_DEFINE_EDGE(edge, trimesh->faces[i][j], trimesh->faces[i][(j + 1) % 3]);

      BfSize k = bfArrayFindSorted(trimesh->edges, edge, (BfCompar)comparEdge);
      BF_ASSERT(k < bfTrimeshGetNumEdges(trimesh));

      BfSize *ef = bfArrayGetPtr(trimesh->ef, k);
      if (!BF_SIZE_OK(ef[0])) {
        ef[0] = i;
      } else {
        BF_ASSERT(!BF_SIZE_OK(ef[1]));
        ef[1] = i;
      }
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void initBoundaryInfo(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  trimesh->isBoundaryEdge = bfMemAllocAndZero(bfTrimeshGetNumEdges(trimesh), sizeof(bool));
  HANDLE_ERROR();

  for (BfSize i = 0; i < trimesh->numFaces; ++i) {
    for (BfSize j = 0; j < 3; ++j) {
      BF_DEFINE_EDGE(edge, trimesh->faces[i][j], trimesh->faces[i][(j + 1) % 3]);

      BfSize k = bfArrayFindSorted(trimesh->edges, edge, (BfCompar)comparEdge);
      BF_ASSERT(BF_SIZE_OK(k));

      trimesh->isBoundaryEdge[k] = !trimesh->isBoundaryEdge[k];
    }
  }

  trimesh->isBoundaryVert = bfMemAllocAndZero(bfTrimeshGetNumVerts(trimesh), sizeof(bool));
  HANDLE_ERROR();

  for (BfSize i = 0; i < bfTrimeshGetNumEdges(trimesh); ++i) {
    if (trimesh->isBoundaryEdge[i]) {
      BfSize2 edge;
      bfArrayGet(trimesh->edges, i, edge);

      for (BfSize j = 0; j < 2; ++j)
        trimesh->isBoundaryVert[edge[j]] = true;
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void initBoundaryEdges(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  trimesh->numBoundaryEdges = 0;
  for (BfSize i = 0; i < bfTrimeshGetNumEdges(trimesh); ++i)
    if (trimesh->isBoundaryEdge[i])
      ++trimesh->numBoundaryEdges;

  trimesh->boundaryEdges = bfMemAlloc(trimesh->numBoundaryEdges, sizeof(BfSize2));
  HANDLE_ERROR();

  BfSize j = 0;
  for (BfSize i = 0; i < bfTrimeshGetNumEdges(trimesh); ++i)
    if (trimesh->isBoundaryEdge[i])
      bfArrayGet(trimesh->edges, i, &trimesh->boundaryEdges[j++]);
  BF_ASSERT(j == trimesh->numBoundaryEdges);

  BF_ERROR_END() {
    BF_DIE();
  }
}

#ifdef BF_EMBREE
static void initEmbree(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  trimesh->device = rtcNewDevice(NULL);
  if (trimesh->device == NULL)
    RAISE_ERROR(BF_ERROR_EMBREE);

  trimesh->geometry = rtcNewGeometry(trimesh->device, RTC_GEOMETRY_TYPE_TRIANGLE);
  if (trimesh->geometry == NULL)
    RAISE_ERROR(BF_ERROR_EMBREE);

  trimesh->scene = rtcNewScene(trimesh->device);
  if (trimesh->scene == NULL)
    RAISE_ERROR(BF_ERROR_EMBREE);

  trimesh->vertexBuffer = rtcSetNewGeometryBuffer(
    /* geometry: */ trimesh->geometry,
    /* type: */ RTC_BUFFER_TYPE_VERTEX,
    /* slot: */ 0,
    /* format: */ RTC_FORMAT_FLOAT3,
    /* byteStride: */ 3*sizeof(float),
    /* numItems: */ bfPoints3GetSize(trimesh->verts));
  if (trimesh->vertexBuffer == NULL)
    RAISE_ERROR(BF_ERROR_EMBREE);

  for (BfSize i = 0; i < bfPoints3GetSize(trimesh->verts); ++i) {
    BfReal const *v = bfPoints3GetPtrConst(trimesh->verts, i);
    for (BfSize j = 0; j < 3; ++j)
      trimesh->vertexBuffer[i][j] = v[j];
  }

  trimesh->indexBuffer = rtcSetNewGeometryBuffer(
    /* geometry: */ trimesh->geometry,
    /* type: */ RTC_BUFFER_TYPE_INDEX,
    /* slot: */ 0,
    /* format: */ RTC_FORMAT_UINT3,
    /* byteStride: */ 3*sizeof(uint32_t),
    /* numItems: */ trimesh->numFaces);
  if (trimesh->indexBuffer == NULL)
    RAISE_ERROR(BF_ERROR_EMBREE);

  for (BfSize i = 0; i < trimesh->numFaces; ++i) {
    BfSize const *f = trimesh->faces[i];
    for (BfSize j = 0; j < 3; ++j)
      trimesh->indexBuffer[i][j] = f[j];
  }

  rtcCommitGeometry(trimesh->geometry);
  rtcAttachGeometry(trimesh->scene, trimesh->geometry);
  rtcReleaseGeometry(trimesh->geometry);
  rtcCommitScene(trimesh->scene);

  BF_ERROR_END() {
    BF_DIE();
  }
}
#endif

static void initFaceCentroids(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  trimesh->faceCentroids = bfMemAllocAndZero(trimesh->numFaces, sizeof(BfPoint3));
  HANDLE_ERROR();

  for (BfSize i = 0; i < trimesh->numFaces; ++i) {
    for (BfSize j = 0; j < 3; ++j) {
      BfReal const *x = bfPoints3GetPtrConst(trimesh->verts, trimesh->faces[i][j]);
      for (BfSize k = 0; k < 3; ++k)
        trimesh->faceCentroids[i][k] += x[k];
    }
    for (BfSize k = 0; k < 3; ++k)
      trimesh->faceCentroids[i][k] /= 3;
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void initFaceAreas(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  trimesh->faceAreas = bfMemAlloc(trimesh->numFaces, sizeof(BfReal));
  HANDLE_ERROR();

  for (BfSize i = 0; i < trimesh->numFaces; ++i) {
    BfReal const *x[3];
    for (BfSize j = 0; j < 3; ++j)
      x[j] = bfPoints3GetPtrConst(trimesh->verts, trimesh->faces[i][j]);
    BfVector3 dx1, dx2, n;
    bfPoint3Sub(x[1], x[0], dx1);
    bfPoint3Sub(x[2], x[0], dx2);
    bfVector3Cross(dx1, dx2, n);
    trimesh->faceAreas[i] = bfVector3Norm(n)/2;
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void initCommon(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  if (bfTrimeshHasDuplicateFaces(trimesh))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  initVf(trimesh);
  HANDLE_ERROR();

  initVv(trimesh);
  HANDLE_ERROR();

  initEdges(trimesh);
  HANDLE_ERROR();

  initEf(trimesh);
  HANDLE_ERROR();

  initBoundaryInfo(trimesh);
  HANDLE_ERROR();

  initBoundaryEdges(trimesh);
  HANDLE_ERROR();

#ifdef BF_EMBREE
  initEmbree(trimesh);
  HANDLE_ERROR();
#endif

  initFaceCentroids(trimesh);
  initFaceAreas(trimesh);

  BF_ERROR_END() {
    BF_DIE();
  }
}

#ifdef BF_EMBREE
static void deinitEmbree(BfTrimesh *trimesh) {
  rtcReleaseScene(trimesh->scene);
  rtcReleaseDevice(trimesh->device);
}
#endif

static void rebuildMesh(BfTrimesh *trimesh) {
  BF_ERROR_BEGIN();

  bfArrayDeinitAndDealloc(&trimesh->edges);
  bfMemFree(trimesh->vfOffset);
  bfMemFree(trimesh->vf);
  bfMemFree(trimesh->vvOffset);
  bfMemFree(trimesh->vv);
  bfArrayDeinitAndDealloc(&trimesh->ef);
  bfMemFree(trimesh->isBoundaryEdge);
  bfMemFree(trimesh->isBoundaryVert);
  bfMemFree(trimesh->boundaryEdges);

#ifdef BF_EMBREE
  deinitEmbree(trimesh);
#endif

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

  bfArrayDeinitAndDealloc(&trimesh->edges);

  bfMemFree(trimesh->vfOffset);
  trimesh->vfOffset = NULL;

  bfMemFree(trimesh->vf);
  trimesh->vf = NULL;

  bfMemFree(trimesh->vvOffset);
  trimesh->vvOffset = NULL;

  bfMemFree(trimesh->vv);
  trimesh->vv = NULL;

  bfArrayDeinitAndDealloc(&trimesh->ef);

  bfMemFree(trimesh->isBoundaryEdge);
  trimesh->isBoundaryEdge = NULL;

  bfMemFree(trimesh->isBoundaryVert);
  trimesh->isBoundaryVert = NULL;

  bfMemFree(trimesh->boundaryEdges);
  trimesh->boundaryEdges = NULL;

  trimesh->numBoundaryEdges = BF_SIZE_BAD_VALUE;

#ifdef BF_EMBREE
  deinitEmbree(trimesh);
#endif
}

void bfTrimeshDealloc(BfTrimesh **trimesh) {
  bfMemFree(*trimesh);
  *trimesh = NULL;
}

void bfTrimeshDeinitAndDealloc(BfTrimesh **trimesh) {
  bfTrimeshDeinit(*trimesh);
  bfTrimeshDealloc(trimesh);
}

BfSize bfTrimeshGetNumVerts(BfTrimesh const *trimesh) {
  return trimesh->verts->size;
}

BfSize bfTrimeshGetNumFaces(BfTrimesh const *trimesh) {
  return trimesh->numFaces;
}

BfPoints3 const *bfTrimeshGetVertsConst(BfTrimesh const *trimesh) {
  return trimesh->verts;
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

  BfSize const *ef = bfArrayGetPtr(trimesh->ef, edgeIndex);
  BF_ASSERT(BF_SIZE_OK(ef[0]) && !BF_SIZE_OK(ef[1]));

  faceIndex = ef[0];

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

  if (edgeInd >= bfTrimeshGetNumEdges(trimesh))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  if (lam <= 0 || 1 <= lam)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize2 edge;
  bfArrayGet(trimesh->edges, edgeInd, edge);

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

  edgeIndex = bfArrayFindSorted(trimesh->edges, edge, (BfCompar)comparEdge);

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
    BfSize2 edge;
    bfArrayGet(trimesh->edges, edgeInd, edge);

    if (faceContainsEdge(trimesh->faces[i], edge)) {
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

BfRealArray *bfTrimeshGetFiedler(BfTrimesh const *trimesh) {
  BF_ERROR_BEGIN();

  /* Get the mask that indicates which vertices of `trimesh` are
   * interior vertices. */

  bool *mask = bfMemAlloc(bfTrimeshGetNumVerts(trimesh), sizeof(bool));
  HANDLE_ERROR();

  for (BfSize i = 0; i < bfTrimeshGetNumVerts(trimesh); ++i)
    mask[i] = !trimesh->isBoundaryVert[i];

  /* Compute the stiffness and mass matrices for the piecewise linear
   * FEM approximation of the LBO and then extract the submatrices
   * corresponding to the interior vertices so we can solve the
   * Dirichlet eigenvalue problem. */

  BfMat *L = NULL, *M = NULL;
  bfTrimeshGetLboFemDiscretization(trimesh, &L, &M);
  HANDLE_ERROR();

  BfMat *LInt = bfMatGetSubmatByMask(L, mask, mask);
  HANDLE_ERROR();

  BfMat *MInt = bfMatGetSubmatByMask(M, mask, mask);
  HANDLE_ERROR();

  /* Solve the shifted eigenvalue problem and get the first nonzero
   * eigenfunction (the "Fiedler vector"). */

  BfMat *PhiTransposeInt = NULL;
  BfVecReal *Lam = NULL;
  bfGetShiftedEigs(LInt, MInt, -0.001, 2, &PhiTransposeInt, &Lam);
  HANDLE_ERROR();

  BfVecReal *phiFiedlerInt = bfVecToVecReal(bfMatGetRowView(PhiTransposeInt, 1));
  HANDLE_ERROR();

  BfVecReal *phiFiedler = bfVecRealNewWithValue(bfTrimeshGetNumVerts(trimesh), 0);
  HANDLE_ERROR();

  bfVecRealSetMask(phiFiedler, mask, bfVecRealToVec(phiFiedlerInt));

  BfRealArray *phiFiedlerArr = bfRealArrayNewFromVecReal(phiFiedler, BF_POLICY_COPY);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMemFree(mask);

  bfMatDelete(&L);
  bfMatDelete(&M);

  bfMatDelete(&LInt);
  bfMatDelete(&MInt);

  bfMatDelete(&PhiTransposeInt);
  bfVecRealDeinitAndDealloc(&Lam);

  bfVecRealDeinitAndDealloc(&phiFiedlerInt);
  bfVecRealDeinitAndDealloc(&phiFiedler);

  return phiFiedlerArr;
}

BfSize bfTrimeshGetNumEdges(BfTrimesh const *trimesh) {
  return bfArrayGetSize(trimesh->edges);
}

BfSizeArray *bfTrimeshGetVertNbs(BfTrimesh const *trimesh, BfSize i) {
  BF_ERROR_BEGIN();

  if (i >= bfTrimeshGetNumVerts(trimesh))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  BfSizeArray *vertNbs = bfSizeArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  for (BfSize j = trimesh->vvOffset[i]; j < trimesh->vvOffset[i + 1]; ++j)
    bfSizeArrayAppend(vertNbs, trimesh->vv[j]);

  BF_ERROR_END() {
    BF_DIE();
  }

  return vertNbs;
}

bool bfTrimeshHasDuplicateFaces(BfTrimesh const *trimesh) {
  BF_ERROR_BEGIN();

  bool hasDuplicateFaces = false;

  BfArray *uniqueFaces = bfArrayNewEmpty(sizeof(BfSize3));
  HANDLE_ERROR();

  for (BfSize i = 0; i < trimesh->numFaces; ++i) {
    BfSize const *j = trimesh->faces[i];

    BF_DEFINE_FACE(face, j[0], j[1], j[2]);

    BfSize pos = bfArrayFindSorted(uniqueFaces, face, (BfCompar)comparFace);
    if (pos < bfArrayGetSize(uniqueFaces)) {
      BfSize3 otherFace;
      bfArrayGet(uniqueFaces, pos, otherFace);
      if (comparFace(face, otherFace, NULL) == 0) {
        hasDuplicateFaces = true;
        break;
      }
    }

    bfArrayInsert(uniqueFaces, pos, face);
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  bfArrayDeinitAndDealloc(&uniqueFaces);

  return hasDuplicateFaces;
}

BfSize bfTrimeshGetNumVertexNeighbors(BfTrimesh const *trimesh, BfSize i) {
  if (i > bfTrimeshGetNumVerts(trimesh))
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
  return trimesh->vvOffset[i + 1] - trimesh->vvOffset[i];
}

BfSize bfTrimeshGetVertexNeighbor(BfTrimesh const *trimesh, BfSize i, BfSize j) {
  if (i > bfTrimeshGetNumVerts(trimesh))
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
  if (j > bfTrimeshGetNumVertexNeighbors(trimesh, i))
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
  return trimesh->vv[trimesh->vvOffset[i] + j];
}

bool bfTrimeshIsBoundaryVertex(BfTrimesh const *trimesh, BfSize i) {
  if (i > bfTrimeshGetNumVerts(trimesh))
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
  return trimesh->isBoundaryVert[i];
}

BfSize bfTrimeshGetNumBoundaryEdges(BfTrimesh const *trimesh) {
  return trimesh->numBoundaryEdges;
}

void bfTrimeshGetBoundaryEdge(BfTrimesh const *trimesh, BfSize i, BfSize2 boundaryEdge) {
  if (i > bfTrimeshGetNumBoundaryEdges(trimesh))
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
  bfMemCopy(trimesh->boundaryEdges[i], 1, sizeof(BfSize2), boundaryEdge);
}

BfSize const *bfTrimeshGetBoundaryEdgeConstPtr(BfTrimesh const *trimesh, BfSize i) {
  if (i > bfTrimeshGetNumBoundaryEdges(trimesh))
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
  return trimesh->boundaryEdges[i];
}

BfSize const *bfTrimeshGetFaceConstPtr(BfTrimesh const *trimesh, BfSize i) {
  if (i >= bfTrimeshGetNumFaces(trimesh))
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
  return trimesh->faces[i];
}

void bfTrimeshGetLboFemDiscretization(BfTrimesh const *trimesh, BfMat **L, BfMat **M) {
  BF_ERROR_BEGIN();

  BfSize *rowptr = NULL;
  BfSize *colind = NULL;
  BfReal *L_data = NULL;
  BfReal *M_data = NULL;

  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);

  BfSize nnz = numVerts + trimesh->vvOffset[numVerts];

  rowptr = bfMemAlloc(numVerts + 1, sizeof(BfSize));
  if (rowptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  colind = bfMemAlloc(nnz, sizeof(BfSize));
  if (colind == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  L_data = bfMemAllocAndZero(nnz, sizeof(BfReal));
  if (L_data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  M_data = bfMemAllocAndZero(nnz, sizeof(BfReal));
  if (M_data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* Fill rowptr: same as vvOffset except we need to add more space
   * for the vertices on the diagonal */
  for (BfSize i = 0; i <= numVerts; ++i)
    rowptr[i] = trimesh->vvOffset[i] + i;

  /* Fill colind */
  for (BfSize i = 0, j = 0; i < numVerts; ++i) {
    BF_ASSERT(j == rowptr[i]);

    /* Find the position of i in vv */
    BfSize kmid = trimesh->vvOffset[i];
    while (kmid < trimesh->vvOffset[i + 1] && trimesh->vv[kmid] < i)
      ++kmid;

    /* Copy over vv, inserting i into the correct sorted order */
    BfSize k = trimesh->vvOffset[i];
    for (; k < kmid; ++k)
      colind[j++] = trimesh->vv[k];
    colind[j++] = i;
    for (; k < trimesh->vvOffset[i + 1]; ++k)
      colind[j++] = trimesh->vv[k];
  }

  BfSize f, j, i, i0, i1, *jptr, k, k0, k1;
  BfPoint3 x, x0, x1, y, y0, y1;
  BfVector3 d, d0, d1, g, g0, g1, n;
  BfReal A;

  /* For each vertex: */
  for (i = 0; i < numVerts; ++i) {
    jptr = &colind[rowptr[i]];

    /* Get current vertex */
    bfTrimeshGetVertex(trimesh, i, x);

    /* Find position of `i` in current block of column indices */
    k = 0;
    while (jptr[k] != i) ++k;

    /* For each incident face: */
    for (j = trimesh->vfOffset[i]; j < trimesh->vfOffset[i + 1]; ++j) {
      f = trimesh->vf[j];

      /* Get the other two vertices in the current face */
      bfTrimeshGetOpFaceVerts(trimesh, f, i, &i0, &i1);
      bfTrimeshGetVertex(trimesh, i0, x0);
      bfTrimeshGetVertex(trimesh, i1, x1);

      /* Find their positions in the current block of column indices */
      k0 = k1 = 0;
      while (jptr[k0] != i0) ++k0;
      while (jptr[k1] != i1) ++k1;

      /* Opposite edge vectors */
      bfPoint3Sub(x1, x0, d);
      bfPoint3Sub(x, x1, d0);
      bfPoint3Sub(x0, x, d1);

      /* Compute orthogonal projection of each vertex onto the
       * opposite face edge */
      bfPoint3GetPointOnRay(x0, d, -bfVector3Dot(d, d1)/bfVector3Dot(d, d), y);
      bfPoint3GetPointOnRay(x1, d0, -bfVector3Dot(d0, d)/bfVector3Dot(d0, d0), y0);
      bfPoint3GetPointOnRay(x, d1, -bfVector3Dot(d1, d0)/bfVector3Dot(d1, d1), y1);

      /* Compute gradients for hat functions centered at each vertex */
      bfPoint3Sub(x, y, g);
      bfPoint3Sub(x0, y0, g0);
      bfPoint3Sub(x1, y1, g1);
      bfVector3Scale(g, 1/bfVector3Dot(g, g));
      bfVector3Scale(g0, 1/bfVector3Dot(g0, g0));
      bfVector3Scale(g1, 1/bfVector3Dot(g1, g1));

      /* Get triangle area */
      bfVector3Cross(d0, d1, n);
      A = bfVector3Norm(n)/2;

      /* Update L values */
      L_data[rowptr[i] + k] += A*bfVector3Dot(g, g);
      L_data[rowptr[i] + k0] += A*bfVector3Dot(g, g0);
      L_data[rowptr[i] + k1] += A*bfVector3Dot(g, g1);

      /* Update M values */
      M_data[rowptr[i] + k] += A/6;
      M_data[rowptr[i] + k0] += A/12;
      M_data[rowptr[i] + k1] += A/12;
    }
  }

  BfMatCsrReal *L_csr = bfMatCsrRealNewFromPtrs(numVerts, numVerts, rowptr, colind, L_data);
  HANDLE_ERROR();

  BfMatCsrReal *M_csr = bfMatCsrRealNewFromPtrs(numVerts, numVerts, rowptr, colind, M_data);
  HANDLE_ERROR();

  *L = bfMatCsrRealToMat(L_csr);
  *M = bfMatCsrRealToMat(M_csr);

  BF_ERROR_END() {
    bfMatCsrRealDeinitAndDealloc(&L_csr);
    bfMatCsrRealDeinitAndDealloc(&M_csr);
  }

  bfMemFree(rowptr);
  bfMemFree(colind);
  bfMemFree(L_data);
  bfMemFree(M_data);
}

#ifdef BF_EMBREE
BfSizeArray *bfTrimeshGetVisibility(BfTrimesh const *trimesh, BfSize srcInd, BfSizeArray const *tgtInds) {
  BF_ERROR_BEGIN();

  BfSizeArray *visTgtInds = bfSizeArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  struct RTCRayHit rayHit;

  struct RTCRay *ray = &rayHit.ray;
  struct RTCHit *hit = &rayHit.hit;

  BfReal const *srcVert = bfTrimeshGetVertPtrConst(trimesh, srcInd);
  ray->org_x = srcVert[0];
  ray->org_y = srcVert[1];
  ray->org_z = srcVert[2];

  for (BfSize i = 0; i < bfSizeArrayGetSize(tgtInds); ++i) {
    BfSize tgtInd = bfSizeArrayGet(tgtInds, i);
    BfReal const *tgtVert = bfTrimeshGetVertPtrConst(trimesh, tgtInd);

    ray->tnear = 0;
    ray->dir_x = tgtVert[0] - srcVert[0];
    ray->dir_y = tgtVert[1] - srcVert[1];
    ray->dir_z = tgtVert[2] - srcVert[2];
    ray->time = 0; /* <-- Unused */
    ray->tfar = BF_INFINITY;
    ray->flags = 0;

    hit->geomID = RTC_INVALID_GEOMETRY_ID;

    rtcIntersect1(trimesh->scene, &rayHit, NULL);

    // If there's a hit, `tfar` should be finite after tracing. Since
    // we're computing pairwise visibilities, there is *always* a hit.
    BF_ASSERT(isfinite(ray->tfar));

    BfSize visTgtInd = hit->primID;

    bfSizeArrayAppend(visTgtInds, visTgtInd);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return visTgtInds;
}
#endif

BfReal const *bfTrimeshGetFaceCentroidConstPtr(BfTrimesh const *trimesh, BfSize i) {
  if (i >= bfTrimeshGetNumFaces(trimesh))
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
  return trimesh->faceCentroids[i];
}

BfReal const *bfTrimeshGetFaceUnitNormalConstPtr(BfTrimesh const *trimesh, BfSize i) {
  if (i >= bfTrimeshGetNumFaces(trimesh))
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
  return trimesh->faceUnitNormals[i];
}

BfReal bfTrimeshGetFaceArea(BfTrimesh const *trimesh, BfSize i) {
  if (i >= bfTrimeshGetNumFaces(trimesh))
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
  return trimesh->faceAreas[i];
}
