#include <bf/trimesh.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/size_array.h>
#include <bf/util.h>
#include <bf/vec_real.h>
#include <bf/vectors.h>

#include "macros.h"

static BfReal getTriArea(BfPoint3 const v0, BfPoint3 const v1, BfPoint3 const v2) {
  BfVector3 dv1, dv2;
  bfPoint3Sub(v1, v0, dv1);
  bfPoint3Sub(v2, v0, dv2);

  BfVector3 n;
  bfVector3Cross(dv1, dv2, n);

  BfReal triArea = bfVector3Norm(n)/2;

  return triArea;
}

static void appendFace(BfPoints3 const *verts, BfPoints3 const *cutVerts,
                       BfSize *numFaces, BfSize *capacity, BfSize3 **faces, BfSize3 const face) {
  BF_ERROR_BEGIN();

  /* Make sure the face isn't degenerate: */

  for (BfSize i = 0; i < 3; ++i) BF_ASSERT(BF_SIZE_OK(face[i]));

  BF_ASSERT(face[0] != face[1] && face[1] != face[2] && face[2] != face[0]);

  BfReal const *v[3] = {NULL, NULL, NULL};
  for (BfSize i = 0; i < 3; ++i)
    v[i] = face[i] < verts->size ?
      verts->data[face[i]] :
      cutVerts->data[face[i] - verts->size];

  BfReal triArea = getTriArea(v[0], v[1], v[2]);
  BF_ASSERT(triArea > 0);

  /* Append the face: */

  if (*numFaces == *capacity) {
    BfSize newCapacity = *capacity;
    newCapacity *= 2;

    BfSize3 *newFaces = bfMemRealloc(*faces, newCapacity, sizeof(BfSize3));
    HANDLE_ERROR();

    *faces = newFaces;
    *capacity = newCapacity;
  }

  bfMemCopy(face, 1, sizeof(BfSize3), *faces + *numFaces);

  ++*numFaces;

  BF_ERROR_END() {
    BF_DIE();
  }
}

static bool appendCutEdge(BfSize *numCutEdges, BfSize *cutEdgesCapacity, BfSize2 **cutEdges, BfSize2 cutEdge) {
  BF_ERROR_BEGIN();

  for (BfSize i = 0; i < *numCutEdges; ++i) {
    BfSize const *otherCutEdge = (*cutEdges)[i];
    if (cutEdge[0] == otherCutEdge[0] && cutEdge[1] == otherCutEdge[1])
      return false;
  }

  if (*numCutEdges == *cutEdgesCapacity) {
    BfSize newCapacity = *cutEdgesCapacity;
    newCapacity *= 2;

    BfSize2 *newCutEdges = bfMemRealloc(*cutEdges, newCapacity, sizeof(BfSize2));
    HANDLE_ERROR();

    *cutEdges = newCutEdges;
    *cutEdgesCapacity = newCapacity;
  }

  bfMemCopy(cutEdge, 1, sizeof(BfSize2), *cutEdges + *numCutEdges);

  ++*numCutEdges;

  BF_ERROR_END() {
    BF_DIE();
  }

  return true;
}

typedef struct {
  BfVector3 dv;
  BfReal const *vInt;
  BfReal const *vZero;
  BfReal const *vt;
} dfdsContext;

static BfReal dfds(BfReal s, dfdsContext const *c) {
  BfPoint3 vs;
  bfPoint3GetPointOnRay(c->vInt, c->dv, s, vs);

  BfVector3 v0_vs;
  bfPoint3Sub(vs, c->vZero, v0_vs);

  BfVector3 vt_vs;
  bfPoint3Sub(vs, c->vt, vt_vs);

  BfReal cos1 = bfVector3Dot(c->dv, v0_vs)/bfVector3Norm(v0_vs);
  BfReal cos2 = bfVector3Dot(c->dv, vt_vs)/bfVector3Norm(vt_vs);

  return cos1 + cos2;
}

/*!
 * Extract a submesh from trimesh using the level set function phi. In
 * this case, phi gives the nodal values of a piecewise linear level
 * set function defined on the triangle mesh. The extracted submesh
 * consists of the subset of trimesh where this level set function is
 * less than or equal to zero.
 *
 * NOTE: if permPtr is used, then afterwards it will contain a pointer
 * to a `BfSizeArray` containing the indices of the vertices in
 * trimesh which have been included in the level set submesh. If `n =
 * bfPtrArrayGetSize(*permPtr)`, then the first `n` vertices of the
 * submesh are vertices from `trimesh`, and `*permPtr` gives the
 * original positions of these vartices in `trimesh`.
 *
 * TODO: if this is slow, it's because it's O(n^2)! easy to fix...
 */
BfTrimesh *bfTrimeshGetLevelSetSubmesh(BfTrimesh const *trimesh, BfVecReal const *phi, bool const *permMask, BfSizeArray **permPtr) {
  BF_ERROR_BEGIN();

  BfTrimesh *submesh = NULL;

  BfSize numFaces = 0;
  BfSize facesCapacity = BF_ARRAY_DEFAULT_CAPACITY;
  BfSize3 *faces = bfMemAlloc(facesCapacity, sizeof(BfSize3));

  BfPoints3 *verts = bfPoints3NewWithDefaultCapacity();
  HANDLE_ERROR();

  BfSize numCutEdges = 0;
  BfSize cutEdgesCapacity = BF_ARRAY_DEFAULT_CAPACITY;
  BfSize2 *cutEdges = bfMemAlloc(cutEdgesCapacity, sizeof(BfSize2));

  BfPoints3 *cutVerts = bfPoints3NewWithDefaultCapacity();
  HANDLE_ERROR();

  /* Check whether caller wants us to fill the permutation array
   * (permPtr != NULL). If it does, then check whether we need to
   * allocate our own array or use one provided by caller. */
  BfSizeArray *perm = NULL;
  if (permPtr != NULL) {
    if (*permPtr == NULL) {
      perm = bfSizeArrayNewWithDefaultCapacity();
      HANDLE_ERROR();
    } else {
      perm = *permPtr;
    }
  }

#if BF_DEBUG
  /* Make sure there are no faces where all three level values are
   * zero. This may not be a problem, but we want to trip this the
   * first time this happens to check for potential problems. */
  for (BfSize i = 0; i < trimesh->numFaces; ++i) {
    BfReal phiFace[3];
    bfVecRealGetValues(phi, 3, trimesh->faces[i], phiFace);
    BF_ASSERT(phiFace[0] != 0 || phiFace[1] != 0 || phiFace[2] != 0);
  }
#endif

  for (BfSize i = 0; i < bfTrimeshGetNumVerts(trimesh); ++i) {
    if (bfVecRealGetElt(phi, i) > 0) continue;

    /* Check if this vertex is isolated: */
    bool isIsolated = true;
    for (BfSize j = trimesh->vvOffset[i]; j < trimesh->vvOffset[i + 1]; ++j) {
      if (bfVecRealGetElt(phi, trimesh->vv[j]) <= 0) {
        isIsolated = false;
        break;
      }
    }
    if (isIsolated) continue; /* (don't add it if it is) */

    BfReal const *v = bfTrimeshGetVertPtrConst(trimesh, i);
    bfPoints3Append(verts, v);

    if (perm != NULL && permMask[i])
      bfSizeArrayAppend(perm, i);
  }

  /* Accumulate the interior faces themselves. We do this now that
   * we've found all of the vertices from trimesh which will be
   * included in submesh, so we can reindex the faces accordingly. */
  for (BfSize i = 0; i < trimesh->numFaces; ++i) {
    BfReal phiFace[3];
    bfVecRealGetValues(phi, 3, trimesh->faces[i], phiFace);

    if (!(phiFace[0] <= 0 && phiFace[1] <= 0 && phiFace[2] <= 0))
      continue;

    BfSize3 newFace;
    for (BfSize j = 0; j < 3; ++j) {
      BfReal const *point = bfPoints3GetPtrConst(trimesh->verts, trimesh->faces[i][j]);
      newFace[j] = bfPoints3Find(verts, point);
      BF_ASSERT(newFace[j] != BF_SIZE_BAD_VALUE);
    }

    appendFace(verts, cutVerts, &numFaces, &facesCapacity, &faces, newFace);
  }

  /* Find the zero level set and accumulate new cut faces: */
  for (BfSize i = 0; i < trimesh->numFaces; ++i) {
    BfReal phiFace[3];
    bfVecRealGetValues(phi, 3, trimesh->faces[i], phiFace);

    /* TODO: weird case which we hope to avoid having to handle: */
    if (phiFace[0] == 0 && phiFace[1] == 0 && phiFace[2] == 0) BF_DIE();

    /* Skip faces that aren't cut by the zero level: */
    if (phiFace[0] <= 0 && phiFace[1] <= 0 && phiFace[2] <= 0) continue;
    if (phiFace[0] >= 0 && phiFace[1] >= 0 && phiFace[2] >= 0) continue;

    /* Figure out which verts are positive, negative, and zero: */
    BfSize numPos = 0, numNeg = 0, numZero = 0;
    BfSize iPos[2] = {BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE};
    BfSize iNeg[2] = {BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE};
    BfSize iZero[2] = {BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE};
    BfReal phiPos[2] = {BF_NAN, BF_NAN};
    BfReal phiNeg[2] = {BF_NAN, BF_NAN};
    for (BfSize j = 0; j < 3; ++j) {
      BF_ASSERT(numPos <= 2 && numNeg <= 2 && numZero <= 1);
      if (phiFace[j] > 0) {
        phiPos[numPos] = phiFace[j];
        iPos[numPos++] = trimesh->faces[i][j];
      } else if (phiFace[j] < 0) {
        phiNeg[numNeg] = phiFace[j];
        iNeg[numNeg++] = trimesh->faces[i][j];
      } else if (phiFace[j] == 0) {
        iZero[numZero++] = trimesh->faces[i][j];
      } else {
        BF_DIE();
      }
    }

    /** Now we'll handle the different cases: */

    if (numPos == 2 && numNeg == 1) {
      BfPoint3 v0, v[2];
      bfTrimeshGetVertex(trimesh, iNeg[0], v0);
      bfTrimeshGetVertex(trimesh, iPos[0], v[0]);
      bfTrimeshGetVertex(trimesh, iPos[1], v[1]);

      /* Compute `t` and `vt`. Afterwards, we'll have determined the
       * coordinates of three face vertices. */
      BfReal t[2];
      BfPoint3 vt[2];
      for (BfSize j = 0; j < 2; ++j) {
        t[j] = -phiNeg[0]/(phiPos[j] - phiNeg[0]);
        BF_ASSERT(0 <= t[j] && t[j] <= 1);

        BfPoint3 dv;
        bfPoint3Sub(v[j], v0, dv);
        bfPoint3GetPointOnRay(v0, dv, t[j], vt[j]);

        /* Due to numerical roundoff, `vt[j]` may exactly equal an existing
         * mesh vertex. We want to avoid duplicate mesh vertices, so
         * we set `t[j]` to `NAN` to signal this and continue. */
        if (bfPoints3Contains(verts, vt[j]))
          t[j] = BF_NAN;
      }

      /* Check if this is a degenerate face: if it is, we need to
       * continue now before adding anything to any of the arrays
       * we're building up. */
      bool coalesced[3] = {
        bfPoint3Dist(v0, vt[0]) <= BF_EPS,
        bfPoint3Dist(v0, vt[1]) <= BF_EPS,
        bfPoint3Dist(vt[0], vt[1]) <= BF_EPS
      };
      bool degenerateFace = coalesced[0] || coalesced[1] || coalesced[2];
      if (degenerateFace)
        continue;

      /* Compute new cut vertices and append them to `cutVerts`: */
      BfSize iCut[2] = {BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE};
      for (BfSize j = 0; j < 2; ++j) {
        if (isnan(t[j])) continue;
        BfSize2 cutEdge = {iNeg[0], iPos[j]};
        SORT2(cutEdge[0], cutEdge[1]);
        if (appendCutEdge(&numCutEdges, &cutEdgesCapacity, &cutEdges, cutEdge)) {
          BF_ASSERT(!bfPoints3Contains(cutVerts, vt[j]));
          iCut[j] = bfPoints3GetSize(cutVerts);
          bfPoints3Append(cutVerts, vt[j]);
        } else {
          BfSize vtIndex = bfPoints3Find(cutVerts, vt[j]);
          BF_ASSERT(vtIndex != BF_SIZE_BAD_VALUE);
          iCut[j] = vtIndex;
        }
      }

      /* Find the "new" index of `v0` (its position in `verts`): */
      BfSize i0New = bfPoints3Find(verts, v0);
      if (i0New == BF_SIZE_BAD_VALUE) {
        i0New = bfPoints3GetSize(verts);
        bfPoints3Append(verts, v0);
      }

      /* Find the "new" indices of the other two verts: */
      BfSize iCutNew[2];
      for (BfSize j = 0; j < 2; ++j) {
        if (isnan(t[j])) {
          iCutNew[j] = bfPoints3Find(verts, vt[j]);
        } else {
          BF_ASSERT(iCut[j] != BF_SIZE_BAD_VALUE);
          iCutNew[j] = verts->size + iCut[j];
        }
      }
      BF_ASSERT(iCutNew[0] != BF_SIZE_BAD_VALUE && iCutNew[1] != BF_SIZE_BAD_VALUE);

      /* Add a new face: */
      BfSize3 newFace = {i0New, iCutNew[0], iCutNew[1]};
      appendFace(verts, cutVerts, &numFaces, &facesCapacity, &faces, newFace);
    }

    else if (numPos == 1 && numNeg == 2) {
      BfPoint3 v0[2], v;
      bfTrimeshGetVertex(trimesh, iNeg[0], v0[0]);
      bfTrimeshGetVertex(trimesh, iNeg[1], v0[1]);
      bfTrimeshGetVertex(trimesh, iPos[0], v);

      /* Compute `t` and `vt`. Afterwards, we'll have determined the
       * coordinates of three face vertices. */
      BfReal t[2];
      BfPoint3 vt[2];
      for (BfSize j = 0; j < 2; ++j) {
        t[j] = -phiNeg[j]/(phiPos[0] - phiNeg[j]);
        BF_ASSERT(0 <= t[j] && t[j] <= 1);

        BfPoint3 dv;
        bfPoint3Sub(v, v0[j], dv);
        bfPoint3GetPointOnRay(v0[j], dv, t[j], vt[j]);

        /* Due to numerical roundoff, `vt[j]` may exactly equal an existing
         * mesh vertex. We want to avoid duplicate mesh vertices, so
         * we set `t[j]` to `NAN` to signal this and continue. */
        if (bfPoints3Contains(verts, vt[j]))
          t[j] = BF_NAN;
      }

      bool coalesced[3] = {
        bfPoint3Dist(v0[0], vt[0]) <= BF_EPS,
        bfPoint3Dist(v0[1], vt[1]) <= BF_EPS,
        bfPoint3Dist(vt[0], vt[1]) <= BF_EPS
      };

      // TODO: quad collapsed to an edge:
      BF_ASSERT(!(coalesced[0] && coalesced[1]));

      /* Compute new cut vertices and append them to `cutVerts`: */
      BfSize iCut[2] = {BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE};
      for (BfSize j = 0; j < 2; ++j) {
        if (coalesced[2] || coalesced[j]) continue;
        if (isnan(t[j])) continue;
        BfSize2 cutEdge = {iNeg[j], iPos[0]};
        SORT2(cutEdge[0], cutEdge[1]);
        if (appendCutEdge(&numCutEdges, &cutEdgesCapacity, &cutEdges, cutEdge)) {
          if (bfPoints3Contains(cutVerts, vt[j])) {
            /* We might have already added this cut vertex---in this
             * case, we should make sure that there's already a cut edge
             * containing the index of the duplicated vertex. */
            BF_ASSERT(t[j] == 0 || t[j] == 1);
            BfSize iAdded = t[j] == 0 ? cutEdge[0] : cutEdge[1];
            bool foundExistingCutVert = false;
            for (BfSize k = 0; k < numCutEdges; ++k)
              if (cutEdges[k][0] == iAdded || cutEdges[k][1] == iAdded)
                foundExistingCutVert = true;
            BF_ASSERT(foundExistingCutVert);

            /* Get the index of `vt[j]` in `cutVerts` since we can
             * conclude now that it's already been added: */
            iCut[j] = bfPoints3Find(cutVerts, vt[j]);
          } else {
            /* If we haven't added `vt[j]` already, then add it now. */
            iCut[j] = bfPoints3GetSize(cutVerts);
            bfPoints3Append(cutVerts, vt[j]);
          }
        } else {
          BfSize vtIndex = bfPoints3Find(cutVerts, vt[j]);
          BF_ASSERT(vtIndex != BF_SIZE_BAD_VALUE);
          iCut[j] = vtIndex;
        }
      }

      /* Find the "new" indices of these two mesh verts: */
      BfSize i0New[2] = {
        bfPoints3Find(verts, v0[0]),
        bfPoints3Find(verts, v0[1]),
      };
      BF_ASSERT(i0New[0] != BF_SIZE_BAD_VALUE && i0New[1] != BF_SIZE_BAD_VALUE);

      /* Find the "new" indices of the cut verts: */
      BfSize iCutNew[2];
      if (coalesced[2]) {
        BF_ASSERT(!coalesced[0] && !coalesced[1]);
        BF_ASSERT(bfPoint3Dist(vt[0], v) <= BF_EPS);
        BF_ASSERT(bfPoint3Dist(vt[1], v) <= BF_EPS);
        BF_ASSERT(!bfPoints3Contains(verts, v));
        for (BfSize j = 0; j < 2; ++j)
          iCutNew[j] = bfPoints3Find(cutVerts, vt[j]);
        BF_ASSERT(BF_SIZE_OK(iCutNew[0]) ^ BF_SIZE_OK(iCutNew[1]));
        if (iCutNew[0] == BF_SIZE_BAD_VALUE)
          iCutNew[0] = iCutNew[1];
        else
          iCutNew[1] = iCutNew[0];
        for (BfSize j = 0; j < 2; ++j) iCutNew[j] += verts->size;
      } else {
        for (BfSize j = 0; j < 2; ++j) {
          if (isnan(t[j])) {
            iCutNew[j] = bfPoints3Find(verts, vt[j]);
          } else if (!coalesced[j]) {
            BF_ASSERT(iCut[j] != BF_SIZE_BAD_VALUE);
            iCutNew[j] = verts->size + iCut[j];
          } else {
            iCutNew[j] = BF_SIZE_BAD_VALUE;
          }
        }
      }

      /** Add new faces. How we do this depends on whether any of the
       ** vertices have coalesced. */

      /* None of the vertices have coalesced: */
      if (!coalesced[0] && !coalesced[1] && !coalesced[2]) {
        BfSize3 newFaces[2] = {
          {i0New[0], iCutNew[0], iCutNew[1]},
          {i0New[0], i0New[1], iCutNew[1]}
        };
        appendFace(verts, cutVerts, &numFaces, &facesCapacity, &faces, newFaces[0]);
        appendFace(verts, cutVerts, &numFaces, &facesCapacity, &faces, newFaces[1]);
      }

      /* The second cut vertex has coalesced: */
      else if (!coalesced[0] && coalesced[1] && !coalesced[2]) {
        BfSize3 newFace = {i0New[0], i0New[1], iCutNew[0]};
        appendFace(verts, cutVerts, &numFaces, &facesCapacity, &faces, newFace);
      }

      /* The two cut vertices have coalesced with each other: */
      else if (!coalesced[0] && !coalesced[1] && coalesced[2]) {
        BF_ASSERT(iCutNew[0] == iCutNew[1]);
        BfSize3 newFace = {i0New[0], i0New[1], iCutNew[0]};
        appendFace(verts, cutVerts, &numFaces, &facesCapacity, &faces, newFace);
      }

      else BF_DIE();
    }

    else if (numPos == 1 && numNeg == 1 && numZero == 1) {
      BfSize iBd = BF_SIZE_BAD_VALUE; /* Index of vertex on boundary */
      BfSize iInt = BF_SIZE_BAD_VALUE; /* Index of vertex in interior */
      BfReal phiInt = BF_NAN; /* Value of phi at interior vertex */

      if (phiPos[0] == BF_EPS && trimesh->isBoundaryVert[iPos[0]]) {
        iBd = iPos[0];
        iInt = iNeg[0];
        phiInt = phiNeg[0];
      } else if (phiNeg[0] == -BF_EPS && trimesh->isBoundaryVert[iNeg[0]]) {
        iBd = iNeg[0];
        iInt = iPos[0];
        phiInt = phiPos[0];
      } else {
        BF_DIE();
      }

      BfSize2 edge = {iBd, iInt};
      SORT2(edge[0], edge[1]);

      BfSize edgeInd = bfTrimeshGetEdgeIndex(trimesh, edge);
      BF_ASSERT(edgeInd != BF_SIZE_BAD_VALUE);

      BfSize2 incFaceInds = {BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE};
      BfSize numIncFaces = bfTrimeshGetFacesIncOnEdge(trimesh, edgeInd, incFaceInds);
      BF_ASSERT(numIncFaces == 2);

      BfSize iOp = BF_SIZE_BAD_VALUE;
      for (BfSize j = 0; j < 2; ++j) {
        BfSize const *f = trimesh->faces[incFaceInds[j]];
        for (BfSize k = 0; k < 3; ++k) {
          if (f[k] != iZero[0] && f[k] != iBd && f[k] != iInt) {
            BF_ASSERT(iOp == BF_SIZE_BAD_VALUE);
            iOp = f[k];
          }
        }
      }
      BF_ASSERT(iOp != BF_SIZE_BAD_VALUE);

      BfReal phiOp = bfVecRealGetElt(phi, iOp);
      BF_ASSERT((phiInt < 0 && phiOp > 0) ^ (phiInt > 0 && phiOp < 0));

      BfReal t = -phiInt/(phiOp - phiInt);
      BF_ASSERT(0 < t && t < 1);

      BfPoint3 vInt, vBd, vZero, vOp;
      bfTrimeshGetVertex(trimesh, iInt, vInt);
      bfTrimeshGetVertex(trimesh, iBd, vBd);
      bfTrimeshGetVertex(trimesh, iZero[0], vZero);
      bfTrimeshGetVertex(trimesh, iOp, vOp);

      /* Find the first cut vertex (`vt`): */

      BfVector3 dv;
      BfPoint3 vt;
      bfPoint3Sub(vOp, vInt, dv);
      bfPoint3GetPointOnRay(vInt, dv, t, vt);

      dfdsContext context = {.vInt = vInt, .vZero = vZero, .vt = vt};
      bfPoint3Sub(vBd, vInt, context.dv);

      /* Find the second cut vertex (`vs`): */

      BfReal s = BF_NAN;
      bool foundZero = bfFindZeroOnInterval(
        (BfReal (*)(BfReal, void *))dfds, 0, 1, &context, &s);
      if (foundZero) {
        BF_ASSERT(0 < s && s < 1);
      } else {
        s = bfPoint3Dist(vInt, vt) < bfPoint3Dist(vBd, vt) ? 0 : 1;
      }

      BfPoint3 vs;
      bfPoint3GetPointOnRay(vInt, context.dv, s, vs);

      if (s == 1) {
        BF_ASSERT(bfPoint3Equal(vs, vBd));

        /* Can't add a new face in this case---skip it */
        continue;
      }

      // TODO: need to figure out what to do here
      BF_ASSERT(0 < s && s < 1);

      /* Add new cut edges corresponding to `vt` and `vs`: */

      // TODO: don't actually need to add `vt` here---when we traverse
      // the opposite face, it will get added anyway...

      BfReal const *vCut[2] = {vt, vs};

      BfSize2 cutEdge[2] = {
        /* vt: */ {iInt, iOp},
        /* vs: */ {iInt, iBd}
      };
      SORT2(cutEdge[0][0], cutEdge[0][1]);
      SORT2(cutEdge[1][0], cutEdge[1][1]);

      bool addedCutEdge[2];

      BfSize iCut[2] = {BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE};
      for (BfSize j = 0; j < 2; ++j) {
        addedCutEdge[j] = appendCutEdge(&numCutEdges,&cutEdgesCapacity,&cutEdges,cutEdge[j]);
        bool addedCutVertAlready = bfPoints3Contains(cutVerts, vCut[j]);
        if (addedCutEdge[j]) {
          BF_ASSERT(!addedCutVertAlready);
          iCut[j] = bfPoints3GetSize(cutVerts);
          bfPoints3Append(cutVerts, vCut[j]);
        } else if (addedCutVertAlready) {
          iCut[j] = bfPoints3Find(cutVerts, vCut[j]);
          BF_ASSERT(iCut[j] != BF_SIZE_BAD_VALUE);
        } else {
          /* We added a cut vertex for this edge previously, but it no
           * longer matches---update it, since we have an improved cut
           * vertex now */
          BF_ASSERT(numCutEdges == cutVerts->size);
          BfSize k = 0;
          for (; k < numCutEdges; ++k)
            if (cutEdges[k][0] == cutEdge[j][0] && cutEdges[k][1] == cutEdge[j][1])
              break;
          BF_ASSERT(k < numCutEdges);
          bfPoints3Set(cutVerts, k, vCut[j]);
          iCut[j] = k;
        }
      }

      /* Add new cut face: */

      BfSize iZeroNew = bfPoints3Find(verts, vZero);
      BF_ASSERT(iZeroNew != BF_SIZE_BAD_VALUE);

      BfSize iNew = bfPoints3Find(verts, phiInt < 0 ? vInt : vBd);
      BF_ASSERT(iNew != BF_SIZE_BAD_VALUE);

      BfSize3 newFace = {iZeroNew, verts->size + iCut[1], iNew};
      appendFace(verts, cutVerts, &numFaces, &facesCapacity, &faces, newFace);
    }

    else BF_DIE();
  }

  BF_ASSERT(bfPoints3AllUnique(verts));
  BF_ASSERT(bfPoints3AllUnique(cutVerts));

  bfPoints3Extend(verts, cutVerts);
  HANDLE_ERROR();

  BF_ASSERT(bfPoints3AllUnique(verts));

  /* Eliminate isolated boundary vertices. These can arise during the
   * splitting process. */

  // TODO: ideally, there should be no isolated vertices!

  bool *isolated = bfMemAlloc(verts->size, sizeof(bool));
  HANDLE_ERROR();

  for (BfSize i = 0; i < verts->size; ++i) isolated[i] = true;

  for (BfSize i = 0; i < numFaces; ++i) {
    for (BfSize j = 0; j < 3; ++j) {
      BfSize k = faces[i][j];
      BF_ASSERT(k < verts->size);
      isolated[k] = false;
    }
  }

  for (BfSize i = 0; i < verts->size; ++i)
    if (isolated[i])
      BF_ASSERT(i >= perm->size);

  for (BfSize i = verts->size; i > 0; --i) {
    if (!isolated[i - 1]) continue;

    bfPoints3Delete(verts, i - 1);

    for (BfSize j = 0; j < numFaces; ++j) {
      for (BfSize k = 0; k < 3; ++k) {
        BF_ASSERT(faces[j][k] != i - 1);
        if (faces[j][k] >= i) --faces[j][k];
      }
    }
  }

  FILE *fp = fopen("verts.bin", "w");
  fwrite(verts->data, sizeof(BfPoint3), verts->size, fp);
  fclose(fp);

  fp = fopen("faces.bin", "w");
  fwrite(faces, sizeof(BfSize3), numFaces, fp);
  fclose(fp);

  submesh = bfTrimeshNewFromVertsAndFaces(verts, numFaces, (BfSize3 const *)faces);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMemFree(isolated);

  if (permPtr != NULL)
    *permPtr = perm;

  return submesh;
}
