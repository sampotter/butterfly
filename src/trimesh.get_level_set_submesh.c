#include <bf/trimesh.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/real_array.h>
#include <bf/size_array.h>
#include <bf/util.h>
#include <bf/vectors.h>

#include "macros.h"

struct NodalDomainBuilder {
  BfTrimesh const *trimesh;
  BfRealArray const *phi;
  BfReal tol;

  /** Variables used for building the new nodal domain mesh: */

  BfSize numFaces;
  BfSize facesCapacity;
  BfSize3 *faces;

  BfPoints3 *verts;

  BfSize numCutEdges;
  BfSize cutEdgesCapacity;
  BfSize2 *cutEdges;

  BfPoints3 *cutVerts;

  BfSizeArray *perm;
  bool const *permMask;

  /** Variables used for each iteration of `addCutFacesAndVerts`: */

  BfSize numPos;
  BfSize numNeg;
  BfSize numZero;
  BfSize iPos[2];
  BfSize iNeg[2];
  BfSize iZero[2];
  BfReal phiPos[2];
  BfReal phiNeg[2];
};

typedef struct NodalDomainBuilder NodalDomainBuilder;

/* Check whether there are any faces for which `values` is equal to
 * zero at all of the face's vertices. */
static bool trimeshHasNodalFaces(BfTrimesh const *trimesh, BfRealArray const *values) {
  for (BfSize i = 0; i < trimesh->numFaces; ++i) {
    BfReal valuesFace[3];
    bfRealArrayGetValues(values, 3, trimesh->faces[i], valuesFace);
    if (valuesFace[0] == 0 && valuesFace[1] == 0 && valuesFace[2] == 0)
      return true;
  }
  return false;
}

static void init(NodalDomainBuilder *builder,
                 BfTrimesh const *trimesh, BfRealArray const *phi, BfReal tol,
                 bool const *permMask, BfSizeArray **permPtr) {
  BF_ERROR_BEGIN();

  builder->trimesh = trimesh;
  builder->phi = phi;
  builder->tol = tol;

  builder->numFaces = 0;
  builder->facesCapacity = BF_ARRAY_DEFAULT_CAPACITY;
  builder->faces = bfMemAlloc(builder->facesCapacity, sizeof(BfSize3));
  HANDLE_ERROR();

  builder->verts = bfPoints3NewWithDefaultCapacity();
  HANDLE_ERROR();

  builder->numCutEdges = 0;
  builder->cutEdgesCapacity = BF_ARRAY_DEFAULT_CAPACITY;
  builder->cutEdges = bfMemAlloc(builder->cutEdgesCapacity, sizeof(BfSize2));
  HANDLE_ERROR();

  builder->cutVerts = bfPoints3NewWithDefaultCapacity();
  HANDLE_ERROR();

  /* Check whether caller wants us to fill the permutation array
   * (permPtr != NULL). If it does, then check whether we need to
   * allocate our own array or use one provided by caller. */
  builder->perm = NULL;
  if (permPtr != NULL) {
    if (*permPtr == NULL) {
      builder->perm = bfSizeArrayNewWithDefaultCapacity();
      HANDLE_ERROR();
    } else {
      builder->perm = *permPtr;
    }
  }

  builder->permMask = permMask;

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void addVertsAndFillPerm(NodalDomainBuilder *builder) {
  BF_ERROR_BEGIN();

  for (BfSize i = 0; i < bfTrimeshGetNumVerts(builder->trimesh); ++i) {
    if (bfRealArrayGetValue(builder->phi, i) > 0) continue;

    /* Check if this vertex is isolated: */
    bool isIsolated = true;
    for (BfSize j = builder->trimesh->vvOffset[i]; j < builder->trimesh->vvOffset[i + 1]; ++j) {
      if (bfRealArrayGetValue(builder->phi, builder->trimesh->vv[j]) <= 0) {
        isIsolated = false;
        break;
      }
    }
    if (isIsolated) continue; /* (don't add it if it is) */

    BfReal const *v = bfTrimeshGetVertPtrConst(builder->trimesh, i);

    bfPoints3Append(builder->verts, v);
    HANDLE_ERROR();

    if (builder->perm != NULL && builder->permMask[i]) {
      bfSizeArrayAppend(builder->perm, i);
      HANDLE_ERROR();
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

static BfReal getTriArea(BfPoint3 const v0, BfPoint3 const v1, BfPoint3 const v2) {
  BfVector3 dv1, dv2;
  bfPoint3Sub(v1, v0, dv1);
  bfPoint3Sub(v2, v0, dv2);

  BfVector3 n;
  bfVector3Cross(dv1, dv2, n);

  BfReal triArea = bfVector3Norm(n)/2;

  return triArea;
}

static void appendFace(NodalDomainBuilder *builder, BfSize3 const face) {
  BF_ERROR_BEGIN();

  /* Make sure the face isn't degenerate: */

  for (BfSize i = 0; i < 3; ++i) BF_ASSERT(BF_SIZE_OK(face[i]));

  BF_ASSERT(face[0] != face[1] && face[1] != face[2] && face[2] != face[0]);

  BfReal const *v[3] = {NULL, NULL, NULL};
  for (BfSize i = 0; i < 3; ++i)
    v[i] = face[i] < builder->verts->size ?
      builder->verts->data[face[i]] :
      builder->cutVerts->data[face[i] - builder->verts->size];

  BfReal triArea = getTriArea(v[0], v[1], v[2]);
  BF_ASSERT(triArea > 0);

  /* Append the face: */

  if (builder->numFaces == builder->facesCapacity) {
    BfSize newCapacity = builder->facesCapacity;
    newCapacity *= 2;

    BfSize3 *newFaces = bfMemRealloc(builder->faces, newCapacity, sizeof(BfSize3));
    HANDLE_ERROR();

    builder->faces = newFaces;
    builder->facesCapacity = newCapacity;
  }

  bfMemCopy(face, 1, sizeof(BfSize3), builder->faces + builder->numFaces);

  ++builder->numFaces;

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void addContainedFaces(NodalDomainBuilder *builder) {
  BF_ERROR_BEGIN();

  /* Accumulate the interior faces themselves. We do this now that
   * we've found all of the vertices from trimesh which will be
   * included in submesh, so we can reindex the faces accordingly. */
  for (BfSize i = 0; i < builder->trimesh->numFaces; ++i) {
    BfReal phiFace[3];
    for (BfSize j = 0; j < 3; ++j)
      phiFace[j] = bfRealArrayGetValue(builder->phi, builder->trimesh->faces[i][j]);

    if (!(phiFace[0] <= 0 && phiFace[1] <= 0 && phiFace[2] <= 0))
      continue;

    BfSize3 newFace;
    for (BfSize j = 0; j < 3; ++j) {
      BfReal const *v = bfPoints3GetPtrConst(builder->trimesh->verts, builder->trimesh->faces[i][j]);
      newFace[j] = bfPoints3FindApprox(builder->verts, v, builder->tol);
      BF_ASSERT(newFace[j] != BF_SIZE_BAD_VALUE);
    }

    appendFace(builder, newFace);
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void invalidateCutFacesAndVertsVars(NodalDomainBuilder *builder) {
  builder->numPos = builder->numNeg = builder->numZero
    = builder->iPos[0] = builder->iPos[1]
    = builder->iNeg[0] = builder->iNeg[1]
    = builder->iZero[0] = builder->iZero[1] = BF_SIZE_BAD_VALUE;

  builder->phiPos[0] = builder->phiPos[1]
    = builder->phiNeg[0] = builder->phiNeg[1] = BF_NAN;
}

static bool prepareForAddCutFacesAndVertsIteration(NodalDomainBuilder *builder, BfSize faceIndex) {
  invalidateCutFacesAndVertsVars(builder);

  BfReal phiFace[3];
  for (BfSize j = 0; j < 3; ++j)
    phiFace[j] = bfRealArrayGetValue(builder->phi, builder->trimesh->faces[faceIndex][j]);

  /* We've already established that there should be no nodal faces: */
  if (phiFace[0] == 0 && phiFace[1] == 0 && phiFace[2] == 0) BF_DIE();

  /* Skip faces that aren't cut by the zero level: */
  if ((phiFace[0] <= 0 && phiFace[1] <= 0 && phiFace[2] <= 0)
    || (phiFace[0] >= 0 && phiFace[1] >= 0 && phiFace[2] >= 0))
    return false;

  /* Figure out which verts are positive, negative, and zero: */
  builder->numPos = builder->numNeg = builder->numZero = 0;
  for (BfSize j = 0; j < 3; ++j) {
    BF_ASSERT(builder->numPos <= 2 && builder->numNeg <= 2 && builder->numZero <= 1);
    if (phiFace[j] > 0) {
      builder->phiPos[builder->numPos] = phiFace[j];
      builder->iPos[builder->numPos++] = builder->trimesh->faces[faceIndex][j];
    } else if (phiFace[j] < 0) {
      builder->phiNeg[builder->numNeg] = phiFace[j];
      builder->iNeg[builder->numNeg++] = builder->trimesh->faces[faceIndex][j];
    } else if (phiFace[j] == 0) {
      builder->iZero[builder->numZero++] = builder->trimesh->faces[faceIndex][j];
    } else {
      BF_DIE();
    }
  }

  return true;
}

static bool appendCutEdge(NodalDomainBuilder *builder, BfSize2 cutEdge) {
  BF_ERROR_BEGIN();

  for (BfSize i = 0; i < builder->numCutEdges; ++i) {
    BfSize const *otherCutEdge = builder->cutEdges[i];
    if (cutEdge[0] == otherCutEdge[0] && cutEdge[1] == otherCutEdge[1])
      return false;
  }

  if (builder->numCutEdges == builder->cutEdgesCapacity) {
    BfSize newCapacity = builder->cutEdgesCapacity;
    newCapacity *= 2;

    BfSize2 *newCutEdges = bfMemRealloc(builder->cutEdges, newCapacity, sizeof(BfSize2));
    HANDLE_ERROR();

    builder->cutEdges = newCutEdges;
    builder->cutEdgesCapacity = newCapacity;
  }

  bfMemCopy(cutEdge, 1, sizeof(BfSize2), builder->cutEdges + builder->numCutEdges);

  ++builder->numCutEdges;

  BF_ERROR_END() {
    BF_DIE();
  }

  return true;
}

static BfSize appendCutVertex(NodalDomainBuilder *builder, BfSize i0, BfSize i1, BfReal t, BfPoint3 const v) {
  BF_DEFINE_EDGE(cutEdge, i0, i1);

  BfSize iCut = BF_SIZE_BAD_VALUE;

  if (appendCutEdge(builder, cutEdge)) {
    if (bfPoints3ContainsApprox(builder->cutVerts, v, builder->tol)) {
      /* We might have already added this cut vertex---in this
       * case, we should make sure that there's already a cut edge
       * containing the index of the duplicated vertex. */
      BF_ASSERT(fabs(t) <= 2.1e2*builder->tol || fabs(1 - t) <= 2.1e2*builder->tol);
      BfSize iAdded = t < 0.5 ? cutEdge[0] : cutEdge[1];
      bool foundExistingCutVert = false;
      for (BfSize k = 0; k < builder->numCutEdges; ++k)
        if (builder->cutEdges[k][0] == iAdded || builder->cutEdges[k][1] == iAdded)
          foundExistingCutVert = true;
      BF_ASSERT(foundExistingCutVert);

      /* Get the index of `v` in `cutVerts` since we can
       * conclude now that it's already been added: */
      iCut = bfPoints3FindApprox(builder->cutVerts, v, builder->tol);
    } else {
      /* If we haven't added `v` already, then add it now. */
      iCut = bfPoints3GetSize(builder->cutVerts);
      bfPoints3Append(builder->cutVerts, v);
    }
  } else {
    BfSize vtIndex = bfPoints3FindApprox(builder->cutVerts, v, builder->tol);
    BF_ASSERT(vtIndex != BF_SIZE_BAD_VALUE);
    iCut = vtIndex;
  }

  return iCut;
}

static void addCutFacesAndVerts_case21(NodalDomainBuilder *builder) {
  BfPoint3 v0, v[2];
  bfTrimeshGetVertex(builder->trimesh, builder->iNeg[0], v0);
  bfTrimeshGetVertex(builder->trimesh, builder->iPos[0], v[0]);
  bfTrimeshGetVertex(builder->trimesh, builder->iPos[1], v[1]);

  /* Compute `t` and `vt`. Afterwards, we'll have determined the
   * coordinates of three face vertices. */
  BfReal t[2];
  BfPoint3 vt[2];
  for (BfSize j = 0; j < 2; ++j) {
    t[j] = -builder->phiNeg[0]/(builder->phiPos[j] - builder->phiNeg[0]);
    BF_ASSERT(0 <= t[j] && t[j] <= 1);

    BfPoint3 dv;
    bfPoint3Sub(v[j], v0, dv);
    bfPoint3GetPointOnRay(v0, dv, t[j], vt[j]);

    /* Due to numerical roundoff, `vt[j]` may exactly equal an existing
     * mesh vertex. We want to avoid duplicate mesh vertices, so
     * we set `t[j]` to `NAN` to signal this and continue. */
    if (bfPoints3ContainsApprox(builder->verts, vt[j], builder->tol))
      t[j] = BF_NAN;
  }

  /* Check if this is a degenerate face: if it is, we return now
   * before adding anything to any of the arrays we're building up. */
  bool coalesced[3] = {
    bfPoint3Dist(v0, vt[0]) <= builder->tol,
    bfPoint3Dist(v0, vt[1]) <= builder->tol,
    bfPoint3Dist(vt[0], vt[1]) <= builder->tol
  };
  if (coalesced[0] || coalesced[1] || coalesced[2])
    return;

  /* Compute new cut vertices and append them to `cutVerts`: */
  BfSize iCut[2] = {BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE};
  for (BfSize j = 0; j < 2; ++j) {
    if (isnan(t[j]))
      continue;
    iCut[j] = appendCutVertex(builder, builder->iNeg[0], builder->iPos[j], t[j], vt[j]);
  }

  /* Find the "new" index of `v0` (its position in `verts`): */
  BfSize i0New = bfPoints3FindApprox(builder->verts, v0, builder->tol);
  if (i0New == BF_SIZE_BAD_VALUE) {
    i0New = bfPoints3GetSize(builder->verts);
    bfPoints3Append(builder->verts, v0);
  }

  /* Find the "new" indices of the other two verts: */
  BfSize iCutNew[2];
  for (BfSize j = 0; j < 2; ++j) {
    if (isnan(t[j])) {
      iCutNew[j] = bfPoints3FindApprox(builder->verts, vt[j], builder->tol);
    } else {
      BF_ASSERT(iCut[j] != BF_SIZE_BAD_VALUE);
      iCutNew[j] = builder->verts->size + iCut[j];
    }
  }
  BF_ASSERT(iCutNew[0] != BF_SIZE_BAD_VALUE && iCutNew[1] != BF_SIZE_BAD_VALUE);

  /* Add a new face: */
  BfSize3 newFace = {i0New, iCutNew[0], iCutNew[1]};
  appendFace(builder, newFace);
}

static void addCutFacesAndVerts_case12(NodalDomainBuilder *builder) {
  BfPoint3 v0[2], v;
  bfTrimeshGetVertex(builder->trimesh, builder->iNeg[0], v0[0]);
  bfTrimeshGetVertex(builder->trimesh, builder->iNeg[1], v0[1]);
  bfTrimeshGetVertex(builder->trimesh, builder->iPos[0], v);

  /* Compute `t` and `vt`. Afterwards, we'll have determined the
   * coordinates of three face vertices. */
  BfReal t[2];
  BfPoint3 vt[2];
  for (BfSize j = 0; j < 2; ++j) {
    t[j] = -builder->phiNeg[j]/(builder->phiPos[0] - builder->phiNeg[j]);
    BF_ASSERT(0 <= t[j] && t[j] <= 1);

    BfPoint3 dv;
    bfPoint3Sub(v, v0[j], dv);
    bfPoint3GetPointOnRay(v0[j], dv, t[j], vt[j]);

    /* Due to numerical roundoff, `vt[j]` may exactly equal an existing
     * mesh vertex. We want to avoid duplicate mesh vertices, so
     * we set `t[j]` to `NAN` to signal this and continue. */
    if (bfPoints3ContainsApprox(builder->verts, vt[j], builder->tol))
      t[j] = BF_NAN;
  }

  bool coalesced[3] = {
    bfPoint3Dist(v0[0], vt[0]) <= builder->tol,
    bfPoint3Dist(v0[1], vt[1]) <= builder->tol,
    bfPoint3Dist(vt[0], vt[1]) <= builder->tol
  };

  // TODO: quad collapsed to an edge:
  BF_ASSERT(!(coalesced[0] && coalesced[1]));

  if (coalesced[2]) BF_ASSERT(!coalesced[0] && !coalesced[1]);

  /* Compute new cut vertices and append them to `cutVerts`: */
  BfSize iCut[2] = {BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE};
  if (coalesced[2]) {
    bool foundVt[2] = {
      bfPoints3ContainsApprox(builder->cutVerts, vt[0], builder->tol),
      bfPoints3ContainsApprox(builder->cutVerts, vt[1], builder->tol),
    };
    BF_ASSERT(!(foundVt[0] ^ foundVt[1]));
    if (!foundVt[0] || !foundVt[1]) {
      iCut[0] = iCut[1] = bfPoints3GetSize(builder->cutVerts);
      bfPoints3Append(builder->cutVerts, vt[0]);
    }
  } else {
    for (BfSize j = 0; j < 2; ++j) {
      if (coalesced[j] || isnan(t[j]))
        continue;
      iCut[j] = appendCutVertex(builder, builder->iNeg[j], builder->iPos[0], t[j], vt[j]);
    }
  }

  /* Find the "new" indices of these two mesh verts: */
  BfSize i0New[2] = {
    bfPoints3FindApprox(builder->verts, v0[0], builder->tol),
    bfPoints3FindApprox(builder->verts, v0[1], builder->tol),
  };
  BF_ASSERT(i0New[0] != BF_SIZE_BAD_VALUE && i0New[1] != BF_SIZE_BAD_VALUE);

  /* Find the "new" indices of the cut verts: */
  BfSize iCutNew[2];
  if (coalesced[2]) {
    BF_ASSERT(!coalesced[0] && !coalesced[1]);
    BF_ASSERT(!bfPoints3ContainsApprox(builder->verts, v, builder->tol));
    for (BfSize j = 0; j < 2; ++j)
      iCutNew[j] = bfPoints3FindApprox(builder->cutVerts, vt[j], builder->tol);
    BF_ASSERT(BF_SIZE_OK(iCutNew[0]) || BF_SIZE_OK(iCutNew[1]));
    if (iCutNew[0] == BF_SIZE_BAD_VALUE) iCutNew[0] = iCutNew[1];
    if (iCutNew[1] == BF_SIZE_BAD_VALUE) iCutNew[1] = iCutNew[0];
    BF_ASSERT(BF_SIZE_OK(iCutNew[0]) && BF_SIZE_OK(iCutNew[1]));
    for (BfSize j = 0; j < 2; ++j) iCutNew[j] += builder->verts->size;
  } else {
    for (BfSize j = 0; j < 2; ++j) {
      if (isnan(t[j])) {
        iCutNew[j] = bfPoints3FindApprox(builder->verts, vt[j], builder->tol);
      } else if (!coalesced[j]) {
        BF_ASSERT(iCut[j] != BF_SIZE_BAD_VALUE);
        iCutNew[j] = builder->verts->size + iCut[j];
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
    appendFace(builder, newFaces[0]);
    appendFace(builder, newFaces[1]);
  }

  /* The second cut vertex has coalesced: */
  else if (!coalesced[0] && coalesced[1] && !coalesced[2]) {
    BfSize3 newFace = {i0New[0], i0New[1], iCutNew[0]};
    appendFace(builder, newFace);
  }

  /* The two cut vertices have coalesced with each other: */
  else if (!coalesced[0] && !coalesced[1] && coalesced[2]) {
    BF_ASSERT(iCutNew[0] == iCutNew[1]);
    BfSize3 newFace = {i0New[0], i0New[1], iCutNew[0]};
    appendFace(builder, newFace);
  }

  else BF_DIE();
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

static void addCutFacesAndVerts_case111(NodalDomainBuilder *builder) {
  BfSize iBd = BF_SIZE_BAD_VALUE; /* Index of vertex on boundary */
  BfSize iInt = BF_SIZE_BAD_VALUE; /* Index of vertex in interior */
  BfReal phiInt = BF_NAN; /* Value of phi at interior vertex */

  /* TODO: instead of checking whether the value of phi is equal to
   * +/- tol, we should use a mask to mark these points (these are the
   * points on the boundary where phi *was* equal to zero, but which
   * have been perturbed to push them into the correct nodal
   * domains)
   *
   * TODO: even better, we should be able to just leave these as
   * zeros without perturbing... */
  if (builder->phiPos[0] == BF_EPS && builder->trimesh->isBoundaryVert[builder->iPos[0]]) {
    iBd = builder->iPos[0];
    iInt = builder->iNeg[0];
    phiInt = builder->phiNeg[0];
  } else if (builder->phiNeg[0] == -BF_EPS && builder->trimesh->isBoundaryVert[builder->iNeg[0]]) {
    iBd = builder->iNeg[0];
    iInt = builder->iPos[0];
    phiInt = builder->phiPos[0];
  } else {
    BF_DIE();
  }

  BfSize2 edge = {iBd, iInt};
  SORT2(edge[0], edge[1]);

  BfSize edgeInd = bfTrimeshGetEdgeIndex(builder->trimesh, edge);
  BF_ASSERT(edgeInd != BF_SIZE_BAD_VALUE);

  BfSize2 incFaceInds = {BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE};
  BfSize numIncFaces = bfTrimeshGetFacesIncOnEdge(builder->trimesh, edgeInd, incFaceInds);
  BF_ASSERT(numIncFaces == 2);

  BfSize iOp = BF_SIZE_BAD_VALUE;
  for (BfSize j = 0; j < 2; ++j) {
    BfSize const *f = builder->trimesh->faces[incFaceInds[j]];
    for (BfSize k = 0; k < 3; ++k) {
      if (f[k] != builder->iZero[0] && f[k] != iBd && f[k] != iInt) {
        BF_ASSERT(iOp == BF_SIZE_BAD_VALUE);
        iOp = f[k];
      }
    }
  }
  BF_ASSERT(iOp != BF_SIZE_BAD_VALUE);

  BfReal phiOp = bfRealArrayGetValue(builder->phi, iOp);
  BF_ASSERT((phiInt < 0 && phiOp > 0) ^ (phiInt > 0 && phiOp < 0));

  BfReal t = -phiInt/(phiOp - phiInt);
  BF_ASSERT(0 < t && t < 1);

  BfPoint3 vInt, vBd, vZero, vOp;
  bfTrimeshGetVertex(builder->trimesh, iInt, vInt);
  bfTrimeshGetVertex(builder->trimesh, iBd, vBd);
  bfTrimeshGetVertex(builder->trimesh, builder->iZero[0], vZero);
  bfTrimeshGetVertex(builder->trimesh, iOp, vOp);

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

    /* Don't need to add a new face in this case */
    return;
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
    addedCutEdge[j] = appendCutEdge(builder, cutEdge[j]);
    bool addedCutVertAlready = bfPoints3ContainsApprox(builder->cutVerts, vCut[j], builder->tol);
    if (addedCutEdge[j]) {
      BF_ASSERT(!addedCutVertAlready);
      iCut[j] = bfPoints3GetSize(builder->cutVerts);
      bfPoints3Append(builder->cutVerts, vCut[j]);
    } else if (addedCutVertAlready) {
      iCut[j] = bfPoints3FindApprox(builder->cutVerts, vCut[j], builder->tol);
      BF_ASSERT(iCut[j] != BF_SIZE_BAD_VALUE);
    } else {
      /* We added a cut vertex for this edge previously, but it no
       * longer matches---update it, since we have an improved cut
       * vertex now */
      BF_ASSERT(builder->numCutEdges == builder->cutVerts->size);
      BfSize k = 0;
      for (; k < builder->numCutEdges; ++k)
        if (builder->cutEdges[k][0] == cutEdge[j][0] && builder->cutEdges[k][1] == cutEdge[j][1])
          break;
      BF_ASSERT(k < builder->numCutEdges);
      bfPoints3Set(builder->cutVerts, k, vCut[j]);
      iCut[j] = k;
    }
  }

  /* Add new cut face: */

  BfSize iZeroNew = bfPoints3FindApprox(builder->verts, vZero, builder->tol);
  BF_ASSERT(iZeroNew != BF_SIZE_BAD_VALUE);

  BfSize iNew = bfPoints3FindApprox(builder->verts, phiInt < 0 ? vInt : vBd, builder->tol);
  BF_ASSERT(iNew != BF_SIZE_BAD_VALUE);

  BfSize3 newFace = {iZeroNew, builder->verts->size + iCut[1], iNew};
  appendFace(builder, newFace);
}

/* Find the zero level set and accumulate new cut faces. */
static void addCutFacesAndVerts(NodalDomainBuilder *builder) {
  BF_ERROR_BEGIN();

  for (BfSize faceIndex = 0; faceIndex < builder->trimesh->numFaces; ++faceIndex) {
    if (!prepareForAddCutFacesAndVertsIteration(builder, faceIndex))
      continue;
    if (builder->numPos == 2 && builder->numNeg == 1)
      addCutFacesAndVerts_case21(builder);
    else if (builder->numPos == 1 && builder->numNeg == 2)
      addCutFacesAndVerts_case12(builder);
    else if (builder->numPos == 1 && builder->numNeg == 1 && builder->numZero == 1)
      addCutFacesAndVerts_case111(builder);
    else BF_DIE();
    HANDLE_ERROR();
  }

  BF_ASSERT(bfPoints3AllUniqueApprox(builder->verts, builder->tol));
  BF_ASSERT(bfPoints3AllUniqueApprox(builder->cutVerts, builder->tol));

  bfPoints3Extend(builder->verts, builder->cutVerts);
  HANDLE_ERROR();

  BF_ASSERT(bfPoints3AllUniqueApprox(builder->verts, builder->tol));

  BF_ERROR_END() {
    BF_DIE();
  }
}

// TODO: ideally, there should be no isolated vertices!
/* Eliminate isolated boundary vertices. These can arise during the
 * splitting process. */
static void eliminateIsolatedVerts(NodalDomainBuilder *builder) {
  BF_ERROR_BEGIN();

  bool *isolated = bfMemAlloc(builder->verts->size, sizeof(bool));
  HANDLE_ERROR();

  for (BfSize i = 0; i < builder->verts->size; ++i) isolated[i] = true;

  for (BfSize i = 0; i < builder->numFaces; ++i) {
    for (BfSize j = 0; j < 3; ++j) {
      BfSize k = builder->faces[i][j];
      BF_ASSERT(k < builder->verts->size);
      isolated[k] = false;
    }
  }

  for (BfSize i = 0; i < builder->verts->size; ++i)
    if (isolated[i])
      BF_ASSERT(i >= builder->perm->size);

  for (BfSize i = builder->verts->size; i > 0; --i) {
    if (!isolated[i - 1]) continue;

    bfPoints3Delete(builder->verts, i - 1);

    for (BfSize j = 0; j < builder->numFaces; ++j) {
      for (BfSize k = 0; k < 3; ++k) {
        BF_ASSERT(builder->faces[j][k] != i - 1);
        if (builder->faces[j][k] >= i)
          --builder->faces[j][k];
      }
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMemFree(isolated);
}

static void deinit(NodalDomainBuilder *builder) {
  builder->trimesh = NULL;
  builder->phi = NULL;
  builder->tol = BF_NAN;

  builder->numFaces = BF_SIZE_BAD_VALUE;
  builder->facesCapacity = BF_SIZE_BAD_VALUE;

  bfMemFree(builder->faces);
  builder->faces = NULL;

  bfPoints3DeinitAndDealloc(&builder->verts);

  builder->numCutEdges = BF_SIZE_BAD_VALUE;
  builder->cutEdgesCapacity = BF_SIZE_BAD_VALUE;

  bfMemFree(builder->cutEdges);
  builder->cutEdges = NULL;

  bfPoints3DeinitAndDealloc(&builder->cutVerts);

  builder->perm = NULL;

  builder->permMask = NULL;

  invalidateCutFacesAndVertsVars(builder);
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
BfTrimesh *bfTrimeshGetLevelSetSubmesh(BfTrimesh const *trimesh, BfRealArray const *phi, BfReal tol, bool const *permMask, BfSizeArray **permPtr) {
  BF_ERROR_BEGIN();

  BfTrimesh *submesh = NULL;

  NodalDomainBuilder builder;

  /* Make sure there are no faces where all three level values are
   * zero. This may not be a problem, but we want to trip this the
   * first time this happens to check for potential problems. */
  if (trimeshHasNodalFaces(trimesh, phi))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  init(&builder, trimesh, phi, tol, permMask, permPtr);
  HANDLE_ERROR();

  addVertsAndFillPerm(&builder);
  HANDLE_ERROR();

  addContainedFaces(&builder);
  HANDLE_ERROR();

  addCutFacesAndVerts(&builder);
  HANDLE_ERROR();

  eliminateIsolatedVerts(&builder);
  HANDLE_ERROR();

  submesh = bfTrimeshNewFromVertsAndFaces(builder.verts, builder.numFaces, (BfSize3 const *)builder.faces);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  if (permPtr != NULL)
    *permPtr = builder.perm;

  deinit(&builder);

  return submesh;
}
