#include <bf/fiedler_tree_node.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/logging.h>
#include <bf/mem.h>
#include <bf/real_array.h>
#include <bf/size_array.h>
#include <bf/trimesh.h>
#include <bf/util.h>
#include <bf/vectors.h>

static BfSize const MAX_NUM_CHILDREN = 2;
static BfSize const LEAF_SIZE_THRESHOLD = 16;

/** Interface: TreeNode */

static BfTreeNodeVtable TreeNodeVtable = {
  .GetType = (__typeof__(&bfTreeNodeGetType))bfFiedlerTreeNodeGetType,
  .Delete = (__typeof__(&bfTreeNodeDelete))bfFiedlerTreeNodeDelete,
};

BfType bfFiedlerTreeNodeGetType(BfFiedlerTreeNode const *treeNode) {
  (void)treeNode;
  return BF_TYPE_FIEDLER_TREE_NODE;
}

void bfFiedlerTreeNodeDelete(BfFiedlerTreeNode **fiedlerTreeNode) {
  bfFiedlerTreeNodeDeinit(*fiedlerTreeNode);
  bfFiedlerTreeNodeDealloc(fiedlerTreeNode);
}

/** Upcasting: FiedlerTreeNode -> TreeNode */

BfTreeNode *bfFiedlerTreeNodeToTreeNode(BfFiedlerTreeNode *node) {
  return &node->super;
}

/** Downcasting: TreeNode -> FiedlerTreeNode */

BfFiedlerTreeNode const *bfTreeNodeConstToFiedlerTreeNodeConst(BfTreeNode const *treeNode) {
  if (!bfTreeNodeInstanceOf(treeNode, BF_TYPE_FIEDLER_TREE_NODE)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfFiedlerTreeNode const *)treeNode;
  }
}

/** Implementation: FiedlerTreeNode */

BfFiedlerTreeNode *bfFiedlerTreeNodeNew() {
  BF_ERROR_BEGIN();

  BfFiedlerTreeNode *node = bfMemAlloc(1, sizeof(BfFiedlerTreeNode));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return node;
}

enum NodalDomainType {UNKNOWN, POSITIVE, NEGATIVE, BOUNDARY};

static BfReal getAveragedKnownPhiValue(BfTrimesh const *trimesh, BfRealArray const *phi,
                                       enum NodalDomainType const *nodalDomainType, BfSize i) {
  BF_ERROR_BEGIN();

  BfSize numNb = trimesh->vvOffset[i + 1] - trimesh->vvOffset[i];

  BfReal *w = bfMemAlloc(numNb, sizeof(BfReal));
  HANDLE_ERROR();

  for (BfSize j = 0; j < numNb; ++j) w[j] = BF_NAN;

  BfReal const *x = bfTrimeshGetVertPtrConst(trimesh, i);
  for (BfSize j = 0; j < numNb; ++j) {
    BfSize iNb = trimesh->vv[trimesh->vvOffset[i] + j];
    if (nodalDomainType[iNb] != UNKNOWN) {
      BfReal const *xNb = bfTrimeshGetVertPtrConst(trimesh, iNb);
      BfReal d = bfPoint3Dist(x, xNb);
      BF_ASSERT(d > 0);
      w[j] = 1/d;
    }
  }

  BfReal averagedPhi = 0;
  for (BfSize j = 0; j < numNb; ++j) {
    BfSize iNb = trimesh->vv[trimesh->vvOffset[i] + j];
    if (nodalDomainType[iNb] != UNKNOWN) {
      BF_ASSERT(isfinite(w[j]));
      BfReal phiNb = bfRealArrayGetValue(phi, iNb);
      averagedPhi += w[j]*phiNb;
    }
  }

  BfReal wSum = 0;
  for (BfSize j = 0; j < numNb; ++j) if(isfinite(w[j])) wSum += w[j];

  averagedPhi /= wSum;

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMemFree(w);

  return averagedPhi;
}

static BfSize getNumUnknown(BfSize numVerts, enum NodalDomainType const *nodalDomainType) {
  BfSize numUnknown = 0;
  for (BfSize i = 0; i < numVerts; ++i)
    if (nodalDomainType[i] == UNKNOWN)
      ++numUnknown;
  return numUnknown;
}

static void fillUnknown(BfTrimesh const *trimesh, BfRealArray const *phiFiedler,
                        enum NodalDomainType *nodalDomainType) {
  BF_ERROR_BEGIN();

  for (BfSize i = 0; i < bfTrimeshGetNumVerts(trimesh); ++i) {
    if (nodalDomainType[i] == UNKNOWN && bfRealArrayGetValue(phiFiedler, i) != 0) {
      BfReal phiNew = getAveragedKnownPhiValue(trimesh, phiFiedler, nodalDomainType, i);
      HANDLE_ERROR();

      if (phiNew == 0)
        continue;

      BfReal phiOld = bfRealArrayGetValue(phiFiedler, i);
      BF_ASSERT((phiOld < 0 && phiNew > 0) ^ (phiOld > 0 && phiNew < 0));

      phiFiedler->data[i] = phiNew;

      nodalDomainType[i] = phiNew > 0 ? POSITIVE : NEGATIVE;
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void dumpType(BfSize numVerts, enum NodalDomainType const *nodalDomainType) {
  FILE *fp = fopen("type.bin", "w");
  for (BfSize i = 0; i < numVerts; ++i) {
    BfSize _ = nodalDomainType[i];
    fwrite(&_, sizeof(BfSize), 1, fp);
  }
  fclose(fp);
}

static void repairTopologyOfNodalDomains(BfTrimesh const *trimesh, BfRealArray *phiFiedler) {
  BF_ERROR_BEGIN();

  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);
  BF_ASSERT(bfRealArrayGetSize(phiFiedler) == numVerts);

  enum NodalDomainType *nodalDomainType = bfMemAlloc(numVerts, sizeof(enum NodalDomainType));
  HANDLE_ERROR();

  for (BfSize i = 0; i < numVerts; ++i)
    nodalDomainType[i] = trimesh->isBoundaryVert[i] ? BOUNDARY : UNKNOWN;

  /* Find the index of the vertex with the largest Fiedler vector value: */
  BfSize positiveSeedIndex = BF_SIZE_BAD_VALUE;
  BfReal phiMax = -BF_INFINITY;
  for (BfSize i = 0; i < numVerts; ++i) {
    BfReal phi = bfRealArrayGetValue(phiFiedler, i);
    if (phi > phiMax) {
      phiMax = phi;
      positiveSeedIndex = i;
    }
  }

  /* ... And find the index of the vertex with the smallest value: */
  BfSize negativeSeedIndex = BF_SIZE_BAD_VALUE;
  BfReal phiMin = BF_INFINITY;
  for (BfSize i = 0; i < numVerts; ++i) {
    BfReal phi = bfRealArrayGetValue(phiFiedler, i);
    if (phi < phiMin) {
      phiMin = phi;
      negativeSeedIndex = i;
    }
  }

  BfSizeArray *queue = bfSizeArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  bfSizeArrayAppend(queue, positiveSeedIndex);
  HANDLE_ERROR();

  while (!bfSizeArrayIsEmpty(queue)) {
    BfSize i = bfSizeArrayGetFirst(queue);
    bfSizeArrayDelete(queue, 0);
    BF_ASSERT(bfRealArrayGetValue(phiFiedler, i) > 0);

    nodalDomainType[i] = POSITIVE;

    for (BfSize j = trimesh->vvOffset[i]; j < trimesh->vvOffset[i + 1]; ++j) {
      BfSize iNb = trimesh->vv[j];
      if (nodalDomainType[iNb] == UNKNOWN
          && bfRealArrayGetValue(phiFiedler, iNb) > 0
          && !bfSizeArrayContains(queue, iNb)) {
        bfSizeArrayAppend(queue, iNb);
        HANDLE_ERROR();
      }
    }
  }

  bfSizeArrayAppend(queue, negativeSeedIndex);
  HANDLE_ERROR();

  while (!bfSizeArrayIsEmpty(queue)) {
    BfSize i = bfSizeArrayGetFirst(queue);
    bfSizeArrayDelete(queue, 0);
    BF_ASSERT(bfRealArrayGetValue(phiFiedler, i) < 0);

    nodalDomainType[i] = NEGATIVE;

    for (BfSize j = trimesh->vvOffset[i]; j < trimesh->vvOffset[i + 1]; ++j) {
      BfSize iNb = trimesh->vv[j];
      if (nodalDomainType[iNb] == UNKNOWN
          && bfRealArrayGetValue(phiFiedler, iNb) < 0
          && !bfSizeArrayContains(queue, iNb)) {
        bfSizeArrayAppend(queue, iNb);
        HANDLE_ERROR();
      }
    }
  }

  dumpType(numVerts, nodalDomainType);

  BfSize numUnknown = getNumUnknown(numVerts, nodalDomainType);
  while (numUnknown > 0) {
    fillUnknown(trimesh, phiFiedler, nodalDomainType);
    BfSize prevNumUnknown = numUnknown;
    numUnknown = getNumUnknown(numVerts, nodalDomainType);
    BF_ASSERT(numUnknown < prevNumUnknown);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  bfSizeArrayDeinitAndDealloc(&queue);
  bfMemFree(nodalDomainType);
}

/* When we do a split in the Fiedler tree, there are two cases: either
 * the zero level is in the interior of the submesh, or it intersects
 * the boundary of the submesh in two places. Splitting along the zero
 * level set in the latter case is complicated by the fact that we use
 * low order---namely, linear---geometry and finite elements. To try
 * to mitigate this problem, we "fix" the submesh before splitting,
 * which is what this function does.
 *
 * Recall that the "Fiedler vector" is a piecewise linear function
 * which is linear over each triangle facet. Consequently, its first
 * derivatives are piecewise constant and may be undefined on the
 * vertices and interior edges of the mesh.
 *
 * We can observe that along the boundary, the normal derivative is
 * positive if a point on the boundary is adjacent to the positive
 * lobe of the "Fiedler vector", negative if it's adjacent to the
 * negative lobe, and lies on the zero level set otherwise. (I think
 * this should be true, anyway... actually, I'm not 100% sure the zero
 * level set should be orthogonal to the boundary... but it seems
 * right...)
 *
 * We first compute the normal derivative on each boundary edge of the
 * submesh. We then approximate the normal derivative on each boundary
 * /vertex/ using weighted linear interpolation. This gives a
 * piecewise linear approximation of the normal derivative over the
 * boundary. We check for the zeros of this piecewise linear
 * function---there will be either 0 or 2 zeros. We insert each of the
 * zeros which are interior to a boundary edge into the triangle mesh.
 *
 * Finally, to give `splitTrimesh` a bit more information, we perturb
 * the value of the "Fiedler vector" at each boundary vertex by +/-
 * EPS, depending on whether the normal derivative is positive or
 * negative (if it's zero, we leave it alone). HOWEVER, doing this
 * using the actual value of the normal derivative can lead to be
 * problems if we've inserted zeros.  */
static void doBoundaryFix(BfTrimesh const *trimesh, BfRealArray const *phiFiedler,
                          BfTrimesh **trimeshFixedPtr, BfRealArray **phiFiedlerFixedPtr) {
  BF_ERROR_BEGIN();

  /* Make a copy of `phiFiedler`---we'll fix it up as we go. */

  BfRealArray *phiFiedlerFixed = bfRealArrayCopy(phiFiedler);
  HANDLE_ERROR();

  repairTopologyOfNodalDomains(trimesh, phiFiedlerFixed);

  /* Estimate the normal derivative of `phiFiedler` at each boundary
   * vertex of `trimesh`. See the comments for
   * `bfTrimeshGetNormalDeriv` for more information about how exactly
   * this works. */

  BfRealArray *normalDeriv = bfTrimeshGetNormalDeriv(trimesh, phiFiedlerFixed, true);
  HANDLE_ERROR();

  /* Set `zeroInds` to the indices of the boundary edges which contain
   * zeros of the normal derivative in their interior: */

  BfSizeArray *interiorZeroInds = bfSizeArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  for (BfSize i = 0; i < trimesh->numBoundaryEdges; ++i) {
    /* Get the normal derivative at the edge endpoints */
    BfReal dudn[2];
    for (BfSize j = 0; j < 2; ++j) {
      BfSize k = trimesh->boundaryEdges[i][j];
      dudn[j] = bfRealArrayGetValue(normalDeriv, k);
    }

    /* Record this index if there's a zero in the interior of the
     * `i`th boundary edge: */
    if ((dudn[0] < 0 && dudn[1] > 0) || (dudn[1] < 0 && dudn[0] > 0))
      bfSizeArrayAppend(interiorZeroInds, i);
  }

  BfTrimesh *trimeshFixed = bfTrimeshCopy(trimesh);
  HANDLE_ERROR();

  /* Insert interior zeros now if we found any: */
  if (!bfSizeArrayIsEmpty(interiorZeroInds)) {
    for (BfSize k = 0; k < bfSizeArrayGetSize(interiorZeroInds); ++k) {
      BfSize boundaryEdgeInd = bfSizeArrayGet(interiorZeroInds, k);
      BfSize const *boundaryEdge = trimesh->boundaryEdges[boundaryEdgeInd];

      BfReal dudn[2];
      for (BfSize j = 0; j < 2; ++j)
        dudn[j] = normalDeriv->data[boundaryEdge[j]];
      BF_ASSERT((dudn[0] < 0 && dudn[1] > 0) ^ (dudn[0] > 0 && dudn[1] < 0));

      BfReal lam = -dudn[0]/(dudn[1] - dudn[0]);
      BF_ASSERT(0 < lam && lam < 1);

      BfSize newVertInd = bfTrimeshGetNumVerts(trimeshFixed);

      BfSize edgeInd = bfTrimeshGetEdgeIndex(trimeshFixed, boundaryEdge);
      bfTrimeshSplitEdge(trimeshFixed, edgeInd, lam);

      bfRealArrayInsert(normalDeriv, newVertInd, BF_NAN);

      bfRealArrayInsert(phiFiedlerFixed, newVertInd, 0);
    }
  }

  bfRealArrayDeinitAndDealloc(&normalDeriv);
  bfSizeArrayDeinitAndDealloc(&interiorZeroInds);

  /* Next, we perturb the components of the Fiedler vector by +/-
   * `BF_EPS`, depending on whether the normal derivative is positive
   * or negative. This should result in two chains of +/- `BF_EPS`
   * values of the same sign, separated by zeros.
   *
   * If we reuse the values of `normalDeriv` above, then we can wind
   * up with +/- `BF_EPS` values with the wrong sign, breaking the
   * chains. One way to avoid this is to compute a new version
   * `normalDeriv` using the fixed versions of `trimesh` and
   * `phiFiedler`. We do that first. */

  BfRealArray *normalDerivFixed
    = bfTrimeshGetNormalDeriv(trimeshFixed, phiFiedlerFixed, true);
  HANDLE_ERROR();

  for (BfSize i = 0; i < bfRealArrayGetSize(normalDerivFixed); ++i) {
    if (trimeshFixed->isBoundaryVert[i]) {
      BF_ASSERT(phiFiedlerFixed->data[i] == 0);
    } else {
      BF_ASSERT(phiFiedlerFixed->data[i] != 0);
      continue;
    }

    BfReal dudn = bfRealArrayGetValue(normalDerivFixed, i);
    if (isnan(dudn)) continue;

    if (dudn > 0)
      phiFiedlerFixed->data[i] += BF_EPS;
    else if (dudn < 0)
      phiFiedlerFixed->data[i] -= BF_EPS;
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  bfRealArrayDeinitAndDealloc(&normalDerivFixed);

#if BF_DEBUG
  for (BfSize i = 0; i < bfTrimeshGetNumFaces(trimesh); ++i)
    for (BfSize j = 0; j < 3; ++j)
      BF_ASSERT(trimesh->faces[i][j] != BF_SIZE_BAD_VALUE);
#endif

  *trimeshFixedPtr = trimeshFixed;
  *phiFiedlerFixedPtr = phiFiedlerFixed;
}

static void splitTrimesh(BfTrimesh const *trimesh, BfRealArray const *phiFiedler,
                         BfReal tol, bool const *permMask,
                         BfTrimesh **submesh1Ptr, BfTrimesh **submesh2Ptr,
                         BfSizeArray **perm1Ptr, BfSizeArray **perm2Ptr) {
  BF_ERROR_BEGIN();

  /* Extract the two nodal domains determined by `phiFiedler`. */

  BfTrimesh *submesh1 = bfTrimeshGetLevelSetSubmesh(trimesh, phiFiedler, tol, permMask, perm1Ptr);
  HANDLE_ERROR();

  BfRealArray *phiFiedlerNeg = bfRealArrayCopy(phiFiedler);
  HANDLE_ERROR();

  bfRealArrayNegate(phiFiedlerNeg);

  BfTrimesh *submesh2 = bfTrimeshGetLevelSetSubmesh(trimesh, phiFiedlerNeg, tol, permMask, perm2Ptr);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  bfRealArrayDeinitAndDealloc(&phiFiedlerNeg);

  *submesh1Ptr = submesh1;
  *submesh2Ptr = submesh2;
}

static void initRecursive(BfFiedlerTreeNode *node,
                          BfFiedlerTree const *tree, BfTrimesh const *trimesh,
                          BfReal tol, bool keepNodeTrimeshes,
                          BfSize i0, BfSize i1, BfSize *perm, BfSize depth) {
  BF_ERROR_BEGIN();

  BF_ASSERT(i0 < i1);

  bfLogDebug("initRecursive(depth = %lu, i0 = %lu, i1 = %lu)\n", depth, i0, i1);

  BfRealArray *phiFiedler = bfTrimeshGetFiedler(trimesh);
  HANDLE_ERROR();

  BfTrimesh *trimeshFixed = NULL;
  BfRealArray *phiFiedlerFixed = NULL;
  doBoundaryFix(trimesh, phiFiedler, &trimeshFixed, &phiFiedlerFixed);
  HANDLE_ERROR();

  /* An assumption we make is that the vertices corresponding to [i0,
   * i1] are the first vertices in `trimeshFixed`. This is a little
   * hacky, but it does make the implementation simpler. We verify
   * this now. */
#if BF_DEBUG
  for (BfSize i = i0; i < i1; ++i)
    BF_ASSERT(bfPoints3Contains(tree->trimesh->verts, trimesh->verts->data[i - i0]));
#endif

  bool *permMask = bfMemAllocAndZero(trimeshFixed->verts->size, sizeof(bool));
  HANDLE_ERROR();

  for (BfSize i = 0; i < i1 - i0; ++i)
    permMask[i] = true;

  /* Split `trimesh` into two pieces by computing the first nonzero
   * eigenfunction of the Dirichlet Laplace-Beltrami eigenvalue
   * problem. This is done using a piecewise linear FEM approximation,
   * so the zero level set separating the two submeshes is also a
   * piecewise linear approximation and doesn't necessarily align with
   * any mesh edges. Afterwards, `submesh` contains the new meshes,
   * which will have new vertices added to them after cutting along
   * the zero level set, and `subperm[i]` gives the position in
   * `trimesh` of the vertices from `trimesh` which were included in
   * `submesh[i]`. */
  BfTrimesh *submesh[2] = {NULL, NULL};
  BfSizeArray *subperm[2] = {NULL, NULL};
  splitTrimesh(trimeshFixed, phiFiedlerFixed, tol, permMask, &submesh[0], &submesh[1], &subperm[0], &subperm[1]);
  HANDLE_ERROR();

  /* Now that we've split the mesh, we can free the fixed triangle
   * mesh and Fiedler vector to save on memory. */
  bfTrimeshDeinitAndDealloc(&trimeshFixed);
  bfRealArrayDeinitAndDealloc(&phiFiedlerFixed);

  /* Update `perm` and compute the offsets for each child */

  BfSize *newPerm = bfMemAlloc(i1 - i0, sizeof(BfSize));
#if BF_DEBUG
  for (BfSize k = 0; k < i1 - i0; ++k) newPerm[k] = BF_SIZE_BAD_VALUE;
#endif

  BfSize i = i0;
  for (BfSize j = 0; j < MAX_NUM_CHILDREN; ++j) {
    node->super.offset[j] = i;
    for (BfSize k = 0; k < bfSizeArrayGetSize(subperm[j]); ++k) {
      BfSize l = i++ - i0;
      BF_ASSERT(l < i1 - i0);
      newPerm[l] = perm[i0 + bfSizeArrayGet(subperm[j], k)];
    }
  }
  BF_ASSERT(i == i1);
  node->super.offset[MAX_NUM_CHILDREN] = i;
  bfMemCopy(newPerm, i1 - i0, sizeof(BfSize), perm + i0);

  /* Free these now that we've permuted the child indices. */
  for (BfSize j = 0; j < 2; ++j)
    bfSizeArrayDeinitAndDealloc(&subperm[j]);

#if BF_DEBUG
  for (BfSize i = i0; i < i1; ++i) BF_ASSERT(perm[i] != BF_SIZE_BAD_VALUE);
  BF_ASSERT(bfSizeIsPerm(tree->super.perm->size, tree->super.perm->index));
#endif

  bfMemFree(newPerm);

  /* Recursively create children */

  for (BfSize j = 0; j < MAX_NUM_CHILDREN; ++j) {
    BfSize i0Child = node->super.offset[j];
    BfSize i1Child = node->super.offset[j + 1];

    /* Don't create a child with no contained points */
    BfSize numChildPoints = i1Child - i0Child;
    if (numChildPoints == 0)
      continue;

    BfFiedlerTreeNode *child = bfFiedlerTreeNodeNew();
    HANDLE_ERROR();

    BfSize childDepth = depth + 1;

    bfTreeNodeInit(&child->super, &TreeNodeVtable, false, (void *)node,
                   MAX_NUM_CHILDREN, j, childDepth);
    HANDLE_ERROR();

    child->trimesh = keepNodeTrimeshes ? submesh[j] : NULL;

    if (numChildPoints > LEAF_SIZE_THRESHOLD) {
      initRecursive(child, tree, submesh[j], tol, keepNodeTrimeshes, i0Child, i1Child, perm, childDepth);
      HANDLE_ERROR();
    }

    node->super.child[j] = bfFiedlerTreeNodeToTreeNode(child);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  bfMemFree(permMask);

  bfRealArrayDeinitAndDealloc(&phiFiedler);
}

void bfFiedlerTreeNodeInitRoot(BfFiedlerTreeNode *node, BfFiedlerTree const *tree, BfReal tol, bool keepNodeTrimeshes) {
  BF_ERROR_BEGIN();

  bfTreeNodeInitRoot(&node->super, &TreeNodeVtable, bfFiedlerTreeConstToTreeConst(tree), MAX_NUM_CHILDREN);
  HANDLE_ERROR();

  if (keepNodeTrimeshes) {
    node->trimesh = bfTrimeshCopy(tree->trimesh);
    HANDLE_ERROR();
  } else {
    node->trimesh = NULL;
  }

  initRecursive(
    node,
    tree,
    tree->trimesh,
    tol,
    keepNodeTrimeshes,
    /* i0: */ 0,
    /* i1: */ bfTrimeshGetNumVerts(tree->trimesh),
    tree->super.perm->index,
    /* currentDepth: */ 0);

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfFiedlerTreeNodeDeinit(BfFiedlerTreeNode *node) {
  if (node->trimesh != NULL)
    bfTrimeshDeinitAndDealloc((BfTrimesh **)&node->trimesh);

  bfTreeNodeDeinit(&node->super);
}

void bfFiedlerTreeNodeDealloc(BfFiedlerTreeNode **node) {
  bfMemFree(*node);
  *node = NULL;
}

void bfFiedlerTreeNodeDeinitAndDealloc(BfFiedlerTreeNode **node) {
  bfFiedlerTreeNodeDeinit(*node);
  bfFiedlerTreeNodeDealloc(node);
}

BfSizeArray *bfTrimeshGetInteriorInds(BfTrimesh const *trimesh) {
  BF_ERROR_BEGIN();

  BfSizeArray *interiorInds = bfSizeArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  for (BfSize i = 0; i < bfTrimeshGetNumVerts(trimesh); ++i) {
    if (!trimesh->isBoundaryVert[i]) {
      bfSizeArrayAppend(interiorInds, i);
      HANDLE_ERROR();
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return interiorInds;
}
