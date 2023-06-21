#include <bf/octree.h>

#include <bf/array.h>
#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/octree_node.h>
#include <bf/points.h>
#include <bf/size_array.h>

typedef struct {
  BfOctree const *octree;
  BfReal const *queryPoint;
  BfArray *pqueue;
} Context;

typedef struct {
  enum {POINT, NODE} type;
  BfReal dist;
  union {
    BfSize pointIndex;
    BfTreeNode const *node;
  };
} PqueueElt;

static bool pqueueEltsAreEqual(PqueueElt const *elt1, PqueueElt const *elt2) {
  if (elt1->type != elt2->type)
    return false;
  else if (elt1->type == POINT)
    return elt1->pointIndex == elt2->pointIndex;
  else
    return elt1->node == elt2->node;
}

static int pqueueEltCmp(PqueueElt const *elt1, PqueueElt const *elt2, void *arg) {
  (void)arg;
  if (elt1->dist < elt2->dist)
    return -1;
  else if (elt1->dist == elt2->dist)
    return 0;
  else
    return 1;
}

static void insertEltIntoPqueue(BfArray *pqueue, PqueueElt const *elt) {
  BF_ERROR_BEGIN();

  BfSize i = bfArrayFindSorted(pqueue, elt, (BfCompar)pqueueEltCmp);

#if BF_DEBUG
  if (i < bfArrayGetSize(pqueue)) {
    PqueueElt const *otherElt = bfArrayGetPtr(pqueue, i);
    BF_ASSERT(!pqueueEltsAreEqual(elt, otherElt));
  }
#endif

  bfArrayInsert(pqueue, i, elt);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void addPoint(Context *c, BfSize pointIndex) {
  BF_ERROR_BEGIN();

  BfReal const *point = bfPoints3GetPtrConst(c->octree->points, pointIndex);
  // TODO: BF_ERROR_OK(); -- add something like this... when we're
  // positive there should be no new errors...

  PqueueElt elt = {
    .type = POINT,
    .dist = bfPoint3Dist(c->queryPoint, point),
    .pointIndex = pointIndex,
  };

  insertEltIntoPqueue(c->pqueue, &elt);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void addInternalNode(Context *c, BfTreeNode const *node) {
  BF_ERROR_BEGIN();

  BF_ASSERT(!bfTreeNodeIsLeaf(node));

  BfOctreeNode const *octreeNode = bfTreeNodeConstToOctreeNodeConst(node);

  PqueueElt elt = {
    .type = NODE,
    .dist = bfBoundingBox3GetDistanceToPoint(&octreeNode->boundingBox, c->queryPoint),
    .node = node,
  };

  insertEltIntoPqueue(c->pqueue, &elt);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void addLeafNode(Context *c, BfTreeNode const *node) {
  BF_ERROR_BEGIN();

  BF_ASSERT(bfTreeNodeIsLeaf(node));

  BfSize i0 = bfTreeNodeGetFirstIndex(node);
  BfSize i1 = bfTreeNodeGetLastIndex(node);

  for (BfSize i = i0; i < i1; ++i) {
    addPoint(c, c->octree->super.perm.index[i]);
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

static void addNode(Context *c, BfTreeNode const *node) {
  if (bfTreeNodeIsLeaf(node))
    addLeafNode(c, node);
  else
    addInternalNode(c, node);
}

BfSizeArray *bfOctreeGetNearestNeighbors(BfOctree const *octree, BfPoint3 const queryPoint, BfSize numNbs) {
  BF_ERROR_BEGIN();

  /* Array containing the indices to the `numNbs` nearest neighbors of
   * `point` in `octree->points`. */
  BfSizeArray *nbInds = bfSizeArrayNewWithCapacity(numNbs);
  HANDLE_ERROR();

  /* A priority queue used to enumerate the internal nodes and points
   * contained in `octree` in order of increasing distance from
   * `point`. */
  BfArray *pqueue = bfArrayNewEmpty(sizeof(PqueueElt));
  HANDLE_ERROR();

  /* Set up context struct: */
  Context c = {.octree = octree, .queryPoint = queryPoint, .pqueue = pqueue};

  /* Initialize `pqueue` with the root node of the octree: */
  addNode(&c, octree->super.root);

  while (!bfArrayIsEmpty(pqueue) && bfSizeArrayGetSize(nbInds) < numNbs) {
    /* Pop the first element from the priority queue. */
    PqueueElt elt;
    bfArrayGet(pqueue, 0, &elt);
    bfArrayRemove(pqueue, 0);

    /* If the element is a point, we've found another nearest
     * neighbor---add it and continue. */
    if (elt.type == POINT) {
      bfSizeArrayAppend(nbInds, elt.pointIndex);
      continue;
    }

    /* If the element is a node, add it to the priority queue (this
     * means adding the node itself if it's an internal node, or all
     * of the contained points if it's a leaf node---handled by
     * `pqueueAddNode`). */
    for (BfSize i = 0; i < elt.node->maxNumChildren; ++i)
      if (elt.node->child[i] != NULL)
        addNode(&c, elt.node->child[i]);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  bfArrayDeinitAndDealloc(&pqueue);

  return nbInds;
}
