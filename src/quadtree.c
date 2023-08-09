#include <bf/quadtree.h>

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <bf/circle.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/points.h>
#include <bf/ptr_array.h>
#include <bf/quadtree_node.h>
#include <bf/tree_level_iter.h>
#include <bf/vectors.h>

static BfSize const NUM_CHILDREN = 4;

/** Interface(Tree, Quadtree) */

static BfTreeVtable TreeVtable = {
  .CopyInto = (__typeof__(&bfTreeCopyInto))bfQuadtreeCopyInto,
  .GetType = (__typeof__(&bfTreeGetType))bfQuadtreeGetType,
};

void bfQuadtreeCopyInto(BfQuadtree *quadtree, BfQuadtree *dstQuadtree) {
  *dstQuadtree = *quadtree;
}

BfType bfQuadtreeGetType(BfQuadtree const *tree) {
  (void)tree;
  return BF_TYPE_QUADTREE;
}

/** Upcasting: Quadtree -> Tree */

BfTree *bfQuadtreeToTree(BfQuadtree *quadtree) {
  return &quadtree->super;
}

BfTree const *bfQuadtreeConstToTreeConst(BfQuadtree const *quadtree) {
  return &quadtree->super;
}

/** Downcasting: Tree -> Quadtree */

BfQuadtree *bfTreeToQuadtree(BfTree *tree) {
  if (!bfTreeInstanceOf(tree, BF_TYPE_QUADTREE)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfQuadtree *)tree;
  }
}

BfQuadtree const *bfTreeConstToQuadtreeConst(BfTree const *tree) {
  if (!bfTreeInstanceOf(tree, BF_TYPE_QUADTREE)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfQuadtree const *)tree;
  }
}

/** Implementation: Quadtree */

BfQuadtree *bfQuadtreeNew() {
  BF_ERROR_BEGIN();

  BfQuadtree *quadtree = bfMemAlloc(1, sizeof(BfQuadtree));
  if (quadtree == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END()
    quadtree = NULL;

  return quadtree;
}

void bfQuadtreeInit(BfQuadtree *tree, BfPoints2 const *points,
                    BfVectors2 const *unitNormals) {
  BF_ERROR_BEGIN()

  BfQuadtreeNode *root = bfQuadtreeNodeNew();
  HANDLE_ERROR();

  bfTreeInit(&tree->super, &TreeVtable, bfQuadtreeNodeToTreeNode(root), points->size);
  HANDLE_ERROR();

  tree->points = points;
  tree->unitNormals = unitNormals;

  bfQuadtreeNodeInitRoot(root, tree);
  HANDLE_ERROR();

  BF_ERROR_END()
    bfQuadtreeDeinit(tree);
}

void bfQuadtreeDeinit(BfQuadtree *tree) {
  bfTreeDeinit(&tree->super);
}

void bfQuadtreeDealloc(BfQuadtree **quadtree) {
  bfMemFree(*quadtree);
  *quadtree = NULL;
}

void bfQuadtreeDeinitAndDealloc(BfQuadtree **quadtree) {
  bfQuadtreeDeinit(*quadtree);
  bfQuadtreeDealloc(quadtree);
}

static void saveBoxesToTextFileRec(BfQuadtreeNode const *node, FILE *fp) {
  BfReal xmin = node->bbox.min[0];
  BfReal xmax = node->bbox.max[0];
  BfReal ymin = node->bbox.min[1];
  BfReal ymax = node->bbox.max[1];

  fprintf(fp, "%g %g %g %g\n", xmin, xmax, ymin, ymax);

  for (BfSize i = 0; i < NUM_CHILDREN; ++i) {
    if (node->super.child[i] == NULL)
      continue;
    saveBoxesToTextFileRec(
      bfTreeNodeConstToQuadtreeNodeConst(node->super.child[i]), fp);
  }
}

void bfQuadtreeSaveBoxesToTextFile(BfQuadtree const *quadtree, char const *path) {
  FILE *fp = fopen(path, "w");

  BfTree const *tree = bfQuadtreeConstToTreeConst(quadtree);

  BfQuadtreeNode const *root =
    bfTreeNodeConstToQuadtreeNodeConst(bfTreeGetRootNodeConst(tree));

  saveBoxesToTextFileRec(root, fp);
  fclose(fp);
}
