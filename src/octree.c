#include <bf/octree.h>

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <bf/circle.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/points.h>
#include <bf/ptr_array.h>
#include <bf/octree_node.h>
#include <bf/tree_level_iter.h>
#include <bf/vectors.h>

/** Interface(Tree, Octree) */

static BfTreeVtable TreeVtable = {
  .GetType = (__typeof__(&bfTreeGetType))bfOctreeGetType,
};

BfType bfOctreeGetType(BfOctree const *tree) {
  (void)tree;
  return BF_TYPE_OCTREE;
}

/** Upcasting: Octree -> Tree */

BfTree *bfOctreeToTree(BfOctree *octree) {
  return &octree->super;
}

BfTree const *bfOctreeConstToTreeConst(BfOctree const *octree) {
  return &octree->super;
}

/** Downcasting: Tree -> Octree */

BfOctree *bfTreeToOctree(BfTree *tree) {
  if (!bfTreeInstanceOf(tree, BF_TYPE_OCTREE)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfOctree *)tree;
  }
}

/** Implementation: Octree */

BfOctree *bfOctreeNew() {
  BF_ERROR_BEGIN();

  BfOctree *octree = bfMemAlloc(1, sizeof(BfOctree));
  if (octree == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    octree = NULL;

  return octree;
}

void bfOctreeInit(BfOctree *tree, BfPoints3 const *points,
                  BfVectors3 const *unitNormals) {
  BF_ERROR_BEGIN()

  BfOctreeNode *root = bfOctreeNodeNew();
  HANDLE_ERROR();

  bfTreeInit(&tree->super, &TreeVtable, bfOctreeNodeToTreeNode(root), points->size);
  HANDLE_ERROR();

  tree->points = points;
  tree->unitNormals = unitNormals;

  bfOctreeNodeInitRoot(root, tree);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfOctreeDeinit(tree);
}

void bfOctreeDeinit(BfOctree *tree) {
  bfTreeDeinit(&tree->super);
}

typedef struct {
  FILE *fp;
} WriteBoxesWkspc;

static void writeBoxes(BfTree const *tree, BfTreeNode const *treeNode, WriteBoxesWkspc *wkspc) {
  (void)tree;

  FILE *fp = wkspc->fp;

  BfSize depth = bfTreeNodeGetDepth(treeNode);

  BfOctreeNode const *octreeNode = bfTreeNodeConstToOctreeNodeConst(treeNode);

  BfBoundingBox3 const *boundingBox = &octreeNode->boundingBox;

  BfReal xmin = boundingBox->min[0];
  BfReal xmax = boundingBox->max[0];

  BfReal ymin = boundingBox->min[1];
  BfReal ymax = boundingBox->max[1];

  BfReal zmin = boundingBox->min[2];
  BfReal zmax = boundingBox->max[2];

  fprintf(fp, "%lu %g %g %g %g %g %g\n", depth, xmin, xmax, ymin, ymax, zmin, zmax);
}

void bfOctreeSaveBoxesToTextFile(BfOctree const *octree, char const *path) {
  BF_ERROR_BEGIN();

  BfTree const *tree = bfOctreeConstToTreeConst(octree);

  FILE *fp = fopen(path, "w");
  if (fp == NULL)
    RAISE_ERROR(BF_ERROR_FILE_ERROR);

  WriteBoxesWkspc wkspc = {.fp = fp};

  bfTreeMapConst(tree, NULL, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, (BfTreeMapConstFunc)writeBoxes, &wkspc);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}

  fclose(fp);
}
