#include <bf/assert.h>
#include <bf/error.h>
#include <bf/fiedler_tree.h>
#include <bf/fiedler_tree_node.h>
#include <bf/util.h>

#include <stdlib.h>

#define DEBUG_OUTPUT 0

void checkPerm(BfTree const *tree, BfTreeNode const *treeNode, void *arg);

int main(int argc, char *argv[]) {
  if (argc < 2 || argc > 3) {
    printf("usage: %s <trimesh.obj> [tol]\n", argv[0]);
    return 1;
  }

  char const *objPath = argv[1];

  BfReal tol = argc >= 3 ? atof(argv[2]) : 1e-15;

  bfToc();

  BfTrimesh *trimesh = bfTrimeshNewFromObjFile(objPath);

  printf("loaded triangle mesh (%lu verts and %lu faces) [%.1fs]\n",
         bfTrimeshGetNumVerts(trimesh), bfTrimeshGetNumFaces(trimesh), bfToc());

  BfFiedlerTree *fiedlerTree = bfFiedlerTreeNewFromTrimesh(trimesh, tol, /* keepNodeTrimeshes: */ true);

  printf("building Fiedler tree (tol = %g) [%.1fs]\n", tol, bfToc());

  BfTree *tree = bfFiedlerTreeToTree(fiedlerTree);

  bfTreeMapConst(tree, NULL, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, (BfTreeMapConstFunc)checkPerm, &tol);

  bfFiedlerTreeDeinitAndDealloc(&fiedlerTree);
  bfTrimeshDeinitAndDealloc(&trimesh);
}

void checkPerm(BfTree const *tree, BfTreeNode const *treeNode, void *arg) {
  BfReal tol = *(BfReal *)arg;

  BfFiedlerTree const *fiedlerTree = bfTreeConstToFiedlerTreeConst(tree);
  BfFiedlerTreeNode const *fiedlerTreeNode = bfTreeNodeConstToFiedlerTreeNodeConst(treeNode);

  BfSize i0 = bfTreeNodeGetFirstIndex(treeNode);
  BfSize i1 = bfTreeNodeGetLastIndex(treeNode);

#if DEBUG_OUTPUT
  BfSize depth = bfTreeNodeGetDepth(treeNode);
  printf("%lu -> [%lu, %lu)\n", depth, i0, i1);
#endif

  for (BfSize i = i0; i < i1; ++i) {
    BfReal const *v = bfPoints3GetPtrConst(
      fiedlerTree->trimesh->verts, tree->perm.index[i]);

    BF_ASSERT(bfPoints3ContainsApprox(fiedlerTreeNode->trimesh->verts, v, tol));
  }
}
