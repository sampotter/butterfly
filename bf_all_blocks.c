#include <stdio.h>

#include "error_macros.h"
#include "fac.h"
#include "quadtree.h"

FILE *fp = NULL;

BfPtrArray getChildrenAsPtrArray(BfQuadtreeNode const *node) {
  BfPtrArray childNodes;
  bfInitPtrArray(&childNodes, 4);
  for (BfSize i = 0; i < 4; ++i)
    if (node->child[i])
      bfPtrArrayAppend(&childNodes, node->child[i]);
  return childNodes;
}

static void makeMlFacRec(BfQuadtree const *tree, BfReal K,
                         BfPtrArray const *srcNodes,
                         BfPtrArray const *tgtNodes,
                         BfSize level) {
  BfSize maxLevel = 3;
  if (level > maxLevel) {
    printf("level = %lu, exceeded maximum level (%lu)\n", level, maxLevel);
    return;
  }

  BEGIN_ERROR_HANDLING();

  BfMatProduct *factorization = NULL;

  BfPtrArray srcChildNodes;
  BfPtrArray tgtChildNodes;

  for (BfSize i = 0; i < bfPtrArraySize(tgtNodes); ++i) {
    BfQuadtreeNode *tgtNode = bfPtrArrayGet(tgtNodes, i);

    for (BfSize j = 0; j < bfPtrArraySize(srcNodes); ++j) {
      BfQuadtreeNode *srcNode = bfPtrArrayGet(srcNodes, j);

      bool separated = bfQuadtreeNodesAreSeparated(srcNode, tgtNode);

      printf("level = %lu, i = %lu, j = %lu", level, i, j);

      if (separated) {
        BfSize numRows = tgtNode->offset[4] - tgtNode->offset[0];
        BfSize numCols = srcNode->offset[4] - srcNode->offset[0];

        printf(", separated: making %lux%lu BF...", numRows, numCols);

        fprintf(fp, "%lu %lu %lu %lu %lu\n",
                tgtNode->offset[0], tgtNode->offset[4],
                srcNode->offset[0], srcNode->offset[4],
                level);

        BfQuadtreeLevelIter srcLevelIter, tgtLevelIter;
        BfSize numFactors = bfFacHelm2Prepare(tree, srcNode, tgtNode, K,
                                              &srcLevelIter, &tgtLevelIter);
        if (numFactors == 0) {
          printf(" FAILED\n");
        } else {
          factorization = bfFacHelm2Make(
            tree, srcNode, tgtNode, K, &srcLevelIter, &tgtLevelIter, numFactors);
          printf(" success\n");
        }

        HANDLE_ERROR();
      } else {
        printf(", not separated: make BFs on next level \n");

        srcChildNodes = getChildrenAsPtrArray(srcNode);
        HANDLE_ERROR();

        tgtChildNodes = getChildrenAsPtrArray(tgtNode);
        HANDLE_ERROR();

        makeMlFacRec(tree, K, &srcChildNodes, &tgtChildNodes, level + 1);
        HANDLE_ERROR();

        bfFreePtrArray(&srcChildNodes);
        bfFreePtrArray(&tgtChildNodes);
      }
    }
  }

  END_ERROR_HANDLING() {
    bfFreePtrArray(&srcChildNodes);
    bfFreePtrArray(&tgtChildNodes);
    bfMatProductDeinitAndDelete(&factorization);
  }
}

void bfMakeMultilevelFactorization(BfQuadtree const *tree, BfReal K) {
  BEGIN_ERROR_HANDLING();

  BfQuadtreeLevelIter levelIter = bfInitQuadtreeLevelIter(
    BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, (BfQuadtreeNode *)tree->root);
  HANDLE_ERROR();

  /* skip to level 2 */
  bfQuadtreeLevelIterNext(&levelIter);
  bfQuadtreeLevelIterNext(&levelIter);

  fp = fopen("blocks.txt", "w");

  BfPtrArray const *levelNodes = &levelIter.levelNodes;

  makeMlFacRec(tree, K, levelNodes, levelNodes, 2);

  fclose(fp);

  END_ERROR_HANDLING() {}

  bfFreeQuadtreeLevelIter(&levelIter);
}

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    printf("usage: %s <points.bin>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  BEGIN_ERROR_HANDLING();

  BfPoints2 points;
  bfReadPoints2FromFile(argv[1], &points);
  HANDLE_ERROR();
  printf("read points from %s\n", argv[1]);

  BfQuadtree tree;
  bfInitQuadtreeFromPoints(&tree, &points);
  HANDLE_ERROR();
  puts("built quadtree");

  BfReal K = 3000;

  bfMakeMultilevelFactorization(&tree, K);

  END_ERROR_HANDLING() {}

  bfFreeQuadtree(&tree);
  bfFreePoints2(&points);
}
