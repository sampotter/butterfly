#pragma once

#include "node_array.h"
#include "tree_node.h"

typedef struct BfFacSpec {
  BfTree *rowTree;
  BfTree *colTree;
  BfSize rowTreeInitDepth;
  BfSize colTreeInitDepth;

  /* The tolerance used to compute truncated SVDs of blocks when
   * streaming the butterfly factorization. */
  BfReal tol;

  /* The minimum number of rows in a block needed before a truncated
   * SVD is attempted. This is used to bottom out when finding
   * epsilon-rank cuts in the row tree. */
  BfSize minNumRows;

  /* The minimum number of columns in a block needed before trying a
   * truncated SVD. This is used to prevent us from starting to low in
   * the column tree. */
  BfSize minNumCols;
} BfFacSpec;

struct BfFac {
  BfTreeNode const *colNode;

  BfConstNodeArray rowNodes;

  BfMat *Psi;

  BfSize numW;
  BfMat **W;
};

void bfFacDelete(BfFac **fac);
BfMat *bfFacGetMat(BfFac const *fac);
BfMatProduct *bfFacGetMatProduct(BfFac const *fac);

BfFac *makeLeafNodePartialFac(BfTreeNode const *colNode, BfPtrArray *PsiBlocks, BfPtrArray *WBlocks);
// void setPartialFacWBlock(PartialFac *partialFac, BfSize k, BfMat *W);
BfSize partialFacGetNumRows(BfFac const *fac);
BfVec *partialFacMulVec(BfFac const *fac, BfVec const *x);
// bool partialFacHasContiguousRowSpan(PartialFac const *partialFac);
// BfPtrArray *getIndexedPsiSubblocksInRowRange(PartialFac const *partialFac, BfSize i0Sel, BfSize i1Sel);
// void getPsiAndW0BlocksByRowNodeForPartialFac(PartialFac const *partialFac, BfTreeNode const *rowNode, BfMat **PsiBlockPtr, BfMat **W0BlockPtr);

bool getPsiAndW(BfFacSpec const *facSpec, BfMat const *mat, BfTreeNode const *rowNode, BfMat **PsiPtr, BfMat **WPtr);

// TODO: all these functions operating on a PtrArray of Facs really
// represent a new type. ... The type is "horizontally concatenated
// BFs". OK, not that sexy. "BfFacRow"? "BfFacs"?

// bool partialFacsHaveContinguousRowSpans(BfPtrArray const *partialFacs);
// bool partialFacsHaveSameRowSpan(BfPtrArray const *partialFacs);
// BfConstPtrArray getFirstRowNodes(BfPtrArray const *partialFacs);
// BfConstPtrArray getLastRowNodes(BfPtrArray const *partialFacs);
// BfConstPtrArray getRowNodesByFirstIndex(BfPtrArray const *partialFacs, BfSize i0);
// BfConstPtrArray getMergeCut(BfPtrArray const *partialFacs);
// void getPsiAndWBlocksByRowNode(BfPtrArray const *currentPartialFacs, BfTreeNode const *rowNode, BfMat **PsiPtr, BfMat **WPtr);
// bool getLowRankApproximation(BfFacStreamer const *facStreamer, BfMat const *PsiStarSubblock, BfMat **PsiSubblockPtr, BfMat **W0SubblockPtr);
// void findEpsilonRankCutAndGetNewBlocks(BfFacStreamer const *facStreamer, BfTreeNode const *rootRowNode, BfMat const *PsiStarBlock, BfConstPtrArray **epsRankCutPtr, BfMat **PsiBlockPtr, BfMat **W0BlockPtr);
BfFac *mergeAndSplit(BfPtrArray const *facs, BfFacSpec const *facSpec);
