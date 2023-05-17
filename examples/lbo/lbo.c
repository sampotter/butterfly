#include "lbo.h"

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/interval_tree.h>
#include <bf/interval_tree_node.h>
#include <bf/linalg.h>
#include <bf/vec_real.h>

#include <math.h>

BfPoints1 *convertEigsToFreqs(BfVecReal const *Lam) {
  BF_ERROR_BEGIN();

  BfPoints1 *freqs = bfPoints1New();
  HANDLE_ERROR();

  bfPoints1InitEmpty(freqs, BF_ARRAY_DEFAULT_CAPACITY);
  HANDLE_ERROR();

  for (BfSize i = 0; i < Lam->super.size; ++i) {
    /* Clamp eigenvalues to be nonnegative if we find any that are
     * negative due to roundoff */
    if (Lam->data[i] < 0)
      Lam->data[i] = 0;

    bfPoints1Append(freqs, sqrt(Lam->data[i]));
    HANDLE_ERROR();
  }

  BF_ERROR_END() {
    BF_ASSERT(false);
  }

  return freqs;
}

void feedFacStreamerNextEigenband(BfFacStreamer *facStreamer, BfPoints1 *freqs,
                                  BfMat const *L, BfMat const *M) {
  BF_ERROR_BEGIN();

  BfMat *Phi = NULL;
  BfVecReal *Lam = NULL;
  BfPoints1 *newFreqs = NULL;

  BfTreeNode *treeNode = bfFacStreamerGetCurrentColumnNode(facStreamer);
  if (!bfTreeNodeIsLeaf(treeNode))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfTree *tree = bfTreeNodeGetTree(treeNode);

  BfIntervalTree *intervalTree = bfTreeToIntervalTree(tree);
  HANDLE_ERROR();

  /* The current node should be empty since we're incrementally
   * constructing the column tree as we stream the band
   * eigenvalues. */
  if (!bfTreeNodeIsEmpty(treeNode))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfIntervalTreeNode *intervalTreeNode = bfTreeNodeToIntervalTreeNode(treeNode);
  HANDLE_ERROR();

  /* The column tree node shouldn't simultaneously be the leftmost and
   * rightmost node (this only happens if it's the root...) */
  if (intervalTreeNode->isLeftmost && intervalTreeNode->isRightmost)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* Set up bracket and shift for computing next eigenband */
  BfReal lam0 = pow(intervalTreeNode->a, 2);
  BfReal lam1 = pow(intervalTreeNode->b, 2);
  BfReal sigma = (lam0 + lam1)/2;
  if (intervalTreeNode->isLeftmost) {
    sigma = lam0;
    lam0 = NAN;
  } else if (intervalTreeNode->isRightmost) {
    sigma = lam1;
    lam1 = NAN;
  }

  /* Compute the next eigenband using Lanczos */
  bfGetEigenband(L, M, lam0, lam1, sigma, &Phi, &Lam);
  HANDLE_ERROR();

  /* Permute the rows of Phi, putting them into row tree order */
  bfMatPermuteRows(Phi, bfFacStreamerGetRowTreeReversePerm(facStreamer));

  /* Convert the new eigenvalues to frequencies: */
  newFreqs = convertEigsToFreqs(Lam);
  HANDLE_ERROR();

  /* Sort the new frequencies into the existing set of frequencies
   *
   * TODO: because of the way we're doing this, we can guarantee that
   * we only need to append the new frequencies to the end of `freqs`,
   * simplifying things a bit. */
  bfPoints1InsertPointsSorted(freqs, newFreqs);
  HANDLE_ERROR();

  /* We set update the tree's point set without rebuilding the
   * tree. This has the effect of adjusting the range of points each
   * tree node points to without actually changing with nodes are in
   * the tree. */
  bfIntervalTreeSetPoints(intervalTree, freqs, /* rebuildTree: */ false);
  HANDLE_ERROR();

  /* The current column tree node should now have exactly the same
   * number of points as the number of newly streamed eigenvalues. */
  if (bfTreeNodeGetNumPoints(treeNode) != Lam->super.size)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* Feed the factorization the streamed band of eigenvectors */
  bfFacStreamerFeed(facStreamer, Phi);
  HANDLE_ERROR();

  BF_ERROR_END() {}

  bfPoints1Deinit(newFreqs);
  bfMatDelete(&Phi);
  bfVecRealDeinitAndDealloc(&Lam);
}
