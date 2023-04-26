#include <math.h>
#include <stdlib.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/fac_streamer.h>
// #include <bf/fiedler_tree.h>
#include <bf/interval_tree.h>
#include <bf/interval_tree_node.h>
#include <bf/lbo.h>
#include <bf/linalg.h>
#include <bf/logging.h>
#include <bf/mat.h>
#include <bf/octree.h>
#include <bf/rand.h>
#include <bf/vec_real.h>

BfPoints1 *convertEigsToFreqs(BfVecReal const *Lam) {
  BEGIN_ERROR_HANDLING();

  bfSetLogLevel(BF_LOG_LEVEL_INFO);

  BfPoints1 *freqs = bfPoints1New();
  HANDLE_ERROR();

  bfPoints1InitEmpty(freqs, BF_ARRAY_DEFAULT_CAPACITY);
  HANDLE_ERROR();

  for (BfSize i = 0; i < Lam->super.size; ++i) {
    /* Clamp eigenvalues to be nonnegative if we find any that are
     * negative due to roundoff */
    if (Lam->data[i] < 0) {
      BF_ASSERT(Lam->data[i] > -1e-13);
      Lam->data[i] = 0;
    }

    bfPoints1Append(freqs, sqrt(Lam->data[i]));
    HANDLE_ERROR();
  }

  END_ERROR_HANDLING() {
    BF_ASSERT(false);
  }

  return freqs;
}

void feedFacStreamerNextEigenband(BfFacStreamer *facStreamer, BfPoints1 *freqs,
                                  BfMat const *L, BfMat const *M) {
  BEGIN_ERROR_HANDLING();

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

  END_ERROR_HANDLING() {}

  bfPoints1Deinit(newFreqs);
  bfMatDelete(&Phi);
  bfVecRealDeinitAndDealloc(&Lam);
}

int main(int argc, char const *argv[]) {
  BEGIN_ERROR_HANDLING() {}

  if (argc != 6) {
    printf("usage: %s objPath tol p rowTreeOffset freqTreeDepth\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  bfSeed(0);

  char const *objPath = argv[1];
  BfReal tol = strtod(argv[2], NULL);
  BfReal p = strtod(argv[3], NULL);
  BfSize rowTreeOffset = strtoull(argv[4], NULL, 10);
  BfSize freqTreeDepth = strtoull(argv[5], NULL, 10);

  /* Load triangle mesh from binary files */
  BfTrimesh trimesh;
  bfTrimeshInitFromObjFile(&trimesh, objPath);
  HANDLE_ERROR();

  BfSize numVerts = bfTrimeshGetNumVerts(&trimesh);
  BfSize numEigs = (BfSize)(p*numVerts);

  printf("- mesh with %lu vertices\n", numVerts);
  printf("- streaming %lu eigenpairs\n", numEigs);

  // TODO: implement
//   /* Build fiedler tree */
//   BfFiedlerTree fiedlerTree;
//   bfFiedlerTreeInitFromTrimesh(&fiedlerTree, trimesh, &L, &M);
//   BfTree *rowTree = bfFiedlerTreeToTree(&fiedlerTree);

  BfOctree octree;
  bfOctreeInit(&octree, &trimesh.verts, NULL);
  HANDLE_ERROR();

  char const *octreeBoxesPath = "octree_boxes.txt";
  bfOctreeSaveBoxesToTextFile(&octree, octreeBoxesPath);
  HANDLE_ERROR();
  printf("- wrote octree cells to %s\n", octreeBoxesPath);

  /* Upcast the spatial tree to get the row tree */
  BfTree *rowTree = bfOctreeToTree(&octree);

  BfSize rowTreeMaxDepth = bfTreeGetMaxDepth(rowTree);
  printf("- row tree has depth %lu\n", rowTreeMaxDepth);

  /* Compute a finite element discretization of the Laplace-Beltrami
   * operator on `trimesh` using linear finite elements. The stiffness
   * matrix is returned in L and the mass matrix is returned in M. The
   * mass matrix isn't diagonal but this isn't too important. */
  BfMat *L, *M;
  bfLboGetFemDiscretization(&trimesh, &L, &M);
  HANDLE_ERROR();

  /* Find the largest eigenvalue. We need this to determine the
   * interval on which we'll build the frequency tree. */
  BfReal lamMax = bfGetMaxEigenvalue(L, M);
  HANDLE_ERROR();
  printf("- lambda_max = %g\n", lamMax);

//   bfMatCsrRealDump(bfMatToMatCsrReal(L), "L_rowptr.bin", "L_colind.bin", "L_data.bin");
//   bfMatCsrRealDump(bfMatToMatCsrReal(M), "M_rowptr.bin", "M_colind.bin", "M_data.bin");

  /* The natural frequency of each eigenvector is the square root of
   * the associated eigenvalue. */
  BfReal freqMax = sqrt(lamMax);

  /* Set up the frequency tree. Note: we build the tree on the
   * frequency scale as opposed to the eigenvalue scale to preserve
   * the time-frequency product in the butterfly factorization. */
  BfIntervalTree *freqTree = bfIntervalTreeNew();
  bfIntervalTreeInitEmpty(freqTree, 0, freqMax, 2, freqTreeDepth);
  HANDLE_ERROR();

  /* Upcast frequency tree to get the column tree */
  BfTree *colTree = bfIntervalTreeToTree(freqTree);

  BfPoints1 *freqs = bfPoints1New();
  HANDLE_ERROR();

  bfPoints1InitEmpty(freqs, BF_ARRAY_DEFAULT_CAPACITY);
  HANDLE_ERROR();

  BfFacSpec spec = {
    .rowTree = rowTree,
    .colTree = colTree,
    .rowTreeInitDepth = rowTreeOffset,
    .colTreeInitDepth = freqTreeDepth, // TODO: this is unused!
    .tol = tol,
    .minNumRows = 20,
    .minNumCols = 20,
  };

  /* Set up the depth-first butterfly factorization streamer. We'll
   * use this below to construct the butterfly factorization
   * incrementally. */
  BfFacStreamer *facStreamer = bfFacStreamerNew();
  bfFacStreamerInit(facStreamer, &spec);
  HANDLE_ERROR();

  while (!bfFacStreamerIsDone(facStreamer)) {
    feedFacStreamerNextEigenband(facStreamer, freqs, L, M);
    HANDLE_ERROR();

    if (freqs->size >= numEigs) break;
  }

  puts("- finished streaming butterfly factorization");

  /* TODO: prune factorization */

  END_ERROR_HANDLING() {}

  /* Clean up */
  bfTreeDelete(&rowTree);
  bfOctreeDeinit(&octree);
  bfMatDelete(&M);
  bfMatDelete(&L);
  bfTrimeshDeinit(&trimesh);
}