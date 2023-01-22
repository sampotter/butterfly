#include <math.h>
#include <stdlib.h>

#include <bf/fac_streamer.h>
// #include <bf/fiedler_tree.h>
#include <bf/interval_tree.h>
#include <bf/interval_tree_node.h>
#include <bf/lbo.h>
#include <bf/linalg.h>
#include <bf/mat.h>
#include <bf/octree.h>

int main(int argc, char const *argv[]) {
  if (argc != 4) {
    printf("usage: %s objPath tol freqTreeDepth\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  char const *objPath = argv[1];
  BfReal tol = strtod(argv[2], NULL);
  BfSize freqTreeDepth = strtoull(argv[3], NULL, 10);

  /* Load triangle mesh from binary files */
  BfTrimesh trimesh;
  bfTrimeshInitFromObjFile(&trimesh, objPath);

  // TODO: implement
//   /* Build fiedler tree */
//   BfFiedlerTree fiedlerTree;
//   bfFiedlerTreeInitFromTrimesh(&fiedlerTree, trimesh, &L, &M);
//   BfTree *rowTree = bfFiedlerTreeToTree(&fiedlerTree);

  BfOctree octree;
  bfOctreeInit(&octree, &trimesh.verts, NULL);

  char const *octreeBoxesPath = "octree_boxes.txt";
  bfOctreeSaveBoxesToTextFile(&octree, octreeBoxesPath);
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

  /* Find the largest eigenvalue. We need this to determine the
   * interval on which we'll build the frequency tree. */
  BfReal lamMax = bfGetMaxEigenvalue(L, M);
  printf("- lambda_max = %g\n", lamMax);

  /* The natural frequency of each eigenvector is the square root of
   * the associated eigenvalue. */
  BfReal freqMax = sqrt(lamMax);

  /* Set up the frequency tree. Note: we build the tree on the
   * frequency scale as opposed to the eigenvalue scale to preserve
   * the time-frequency product in the butterfly factorization. */
  BfIntervalTree *freqTree = bfIntervalTreeNew();
  bfIntervalTreeInitEmpty(freqTree, 0, freqMax, 2, rowTreeMaxDepth);

  /* Upcast frequency tree to get the column tree */
  BfTree *colTree = bfIntervalTreeToTree(freqTree);
  BfSize colTreeMaxDepth = bfTreeGetMaxDepth(colTree);
  assert(colTreeMaxDepth == rowTreeMaxDepth);

  /* Set up the depth-first butterfly factorization streamer. We'll
   * use this below to construct the butterfly factorization
   * incrementally. */
  BfFacStreamer *facStreamer = bfFacStreamerNew();
  bfFacStreamerInit(facStreamer, rowTree, colTree, 1, freqTreeDepth, tol);

  /* Feed eigenvalues until done */
  BfPtrArray colTreeLeafNodes = bfTreeGetLevelPtrArray(colTree, colTreeMaxDepth);
  for (BfSize i = 0; i < bfPtrArraySize(&colTreeLeafNodes); ++i) {
    BfIntervalTreeNode const *node = bfPtrArrayGet(&colTreeLeafNodes, i);
    BfReal lam0 = pow(node->a, 2);
    BfReal lam1 = pow(node->b, 2);
    BfMat *Phi, *Lam;
    bfGetEigenband(L, M, lam0, lam1, &Phi, &Lam);
    bfFacStreamerFeed(facStreamer, Phi);
    bfMatDelete(&Phi);
    bfMatDelete(&Lam);
  }

  BfMat *fac = bfFacStreamerGetFac(facStreamer);

  /* Do numerical tests */

  /* Clean up */

  bfMatDelete(&fac);
  bfPtrArrayDeinit(&colTreeLeafNodes);
  bfTreeDelete(&rowTree);
  bfOctreeDeinit(&octree);
  bfMatDelete(&M);
  bfMatDelete(&L);
  bfTrimeshDeinit(&trimesh);
}
