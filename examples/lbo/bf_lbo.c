#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/fac_streamer.h>
#include <bf/interval_tree.h>
#include <bf/lbo.h>
#include <bf/linalg.h>
#include <bf/logging.h>
#include <bf/octree.h>
#include <bf/rand.h>

#include <math.h>
#include <stdlib.h>

#include "lbo.h"

int main(int argc, char const *argv[]) {
  BEGIN_ERROR_HANDLING() {}

  if (argc != 6) {
    printf("usage: %s objPath tol p rowTreeOffset freqTreeDepth\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  bfSeed(0);
  bfSetLogLevel(BF_LOG_LEVEL_INFO);

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
