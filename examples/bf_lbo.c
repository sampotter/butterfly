#include <math.h>
#include <stdlib.h>

#include <bf/eig.h>
#include <bf/fac_streamer.h>
// #include <bf/fiedler_tree.h>
#include <bf/interval_tree.h>
#include <bf/lbo.h>
#include <bf/octree.h>
#include <bf/spmat.h>

int main(int argc, char const *argv[]) {
  if (argc != 5) {
    printf("usage: %s vertsPath facesPath tol initDepth\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  char const *vertsPath = argv[1];
  char const *facesPath = argv[2];
  BfReal tol = strtod(argv[3], NULL);
  BfSize colTreeInitDepth = strtoull(argv[4], NULL, 10);

  /* Load triangle mesh from binary files */
  BfTrimesh *trimesh = bfTrimeshNew();
  bfTrimeshInitFromBinaryFiles(trimesh, vertsPath, facesPath);

  /* Stiffness and mass matrices */
  BfSpmat *L, *M;

  // TODO: implement
//   /* Build fiedler tree */
//   BfFiedlerTree fiedlerTree;
//   bfFiedlerTreeInitFromTrimesh(&fiedlerTree, trimesh, &L, &M);
//   BfTree *rowTree = bfFiedlerTreeToTree(&fiedlerTree);

  BfOctree octree;
  bfOctreeInitFromTrimesh(&octree, trimesh);
  BfTree *rowTree = bfOctreeToTree(&octree);

  bfLboGetFemDiscretization(trimesh, &L, &M);

  /* Find largest eigenvalue */
  BfReal lamMax = bfGetEigMax(bfSpmatToMat(L), bfSpmatToMat(M));

  /* Set up frequency tree */
  BfReal freqMax = sqrt(lamMax);
  BfIntervalTree *intervalTree = bfIntervalTreeNew();
  bfIntervalTreeInit(intervalTree, 0, freqMax, colTreeInitDepth, 2);
  BfTree *colTree = bfIntervalTreeToTree(intervalTree);

  /* Set up fac streamer */
  BfFacStreamer *facStreamer = bfFacStreamerNew();
  bfFacStreamerInit(facStreamer, rowTree, colTree,
                    /* rowTreeInitDepth */ 1, colTreeInitDepth, tol);

  /* Feed eigenvalues until done */
  BfSize n = 1 << colTreeInitDepth;
  for (BfSize j = 0; j < n; ++j) {
    BfReal j0 = j, j1 = j + 1;
    BfReal lam0 = lamMax*j0/n, lam1 = lamMax*j1/n;
    BfMat *Phi;
    BfMat *Lam;
    bfGetEigBand(L, M, lam0, lam1, &Phi, &Lam);
    bfFacStreamerFeed(facStreamer, Phi);
    bfMatDelete(&Phi);
    bfMatDelete(&Lam);
  }

  BfMat *fac = bfFacStreamerGetFac(facStreamer);

  /* Do numerical tests */

  /* Clean up */

  bfMatDelete(&fac);
}
