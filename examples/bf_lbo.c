#include <math.h>
#include <stdlib.h>

#include <bf/fac_streamer.h>
// #include <bf/fiedler_tree.h>
#include <bf/interval_tree.h>
#include <bf/lbo.h>
#include <bf/mat.h>
#include <bf/octree.h>

int main(int argc, char const *argv[]) {
  if (argc != 4) {
    printf("usage: %s objPath tol initDepth\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  char const *objPath = argv[1];
  BfReal tol = strtod(argv[2], NULL);
  BfSize colTreeInitDepth = strtoull(argv[3], NULL, 10);

  /* Load triangle mesh from binary files */
  BfTrimesh trimesh;
  bfTrimeshInitFromObjFile(&trimesh, objPath);

  /* Stiffness and mass matrices */
  BfMatCsrReal *L, *M;

  // TODO: implement
//   /* Build fiedler tree */
//   BfFiedlerTree fiedlerTree;
//   bfFiedlerTreeInitFromTrimesh(&fiedlerTree, trimesh, &L, &M);
//   BfTree *rowTree = bfFiedlerTreeToTree(&fiedlerTree);

  BfOctree octree;
  bfOctreeInitFromPoints(&octree, &trimesh.verts, NULL);
  BfTree *rowTree = bfOctreeToTree(&octree);

  bfLboGetFemDiscretization(&trimesh, &L, &M);

  bfMatCsrRealDump(L, "L_rowptr.bin", "L_colind.bin", "L_data.bin");
  bfMatCsrRealDump(M, "M_rowptr.bin", "M_colind.bin", "M_data.bin");

  /* Find largest eigenvalue */
  BfReal lamMax = bfMatGetEigMaxGen(bfMatCsrRealToMat(L), bfMatCsrRealToMat(M));

  /* Set up frequency tree */
  BfReal freqMax = sqrt(lamMax);
  BfIntervalTree *intervalTree = bfIntervalTreeNew();
  bfIntervalTreeInitEmpty(intervalTree, 0, freqMax, 2);
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
    bfMatGetEigBandGen(bfMatCsrRealToMat(L), bfMatCsrRealToMat(M), lam0, lam1, &Phi, &Lam);
    bfFacStreamerFeed(facStreamer, Phi);
    bfMatDelete(&Phi);
    bfMatDelete(&Lam);
  }

  BfMat *fac = bfFacStreamerGetFac(facStreamer);

  /* Do numerical tests */

  /* Clean up */

  bfMatDelete(&fac);
}
