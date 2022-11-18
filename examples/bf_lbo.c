#include <bf/trimesh.h>

int main(int argc, char const *argv[]) {
  char const *vertsPath = argv[1];
  char const *facesPath = argv[2];

  // load triangle mesh from binary files
  BfTrimesh *trimesh = bfTrimeshNew();
  bfTrimeshInitFromBinaryFiles(trimesh, vertsPath, facesPath);

  // build fiedler tree
  BfTree *rowTree = bfFiedlerTreeNew();
  bfFiedlerTreeInitFromTrimesh(rowTree, trimesh);

  // set up stiffness and mass matrices
  BfSpmat *L, *M;
  bfLboGetFemDiscretization(trimesh, &L, &M);

  // find largest eigenvalue
  BfReal lamMax = bfGetMaxEig(L, M);

  // set up frequency tree
  BfReal freqMax = sqrt(lamMax);
  BfTree *colTree = bfIntervalTreeNew();
  bfIntervalTreeInit(0, freqMax, colTreeInitDepth);

  // set up fac streamer
  BfFacStreamer *facStreamer = bfFacStreamerNew();
  bfFacStreamerInit(facStreamer, rowTree, colTree, BF_STREAM_MODE_POST_ORDER,
                    /* rowTreeInitDepth */ 1, colTreeInitDepth, tol);

  // feed eigenvalues until done
  BfSize n = 1 << colTreeInitDepth;
  for (BfSize j = 0; j < n; ++j) {
    BfReal j0 = j, j1 = j + 1;
    BfReal lam0 = lamMax*j0/n, lam1 = lamMax*j1/n;
    BfMat *Phi;
    BfMatDiagReal *Lam;
    bfGetEigBand(L, M, lam0, lam1, &Phi, &Lam);
    bfFacStreamerFeed(facStreamer, Phi);
    bfMatDelete(&Phi);
    bfMatDiagRealDeinitAndDealloc(&Lam);
  }

  bfFacStreamerFinish(facStreamer);

  BfMat *fac = bfFacStreamerGetFac(facStreamer);

  // do numerical tests
}
