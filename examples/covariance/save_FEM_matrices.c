#include <bf/assert.h>
#include <bf/const.h>
#include <bf/fac_span.h>
#include <bf/fac_streamer.h>
#include <bf/interval_tree.h>
#include <bf/linalg.h>
#include <bf/logging.h>
#include <bf/mat_csr_real.h>
#include <bf/mat_dense_real.h>
#include <bf/octree.h>
#include <bf/rand.h>
#include <bf/util.h>
#include <bf/vec_real.h>
#include <bf/trimesh.h>
#include <bf/lbo.h>

#include <time.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char const *argv[]) {
  if (argc < 2) {
    printf("usage: %s mesh.obj \n", argv[0]);
    exit(EXIT_FAILURE);
  }

  bfSeed(0);
  bfSetLogLevel(BF_LOG_LEVEL_INFO);

  char const *objPath = argv[1];
  BfTrimesh *trimesh = bfTrimeshNewFromObjFile(objPath);
  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);
  printf("triangle mesh with %lu verts\n", numVerts);

  BfMat *L, *M;
  bfTrimeshGetLboFemDiscretization(trimesh, &L, &M);
  printf("set up FEM discretization [%0.1fs]\n", bfToc());

  bfMatCsrRealDump(bfMatToMatCsrReal(L), "L_rowptr.bin", "L_colind.bin", "L_data.bin");
  bfMatCsrRealDump(bfMatToMatCsrReal(M), "M_rowptr.bin", "M_colind.bin", "M_data.bin");

  // // nullspace of LBO?
  // int nev = numVerts-2;
  // BfMat *Phi = NULL;
  // BfVecReal *Lam = NULL;
  // bfGetShiftedEigs(L, M, -0.001, nev, &Phi, &Lam);
  // printf("Computed dense eigendecomposition [%0.1fs]\n", bfToc());

  // bfMatDenseRealSave(bfMatToMatDenseReal(Phi), "Phi.bin");
  // bfVecRealSave(bfVecToVecReal(Lam), "Lam.bin");
  // printf("Saved dense eigendecomposition\n");

  /* Clean up */
  // bfMatDelete(&M);
  // bfMatDelete(&L);
  // bfMatDelete(&Phi);
  // bfVecDelete(&Lam);
  // bfTrimeshDeinit(&trimesh);
}
