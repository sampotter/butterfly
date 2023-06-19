#include <bf/assert.h>
#include <bf/error.h>
#include <bf/fiedler_tree.h>
#include <bf/util.h>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("usage: %s <trimesh.obj>\n", argv[0]);
    return 1;
  }

  char const *objPath = argv[1];

  BfReal tol = 1e-15;

  bfToc();

  BfTrimesh *trimesh = bfTrimeshNewFromObjFile(objPath);

  printf("loaded triangle mesh (%lu verts and %lu faces) [%.1fs]\n",
         bfTrimeshGetNumVerts(trimesh), bfTrimeshGetNumFaces(trimesh), bfToc());

  BfFiedlerTree *fiedlerTree = bfFiedlerTreeNewFromTrimesh(trimesh, tol, /* keepNodeTrimeshes: */ true);

  printf("built Fiedler tree [%.1fs]\n", bfToc());

  bfFiedlerTreeDeinitAndDealloc(&fiedlerTree);
  bfTrimeshDeinitAndDealloc(&trimesh);
}
