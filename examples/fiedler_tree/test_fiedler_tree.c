#include <bf/assert.h>
#include <bf/error.h>
#include <bf/fiedler_tree.h>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("usage: %s <trimesh.obj>\n", argv[0]);
    return 1;
  }

  char const *objPath = argv[1];

  BfTrimesh *trimesh = bfTrimeshNewFromObjFile(objPath);

  BfReal tol = 1e-15;

  printf("\nbuilding Fielder tree:\n");
  BfFiedlerTree *fiedlerTree = bfFiedlerTreeNewFromTrimesh(trimesh, tol, /* keepNodeTrimeshes: */ true);

  bfFiedlerTreeDeinitAndDealloc(&fiedlerTree);
  bfTrimeshDeinitAndDealloc(&trimesh);
}
