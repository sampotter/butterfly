#include <bf/fiedler_tree.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/lbo.h>
#include <bf/spmat.h>

void bfFiedlerTreeInitFromTrimesh(BfFiedlerTree *fiedlerTree,
                                  BfTrimesh const *trimesh) {
  BEGIN_ERROR_HANDLING();

  (void)fiedlerTree;
  (void)trimesh;

  RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  END_ERROR_HANDLING() {}
}
