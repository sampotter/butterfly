#include <bf/fiedler_tree.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/lbo.h>

void bfFiedlerTreeInitFromTrimesh(BfFiedlerTree *fiedlerTree,
                                  BfTrimesh const *trimesh) {
  BF_ERROR_BEGIN();

  (void)fiedlerTree;
  (void)trimesh;

  RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BF_ERROR_END() {}
}
