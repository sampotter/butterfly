#include <bf/fac_span.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_block_dense.h>
#include <bf/mem.h>

BfFacSpan *bfFacSpanNewFromPtrArray(BfPtrArray const *facs) {
  BF_ERROR_BEGIN();

  BfFacSpan *facSpan = bfMemAlloc(1, sizeof(BfFacSpan));
  HANDLE_ERROR();

  bfFacSpanInitFromPtrArray(facSpan, facs);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    BF_DIE();
  }

  return facSpan;
}

static bool facsHaveSameRowSpan(BfPtrArray const *facs) {
  BF_ASSERT(!bfPtrArrayIsEmpty(facs));

  BF_ERROR_BEGIN();

  /* Get the first factorization and compare all the others' row spans
   * to this one's row span. */
  BfFac const *fac = bfPtrArrayGet(facs, 0);

  bool sameRowSpan = true;
  for (BfSize i = 1; i < bfPtrArraySize(facs); ++i) {
    BfFac const *otherFac = bfPtrArrayGet(facs, i);

    bool same = bfConstNodeArrayIsSameSpan(&fac->rowNodes, &otherFac->rowNodes);
    HANDLE_ERROR();

    if (!same) {
      sameRowSpan = false;
      break;
    }
  }

  END_ERROR_HANDLING() {
    BF_DIE();
  }

  return sameRowSpan;
}

static bool facsDetermineContiguousColNodeSpan(BfPtrArray const *facs) {
  BF_ASSERT(!bfPtrArrayIsEmpty(facs));

  BF_ERROR_BEGIN();

  bool contiguous = false;

  BfConstNodeArray colNodes;
  bfConstNodeArrayInitWithDefaultCapacity(&colNodes);
  HANDLE_ERROR();

  for (BfSize i = 0; i < bfPtrArraySize(facs); ++i) {
    BfFac const *fac = bfPtrArrayGet(facs, i);

    bfConstNodeArrayAppend(&colNodes, fac->colNode);
    HANDLE_ERROR();
  }

  contiguous = bfConstNodeArrayIsContiguous(&colNodes);

  END_ERROR_HANDLING() {
    BF_DIE();
  }

  bfConstNodeArrayDeinit(&colNodes);

  return contiguous;
}

void bfFacSpanInitFromPtrArray(BfFacSpan *facSpan, BfPtrArray const *facs) {
  BF_ERROR_BEGIN();

  if (bfPtrArrayIsEmpty(facs))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (!facsHaveSameRowSpan(facs))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (!facsDetermineContiguousColNodeSpan(facs))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  facSpan->numFacs = bfPtrArraySize(facs);

  facSpan->fac = bfMemAllocAndZero(facSpan->numFacs, sizeof(BfFac *));
  HANDLE_ERROR();

  for (BfSize i = 0; i < facSpan->numFacs; ++i)
    facSpan->fac[i] = bfPtrArrayGet(facs, i);

  END_ERROR_HANDLING() {
    BF_DIE();
  }
}

BfMat *bfFacSpanGetMat(BfFacSpan const *facSpan) {
  BF_ERROR_BEGIN();

  BfMatBlockDense *matBlockDense = NULL;

  BfPtrArray *rowBlocks = bfPtrArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  for (BfSize i = 0; i < facSpan->numFacs; ++i) {
    BfMat *rowBlock = bfFacGetMat(facSpan->fac[i]);
    HANDLE_ERROR();

    bfPtrArrayAppend(rowBlocks, rowBlock);
    HANDLE_ERROR();
  }

  matBlockDense = bfMatBlockDenseNewRowFromBlocks(rowBlocks, BF_POLICY_STEAL);

  END_ERROR_HANDLING() {
    BF_DIE();
  }

  bfPtrArrayDelete(&rowBlocks);

  return bfMatBlockDenseToMat(matBlockDense);
}
