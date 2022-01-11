#include "mat_block_dense.h"

BfMatBlockDense *bfMatBlockDenseNew() {
  return malloc(sizeof(BfMatBlockDense));
}

void BfMatBlockInit(

void bfMatBlockDenseDeinit(BfMatBlockDense *mat) {

}

void bfMatBlockDenseDelete(BfMatBlockDense **mat) {
}

void bfMatBlockDenseDeinitAndDelete(BfMatBlockDense **mat) {
}

BfMat *bfMatBlockDenseGetMatPtr(BfMatBlockDense *mat) {
}

void bfMatBlockDenseSetBlock(BfMatBlockDense *mat, BfSize i, BfSize j, BfMat *block) {
}
