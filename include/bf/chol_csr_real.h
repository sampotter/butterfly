#pragma once

#include "def.h"
#include "chol.h"

/** Interface: Chol */

BfMat *bfCholCsrRealSolve(BfCholCsrReal const *cholCsrReal, BfMat const *B);
BfVec *bfCholCsrRealSolveVec(BfCholCsrReal const *cholCsrReal, BfVec const *b);

/** Upcasting: CholCsrReal -> Chol */

/** Downcasting: Chol -> CholCsrReal */

/** Implementation: CholCsrReal */

typedef struct BfCholCsrRealImpl BfCholCsrRealImpl;

struct BfCholCsrReal {
  BfChol super;
  BfCholCsrRealImpl *impl;
};

BfCholCsrReal *bfCholCsrRealNew(void);
void bfCholCsrRealInit(BfCholCsrReal *cholCsrReal, BfMat const *mat);
void bfCholCsrRealDeinit(BfCholCsrReal *cholCsrReal);
void bfCholCsrRealDealloc(BfCholCsrReal **cholCsrReal);
void bfCholCsrRealDeinitAndDealloc(BfCholCsrReal **cholCsrReal);
