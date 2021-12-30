#pragma once

#include "error.h"

enum BfDecomp {
  BF_DECOMP_SVD
};

enum BfDomain {
  BF_DOMAIN_R1,
  BF_DOMAIN_R2,
  BF_DOMAIN_R3
};

enum BfAccel {
  BF_ACCEL_QUADTREE
};

enum BfMesh {
  BF_MESH_POLY_LOOP
};

enum BfKernel {
  BF_KERNEL_HELMHOLTZ
};

enum BfQuad {
  BF_QUAD_TRAPEZOID
};

typedef struct {
  enum BfDecomp decomp;
  enum BfDomain domain;
  enum BfAccel accel;
  enum BfMesh mesh;
  enum BfKernel kernel;
  enum BfQuad quad;
} BfSpec;
