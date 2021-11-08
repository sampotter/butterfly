#pragma once

#include "error.h"

enum BfDecomp {
  BF_DECOMP_SVD
};

enum BfDim {
  BF_DIM_1D,
  BF_DIM_2D,
  BF_DIM_3D
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
  BF_QUAD_TRAPEZOIDAL
};

struct BfSpec {
  enum BfDecomp decomp;
  enum BfDim dim;
  enum BfAccel accel;
  enum BfMesh mesh;
  enum BfKernel kernel;
  enum BfQuad quad;
};

BfError bfSpecValidate(BfSpec const *spec);
