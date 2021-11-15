#pragma once

#include <stdlib.h>

#include "error.h"

typedef struct BfQuadtree BfQuadtree;
typedef struct BfQuadtreeNode BfQuadtreeNode;

struct BfQuadtreeNode {
  size_t offset[5]; // offset[0] == 0 and offset[4] == size are
                    // sentinel values
  BfQuadtreeNode *child[4];
};

struct BfQuadtree {
  size_t num_points;
  double const (* points)[2];

  size_t *perm;

  BfQuadtreeNode *root;
};

enum BfError
bfInitQuadtreeFromPoints(BfQuadtree *tree,
                         size_t num_points, double const (*points)[2]);
