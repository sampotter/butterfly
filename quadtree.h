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

  size_t depth;
};

enum BfError
bfInitQuadtreeFromPoints(BfQuadtree *tree,
                         size_t num_points, double const (*points)[2]);

enum BfError
bfGetQuadtreeIndexRange(BfQuadtree const *tree,
                        size_t level, size_t i,
                        size_t *num_indices, size_t **indices);
