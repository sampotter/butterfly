#pragma once

#include <stdlib.h>

typedef double point2[2];

typedef struct bfQuadtreeNode bfQuadtreeNode;

struct bfQuadtreeNode {
  quadtree_s const *root;
  size_t offset[5]; // offset[0] == 0 and offset[4] == size are sentinel values
  bfQuadtreeNode *child[4];
};

typedef struct {
  point2 const *points;
  size_t num_points;
  size_t *perm;
  bfQuadtreeNode *root;
} bfQuadtree;
