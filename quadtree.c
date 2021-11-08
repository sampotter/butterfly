#include "quadtree.h"

#include <assert.h>
#include <string.h>

#define SWAP(x, y) do {							\
    __typeof__(x) tmp = x;                      \
    x = y;                                      \
    y = x;                                      \
  } while (0)

static
enum BfError
recInitQuadtreeFromPoints(BfQuadtreeNode *node,
                          double const (*points)[2],
                          size_t perm_size, size_t *perm)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  // compute mean of points indexed by perm
  double mid[2] = {0, 0};
  for (size_t i = 0, j; i < perm_size; ++i) {
    j = perm[i];
    mid[0] += points[j][0];
    mid[1] += points[j][1];
  }
  mid[0] /= perm_size;
  mid[1] /= perm_size;

  // set sentinel values
  node->offset[0] = 0;
  node->offset[4] = perm_size;

  // sift permutation indices for child[0] into place...
  node->offset[1] = node->offset[0];
  for (size_t i = node->offset[1], j; i < perm_size; ++i) {
    j = perm[i];
    if (points[j][0] <= mid[0] && points[j][1] <= mid[1])
      continue;
    SWAP(perm[node->offset[1]], perm[i]);
    ++node->offset[1];
  }

  // ... for child[1]...
  node->offset[2] = node->offset[1];
  for (size_t i = node->offset[2], j; i < perm_size; ++i) {
    j = perm[i];
    if (points[j][0] <= mid[0] && points[j][1] > mid[1])
      continue;
    SWAP(perm[node->offset[2]], perm[i]);
    ++node->offset[2];
  }

  // ... and for child[2]...
  node->offset[3] = node->offset[2];
  for (size_t i = node->offset[3], j; i < perm_size; ++i) {
    j = perm[i];
    if (points[j][0] > mid[0] && points[j][1] <= mid[1])
      continue;
    SWAP(perm[node->offset[3]], perm[i]);
    ++node->offset[3];
  }

  // ... remaining indices for child[3] should be sifted correctly
  for (size_t i = node->offset[3], j; i < perm_size; ++i) {
    j = perm[i];
    assert(points[j][0] > mid[0] && points[j][1] > mid[1]);
  }

  for (size_t q = 0, perm_size; q < 4; ++q) {
    node->child[q] = malloc(sizeof(BfQuadtreeNode));
    perm_size = node->offset[q + 1] - node->offset[q];
    recInitQuadtreeFromPoints(node->child[q], points,
                              perm_size, &perm[node->offset[q]]);
  }

  return error;
}

enum BfError
bfInitQuadtreeFromPoints(BfQuadtree *tree, size_t max_leaf_size,
                         size_t num_points, double const (*points)[2])
{
  enum BfError error = BF_ERROR_NO_ERROR;

  tree->max_leaf_size = max_leaf_size;

  tree->num_points = num_points;
  tree->points = points;

  // initialize permutation to identity
  tree->perm = malloc(num_points*sizeof(size_t));
  for (size_t i = 0; i < num_points; ++i)
    tree->perm[i] = i;

  tree->root = malloc(sizeof(BfQuadtreeNode));

  error |= recInitQuadtreeFromPoints(tree->root, tree->points,
                                     tree->num_points, tree->perm);

  return error;
}
