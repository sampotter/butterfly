#include "quadtree.h"

#include <assert.h>
#include <string.h>

#define SWAP(x, y) do {							\
    __typeof__(x) tmp = x;                      \
    x = y;                                      \
    y = x;                                      \
  } while (0)

void quadtreeNode_recursiveBuildFromPoints(quadtreeNode_s *node,
										   point2 const *points,
										   size_t *perm, size_t perm_size) {
  // compute mean of points indexed by perm
  point2 mid = {0, 0};
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

  for (size_t q = 0; q < 4; ++q) {
    node->child[q] = malloc(sizeof(quadtreeNode_s));
    quadtreeNode_recursiveBuildFromPoints(
      node->child[q], points, &perm[node->offset[q]],
      node->offset[q + 1] - node->offset[q]);
  }
}

void quadtree_initFromPoints2(quadtree_s *qt, point2 const *points, size_t num_points) {
  qt->points = points;
  qt->num_points = num_points;

  // initialize permutation to identity
  qt->perm = malloc(num_points*sizeof(size_t));
  for (size_t i = 0; i < num_points; ++i)
    qt->perm[i] = i;

  qt->root = malloc(sizeof(quadtreeNode_s));

  quadtreeNode_recursiveBuildFromPoints(qt->root, qt->points, qt->perm,
                                        qt->num_points);
}
