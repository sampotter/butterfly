#include "quadtree.h"

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#define SWAP(x, y) do {							\
    __typeof__(x) tmp = x;                      \
    x = y;                                      \
    y = tmp;                                    \
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

  size_t i, j;

  // sift permutation indices for child[0] into place...
  i = node->offset[0];
  while (i < perm_size &&
         (points[perm[i]][0] <= mid[0] && points[perm[i]][1] <= mid[1]))
    ++i;
  j = i == perm_size ? i : i + 1;
  while (j < perm_size) {
    if (points[perm[j]][0] <= mid[0] && points[perm[j]][1] <= mid[1] &&
        (points[perm[i]][0] > mid[0] || points[perm[i]][1] > mid[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  node->offset[1] = i;

  // (verify that child[0]'s indices are correctly sifted)
  for (size_t k = 0; k < perm_size; ++k) {
    bool in_quadrant =
      points[perm[k]][0] <= mid[0] && points[perm[k]][1] <= mid[1];
    assert(node->offset[0] <= k && k < node->offset[1] ?
           in_quadrant : !in_quadrant);
  }

  // ... for child[1]...
  i = node->offset[1];
  while (i < perm_size &&
         (points[perm[i]][0] <= mid[0] && points[perm[i]][1] > mid[1]))
    ++i;
  j = i == perm_size ? i : i + 1;
  while (j < perm_size) {
    if (points[perm[j]][0] <= mid[0] && points[perm[j]][1] > mid[1] &&
        (points[perm[i]][0] > mid[0] || points[perm[i]][1] <= mid[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  node->offset[2] = i;

  // (verify that child[1]'s indices are correctly sifted)
  for (size_t k = 0; k < perm_size; ++k) {
    bool in_quadrant =
      points[perm[k]][0] <= mid[0] && points[perm[k]][1] > mid[1];
    assert(node->offset[1] <= k && k < node->offset[2] ?
           in_quadrant : !in_quadrant);
  }

  // ... for child[2]...
  i = node->offset[2];
  while (i < perm_size &&
         (points[perm[i]][0] > mid[0] && points[perm[i]][1] <= mid[1]))
    ++i;
  j = i == perm_size ? i : i + 1;
  while (j < perm_size) {
    if (points[perm[j]][0] > mid[0] && points[perm[j]][1] <= mid[1] &&
        (points[perm[i]][0] <= mid[0] || points[perm[i]][1] > mid[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  node->offset[3] = i;

  // (verify that child[2]'s indices are correctly sifted)
  for (size_t k = 0; k < perm_size; ++k) {
    bool in_quadrant =
      points[perm[k]][0] > mid[0] && points[perm[k]][1] <= mid[1];
    assert(node->offset[2] <= k && k < node->offset[3] ?
           in_quadrant : !in_quadrant);
  }

  // we don't have to do any sifting for the last child! nice.
  // ... but let's verify that things are as they should be.

  assert(j == perm_size);
  assert(j == node->offset[4]);

  // (verify that child[2]'s indices are correctly sifted)
  for (size_t k = 0; k < perm_size; ++k) {
    bool in_quadrant =
      points[perm[k]][0] > mid[0] && points[perm[k]][1] > mid[1];
    assert(node->offset[3] <= k && k < node->offset[4] ?
           in_quadrant : !in_quadrant);
  }

  // (verify that child[3]'s indices are correctly sifted)
  for (size_t i = node->offset[3], j; i < perm_size; ++i) {
    j = perm[i];
    assert(points[j][0] > mid[0] && points[j][1] > mid[1]);
  }

  for (size_t q = 0, perm_size; q < 4; ++q) {
    perm_size = node->offset[q + 1] - node->offset[q];

    // TODO: for now, we build the quadtree out to the max depth
    if (perm_size <= 1) {
      node->child[q] = NULL;
    } else {
      node->child[q] = malloc(sizeof(BfQuadtreeNode));
      recInitQuadtreeFromPoints(
        node->child[q], points, perm_size, &perm[node->offset[q]]);
    }
  }

  return error;
}

enum BfError
bfInitQuadtreeFromPoints(BfQuadtree *tree,
                         size_t num_points, double const (*points)[2])
{
  enum BfError error = BF_ERROR_NO_ERROR;

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
