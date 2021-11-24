#include "quadtree.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define SWAP(x, y) do {                         \
    __typeof__(x) tmp = x;                      \
    x = y;                                      \
    y = tmp;                                    \
  } while (0)

enum BfQuadtreeNodeFlags child_flags[4] = {
  BF_QUADTREE_NODE_FLAG_CHILD_0,
  BF_QUADTREE_NODE_FLAG_CHILD_1,
  BF_QUADTREE_NODE_FLAG_CHILD_2,
  BF_QUADTREE_NODE_FLAG_CHILD_3
};

enum BfQuadtreeNodeFlags child_mask =
  BF_QUADTREE_NODE_FLAG_CHILD_0 |
  BF_QUADTREE_NODE_FLAG_CHILD_1 |
  BF_QUADTREE_NODE_FLAG_CHILD_2 |
  BF_QUADTREE_NODE_FLAG_CHILD_3;

size_t child_flag_to_index[] = {
  [BF_QUADTREE_NODE_FLAG_CHILD_0] = 0,
  [BF_QUADTREE_NODE_FLAG_CHILD_1] = 1,
  [BF_QUADTREE_NODE_FLAG_CHILD_2] = 2,
  [BF_QUADTREE_NODE_FLAG_CHILD_3] = 3
};

static
enum BfError
recInitQuadtreeFromPoints(BfQuadtreeNode *node,
                          BfPoint2 const *points, BfBbox2 bbox,
                          size_t perm_size, size_t *perm)
{
  enum BfError error = BF_ERROR_NO_ERROR;
  size_t i, j;

  BfReal const *split = node->split;

  // set sentinel values
  node->offset[0] = 0;
  node->offset[4] = perm_size;

  // sift permutation indices for child[0] into place...
  i = node->offset[0];
  while (i < perm_size &&
         (points[perm[i]][0] <= split[0] && points[perm[i]][1] <= split[1]))
    ++i;
  j = i == perm_size ? i : i + 1;
  while (j < perm_size) {
    if (points[perm[j]][0] <= split[0] && points[perm[j]][1] <= split[1] &&
        (points[perm[i]][0] > split[0] || points[perm[i]][1] > split[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  node->offset[1] = i;

#ifdef BF_DEBUG
  // (verify that child[0]'s indices are correctly sifted)
  for (size_t k = 0; k < perm_size; ++k) {
    bool in_quadrant =
      points[perm[k]][0] <= split[0] && points[perm[k]][1] <= split[1];
    assert(node->offset[0] <= k && k < node->offset[1] ?
           in_quadrant : !in_quadrant);
  }
#endif

  // ... for child[1]...
  i = node->offset[1];
  while (i < perm_size &&
         (points[perm[i]][0] <= split[0] && points[perm[i]][1] > split[1]))
    ++i;
  j = i == perm_size ? i : i + 1;
  while (j < perm_size) {
    if (points[perm[j]][0] <= split[0] && points[perm[j]][1] > split[1] &&
        (points[perm[i]][0] > split[0] || points[perm[i]][1] <= split[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  node->offset[2] = i;

#ifdef BF_DEBUG
  // (verify that child[1]'s indices are correctly sifted)
  for (size_t k = 0; k < perm_size; ++k) {
    bool in_quadrant =
      points[perm[k]][0] <= split[0] && points[perm[k]][1] > split[1];
    assert(node->offset[1] <= k && k < node->offset[2] ?
           in_quadrant : !in_quadrant);
  }
#endif

  // ... for child[2]...
  i = node->offset[2];
  while (i < perm_size &&
         (points[perm[i]][0] > split[0] && points[perm[i]][1] <= split[1]))
    ++i;
  j = i == perm_size ? i : i + 1;
  while (j < perm_size) {
    if (points[perm[j]][0] > split[0] && points[perm[j]][1] <= split[1] &&
        (points[perm[i]][0] <= split[0] || points[perm[i]][1] > split[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  node->offset[3] = i;

#ifdef BF_DEBUG
  // (verify that child[2]'s indices are correctly sifted)
  for (size_t k = 0; k < perm_size; ++k) {
    bool in_quadrant =
      points[perm[k]][0] > split[0] && points[perm[k]][1] <= split[1];
    assert(node->offset[2] <= k && k < node->offset[3] ?
           in_quadrant : !in_quadrant);
  }
#endif

  // we don't have to do any sifting for the last child! nice.
  // ... but let's verify that things are as they should be.

#ifdef BF_DEBUG
  assert(j == perm_size);
  assert(j == node->offset[4]);
#endif

#ifdef BF_DEBUG
  // (verify that child[2]'s indices are correctly sifted)
  for (size_t k = 0; k < perm_size; ++k) {
    bool in_quadrant =
      points[perm[k]][0] > split[0] && points[perm[k]][1] > split[1];
    assert(node->offset[3] <= k && k < node->offset[4] ?
           in_quadrant : !in_quadrant);
  }
#endif

#ifdef BF_DEBUG
  // (verify that child[3]'s indices are correctly sifted)
  for (size_t i = node->offset[3], j; i < perm_size; ++i) {
    j = perm[i];
    assert(points[j][0] > split[0] && points[j][1] > split[1]);
  }
#endif

  // compute the bounding boxes each child node
  BfBbox2 child_bbox[4] = {
    [0] = {
      .min = {bbox.min[0], bbox.min[1]},
      .max = {split[0], split[1]}
    },
    [1] = {
      .min = {bbox.min[0], split[1]},
      .max = {split[0], bbox.max[1]}
    },
    [2] = {
      .min = {split[0], bbox.min[1]},
      .max = {bbox.max[0], split[1]}
    },
    [3] = {
      .min = {split[0], split[1]},
      .max = {bbox.max[0], bbox.max[1]}
    }
  };

  for (size_t q = 0, perm_size; q < 4; ++q) {
    perm_size = node->offset[q + 1] - node->offset[q];

    // TODO: for now, we build the quadtree out to the maximum depth possible
    if (perm_size <= 1) {
      node->child[q] = NULL;
      continue;
    }

    node->child[q] = malloc(sizeof(BfQuadtreeNode));
    node->child[q]->flags = child_flags[q];
    node->child[q]->parent = node;
    node->child[q]->bbox = child_bbox[q];

    // compute the split for the `q`th child node
    node->child[q]->split[0] = (child_bbox[q].min[0] + child_bbox[q].max[0])/2;
    node->child[q]->split[1] = (child_bbox[q].min[1] + child_bbox[q].max[1])/2;

    error |= recInitQuadtreeFromPoints(
      node->child[q], points, child_bbox[q], perm_size, &perm[node->offset[q]]);

    if (error)
      break;
  }

  if (error)
    for (size_t q = 0; q < 4; ++q)
      free(node->child[q]);

  return error;
}

static
size_t getQuadtreeNodeDepth(BfQuadtreeNode const *node,
                            size_t current_depth) {
  size_t next_depth = current_depth;
  for (size_t i = 0; i < 4; ++i) {
    size_t child_depth = node->child[i] == NULL ?
      current_depth :
      getQuadtreeNodeDepth(node->child[i], current_depth + 1);
    next_depth = child_depth > next_depth ? child_depth : next_depth;
  }
  return next_depth;
}

enum BfError
bfInitQuadtreeFromPoints(BfQuadtree *tree,
                         size_t num_points, BfPoint2 const *points)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  tree->num_points = num_points;
  tree->points = points;

  // initialize permutation to identity
  tree->perm = malloc(num_points*sizeof(size_t));
  for (size_t i = 0; i < num_points; ++i)
    tree->perm[i] = i;

  tree->root = malloc(sizeof(BfQuadtreeNode));
  tree->root->flags = BF_QUADTREE_NODE_FLAG_ROOT;
  tree->root->parent = tree;

  // compute the bounding box for the entire quadtree
  BfBbox2 bbox = {
    .min = {INFINITY, INFINITY},
    .max = {-INFINITY, -INFINITY}
  };
  for (size_t i = 0; i < num_points; ++i) {
    bbox.min[0] = fmin(bbox.min[0], points[i][0]);
    bbox.max[0] = fmax(bbox.max[0], points[i][0]);
    bbox.min[1] = fmin(bbox.min[1], points[i][1]);
    bbox.max[1] = fmax(bbox.max[1], points[i][1]);
  }

  assert(bbox.min[0] < bbox.max[0]);
  assert(bbox.min[1] < bbox.max[1]);

  // rescale the bounding box so that it's square
  BfReal w = bbox.max[0] - bbox.min[0];
  BfReal h = bbox.max[1] - bbox.min[1];
  if (w > h) {
    BfReal c = (bbox.min[1] + bbox.max[1])/2;
    bbox.min[1] = w*(bbox.min[1] - c)/h + c;
    bbox.max[1] = w*(bbox.max[1] - c)/h + c;
  } else {
    BfReal c = (bbox.min[0] + bbox.max[1])/2;
    bbox.min[0] = h*(bbox.min[0] - c)/h + c;
    bbox.max[0] = h*(bbox.max[0] - c)/h + c;
  }

  tree->root->bbox = bbox;

  // compute the split for the root node
  tree->root->split[0] = (bbox.min[0] + bbox.max[0])/2;
  tree->root->split[1] = (bbox.min[1] + bbox.max[1])/2;

  error |= recInitQuadtreeFromPoints(tree->root, tree->points, bbox,
                                     tree->num_points, tree->perm);

  tree->depth = getQuadtreeNodeDepth(tree->root, 0);

  return error;
}

enum BfError
bfGetQuadtreeNode(BfQuadtree const *tree, size_t depth, size_t node_index,
                  BfQuadtreeNode **node)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  size_t nodes_at_depth = pow(4.0, depth);
  if (node_index >= nodes_at_depth)
    return error | BF_ERROR_INVALID_ARGUMENTS;

  size_t i, r = node_index;

  BfQuadtreeNode *current_node = tree->root;
  do {
    nodes_at_depth /= 4;
    i = r/nodes_at_depth;
    r = r % nodes_at_depth;
    current_node = current_node->child[i];
    if (!current_node)
      return error | BF_ERROR_INVALID_ARGUMENTS | BF_ERROR_RUNTIME_ERROR;
  } while (--depth > 0);

  *node = current_node;

  return error;
}

enum BfError
bfGetQuadtreeNodeIndices(BfQuadtreeNode const *node,
                         size_t *num_indices, size_t **indices)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  assert(node->offset[0] == 0); // I think?

  *num_indices = node->offset[4];

  BfQuadtreeNode const *parent_node;
  size_t q, offset = 0;

  while (node->flags & child_mask) {
    q = child_flag_to_index[node->flags];
    parent_node = node->parent;
    offset += parent_node->offset[q];
    node = parent_node;
  }

  BfQuadtree const *tree = node->parent;

  *indices = tree->perm + offset;

  return error;
}
