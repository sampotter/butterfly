#include <assert.h>
#include <stdio.h>

#include "helm2.h"
#include "mat.h"
#include "quadtree.h"

static void bf_one_block(BfQuadtree const *tree, double k) {
  enum BfError error;

  /* Get source and target nodes from quadtree and check that their
   * indices are OK */

  size_t src_depth = 3, src_node_index = 20, src_num_inds, *src_inds;
  size_t tgt_depth = 3, tgt_node_index = 16*3 + 4*3 + 2, tgt_num_inds, *tgt_inds;

  BfQuadtreeNode *src_node;
  bfGetQuadtreeNode(tree, src_depth, src_node_index, &src_node);
  bfGetQuadtreeNodeIndices(src_node, &src_num_inds, &src_inds);

  printf("source node:\n");
  printf("- depth: %lu\n", src_depth);
  printf("- node index: %lu\n", src_node_index);
  printf("- no. of indices: %lu\n", src_num_inds);
  printf("- bbox: [%g, %g] x [%g, %g]\n",
         src_node->bbox.min[0], src_node->bbox.max[0],
         src_node->bbox.min[1], src_node->bbox.max[1]);
  printf("- indices:");
  if (src_num_inds > 10) {
    for (size_t i = 0; i < 10; ++i)
      printf(" %lu,", src_inds[i]);
    printf(" ...\n");
  } else {
    for (size_t i = 0; i < src_num_inds - 1; ++i)
      printf(" %lu,", src_inds[i]);
    printf(" %lu\n", src_inds[src_num_inds - 1]);
  }

  for (size_t i = 0; i < src_num_inds; ++i)
    assert(bfBbox2ContainsPoint(&src_node->bbox, tree->points[src_inds[i]]));

  BfQuadtreeNode *tgt_node;
  bfGetQuadtreeNode(tree, tgt_depth, tgt_node_index, &tgt_node);
  bfGetQuadtreeNodeIndices(tgt_node, &tgt_num_inds, &tgt_inds);

  printf("target node:\n");
  printf("- depth: %lu\n", tgt_depth);
  printf("- node index: %lu\n", tgt_node_index);
  printf("- no. of indices: %lu\n", tgt_num_inds);
  printf("- bbox: [%g, %g] x [%g, %g]\n",
         tgt_node->bbox.min[0], tgt_node->bbox.max[0],
         tgt_node->bbox.min[1], tgt_node->bbox.max[1]);
  printf("- indices:");
  if (tgt_num_inds > 10) {
    for (size_t i = 0; i < 10; ++i)
      printf(" %lu,", tgt_inds[i]);
    printf(" ...\n");
  } else {
    for (size_t i = 0; i < tgt_num_inds - 1; ++i)
      printf(" %lu,", tgt_inds[i]);
    printf(" %lu\n", tgt_inds[tgt_num_inds - 1]);
  }

  for (size_t i = 0; i < tgt_num_inds; ++i)
    assert(bfBbox2ContainsPoint(&tgt_node->bbox, tree->points[tgt_inds[i]]));

  /* Compute the groundtruth subblock of the kernel matrix induced by
   * the source and target nodes */

  BfMat mat;
  BfSize shape[] = {tgt_num_inds, src_num_inds};
  bfMakeEmptyMat(&mat, BF_DTYPE_COMPLEX, 2, shape);

  BfComplex *row;
  for (size_t i = 0; i < mat.shape[0]; ++i) {
    row = (BfComplex *)mat.data + i*mat.shape[1];
    for (size_t j = 0; j < mat.shape[1]; ++j) {
      row[j] = bfHelm2GetKernelValue(
        tree->points[tgt_inds[i]], tree->points[src_inds[j]], k);
    }
  }

  BfSize num_bytes;
  bfMatNumBytes(&mat, &num_bytes);

  printf("computed groundtruth subblock of kernel matrix:\n");
  printf("- rows: %lu\n", mat.shape[0]);
  printf("- columns: %lu\n", mat.shape[1]);
  printf("- size: %1.2f MB\n", ((double)num_bytes)/(1024*1024));

  /* Traverse leaves below source node at finest level */

  // TODO: what do we do when the leaves aren't all at the same
  // height?

  BfCircle2 tgt_circ = bfGetQuadtreeNodeBoundingCircle(tgt_node);

  enum BfError print(BfQuadtreeNode const *node, void *arg) {
    (void)arg;

    BfCircle2 src_leaf_circ = bfGetQuadtreeNodeBoundingCircle(node);

    BfReal rank_estimate = bfHelm2RankEstForTwoCircles(
      tgt_circ, src_leaf_circ, k, 1, 1e-15);

    printf("[%p] depth: %lu, rank: %lu\n",
           node, bfQuadtreeNodeDepth(node), rank_estimate);

    return BF_ERROR_NO_ERROR;
  }

  printf("src nodes:\n");
  bfMapQuadtreeNodeLeaves(src_node, print, NULL);

  /* Clean up */

  bfFreeMat(&mat);
}

int main(int argc, char const *argv[]) {
  if (argc != 2) {
    printf("usage: %s <points.bin>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  FILE *fp = fopen(argv[1], "r");

  // Get size of the file
  fseek(fp, 0, SEEK_END);
  size_t size = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  printf("size: %lu\n", size);

  // Calculate the number of points stored in points.bin
  size_t num_points = size/sizeof(double[2]);

  printf("num_points: %lu\n", num_points);

  // Allocate space for the points contained in points.bin and read
  // them in
  double (*points)[2] = malloc(size);
  fread(points, sizeof(double[2]), num_points, fp);

  printf("building quadtree...\n");

  BfQuadtree tree;
  bfInitQuadtreeFromPoints(&tree, num_points, points);

  BfReal k = 3000;

  bf_one_block(&tree, k);

  // Clean up
  free(points);
  fclose(fp);
}
