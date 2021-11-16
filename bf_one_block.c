#include <assert.h>
#include <stdio.h>

#include "quadtree.h"

static void bf_one_block(BfQuadtree const *tree) {
  enum BfError error;

  size_t src_depth = 3, src_node_index = 20;
  size_t tgt_depth = 3, tgt_node_index = 16*3 + 4*3 + 2;

  size_t src_num_inds, *src_inds;
  error = bfGetQuadtreeIndexRange(
    tree, src_depth, src_node_index, &src_num_inds, &src_inds);
  assert(error == BF_ERROR_NO_ERROR);

  printf("source node:\n");
  printf("- depth: %lu\n", src_depth);
  printf("- node index: %lu\n", src_node_index);
  printf("- no. of indices: %lu\n", src_num_inds);

  size_t tgt_num_inds, *tgt_inds;
  error = bfGetQuadtreeIndexRange(
    tree, tgt_depth, tgt_node_index, &tgt_num_inds, &tgt_inds);
  assert(error == BF_ERROR_NO_ERROR);

  printf("target node:\n");
  printf("- depth: %lu\n", tgt_depth);
  printf("- node index: %lu\n", tgt_node_index);
  printf("- no. of indices: %lu\n", tgt_num_inds);
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

  bf_one_block(&tree);

  // Clean up
  free(points);
  fclose(fp);
}
