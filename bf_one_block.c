#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "error.h"
#include "helm2.h"
#include "mat.h"
#include "rand.h"
#include "quadtree.h"

static void bf_one_block(BfQuadtree const *tree, double k) {
  enum BfError error = BF_ERROR_NO_ERROR;

  /* Seed the random number generators */

  bfSeed(1234);

  /* Get source and target nodes from quadtree and check that their
   * indices are OK */

  size_t src_depth = 3, src_node_index = 20, src_num_inds, *src_inds;
  size_t tgt_depth = 3, tgt_node_index = 16*3 + 4*3 + 2, tgt_num_inds, *tgt_inds;

  BfQuadtreeNode *src_node;
  bfGetQuadtreeNode(tree, src_depth, src_node_index, &src_node);
  bfGetQuadtreeNodeIndices(src_node, &src_num_inds, &src_inds);

  for (size_t i = 0; i < src_num_inds; ++i)
    assert(bfBbox2ContainsPoint(&src_node->bbox, tree->points[src_inds[i]]));

  BfQuadtreeNode *tgt_node;
  bfGetQuadtreeNode(tree, tgt_depth, tgt_node_index, &tgt_node);
  bfGetQuadtreeNodeIndices(tgt_node, &tgt_num_inds, &tgt_inds);

  for (size_t i = 0; i < tgt_num_inds; ++i)
    assert(bfBbox2ContainsPoint(&tgt_node->bbox, tree->points[tgt_inds[i]]));

  /* Compute the groundtruth subblock of the kernel matrix induced by
   * the source and target nodes */

  BfSize m = bfQuadtreeNodeNumPoints(tgt_node);
  BfSize n = bfQuadtreeNodeNumPoints(src_node);

  BfMat tgt_pts;
  bfInitEmptyMat(&tgt_pts, BF_DTYPE_REAL, BF_MAT_PROP_NONE, (BfSize[]) {m, 2});
  bfGetQuadtreeNodePoints(tgt_node, &tgt_pts);

  bfSaveMat(&tgt_pts, "tgt_pts.bin");
  puts("wrote tgt_pts.bin");

  BfMat src_pts;
  bfInitEmptyMat(&src_pts, BF_DTYPE_REAL, BF_MAT_PROP_NONE, (BfSize[]) {n, 2});
  bfGetQuadtreeNodePoints(src_node, &src_pts);

  bfSaveMat(&src_pts, "src_pts.bin");
  puts("wrote src_pts.bin");

  BfMat Z_gt;
  bfInitEmptyMat(&Z_gt, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, (BfSize[]) {m, n});
  bfHelm2KernelMatrixFromPoints(&Z_gt, &src_pts, &tgt_pts, k);

  BfSize num_bytes;
  bfMatNumBytes(&Z_gt, &num_bytes);

  printf("computed groundtruth subblock of kernel matrix:\n");
  printf("- rows: %lu\n", Z_gt.shape[0]);
  printf("- columns: %lu\n", Z_gt.shape[1]);
  printf("- size: %1.2f MB\n", ((double)num_bytes)/(1024*1024));

  bfSaveMat(&Z_gt, "Z_gt.bin");
  printf("wrote Z_gt.bin\n");

  BfMat U_gt, S_gt, Vt_gt;
  bfInitEmptySvdMats(&Z_gt, &U_gt, &S_gt, &Vt_gt);
  bfComputeMatSvd(&Z_gt, &U_gt, &S_gt, &Vt_gt);
  printf("computed SVD of groundtruth subblock\n");

  bfSaveMat(&U_gt, "U_gt.bin");
  bfSaveMat(&S_gt, "S_gt.bin");
  bfSaveMat(&Vt_gt, "Vt_gt.bin");
  printf("wrote U_gt.bin, S_gt.bin, and Vt_gt.bin\n");

  /* Do a single stage of the bf fac and check the error */

  BfCircle2 src_circ = bfGetQuadtreeNodeBoundingCircle(src_node);
  BfCircle2 tgt_circ = bfGetQuadtreeNodeBoundingCircle(tgt_node);

  BfReal p_hat;
  error = bfHelm2RankEstForTwoCircles(src_circ, tgt_circ, k, 1, 1e-15, &p_hat);
  assert(!error);

  BfSize p = ceil(p_hat);

  printf("rank estimate for single stage: p = %lu\n", p);

  BfMat src_circ_pts;
  bfInitEmptyMat(&src_circ_pts,BF_DTYPE_REAL,BF_MAT_PROP_NONE,(BfSize[]){p,2});
  bfSamplePointsOnCircle2(&src_circ, &src_circ_pts);

  bfSaveMat(&src_circ_pts, "src_circ_pts.bin");
  puts("wrote src_circ_pts.bin");

  BfMat tgt_circ_pts;
  bfInitEmptyMat(&tgt_circ_pts,BF_DTYPE_REAL,BF_MAT_PROP_NONE,(BfSize[]){p,2});
  bfSamplePointsOnCircle2(&tgt_circ, &tgt_circ_pts);

  bfSaveMat(&tgt_circ_pts, "tgt_circ_pts.bin");
  puts("wrote tgt_circ_pts.bin");

  printf("computed Z1\n");

  BfMat Z1;
  bfInitEmptyMat(&Z1, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, (BfSize[]) {p, n});
  error = bfHelm2KernelMatrixFromPoints(&Z1, &src_pts, &tgt_circ_pts, k);
  assert(!error);

  printf("computed Z2\n");

  BfMat Z2;
  bfInitEmptyMat(&Z2, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, (BfSize[]) {p, p});
  error = bfHelm2KernelMatrixFromPoints(&Z2, &src_circ_pts, &tgt_circ_pts, k);
  assert(!error);

  printf("computed Z3\n");

  BfMat Z3;
  bfInitEmptyMat(&Z3, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, (BfSize[]) {m, p});
  bfHelm2KernelMatrixFromPoints(&Z3, &src_circ_pts, &tgt_pts, k);

  bfSaveMat(&Z1, "Z1.bin");
  bfSaveMat(&Z2, "Z2.bin");
  bfSaveMat(&Z3, "Z3.bin");

  printf("wrote Z1.bin, Z2.bin, and Z3.bin\n");

  // compute the SVD of Z2:

  BfMat U2, S2, Vt2;
  bfInitEmptySvdMats(&Z2, &U2, &S2, &Vt2);
  bfComputeMatSvd(&Z2, &U2, &S2, &Vt2);

  bfSaveMat(&U2, "U2.bin");
  bfSaveMat(&S2, "S2.bin");
  bfSaveMat(&Vt2, "Vt2.bin");
  puts("wrote U2.bin, S2.bin, and Vt2.bin");

  // check that b_gt = Z_gt*x and b = Z3*(Z2\(Z1*x) are close:

  BfSize num_trials = 10;

  BfMat x;
  bfInitEmptyMat(&x,BF_DTYPE_COMPLEX,BF_MAT_PROP_NONE,(BfSize[]){n,num_trials});
  bfFillMatRandn(&x);

  bfSaveMat(&x, "x.bin");
  puts("wrote x.bin");

  BfMat b_gt;
  bfInitEmptyMat(&b_gt,BF_DTYPE_COMPLEX,BF_MAT_PROP_NONE,(BfSize[]){m,num_trials});
  bfMatMul(&Z_gt, &x, &b_gt);

  bfSaveMat(&b_gt, "b_gt.bin");

  puts("wrote b_gt.bin");

  BfMat b1, b2, b;
  bfInitEmptyMat(&b1,BF_DTYPE_COMPLEX,BF_MAT_PROP_NONE,(BfSize[]){p,num_trials});
  bfInitEmptyMat(&b2,BF_DTYPE_COMPLEX,BF_MAT_PROP_NONE,(BfSize[]){p,num_trials});
  bfInitEmptyMat(&b,BF_DTYPE_COMPLEX,BF_MAT_PROP_NONE,(BfSize[]){m,num_trials});

  error |= bfMatMul(&Z1, &x, &b1); // b1 = Z1*x
  error |= bfMatSolve(&U2, &b1, &b2); // b2 = U2\b1
  error |= bfMatSolve(&S2, &b2, &b1); // b1 = S2\b2
  error |= bfMatSolve(&Vt2, &b1, &b2); // b2 = Vt2\b1
  error |= bfMatMul(&Z3, &b2, &b); // b = Z3*b2
  assert(!error);

  bfSaveMat(&b, "b.bin");
  puts("wrote b.bin");

  /* Check that LR level order works correctly */

  enum BfError printNode(BfQuadtreeNode const *node, void *arg) {
    (void)arg;

    BfSize depth = bfQuadtreeNodeDepth(node);

    if (node->flags & BF_QUADTREE_NODE_FLAG_ROOT) {
      printf("%lu root\n", depth);
    } else if (node->flags & BF_QUADTREE_NODE_FLAG_CHILD_0) {
      printf("%lu 0\n", depth);
    } else if (node->flags & BF_QUADTREE_NODE_FLAG_CHILD_1) {
      printf("%lu 1\n", depth);
    } else if (node->flags & BF_QUADTREE_NODE_FLAG_CHILD_2) {
      printf("%lu 2\n", depth);
    } else if (node->flags & BF_QUADTREE_NODE_FLAG_CHILD_3) {
      printf("%lu 3\n", depth);
    }

    return BF_ERROR_NO_ERROR;
  }

  bfMapQuadtreeNodes(src_node,BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER,printNode,NULL);

  /* Figure out the common maximum depth */

  BfSize src_max_depth = 0, tgt_max_depth = 0;

  enum BfError findMaxDepth(BfQuadtreeNode const *node, void *arg) {
    BfSize *max_depth = arg;

    BfSize depth = bfQuadtreeNodeDepth(node);

    if (depth > *max_depth)
      *max_depth = depth;

    return BF_ERROR_NO_ERROR;
  }

  bfMapQuadtreeNodes(src_node, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, findMaxDepth, &src_max_depth);
  bfMapQuadtreeNodes(src_node, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, findMaxDepth, &tgt_max_depth);

  printf("src_max_depth: %lu\n", src_max_depth);
  printf("tgt_max_depth: %lu\n", tgt_max_depth);

  /* Clean up */

  bfFreeMat(&tgt_pts);
  bfFreeMat(&src_pts);
  bfFreeMat(&Z_gt);
  bfFreeMat(&U_gt);
  bfFreeMat(&S_gt);
  bfFreeMat(&Vt_gt);
  bfFreeMat(&src_circ_pts);
  bfFreeMat(&tgt_circ_pts);
  bfFreeMat(&Z1);
  bfFreeMat(&Z2);
  bfFreeMat(&Z3);
  bfFreeMat(&U2);
  bfFreeMat(&S2);
  bfFreeMat(&Vt2);
  bfFreeMat(&x);
  bfFreeMat(&b);
  bfFreeMat(&b1);
  bfFreeMat(&b2);
  bfFreeMat(&b_gt);
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
