#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "block_mat.h"
#include "error.h"
#include "helm2.h"
#include "mat.h"
#include "rand.h"
#include "quadtree.h"

static BfReal const PINV_ATOL = 1e-13;
static BfReal const PINV_RTOL = 1e-13;

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

//   enum BfError printNode(BfQuadtreeNode const *node, void *arg) {
//     (void)arg;

//     BfSize depth = bfQuadtreeNodeDepth(node);

//     if (node->flags & BF_QUADTREE_NODE_FLAG_ROOT) {
//       printf("%lu root\n", depth);
//     } else if (node->flags & BF_QUADTREE_NODE_FLAG_CHILD_0) {
//       printf("%lu 0\n", depth);
//     } else if (node->flags & BF_QUADTREE_NODE_FLAG_CHILD_1) {
//       printf("%lu 1\n", depth);
//     } else if (node->flags & BF_QUADTREE_NODE_FLAG_CHILD_2) {
//       printf("%lu 2\n", depth);
//     } else if (node->flags & BF_QUADTREE_NODE_FLAG_CHILD_3) {
//       printf("%lu 3\n", depth);
//     }

//     return BF_ERROR_NO_ERROR;
//   }

//   bfMapQuadtreeNodes(src_node,BF_TREE_TRAVERSAL_LR_LEVEL_ORDER,printNode,NULL);

  /* Try level-by-level traversal */

//   BfQuadtreeLevelIter iter;
//   bfInitQuadtreeLevelIter(&iter, BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER, src_node);
//   while (!bfQuadtreeLevelIterIsDone(&iter)) {
//     bfMapPtrArray(&iter.level_nodes, (BfPtrFunc)printNode, NULL);
//     bfQuadtreeLevelIterNext(&iter);
//   }

  /* Figure out the common maximum depth */

  BfSize src_max_depth = bfGetMaxDepthBelowQuadtreeNode(src_node);
  BfSize tgt_max_depth = bfGetMaxDepthBelowQuadtreeNode(tgt_node);

  printf("src_max_depth: %lu\n", src_max_depth);
  printf("tgt_max_depth: %lu\n", tgt_max_depth);

  /* Do a bit of clean up before moving on */

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

  /* Let's try manually doing a couple levels */

  BfQuadtreeLevelIter src_level_iter;
  bfInitQuadtreeLevelIter(
    &src_level_iter, BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER, src_node);

  BfQuadtreeLevelIter tgt_level_iter;
  bfInitQuadtreeLevelIter(
    &tgt_level_iter, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, tgt_node);

  // skip a few levels on the source side
  for (BfSize i = 0; i < 6; ++i)
    bfQuadtreeLevelIterNext(&src_level_iter);

  BfSize current_src_depth;
  bfQuadtreeLevelIterCurrentDepth(&src_level_iter, &current_src_depth);
  printf("starting source depth: %lu\n", current_src_depth);

  BfSize current_tgt_depth;
  bfQuadtreeLevelIterCurrentDepth(&tgt_level_iter, &current_tgt_depth);
  printf("starting target depth: %lu\n", current_tgt_depth);

  BfQuadtreeNode *current_src_node, *current_tgt_node;
  BfCircle2 current_src_circ, current_tgt_circ;

  BfMat src_node_pts;

  BfSize Z1_num_blocks = bfPtrArraySize(&src_level_iter.level_nodes);

  BfSize *Z1_row_block_sizes = malloc(Z1_num_blocks*sizeof(BfSize));
  BfSize *Z1_col_block_sizes = malloc(Z1_num_blocks*sizeof(BfSize));

  /* get the current target node and its bounding circle (there's only
   * one target node on the first level/for the first stage) */
  bfPtrArrayGetFirst(&tgt_level_iter.level_nodes, (BfPtr *)&current_tgt_node);
  current_tgt_circ = bfGetQuadtreeNodeBoundingCircle(current_tgt_node);

  printf("initial target node: %p (%lu points)\n",
         current_tgt_node, bfQuadtreeNodeNumPoints(current_tgt_node));

  /* iterate over each source block and use the a priori rank estimate
   * to determine the number of equivalent sources to place on each
   * pair of circles */
  for (BfSize i = 0; i < Z1_num_blocks; ++i) {
    bfPtrArrayGet(&src_level_iter.level_nodes, i, (BfPtr *)&current_src_node);

    current_src_circ = bfGetQuadtreeNodeBoundingCircle(current_src_node);

    bfHelm2RankEstForTwoCircles(
      current_src_circ, current_tgt_circ, k, 1, 1e-15, &p_hat);

    Z1_row_block_sizes[i] = (BfSize)ceil(p_hat);
    Z1_col_block_sizes[i] = bfQuadtreeNodeNumPoints(current_src_node);

    printf("* i = %lu: %p (block: %lu x %lu)\n",
           i, current_src_node, Z1_row_block_sizes[i], Z1_col_block_sizes[i]);
  }

  /* allocate block diagonal matrix which will hold the kernel
   * matrices mapping the charges at the original points contained in
   * each node to potentials at the target circle points */
  BfBlockMat Z1_block_mat;
  bfInitEmptyBlockMat(&Z1_block_mat,
                      BF_DTYPE_COMPLEX, BF_BLOCK_MAT_PROP_DIAGONAL,
                      Z1_num_blocks, Z1_row_block_sizes,
                      Z1_num_blocks, Z1_col_block_sizes);

  for (BfSize i = 0; i < Z1_num_blocks; ++i) {
    bfPtrArrayGet(&src_level_iter.level_nodes, i, (BfPtr *)&current_src_node);

    /* sample points on the target circle */
    bfInitEmptyMat(&tgt_circ_pts, BF_DTYPE_REAL, BF_MAT_PROP_NONE,
                   (BfSize[]) {Z1_row_block_sizes[i], 2});
    bfSamplePointsOnCircle2(&current_tgt_circ, &tgt_circ_pts);

    /* get the original points in the source box */
    bfInitEmptyMat(&src_node_pts, BF_DTYPE_REAL, BF_MAT_PROP_NONE,
                   (BfSize[]) {Z1_col_block_sizes[i], 2});
    bfGetQuadtreeNodePoints(current_src_node, &src_node_pts);

    /* set the current block to the kernel matrix */
    BfMat *Z1_block = bfGetBlock(&Z1_block_mat, i, i);
    bfHelm2KernelMatrixFromPoints(Z1_block, &src_node_pts, &tgt_circ_pts, k);

    bfFreeMat(&tgt_circ_pts); // TODO: better resize instead of re-allocating
    bfFreeMat(&src_node_pts); // and re-freeing over and over again
  }

  /* probably makes more sense to build the array of matrices and then
   * construct the block diagonal matrix as a view of them... */

  /* also, we can infer the row and col block sizes from the matrices
   * stored in the BfBlockMat themselves, so no need to additionally
   * hold that pointer */

  /* Z1 and Z2 have the same number of blocks */
  BfSize Z2_num_blocks = Z1_num_blocks;
  BfMat *Z2_blocks = malloc(Z1_num_blocks*sizeof(BfMat));
  BfMat kernel;

  for (BfSize i = 0; i < Z2_num_blocks; ++i) {
    bfPtrArrayGet(&src_level_iter.level_nodes, i, (BfPtr *)&current_src_node);

    p = Z1_row_block_sizes[i];
    BfSize pts_shape[] = {p, 2};

    /* sample points on the target circle */
    bfInitEmptyMat(&tgt_circ_pts, BF_DTYPE_REAL, BF_MAT_PROP_NONE, pts_shape);
    bfSamplePointsOnCircle2(&current_tgt_circ, &tgt_circ_pts);

    /* sample points on the source circle */
    bfInitEmptyMat(&src_circ_pts, BF_DTYPE_REAL, BF_MAT_PROP_NONE, pts_shape);
    bfSamplePointsOnCircle2(&current_src_circ, &src_circ_pts);

    /* compute the kernel matrix for this set of interactions */
    bfInitEmptyMat(&kernel, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE,
                            (BfSize[]) {p, p});
    bfHelm2KernelMatrixFromPoints(&kernel, &src_circ_pts, &tgt_circ_pts, k);

    bfFreeMat(&tgt_circ_pts);
    bfFreeMat(&src_circ_pts);

    /* compute its pseudoinverse and store it in the current block */
    bfComputePinv(&kernel, PINV_ATOL, PINV_RTOL, &Z2_blocks[i]);

    bfFreeMat(&kernel);
  }

  /* next, we want to create the first layer of butterfly matrices */

  BfSize getChildren(BfQuadtreeNode const *node, BfQuadtreeNode const *ch[4]) {
    for (BfSize i = 0; i < 4; ++i)
      ch[i] = NULL;
    BfSize num_children = 0;
    for (BfSize i = 0; i < 4; ++i)
      if (node->child[i] != NULL)
        ch[num_children++] = node->child[i];
    return num_children;
  }

  void makeBfacLayer() {
    printf("makeBfacLayer()\n");

    bfQuadtreeLevelIterNext(&src_level_iter); // go up a level on the source tree

    BfSize depth;

    bfQuadtreeLevelIterCurrentDepth(&src_level_iter, &depth);
    printf("* src_level_iter depth: %lu\n", depth);

    bfQuadtreeLevelIterCurrentDepth(&tgt_level_iter, &depth);
    printf("* tgt_level_iter depth: %lu\n", depth);

    BfSize num_tgt_ch, num_src_ch;
    BfQuadtreeNode const *tgt_child[4], *src_child[4];
    BfCircle2 tgt_child_circ[4], src_child_circ[4];

    for (BfSize _ = 0; _ < bfPtrArraySize(&tgt_level_iter.level_nodes); ++_) {
      bfPtrArrayGet(&tgt_level_iter.level_nodes, _, (BfPtr *)&current_tgt_node);

      printf("  - current_tgt_node #%lu (%p)\n", _, current_tgt_node);

      num_tgt_ch = getChildren(current_tgt_node, tgt_child);
      for (BfSize i = 0; i < num_tgt_ch; ++i)
        tgt_child_circ[i] = bfGetQuadtreeNodeBoundingCircle(tgt_child[i]);

      for (BfSize i = 0; i < bfPtrArraySize(&src_level_iter.level_nodes); ++i) {
        bfPtrArrayGet(&src_level_iter.level_nodes, i, (BfPtr *)&current_src_node);

        printf("    - current_src_node #%lu (%p)\n", i, current_src_node);

        /* get the current source circle---we should have this handy
         * already so that we don't need to recompute or get the points
         * again */
        current_src_circ = bfGetQuadtreeNodeBoundingCircle(current_src_node);

        num_src_ch = getChildren(current_src_node, src_child);
        for (BfSize j = 0; j < num_src_ch; ++j)
          src_child_circ[j] = bfGetQuadtreeNodeBoundingCircle(src_child[j]);

        for (BfSize j = 0; j < num_tgt_ch; ++j) {
          bfHelm2RankEstForTwoCircles(
            current_src_circ, tgt_child_circ[j], k, 1, 1e-15, &p_hat);

          p = (BfSize)ceil(p_hat);

          for (BfSize l = 0; l < num_src_ch; ++l) {

            printf("      - src_child[%lu] -> tgt_child[%lu]\n", l, j);

            BfSize q = Z1_row_block_sizes[i];

            BfMat Y_child;
            bfInitEmptyMat(&Y_child, BF_DTYPE_REAL, BF_MAT_PROP_NONE, (BfSize[]) {p, 2});
            bfSamplePointsOnCircle2(&tgt_child_circ[j], &Y_child);

            BfMat X_child;
            bfInitEmptyMat(&X_child, BF_DTYPE_REAL, BF_MAT_PROP_NONE, (BfSize[]) {q, 2});
            bfSamplePointsOnCircle2(&src_child_circ[l], &X_child);

            BfMat Z_eval;
            bfInitEmptyMat(&Z_eval, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, (BfSize[]) {q, p});
            bfHelm2KernelMatrixFromPoints(&Z_eval, &X_child, &Y_child, k);

            BfMat X;
            bfInitEmptyMat(&X, BF_DTYPE_REAL, BF_MAT_PROP_NONE, (BfSize[]) {p, 2});
            bfSamplePointsOnCircle2(&current_src_circ, &X);

            BfMat Z_eq;
            bfInitEmptyMat(&Z_eq, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE, (BfSize[]) {p, p});
            bfHelm2KernelMatrixFromPoints(&Z_eq, &X, &Y_child, k);

            BfMat Z_eq_pinv;
            bfInitEmptyMat(&Z_eq_pinv, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE,
                           (BfSize[]) {p, p});
            bfComputePinv(&Z_eq, PINV_ATOL, PINV_RTOL, &Z_eq_pinv);

            // set block to Z_eq\Z_eval

            BfMat block;
            bfInitEmptyMat(&block, BF_DTYPE_COMPLEX, BF_MAT_PROP_NONE,
                           (BfSize[]) {p, q});
            bfMatMul(&Z_eq_pinv, &Z_eq, &block);

            bfFreeMat(&Y_child);
            bfFreeMat(&X_child);
            bfFreeMat(&Z_eval);
            bfFreeMat(&X);
            bfFreeMat(&Z_eq);
            bfFreeMat(&Z_eq_pinv);
            bfFreeMat(&block);
          }
        }
      }
    }

    bfQuadtreeLevelIterNext(&tgt_level_iter);
  }

  makeBfacLayer();
  makeBfacLayer();
  makeBfacLayer();
  makeBfacLayer();

  /* clean up */

  free(Z2_blocks);
  free(Z1_row_block_sizes);
  free(Z1_col_block_sizes);

  // TODO: next we'll iterate over the edges connecting each level

  // TODO: finally, we'll set up the block diagonal matrix used to evaluate

  /* Clean up */

  bfFreeMat(&tgt_pts);
  bfFreeMat(&src_pts);
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
