#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/fac_span.h>
#include <bf/fac_streamer.h>
#include <bf/interval_tree.h>
#include <bf/lbo.h>
#include <bf/linalg.h>
#include <bf/logging.h>
#include <bf/octree.h>
#include <bf/rand.h>
#include <bf/util.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "argtable3.h"

#include "lbo.h"

int const MAX_NUM_ARG_ERRORS = 20;

typedef struct {
  BfSize seed;
  BfLogLevel logLevel;

  char const *objPath;

  bool useOctree;
  bool useFiedlerTree;

  BfReal tol;
  BfSize numLeafNodes;
  BfSize rowTreeOffset;
  BfSize freqTreeOffset;

  bool compareRelativeErrors;
} Opts;

bool parseArgs(Opts *opts, int argc, char *argv[]) {
  bool success = true;

  struct arg_lit *help;

  struct arg_int *seed;
  struct arg_str *logLevel;

  struct arg_str *objPath;

  struct arg_lit *useOctree;
  struct arg_lit *useFiedlerTree;

  struct arg_dbl *tol;
  struct arg_int *numLeafNodes;
  struct arg_int *rowTreeOffset;
  struct arg_int *freqTreeOffset;

  struct arg_lit *compareRelativeErrors;

  struct arg_end *end;

  void *argtable[] = {
    help = arg_lit0(NULL, "help", "Display help and exit"),

    seed = arg_int0(NULL, "seed", NULL, "Seed for random number generator (default: 0)"),
    logLevel = arg_str0(NULL, "logLevel", NULL, NULL),

    objPath = arg_str0(NULL, "objPath", NULL, NULL),

    useOctree = arg_lit0(NULL, "useOctree", NULL),
    useFiedlerTree = arg_lit0(NULL, "useFiedlerTree", NULL),

    tol = arg_dbl0(NULL, "tol", NULL, NULL),
    numLeafNodes = arg_int0(NULL, "numLeafNodes", NULL, NULL),
    rowTreeOffset = arg_int0(NULL, "rowTreeOffset", NULL, NULL),
    freqTreeOffset = arg_int0(NULL, "freqTreeOffset", NULL, NULL),

    compareRelativeErrors = arg_lit0(NULL, "compareRelativeErrors", NULL),

    end = arg_end(MAX_NUM_ARG_ERRORS)
  };

  *seed->ival = 0;
  *logLevel->sval = "error";

  useOctree->count = 1;
  useFiedlerTree->count = 0;

  *tol->dval = 1e-3;
  *numLeafNodes->ival = -1;
  *rowTreeOffset->ival = 0;
  *freqTreeOffset->ival = 2;

  compareRelativeErrors->count = 0;

  BfSize numErrors = arg_parse(argc, argv, argtable);

  if (help->count > 0) {
    printf("Usage: %s\n", argv[0]);
    arg_print_syntax(stdout, argtable, "\n");
    printf("Test driver for butterfly compression of LBO\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");
    success = false;
    goto cleanup;
  }

  if (numErrors > 0) {
    arg_print_errors(stdout, end, argv[0]);
    printf("Try '%s --help' for more information.\n", argv[0]);
    success = false;
    goto cleanup;
  }

  opts->seed = *seed->ival;

  if (!strcmp(*logLevel->sval, "todo")) {
    opts->logLevel = BF_LOG_LEVEL_TODO;
  } else if (!strcmp(*logLevel->sval, "debug")) {
    opts->logLevel = BF_LOG_LEVEL_DEBUG;
  } else if (!strcmp(*logLevel->sval, "info")) {
    opts->logLevel = BF_LOG_LEVEL_INFO;
  } else if (!strcmp(*logLevel->sval, "warn")) {
    opts->logLevel = BF_LOG_LEVEL_WARN;
  } else if (!strcmp(*logLevel->sval, "error")) {
    opts->logLevel = BF_LOG_LEVEL_ERROR;
  } else {
    printf("--logLevel must be one of: \"todo\", \"debug\", \"info\", \"warn\", \"error\"\n");
    success = false;
    goto cleanup;
  }

  opts->objPath = *objPath->sval;

  opts->useOctree = useOctree->count > 0;
  opts->useFiedlerTree = useFiedlerTree->count > 0;

  opts->tol = *tol->dval;
  opts->numLeafNodes = *numLeafNodes->ival;
  opts->rowTreeOffset = *rowTreeOffset->ival;
  opts->freqTreeOffset = *freqTreeOffset->ival;

  opts->compareRelativeErrors = compareRelativeErrors->count > 0;

cleanup:
  arg_freetable(argtable, sizeof(argtable)/sizeof(argtable[0]));

  return success;
}

int main(int argc, char *argv[]) {
  Opts opts;
  if (!parseArgs(&opts, argc, argv))
    return EXIT_SUCCESS;

  bfSeed(opts.seed);
  bfSetLogLevel(opts.logLevel);

  BF_ERROR_BEGIN() {}

  BfTrimesh trimesh;
  bfTrimeshInitFromObjFile(&trimesh, opts.objPath);
  HANDLE_ERROR();

  BfSize numVerts = bfTrimeshGetNumVerts(&trimesh);

  printf("loaded triangle mesh from %s (%lu verts)\n", opts.objPath, numVerts);

  // TODO: implement
//   /* Build fiedler tree */
//   BfFiedlerTree fiedlerTree;
//   bfFiedlerTreeInitFromTrimesh(&fiedlerTree, trimesh, &L, &M);
//   BfTree *rowTree = bfFiedlerTreeToTree(&fiedlerTree);

  BfOctree octree;
  bfOctreeInit(&octree, &trimesh.verts, NULL);
  HANDLE_ERROR();

  char const *octreeBoxesPath = "octree_boxes.txt";
  bfOctreeSaveBoxesToTextFile(&octree, octreeBoxesPath);
  HANDLE_ERROR();
  printf("wrote octree cells to %s\n", octreeBoxesPath);

  /* Upcast the spatial tree to get the row tree */
  BfTree *rowTree = bfOctreeToTree(&octree);

  BfSize rowTreeDepth = bfTreeGetMaxDepth(rowTree);
  printf("row tree has depth %lu\n", rowTreeDepth);

  /* Compute a finite element discretization of the Laplace-Beltrami
   * operator on `trimesh` using linear finite elements. The stiffness
   * matrix is returned in L and the mass matrix is returned in M. The
   * mass matrix isn't diagonal but this isn't too important. */
  BfMat *L, *M;
  bfLboGetFemDiscretization(&trimesh, &L, &M);
  HANDLE_ERROR();

  /* Find the largest eigenvalue. We need this to determine the
   * interval on which we'll build the frequency tree. */
  BfReal lamMax = bfGetMaxEigenvalue(L, M);
  HANDLE_ERROR();
  printf("lambda_max = %g\n", lamMax);

  /* The natural frequency of each eigenvector is the square root of
   * the associated eigenvalue. */
  BfReal freqMax = sqrt(lamMax);

  BfSize freqTreeDepth = rowTreeDepth - opts.freqTreeOffset;
  printf("building frequency tree with depth %lu\n", freqTreeDepth);

  /* Set up the frequency tree. Note: we build the tree on the
   * frequency scale as opposed to the eigenvalue scale to preserve
   * the time-frequency product in the butterfly factorization. */
  BfIntervalTree *freqTree = bfIntervalTreeNew();
  bfIntervalTreeInitEmpty(freqTree, 0, freqMax, 2, freqTreeDepth);
  HANDLE_ERROR();

  /* Upcast frequency tree to get the column tree */
  BfTree *colTree = bfIntervalTreeToTree(freqTree);

  BfPoints1 *freqs = bfPoints1New();
  HANDLE_ERROR();

  bfPoints1InitEmpty(freqs, BF_ARRAY_DEFAULT_CAPACITY);
  HANDLE_ERROR();

  BfFacSpec spec = {
    .rowTree = rowTree,
    .colTree = colTree,
    .rowTreeInitDepth = opts.rowTreeOffset,
    .colTreeInitDepth = freqTreeDepth, // TODO: this is unused!
    .tol = opts.tol,
    .minNumRows = 20,
    .minNumCols = 1,
    .compareRelativeErrors = opts.compareRelativeErrors,
  };

  /* Set up the depth-first butterfly factorization streamer. We'll
   * use this below to construct the butterfly factorization
   * incrementally. */
  BfFacStreamer *facStreamer = bfFacStreamerNew();
  bfFacStreamerInit(facStreamer, &spec);
  HANDLE_ERROR();

  if (opts.numLeafNodes != BF_SIZE_BAD_VALUE)
    printf("streaming %lu leaf nodes\n", opts.numLeafNodes);
  else
    printf("streaming full eigenvector matrix\n");

  bfToc();

  BfSize numStreamed = 0;
  while (!bfFacStreamerIsDone(facStreamer)) {
    feedFacStreamerNextEigenband(facStreamer, freqs, L, M);
    HANDLE_ERROR();

    if (++numStreamed >= opts.numLeafNodes) break;
  }

  printf("finished streaming butterfly factorization [%.2fs]\n", bfToc());

  BfFacSpan *facSpan = bfFacStreamerGetFacSpan(facStreamer);
  BfMat *mat = bfFacSpanGetMat(facSpan);

  BfSize numEigs = bfMatGetNumCols(mat);
  printf("- streamed %lu eigs (%1.1f%% of total)\n", numEigs, (100.0*numEigs)/numVerts);

  BfReal numBytesUncompressed = sizeof(BfReal)*numVerts*numEigs;
  BfReal numBytesCompressed = bfMatNumBytes(mat);

  printf("- compressed size: %.1f MB\n", numBytesCompressed/pow(1024, 2));
  printf("- uncompressed size: %.1f MB\n", numBytesUncompressed/pow(1024, 2));
  printf("- compression rate: %.1f\n", numBytesUncompressed/numBytesCompressed);

  END_ERROR_HANDLING() {}

  /* Clean up */
  bfTreeDelete(&rowTree);
  bfOctreeDeinit(&octree);
  bfMatDelete(&M);
  bfMatDelete(&L);
  bfTrimeshDeinit(&trimesh);
}
