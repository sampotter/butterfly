#include <bf/lbo.h>

#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/vectors.h>

void bfLboGetFemDiscretization(BfTrimesh const *trimesh, BfMatCsrReal **L,
                               BfMatCsrReal **M) {
  BEGIN_ERROR_HANDLING();

  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);

  BfSize nnz = numVerts + trimesh->vvOffset[numVerts];

  BfSize *rowptr = malloc((numVerts + 1)*sizeof(BfSize));
  if (rowptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfSize *colind = malloc(nnz*sizeof(BfSize));
  if (colind == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *L_data = calloc(nnz, sizeof(BfReal));
  if (L_data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BfReal *M_data = calloc(nnz, sizeof(BfReal));
  if (M_data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* Fill rowptr: same as vvOffset except we need to add more space
   * for the vertices on the diagonal */
  for (BfSize i = 0; i <= numVerts; ++i)
    rowptr[i] = trimesh->vvOffset[i] + i;

  /* Fill colind */
  for (BfSize i = 0; i < numVerts; ++i) {
    BfSize j = rowptr[i];
    BfSize jend = rowptr[i + 1] - 1;
    BfSize *vv = &trimesh->vv[trimesh->vvOffset[i]];
    while (j < jend && *vv < i)
      colind[j++] = *vv++;
    colind[j++] = i;
    while (j < jend)
      colind[j++] = *vv++;
  }

  BfSize f, j, i, i0, i1, *jptr, k, k0, k1;
  BfPoint3 x, x0, x1, y, y0, y1;
  BfVector3 d, d0, d1, g, g0, g1, n;
  BfReal A;

  /* For each vertex: */
  for (i = 0; i < numVerts; ++i) {
    jptr = &colind[rowptr[i]];

    /* Get current vertex */
    bfTrimeshGetVertex(trimesh, i, x);

    /* Find position of `i` in current block of column indices */
    k = 0;
    while (jptr[k] != i) ++k;

    /* For each incident face: */
    for (j = trimesh->vfOffset[i]; j < trimesh->vfOffset[i + 1]; ++j) {
      f = trimesh->vf[j];

      /* Get the other two vertices in the current face */
      bfTrimeshGetOpFaceVerts(trimesh, f, i, &i0, &i1);
      bfTrimeshGetVertex(trimesh, i0, x0);
      bfTrimeshGetVertex(trimesh, i1, x1);

      /* Find their positions in the current block of column indices */
      k0 = k1 = 0;
      while (jptr[k0] != i0) ++k0;
      while (jptr[k1] != i1) ++k1;

      /* Opposite edge vectors */
      bfPoint3Sub(x1, x0, d);
      bfPoint3Sub(x, x1, d0);
      bfPoint3Sub(x0, x, d1);

      /* Compute orthogonal projection of each vertex onto the
       * opposite face edge */
      bfPoint3GetPointOnRay(x0, d, -bfVector3Dot(d, d1)/bfVector3Dot(d, d), y);
      bfPoint3GetPointOnRay(x1, d0, -bfVector3Dot(d0, d)/bfVector3Dot(d0, d0), y0);
      bfPoint3GetPointOnRay(x, d1, -bfVector3Dot(d1, d0)/bfVector3Dot(d1, d1), y1);

      /* Compute gradients for hat functions centered at each vertex */
      bfPoint3Sub(x, y, g);
      bfPoint3Sub(x0, y0, g0);
      bfPoint3Sub(x1, y1, g1);
      bfVector3Scale(g, 1/bfVector3Dot(g, g));
      bfVector3Scale(g0, 1/bfVector3Dot(g0, g0));
      bfVector3Scale(g1, 1/bfVector3Dot(g1, g1));

      /* Get triangle area */
      bfVector3Cross(d0, d1, n);
      A = bfVector3Norm(n)/2;

      /* Update L values */
      L_data[rowptr[i] + k] += A*bfVector3Dot(g, g);
      L_data[rowptr[i] + k0] += A*bfVector3Dot(g, g0);
      L_data[rowptr[i] + k1] += A*bfVector3Dot(g, g1);

      /* Update M values */
      M_data[rowptr[i] + k] += A/6;
      M_data[rowptr[i] + k0] += A/12;
      M_data[rowptr[i] + k1] += A/12;
    }
  }

  *L = bfMatCsrRealNew();
  HANDLE_ERROR();

  bfMatCsrRealInit(*L, numVerts, numVerts, rowptr, colind, L_data);
  HANDLE_ERROR();

  *M = bfMatCsrRealNew();
  HANDLE_ERROR();

  bfMatCsrRealInit(*M, numVerts, numVerts, rowptr, colind, M_data);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfMatCsrRealDeinitAndDealloc(L);
    bfMatCsrRealDeinitAndDealloc(M);
  }

  free(rowptr);
  free(colind);
  free(L_data);
  free(M_data);
}
