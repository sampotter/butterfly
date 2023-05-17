#include <bf/lbo.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/vectors.h>

void bfLboGetFemDiscretization(BfTrimesh const *trimesh, BfMat **L, BfMat **M) {
  BF_ERROR_BEGIN();

  BfSize *rowptr = NULL;
  BfSize *colind = NULL;
  BfReal *L_data = NULL;
  BfReal *M_data = NULL;

  BfSize numVerts = bfTrimeshGetNumVerts(trimesh);

  BfSize nnz = numVerts + trimesh->vvOffset[numVerts];

  rowptr = bfMemAlloc(numVerts + 1, sizeof(BfSize));
  if (rowptr == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  colind = bfMemAlloc(nnz, sizeof(BfSize));
  if (colind == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  L_data = bfMemAllocAndZero(nnz, sizeof(BfReal));
  if (L_data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  M_data = bfMemAllocAndZero(nnz, sizeof(BfReal));
  if (M_data == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* Fill rowptr: same as vvOffset except we need to add more space
   * for the vertices on the diagonal */
  for (BfSize i = 0; i <= numVerts; ++i)
    rowptr[i] = trimesh->vvOffset[i] + i;

  /* Fill colind */
  for (BfSize i = 0, j = 0; i < numVerts; ++i) {
    BF_ASSERT(j == rowptr[i]);

    /* Find the position of i in vv */
    BfSize kmid = trimesh->vvOffset[i];
    while (kmid < trimesh->vvOffset[i + 1] && trimesh->vv[kmid] < i)
      ++kmid;

    /* Copy over vv, inserting i into the correct sorted order */
    BfSize k = trimesh->vvOffset[i];
    for (; k < kmid; ++k)
      colind[j++] = trimesh->vv[k];
    colind[j++] = i;
    for (; k < trimesh->vvOffset[i + 1]; ++k)
      colind[j++] = trimesh->vv[k];
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

  BfMatCsrReal *L_csr = bfMatCsrRealNew();
  HANDLE_ERROR();

  bfMatCsrRealInit(L_csr, numVerts, numVerts, rowptr, colind, L_data);
  HANDLE_ERROR();

  BfMatCsrReal *M_csr = bfMatCsrRealNew();
  HANDLE_ERROR();

  bfMatCsrRealInit(M_csr, numVerts, numVerts, rowptr, colind, M_data);
  HANDLE_ERROR();

  *L = bfMatCsrRealToMat(L_csr);
  *M = bfMatCsrRealToMat(M_csr);

  END_ERROR_HANDLING() {
    bfMatCsrRealDeinitAndDealloc(&L_csr);
    bfMatCsrRealDeinitAndDealloc(&M_csr);
  }

  bfMemFree(rowptr);
  bfMemFree(colind);
  bfMemFree(L_data);
  bfMemFree(M_data);
}
