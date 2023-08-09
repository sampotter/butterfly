#include <bf/mat_python.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/macros.h>
#include <bf/mat_dense_complex.h>
#include <bf/mem.h>

#define NO_IMPORT_ARRAY
#include "numpy.h"

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .GetView = (__typeof__(&bfMatGetView))bfMatPythonGetView,
  .GetType = (__typeof__(&bfMatGetType))bfMatPythonGetType,
  .GetNumRows = (__typeof__(&bfMatGetNumRows))bfMatPythonGetNumRows,
  .GetNumCols = (__typeof__(&bfMatGetNumCols))bfMatPythonGetNumCols,
  .Mul = (__typeof__(&bfMatMul))bfMatPythonMul,
};

BfMat *bfMatPythonGetView(BfMatPython *matPython) {
  BF_ERROR_BEGIN();

  BfMat *matView = NULL;

  BfMatPython *matPythonView = bfMatPythonAlloc();
  HANDLE_ERROR();

  *matPythonView = *matPython;

  matView = bfMatPythonToMat(matPythonView);
  matView->props |= BF_MAT_PROPS_VIEW;

  BF_ERROR_END() {
    BF_DIE();
  }

  return matView;
}

BfType bfMatPythonGetType(BfMatPython const *matPython) {
  (void)matPython;
  return BF_TYPE_MAT_PYTHON;
}

BfSize bfMatPythonGetNumRows(BfMatPython const *matPython) {
  BfMat const *mat = bfMatPythonConstToMatConst(matPython);
  return bfMatIsTransposed(mat) ? mat->numCols : mat->numRows;
}

BfSize bfMatPythonGetNumCols(BfMatPython const *matPython) {
  BfMat const *mat = bfMatPythonConstToMatConst(matPython);
  return bfMatIsTransposed(mat) ? mat->numRows : mat->numCols;
}

static BfMat *mul_matDenseComplex(BfMatPython const *matPython, BfMatDenseComplex const *matDenseComplex) {
  BF_ERROR_BEGIN();

  BfMat *res = NULL;

  if (!PyObject_HasAttrString(matPython->obj, "_Mul"))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  PyObject *_MulObj = PyObject_GetAttrString(matPython->obj, "_Mul");

  int nd = 2;
  npy_intp dims[2] = {
    bfMatDenseComplexGetNumRows(matDenseComplex),
    bfMatDenseComplexGetNumCols(matDenseComplex)
  };
  int typenum = BF_COMPLEX_TYPENUM;
  void *data = matDenseComplex->data;
  PyObject *matDenseComplexObj = PyArray_SimpleNewFromData(nd, dims, typenum, data);

  PyObject *argsObj = PyTuple_New(1);
  PyTuple_SetItem(argsObj, 0, matDenseComplexObj);

  PyObject *resObj = PyObject_Call(_MulObj, argsObj, NULL);
  if (!PyArray_Check(resObj))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  PyArrayObject *resArrayObj = (PyArrayObject *)resObj;
  BF_ASSERT(PyArray_NDIM(resArrayObj) == 2);

  res = BF_TO_MAT(bfMatDenseComplexNewViewFromPyArray((BfPtr)resObj));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  // Py_DECREF(_MulObj);
  // Py_DECREF(matDenseComplexObj);
  // Py_DECREF(argsObj);
  // Py_DECREF(resObj);

  return res;
}

BfMat *bfMatPythonMul(BfMatPython const *matPython, BfMat const *otherMat) {
  BF_ERROR_BEGIN();

  if (bfMatPythonGetNumCols(matPython) != bfMatGetNumRows(otherMat))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfMat *result = NULL;
  switch (bfMatGetType(otherMat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    result = mul_matDenseComplex(matPython, bfMatConstToMatDenseComplexConst(otherMat));
    break;
  default:
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);
  }

  BF_ERROR_END() {
    BF_DIE();
  }

  return result;
}

/** Upcasting: MatPython -> Mat */

BfMat *bfMatPythonToMat(BfMatPython *matPython) {
  return &matPython->super;
}

BfMat const *bfMatPythonConstToMatConst(BfMatPython const *matPython) {
  return &matPython->super;
}

/** Implementation: MatPython */

BfMatPython *bfMatPythonAlloc(void) {
  BF_ERROR_BEGIN();

  BfMatPython *matPython = bfMemAlloc(1, sizeof(BfMatPython));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return matPython;
}

BfMatPython *bfMatPythonNewFromPyObject(PyObject *obj, BfSize numRows, BfSize numCols) {
  BF_ERROR_BEGIN();

  BfMatPython *matPython = bfMatPythonAlloc();
  HANDLE_ERROR();

  bfMatPythonInitFromPyObject(matPython, obj, numRows, numCols);

  BF_ERROR_END() {
    BF_DIE();
  }

  return matPython;
}

void bfMatPythonInitFromPyObject(BfMatPython *matPython, PyObject *obj, BfSize numRows, BfSize numCols) {
  bfMatInit(&matPython->super, &MAT_VTABLE, numRows, numCols);

  matPython->obj = obj;
}
