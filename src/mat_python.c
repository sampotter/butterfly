#include <bf/mat_python.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/macros.h>
#include <bf/mat_dense_complex.h>
#include <bf/mem.h>

#ifndef BF_PYTHON
#  error "Building mat_python.c without python Meson feature enabled"
#endif

#define NO_IMPORT_ARRAY
#include "numpy.h"

/** Interface: Mat */

static BfMatVtable MAT_VTABLE = {
  .GetView = (__typeof__(&bfMatGetView))bfMatPythonGetView,
  .GetType = (__typeof__(&bfMatGetType))bfMatPythonGetType,
  .GetNumRows = (__typeof__(&bfMatGetNumRows))bfMatPythonGetNumRows,
  .GetNumCols = (__typeof__(&bfMatGetNumCols))bfMatPythonGetNumCols,
  .Mul = (__typeof__(&bfMatMul))bfMatPythonMul,
  .Rmul = (__typeof__(&bfMatRmul))bfMatPythonRmul,
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

  if (bfMatPythonGetNumCols(matPython) != bfMatDenseComplexGetNumRows(matDenseComplex))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (!PyObject_HasAttrString(matPython->obj, "_Mul"))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  PyObject *_MulObj = PyObject_GetAttrString(matPython->obj, "_Mul");

  /* Create PyArray w/ a view of the contents of `matDenseComplex`: */
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
  switch (bfMatGetType(otherMat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    return mul_matDenseComplex(matPython, bfMatConstToMatDenseComplexConst(otherMat));
    break;
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

static BfMat *rmul_matDenseComplex(BfMatPython const *matPython, BfMatDenseComplex const *matDenseComplex) {
  BF_ERROR_BEGIN();

  BfMat *result = NULL;

  if (bfMatDenseComplexGetNumCols(matDenseComplex) != bfMatPythonGetNumRows(matPython))
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (!PyObject_HasAttrString(matPython->obj, "_Rmul"))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* Get the `_Rmul` method of the `bf.MatPython` instance pointed to
   * by `matPython->obj`: */
  PyObject *_RmulObj = PyObject_GetAttrString(matPython->obj, "_Rmul");

  /* Create a `PyArray` w/ a view of the contents of `matDenseComplex`: */
  int nd = 2;
  npy_intp dims[2] = {
    bfMatDenseComplexGetNumRows(matDenseComplex),
    bfMatDenseComplexGetNumCols(matDenseComplex)
  };
  int typenum = BF_COMPLEX_TYPENUM;
  void *data = matDenseComplex->data;
  PyObject *matDenseComplexObj = PyArray_SimpleNewFromData(nd, dims, typenum, data);

  /* Prepare the arguments for the call to `_Rmul`: */
  PyObject *argsObj = PyTuple_New(1);
  PyTuple_SetItem(argsObj, 0, matDenseComplexObj);

  /* Call `_Rmul`: */
  PyObject *resObj = PyObject_Call(_RmulObj, argsObj, NULL);
  if (!PyArray_Check(resObj))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* Make sure the result has two dimensions: */
  PyArrayObject *resArrayObj = (PyArrayObject *)resObj;
  BF_ASSERT(PyArray_NDIM(resArrayObj) == 2);

  /* Create a `BfMatDenseComplex` with a view of the the result: */
  result = BF_TO_MAT(bfMatDenseComplexNewViewFromPyArray((BfPtr)resObj));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  // Py_DECREF(_RmulObj);
  // Py_DECREF(matDenseComplexObj);
  // Py_DECREF(argsObj);
  // Py_DECREF(resObj);

  return result;
}

BfMat *bfMatPythonRmul(BfMatPython const *matPython, BfMat const *otherMat) {
  switch (bfMatGetType(otherMat)) {
  case BF_TYPE_MAT_DENSE_COMPLEX:
    return rmul_matDenseComplex(matPython, bfMatConstToMatDenseComplexConst(otherMat));
    break;
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

/** Upcasting: MatPython -> Mat */

BfMat *bfMatPythonToMat(BfMatPython *matPython) {
  return &matPython->super;
}

BfMat const *bfMatPythonConstToMatConst(BfMatPython const *matPython) {
  return &matPython->super;
}

/** Downcasting: Mat -> MatPython */

BfMatPython *bfMatToMatPython(BfMat *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_PYTHON)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatPython *)mat;
  }
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
