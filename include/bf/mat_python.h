#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "mat.h"

/** Interface: Mat */

BfMat *bfMatPythonGetView(BfMatPython *matPython);
BfType bfMatPythonGetType(BfMatPython const *matPython);
BfSize bfMatPythonGetNumRows(BfMatPython const *matPython);
BfSize bfMatPythonGetNumCols(BfMatPython const *matPython);
BfMat *bfMatPythonMul(BfMatPython const *matPython, BfMat const *otherMat);
BfMat *bfMatPythonRmul(BfMatPython const *matPython, BfMat const *otherMat);

/** Upcasting: MatPython -> Mat */

BfMat *bfMatPythonToMat(BfMatPython *matPython);
BfMat const *bfMatPythonConstToMatConst(BfMatPython const *matPython);

/** Implementation: MatPython */

struct BfMatPython {
  BfMat super;

  /* A pointer to the Python extension class backing this instance. */
  PyObject *obj;
};

BfMatPython *bfMatPythonAlloc(void);
BfMatPython *bfMatPythonNewFromPyObject(PyObject *obj, BfSize numRows, BfSize numCols);
void bfMatPythonInitFromPyObject(BfMatPython *matPython, PyObject *obj, BfSize numRows, BfSize numCols);
