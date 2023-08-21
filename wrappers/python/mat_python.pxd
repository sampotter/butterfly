from defs cimport BfSize
from types cimport BfMat, BfMatPython

cdef extern from "bf/mat_python.h":
    BfMat *bfMatPythonToMat(BfMatPython *matPython)
    BfMatPython *bfMatToMatPython(BfMat *mat)
    BfMatPython *bfMatPythonNewFromPyObject(object obj, BfSize numRows, BfSize numCols)
