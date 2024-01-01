from defs cimport BfSize
from types cimport BfTrimesh
from vectors cimport BfVectors3

cdef extern from "bf/trimesh.h":
    BfTrimesh *bfTrimeshNewFromObjFile(const char *objPath)
    BfSize bfTrimeshGetNumFaces(const BfTrimesh *trimesh)
    bint bfTrimeshHasFaceNormals(const BfTrimesh *trimesh)
    bint bfTrimeshHasVertexNormals(const BfTrimesh *trimesh)
    BfVectors3 *bfTrimeshGetFaceNormalsPtr(BfTrimesh *trimesh)
    BfVectors3 *bfTrimeshGetVertexNormalsPtr(BfTrimesh *trimesh)
    void bfTrimeshComputeFaceNormalsMatchingVertexNormals(BfTrimesh *trimesh)
