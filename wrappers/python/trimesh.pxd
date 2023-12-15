from defs cimport BfSize
from types cimport BfTrimesh

cdef extern from "bf/trimesh.h":
    BfTrimesh *bfTrimeshNewFromObjFile(const char *objPath)
    BfSize bfTrimeshGetNumFaces(const BfTrimesh *trimesh)
