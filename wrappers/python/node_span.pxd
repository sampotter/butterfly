from defs cimport BfPolicy
from types cimport BfPtrArray, BfNodeSpan

cdef extern from "bf/node_span.h":
    BfNodeSpan *bfNodeSpanNewFromPtrArray(BfPtrArray *ptrArray, BfPolicy policy)
