cdef extern from "bf/tree_traversals.h":
    cdef enum BfTreeTraversal:
        BF_TREE_TRAVERSAL_UNKNOWN
        BF_TREE_TRAVERSAL_LR_LEVEL_ORDER
        BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER
        BF_TREE_TRAVERSAL_POST_ORDER
