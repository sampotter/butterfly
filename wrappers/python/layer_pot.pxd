cdef extern from "bf/layer_pot.h":
    cdef enum BfLayerPotentials:
        BF_LAYER_POTENTIAL_UNKNOWN
        BF_LAYER_POTENTIAL_SINGLE
        BF_LAYER_POTENTIAL_PV_DOUBLE
        BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE
        BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_DOUBLE
        BF_LAYER_POTENTIAL_COMBINED_FIELD
        BF_LAYER_POTENTIAL_COUNT

    ctypedef BfLayerPotentials BfLayerPotential
