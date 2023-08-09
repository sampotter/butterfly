from fac cimport BfFacSpec
from types cimport BfMat

cdef extern from "bf/fac_streamer.h":
    ctypedef struct BfFacStreamer:
        pass

    BfFacStreamer *bfFacStreamerNew()
    void bfFacStreamerInit(BfFacStreamer *facStreamer, const BfFacSpec *facSpec)
    void bfFacStreamerFeed(BfFacStreamer *facStreamer, BfMat *Phi)
    bint bfFacStreamerIsDone(const BfFacStreamer *facStreamer)
