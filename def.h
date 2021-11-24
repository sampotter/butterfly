#pragma once

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef BF_DOUBLE
#  ifdef BF_SINGLE
#    error "both double and single precision requested"
#  endif
typedef double BfReal;
typedef double _Complex BfComplex;
typedef double BfPoint2[2];
typedef size_t BfSize;
#elif BF_SINGLE
#  error "single precision not yet implemented"
#else
#  error "neither single nor double precision specified"
#endif
