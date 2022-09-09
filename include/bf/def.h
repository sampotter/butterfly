#pragma once

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef BF_DOUBLE
#  ifdef BF_SINGLE
#    error "both double and single precision requested"
#  endif
typedef double BfReal;
typedef double _Complex BfComplex;
#define BF_EPS_MACH 2.220446049250313e-16
#elif BF_SINGLE
#  error "single precision not yet implemented"
#else
#  error "neither single nor double precision specified"
#endif

typedef uint8_t BfByte;
typedef size_t BfSize;

#define BF_SIZE_BAD_VALUE (BfSize)(-1)

typedef void *BfPtr;
