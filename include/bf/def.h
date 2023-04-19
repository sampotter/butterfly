#pragma once

/* For getline */
#define _POSIX_C_SOURCE 200809L

/* For j0, j1, y0, and y1: */
#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE 1
#endif

/* For qsort_r: */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

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

/* TODO: only define if we're using LAPACK? make sure it's right,
 * anyway... should get this info from Meson */
typedef int BfLapackInt;

typedef BfSize BfSize3[3];

static BfSize const BF_SIZE_BAD_VALUE = -1;

typedef void *BfPtr;
typedef void const *BfConstPtr;

static BfSize const BF_DEFAULT_STRIDE = 1;

static BfSize const BF_ARRAY_DEFAULT_CAPACITY = 128;

typedef int (*BfCompar)(const void *, const void *, void *);
