#pragma once

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

/* NOTE: in all translation units except bf.c, the correct way to
 * include this file is:
 *
 *   #define NO_IMPORT_ARRAY
 *   #include "numpy.h"
 *
 * see the link to the NumPy docs in the comment below for more
 * explanation about NumPy's special needs. */

/* NumPy has some wacky requirements to be initialized correctly. See
 * this link:
 *
 *   https://numpy.org/devdocs/reference/c-api/array.html#importing-the-api
 *
 * for an explanation... */
#define PY_ARRAY_UNIQUE_SYMBOL bf_ARRAY_API
#include <numpy/arrayobject.h>

#ifdef BF_DOUBLE
static int const BF_COMPLEX_TYPENUM = NPY_COMPLEX128;
static int const BF_REAL_TYPENUM = NPY_FLOAT64;
#else
#  error "Not implemented yet!"
#endif
