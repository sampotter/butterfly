#pragma once

#if defined BF_LINUX
#  include <openblas/cblas.h>
#  include <openblas/lapacke.h>
#elif defined BF_DARWIN
#  include <cblas.h>
#  include <lapacke.h>
#endif
