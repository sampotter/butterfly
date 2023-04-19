#pragma once

#ifdef BF_DEBUG
#if __GNUC__
#define BF_ASSERT(c) if (!(c)) __builtin_trap()
#elif _MSC_VER
#define BF_ASSERT(c) if (!(c)) __debugbreak()
#else
#define BF_ASSERT(c) if (!(c)) *(volatile int *)0 = 0
#endif
#else
#include <assert.h>
#define BF_ASSERT(c) assert(c)
#endif

#define BF_DIE() BF_ASSERT(false)
