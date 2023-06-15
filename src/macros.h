#pragma once

#define SWAP(x, y) do {                         \
    __typeof__(x) tmp = x;                      \
    x = y;                                      \
    y = tmp;                                    \
  } while (0)

#define SORT2(x, y) do {                        \
    __typeof__(x) tmp = x < y ? x : y;          \
    y = x < y ? y : x;                          \
    x = tmp;                                    \
  } while (0)
