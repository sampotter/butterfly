#pragma once

#define SWAP(x, y) do {                         \
    __typeof__(x) tmp = x;                      \
    x = y;                                      \
    y = tmp;                                    \
  } while (0)
