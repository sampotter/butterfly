#pragma once

#define BF_DEFINE_EDGE(VAR_NAME, I0, I1) \
  BfSize2 VAR_NAME = {I0, I1};           \
  SORT2(VAR_NAME[0], VAR_NAME[1]);

#define BF_DEFINE_FACE(VAR_NAME, I0, I1, I2)    \
  BfSize3 VAR_NAME = {I0, I1, I2};              \
  SORT3(VAR_NAME[0], VAR_NAME[1], VAR_NAME[2]);

#define BF_SIZE_OK(s) (s != BF_SIZE_BAD_VALUE)

#define SWAP(x, y) do {                         \
    __typeof__(x) tmp = x;                      \
    x = y;                                      \
    y = tmp;                                    \
  } while (0)

#define BF_SWAP SWAP

#define SORT2(x, y) do {                        \
    __typeof__(x) tmp = x < y ? x : y;          \
    y = x < y ? y : x;                          \
    x = tmp;                                    \
  } while (0)

#define SORT3(x, y, z)                          \
  SORT2(x, y);                                  \
  SORT2(x, z);                                  \
  SORT2(y, z);
