#pragma once

#define BF_ERROR_BEGIN()                        \
  enum BfError __bf_error__;                    \
  bool __bf_erred__ = false;                    \
  (void)__bf_error__;                           \
  (void)__bf_erred__;

#define RAISE_ERROR(error) do { \
    bfSetError(error);          \
    __bf_erred__ = true;        \
    goto __bf_cleanup__;        \
  } while (0);

#define HANDLE_ERROR() do {                     \
    __bf_error__ = bfGetError();                \
    if (__bf_error__) {                         \
      bfSetError(__bf_error__);                 \
      __bf_erred__ = true;                      \
      goto __bf_cleanup__;                      \
    }                                           \
  } while (0);

#define END_ERROR_HANDLING()                    \
__bf_cleanup__:                                 \
  if (__bf_erred__)
