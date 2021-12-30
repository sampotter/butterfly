#pragma once

#define RAISE_ERROR(error) do { \
    bfSetError(error);          \
    erred = true;               \
    goto cleanup;               \
  } while (0);

#define HANDLE_ERROR() do {                     \
    error = bfGetError();                       \
    if (error) {                                \
      bfSetError(error);                        \
      erred = true;                             \
      goto cleanup;                             \
    }                                           \
  } while (0);
