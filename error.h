#pragma once

enum BfError {
  BF_ERROR_NO_ERROR = 0,
  BF_ERROR_INVALID_ARGUMENTS = (1 << 0),
  BF_ERROR_RUNTIME_ERROR = (1 << 1)
};
