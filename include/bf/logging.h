#pragma once

#include <stdio.h>

typedef enum BfLogLevels {
  BF_LOG_LEVEL_TODO = 0,
  BF_LOG_LEVEL_DEBUG = 1,
  BF_LOG_LEVEL_INFO = 2,
  BF_LOG_LEVEL_WARN = 3,
  BF_LOG_LEVEL_ERROR = 4,
} BfLogLevel;

void bfSetLogLevel(BfLogLevel logLevel);
void bfSetLogStream(FILE *restrict logStream);
void bfLog(BfLogLevel logLevel, char const *restrict format, ...);
void bfLogTodo(char const *restrict format, ...);
void bfLogDebug(char const *restrict format, ...);
void bfLogInfo(char const *restrict format, ...);
