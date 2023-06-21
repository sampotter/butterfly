#include <bf/logging.h>

#include <stdarg.h>

static BfLogLevel LOG_LEVEL = BF_LOG_LEVEL_ERROR;

static char const *LOG_LEVEL_STRING[] = {
  [BF_LOG_LEVEL_TODO] = "TODO",
  [BF_LOG_LEVEL_DEBUG] = "DEBUG",
  [BF_LOG_LEVEL_INFO] = "INFO",
  [BF_LOG_LEVEL_WARN] = "WARN",
  [BF_LOG_LEVEL_ERROR] = "ERROR",
};

static FILE *LOG_STREAM = NULL;

void bfSetLogLevel(BfLogLevel logLevel) {
  LOG_LEVEL = logLevel;
}

void bfSetLogStream(FILE *restrict logStream) {
  LOG_STREAM = logStream;
}

FILE *bfGetLogStream() {
  return LOG_STREAM == NULL ? stdout : LOG_STREAM;
}

void bfLog(BfLogLevel logLevel, char const *restrict format, ...) {
  if (logLevel < LOG_LEVEL)
    return;
  FILE *logStream = bfGetLogStream();
  fprintf(logStream, "bf: %s: ", LOG_LEVEL_STRING[logLevel]);
  va_list argp;
  va_start(argp, format);
  vfprintf(logStream, format, argp);
  va_end(argp);
}

void bfLogTodo(char const *restrict format, ...) {
  BfLogLevel logLevel = BF_LOG_LEVEL_TODO;
  if (logLevel < LOG_LEVEL)
    return;
  FILE *logStream = bfGetLogStream();
  fprintf(logStream, "bf: %s: ", LOG_LEVEL_STRING[logLevel]);
  va_list argp;
  va_start(argp, format);
  vfprintf(logStream, format, argp);
  va_end(argp);
}

void bfLogDebug(char const *restrict format, ...) {
  BfLogLevel logLevel = BF_LOG_LEVEL_DEBUG;
  if (logLevel < LOG_LEVEL)
    return;
  FILE *logStream = bfGetLogStream();
  fprintf(logStream, "bf: %s: ", LOG_LEVEL_STRING[logLevel]);
  va_list argp;
  va_start(argp, format);
  vfprintf(logStream, format, argp);
  va_end(argp);
}

void bfLogInfo(char const *restrict format, ...) {
  BfLogLevel logLevel = BF_LOG_LEVEL_INFO;
  if (logLevel < LOG_LEVEL)
    return;
  FILE *logStream = bfGetLogStream();
  fprintf(logStream, "bf: %s: ", LOG_LEVEL_STRING[logLevel]);
  va_list argp;
  va_start(argp, format);
  vfprintf(logStream, format, argp);
  va_end(argp);
}
