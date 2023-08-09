#include <bf/bf.h>

#include <stdio.h>

#include "numpy.h"

void bfInit(void) {
  printf("bfInit()\n");

  _import_array();
}
