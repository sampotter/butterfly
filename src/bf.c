#include <bf/bf.h>

#ifdef BF_PYTHON
#  include "numpy.h"
#endif

void bfInit(void) {
#ifdef BF_PYTHON
  _import_array();
#endif
}
