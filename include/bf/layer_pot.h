#pragma once

#include "def.h"

/* An enumeration defining the different layer potential operators
 * commonly encountered when working with boundary integral equations.
 *
 * If the kernel or Green's function of interest is G(x - y), then the
 * single-layer potential is:
 *
 *   [Sσ](x) = ∫G(x, y)σ(y)dy,
 *
 * the double-layer potential is:
 *
 *   [Dσ](x) = ∫ν(y)∙∇G(x - y)σ(y)dy,
 *
 * and the normal derivative of the single-layer potential is:
 *
 *   [S'σ](x) = ∫ν(x)∙∇G(x - y)σ(y)dy,
 *
 * where D and S' are both evaluated in the principal value sense.
 *
 * The combined field layer potential is:
 *
 *   αS + βD
 *
 * where α and β are scalars. */
typedef enum BfLayerPotentials {
  /* A default, unknown layer potential: */
  BF_LAYER_POTENTIAL_UNKNOWN = 0,

  /* The layer potentials: */
  BF_LAYER_POTENTIAL_SINGLE,
  BF_LAYER_POTENTIAL_PV_DOUBLE,
  BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE,
  BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_DOUBLE,
  BF_LAYER_POTENTIAL_COMBINED_FIELD,

  /* The number of different layer potentials: */
  BF_LAYER_POTENTIAL_COUNT,
} BfLayerPotential;

static bool const BF_LAYER_POT_USES_SRC_NORMALS[BF_LAYER_POTENTIAL_COUNT] = {
  [BF_LAYER_POTENTIAL_PV_DOUBLE] = true,
  [BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_DOUBLE] = true,
  [BF_LAYER_POTENTIAL_COMBINED_FIELD] = true,
};

static bool const BF_LAYER_POT_USES_TGT_NORMALS[BF_LAYER_POTENTIAL_COUNT] = {
  [BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE] = true,
  [BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_DOUBLE] = true,
};

/* A lookup table consisting of which layer potentials to use for
 * reexpansion when constructing proxy sets.
 *
 * Note that this map basically just has the effect of stripping off a
 * normal derivative on the target points.
 *
 * NOTE: anything uninitialized here will map to
 * `BF_LAYER_POTENTIAL_UNKNOWN` by default because of zero
 * initialization. */
static BfLayerPotential const BF_PROXY_LAYER_POT[BF_LAYER_POTENTIAL_COUNT] = {
  [BF_LAYER_POTENTIAL_SINGLE] = BF_LAYER_POTENTIAL_SINGLE,
  [BF_LAYER_POTENTIAL_PV_DOUBLE] = BF_LAYER_POTENTIAL_PV_DOUBLE,
  [BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE] = BF_LAYER_POTENTIAL_SINGLE,
  [BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_DOUBLE] = BF_LAYER_POTENTIAL_PV_DOUBLE,
  [BF_LAYER_POTENTIAL_COMBINED_FIELD] = BF_LAYER_POTENTIAL_COMBINED_FIELD
};
