#pragma once

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
 * where D and S' are both evaluated in the principal value sense. */
typedef enum BfLayerPotentials {
  BF_LAYER_POTENTIAL_SINGLE,
  BF_LAYER_POTENTIAL_PV_DOUBLE,
  BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE,

  /* The number of different layer potentials */
  BF_LAYER_POTENTIAL_COUNT,

  /* Unknown layer potential value */
  BF_LAYER_POTENTIAL_UNKNOWN
} BfLayerPotential;

/* A lookup table consisting of which layer potentials to use for
 * reexpansion when constructing proxy sets. */
static BfLayerPotential const BF_PROXY_LAYER_POT[BF_LAYER_POTENTIAL_COUNT] = {
  [BF_LAYER_POTENTIAL_SINGLE] = BF_LAYER_POTENTIAL_SINGLE,
  [BF_LAYER_POTENTIAL_PV_DOUBLE] = BF_LAYER_POTENTIAL_PV_DOUBLE,
  [BF_LAYER_POTENTIAL_PV_NORMAL_DERIV_SINGLE] = BF_LAYER_POTENTIAL_SINGLE
};
