#include <bf/rand.h>

#include <math.h>
#include <stdint.h>

#include <bf/const.h>

/* Declaration of PRNGs. Implementations are in xoshiro256plus.c and
 * splitmix64.c. These PRNGs come from https://prng.di.unimi.it/ and
 * supposedly overcome some well-known defects belonging to the
 * Mersenne Twister.  */

void splitmix64_seed(uint64_t seed);
uint64_t splitmix64_next();

void xoshiro256plus_seed(uint64_t const seed[4]);
uint64_t xoshiro256plus_next();

void bfSeed(BfSize seed) {
  /* seed splitmix64 using seed... */
  splitmix64_seed(seed);

  /* ... generate seeds for xoshiro256 using splitmix64... */
  uint64_t seeds[4];
  seeds[0] = splitmix64_next();
  seeds[1] = splitmix64_next();
  seeds[2] = splitmix64_next();
  seeds[3] = splitmix64_next();

  /* ... and seed xoshiro256 */
  xoshiro256plus_seed(seeds);
}

/* Compute a random BfSize in the half-open interval [low, high). */
BfSize bfSizeUniform1(BfSize low, BfSize high) {
  double x = bfRealUniform1();
  x *= high - low;
  x += low;
  return (BfSize)floor(x);
}

BfReal bfRealUniform1() {
  /* For an explanation, see https://prng.di.unimi.it/, section
   * "Generating uniform doubles in the unit interval". */
  return (xoshiro256plus_next() >> 11) * 0x1.0p-53;
}

void bfRealUniform(BfSize n, BfReal *x) {
  for (BfSize i = 0; i < n; ++i) x[i] = bfRealUniform1();
}

void bfRealRandn(BfSize n, BfReal *x) {
  /* fill the first 2*floor(n/2) entries of x with uniform deviates
   * from [0, 1] */
  bfRealUniform(2*(n/2), x);

  /* use Box-Muller to compute 2*(n/2) independent N(0, 1) deviates */
  BfReal mag, theta;
  for (BfSize i = 0; i < n/2; ++i) {
    mag = sqrt(-2*log(x[2*i]));
    theta = BF_TWO_PI*x[2*i + 1];
    x[2*i] = mag*cos(theta);
    x[2*i + 1] = mag*sin(theta);
  }

  /* if x has an odd length, sample N(0, 1) one more time */
  if (n % 2 == 1) {
    BfReal u[2];
    bfRealUniform(2, u);
    x[n - 1] = sqrt(-2*log(u[0]))*cos(BF_TWO_PI*u[1]);
  }
}

void bfComplexRandn(BfSize n, BfComplex *x) {
  bfRealRandn(2*n, (BfReal *)x);
}

void bfSampleRandomUnitVector2(BfVector2 u) {
  bfRealRandn(2, u);
  BfReal uMag = hypot(u[0], u[1]);
  u[0] /= uMag;
  u[1] /= uMag;
}
