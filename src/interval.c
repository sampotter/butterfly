#include <bf/interval.h>

#include <bf/assert.h>
#include <bf/const.h>

void bfIntervalDeinit(BfInterval *interval) {
  interval->endpoint[0] = interval->endpoint[1] = BF_NAN;
  interval->closed[0] = interval->closed[1] = false;
}

bool bfIntervalLeftOf(BfInterval const *I1, BfInterval const *I2) {
  return I1->closed[1] && I2->closed[0] ?
    I1->endpoint[1] < I2->endpoint[0] :
    I1->endpoint[1] <= I2->endpoint[0];
}

bool bfIntervalRightOf(BfInterval const *I1, BfInterval const *I2) {
  return bfIntervalLeftOf(I2, I1);
}

bool bfIntervalOverlaps(BfInterval const *I1, BfInterval const *I2) {
  return !bfIntervalLeftOf(I1, I2) && !bfIntervalRightOf(I1, I2);
}

bool bfIntervalContains(BfInterval const *I1, BfInterval const *I2) {
  bool leftContained = I1->closed[0] ?
    I1->endpoint[0] <= I2->endpoint[0] : I1->endpoint[0] < I2->endpoint[0];
  bool rightContained = I1->closed[1] ?
    I2->endpoint[1] <= I1->endpoint[1] : I2->endpoint[1] < I1->endpoint[1];
  return leftContained && rightContained;
}

BfSize bfIntervalDifference(BfInterval const *I1, BfInterval const *I2, BfInterval *Idiff) {
  if (!bfIntervalOverlaps(I1, I2))
    return 0;

  if (bfIntervalContains(I1, I2)) {
    BfSize i = 0;

    Idiff[i].endpoint[0] = I1->endpoint[0];
    Idiff[i].endpoint[1] = I2->endpoint[0];

    Idiff[i].closed[0] = I1->closed[0];
    Idiff[i].closed[1] = !I2->closed[0];

    BF_ASSERT(Idiff[i].endpoint[0] <= Idiff[i].endpoint[1]);
    if (Idiff[i].endpoint[0] < Idiff[i].endpoint[1] ||
        (Idiff[i].closed[0] && Idiff[i].closed[1]))
      ++i;

    Idiff[i].endpoint[0] = I2->endpoint[1];
    Idiff[i].endpoint[1] = I1->endpoint[1];

    Idiff[i].closed[0] = !I2->closed[1];
    Idiff[i].closed[1] = I1->closed[1];

    BF_ASSERT(Idiff[i].endpoint[0] <= Idiff[i].endpoint[1]);
    if (Idiff[i].endpoint[0] < Idiff[i].endpoint[1] ||
        (Idiff[i].closed[0] && Idiff[i].closed[1]))
      ++i;

    for (BfSize j = i; j < 2; ++j)
      bfIntervalDeinit(&Idiff[j]);

    return i;
  }

  if (I2->endpoint[0] <= I1->endpoint[0] &&
      I1->endpoint[0] <= I2->endpoint[1] &&
      I2->endpoint[1] <= I1->endpoint[1]) {
    Idiff[0].endpoint[0] = I2->endpoint[1];
    Idiff[0].endpoint[1] = I1->endpoint[1];

    Idiff[0].closed[0] = !I2->closed[1];
    Idiff[0].closed[1] = I2->closed[1];

    return 1;
  }

  if (I1->endpoint[0] <= I2->endpoint[0] &&
      I2->endpoint[0] <= I1->endpoint[1] &&
      I1->endpoint[1] <= I2->endpoint[1]) {
    Idiff[0].endpoint[0] = I1->endpoint[0];
    Idiff[0].endpoint[1] = I2->endpoint[0];

    Idiff[0].closed[0] = I1->closed[0];
    Idiff[0].closed[1] = !I2->closed[0];

    return 1;
  }

  BF_DIE();
}
