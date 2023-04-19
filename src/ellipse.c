#include <bf/ellipse.h>

#include <assert.h>

#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/points.h>
#include <bf/vec_real.h>
#include <bf/vectors.h>

BfReal bfEllipseGetPerimeter(BfEllipse const *ellipse) {
  BfReal a = ellipse->semiMajorAxis;
  BfReal b = ellipse->semiMinorAxis;
  BfReal h = (a - b)/(a + b);

  /* Compute ellipse perimeter using Gauss-Kummer series:
   * https://mathworld.wolfram.com/Ellipse.html */
  BfSize m = 1;
  BfReal sum = 0, term = 1;
  while (fabs(term) > 1e-15) {
    sum += term;
    term = pow(tgamma(1 + 0.5)/(tgamma(1 + 0.5 - m)*tgamma(1 + m)), 2);
    term *= pow(h, 2*m);
    ++m;
  }

  BfReal p = BF_PI*(a + b)*sum;

  return p;
}

void bfEllipseSampleLinspaced(BfEllipse const *ellipse, BfSize numPoints, BfPoints2 *points, BfVectors2 *unitTangents, BfVectors2 *unitNormals, BfRealArray *weights) {
  BfReal a = ellipse->semiMajorAxis;
  BfReal b = ellipse->semiMinorAxis;

  BfReal h = BF_TWO_PI/numPoints;

  for (BfSize i = 0; i < numPoints; ++i) {
    BfReal theta = h*i;

    /* Compute point: */
    BfReal px = a*cos(theta);
    BfReal py = b*sin(theta);
    BfPoint2 p = {px, py};
    bfPoint2RotateAboutOrigin(p, ellipse->theta);
    bfPoint2Translate(p, ellipse->center);
    if (points) bfPoints2Append(points, p);

    /* Compute unit tangent vector: */
    BfReal tx = -a*sin(theta);
    BfReal ty = b*cos(theta);
    BfVector2 t = {tx, ty};

    /* If we're computing weights, get the current one now: */
    if (weights != NULL)
      bfRealArrayAppend(weights, h*hypot(tx, ty));

    bfVector2Normalize(t);

    /* Compute unit normal vector: */
    BfReal nx = -a*cos(theta);
    BfReal ny = -b*sin(theta);
    BfVector2 n = {nx, ny};
    bfVector2Reject(n, t);
    bfVector2Normalize(n);
    bfVector2Rotate(n, ellipse->theta);
    bfVector2Negate(n); // get outward-facing normal
    if (unitNormals) bfVectors2Append(unitNormals, n);

    bfVector2Rotate(t, ellipse->theta);
    if (unitTangents) bfVectors2Append(unitTangents, t);
  }
}

void bfEllipseSampleEquispaced(BfEllipse const *ellipse, BfSize numPoints,
                               BfPoints2 *points, BfVectors2 *unitTangents,
                               BfVectors2 *unitNormals) {
  BfReal a = ellipse->semiMajorAxis;
  BfReal b = ellipse->semiMinorAxis;

  BfReal dtheta = BF_TWO_PI/numPoints;

  BfReal *D = bfMemAlloc((numPoints + 1), sizeof(BfReal));
  D[0] = 0;
  for (BfSize i = 1; i <= numPoints; ++i) {
    BfReal dx = a*cos(i*dtheta) - a*cos((i - 1)*dtheta);
    BfReal dy = b*sin(i*dtheta) - b*sin((i - 1)*dtheta);
    D[i] = D[i - 1] + sqrt(dx*dx + dy*dy);
  }

  BfReal Dmax = D[numPoints];

  /* TODO: bad O(n^2) algorithm here... fix */
  for (BfSize i = 0; i < numPoints; ++i) {
    BfReal d = (Dmax/numPoints)*i;

    BfSize j = 0;
    while (j < numPoints && !(D[j] <= d && d < D[j + 1])) ++j;

    BfReal lam = (d - D[j])/(D[j + 1] - D[j]);
    BfReal theta = (j + lam)*dtheta;

    BfReal px = a*cos(theta);
    BfReal py = b*sin(theta);
    BfPoint2 p = {px, py};
    bfPoint2RotateAboutOrigin(p, ellipse->theta);
    bfPoint2Translate(p, ellipse->center);
    if (points) bfPoints2Append(points, p);

    /* Compute unit tangent vector: */
    BfReal tx = -a*sin(theta);
    BfReal ty = b*cos(theta);
    BfVector2 t = {tx, ty};
    bfVector2Normalize(t);

    /* Compute unit normal vector: */
    BfReal nx = -a*cos(theta);
    BfReal ny = -b*sin(theta);
    BfVector2 n = {nx, ny};
    bfVector2Reject(n, t);
    bfVector2Normalize(n);
    bfVector2Rotate(n, ellipse->theta);
    bfVector2Negate(n); // get outward-facing normal
    if (unitNormals) bfVectors2Append(unitNormals, n);

    bfVector2Rotate(t, ellipse->theta);
    if (unitTangents) bfVectors2Append(unitTangents, t);
  }

  bfMemFree(D);
}

void bfEllipseSampleWithInverseCurvatureSpacing(BfEllipse const *ellipse, BfSize numPoints,
                                                BfPoints2 *points, BfVectors2 *unitNormals) {
  BfReal a = ellipse->semiMajorAxis;
  BfReal b = ellipse->semiMinorAxis;

  BfReal dt = BF_TWO_PI/numPoints;

  BfReal *D = bfMemAlloc((numPoints + 1), sizeof(BfReal));
  D[0] = 0;
  for (BfSize i = 1; i <= numPoints; ++i) {
    BfReal dx = a*cos(i*dt) - a*cos((i - 1)*dt);
    BfReal dy = b*sin(i*dt) - b*sin((i - 1)*dt);
    D[i] = D[i - 1] + sqrt(dx*dx + dy*dy);
  }
  BfReal Dmax = D[numPoints];

  BfReal *S = bfMemAlloc((numPoints + 1), sizeof(BfReal));
  S[0] = 0;
  for (BfSize i = 1; i <= numPoints; ++i) {
    BfReal t = i*dt;
    BfReal nx = -a*cos(t);
    BfReal ny = -b*sin(t);
    BfReal kappa = sqrt(nx*nx + ny*ny);
    BfReal rho = 1/kappa;
    S[i] = S[i - 1] + rho;
  }
  BfReal Smax = S[numPoints];

  BfReal c = cos(ellipse->theta);
  BfReal s = sin(ellipse->theta);

  /* TODO: bad O(n^2) algorithm here... fix */
  for (BfSize i = 0; i < numPoints; ++i) {
    BfReal d = Dmax*S[i]/Smax;

    BfSize j = 0;
    while (j < numPoints && !(D[j] <= d && d < D[j + 1])) ++j;

    BfReal lam = (d - D[j])/(D[j + 1] - D[j]);
    BfReal t = (j + lam)*dt;

    BfReal x = a*cos(t);
    BfReal y = b*sin(t);
    BfPoint2 p = {c*x - s*y + ellipse->center[0], s*x + c*y + ellipse->center[1]};
    bfPoints2Append(points, p);

    BfReal nx = -b*cos(t);
    BfReal ny = -a*sin(t);
    BfVector2 n = {c*nx - s*ny, s*nx + c*ny};
    bfVector2Normalize(n);
    bfVectors2Append(unitNormals, n);
  }

  bfMemFree(D);
  bfMemFree(S);
}
