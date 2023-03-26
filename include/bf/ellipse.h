#pragma once

#include "geom.h"

struct BfEllipse {
  BfReal semiMajorAxis;
  BfReal semiMinorAxis;
  BfPoint2 center;
  BfReal theta; /* rotation angle */
};

BfReal bfEllipseGetPerimeter(BfEllipse const *ellipse);
void bfEllipseSampleEquispaced(BfEllipse const *ellipse, BfSize numPoints, BfPoints2 *points, BfVectors2 *unitTangents, BfVectors2 *unitNormals);
void bfEllipseSampleWithInverseCurvatureSpacing(BfEllipse const *ellipse, BfSize numPoints, BfPoints2 *points, BfVectors2 *unitNormals);
