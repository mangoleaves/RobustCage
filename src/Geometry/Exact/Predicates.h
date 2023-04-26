#pragma once

#include "Geometry/Basic/Types.h"
#include "Geometry/Basic/Triangle.h"
#include "Geometry/Exact/ExactPoint.h"
#include "indirect_predicates.h"
#include "ip_filtered_ex.h"

namespace Cage
{
namespace Geometry
{
/*********************************/
/********* Orientation ***********/
/*********************************/

int orient3d_wrapped(const Point& p, const Point& q, const Point& r, const Point& query);

int orient3d_wrapped(const Triangle& tri, const Point& point);

bool are_points_colinear(const Point& p0, const Point& p1, const Point& p2);

bool are_points_colinear(
  const Point& p0, const Point& p1, const Point& p2,
  const ExactPoint* ep0, const ExactPoint* ep1, const ExactPoint* ep2);

int coplanar_orient2d(const Point& p0, const Point& p1, const Point& p2);

/*********************************/
/************* Sign **************/
/*********************************/

inline int dot_sign(const Vec3d vec1, const Vec3d vec2)
{
  return dot_product_sign(
    vec1.x(), vec1.y(), vec1.z(),
    vec2.x(), vec2.y(), vec2.z());
}

/*********************************/
/********** Degenerate ***********/
/*********************************/

enum class TriDeType
{
  NOT_DE,
  SEG_DE,
  POINT_DE
};

TriDeType triangle_degenerate_type(const Triangle& triangle);

/*********************************/
/*********** Accurate ************/
/*********************************/

Vec3d accurate_triangle_normal(const Triangle& triangle);

Vec3d accurate_cross_normalized(const Vec3d firstBegin, const Vec3d firstEnd, const Vec3d secondBegin, const Vec3d secondEnd);
}// namespace Geometry
}// namespace Cage