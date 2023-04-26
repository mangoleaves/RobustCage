#include "Predicates.h"

namespace Cage
{
namespace Geometry
{

/*********************************/
/********* Orientation ***********/
/*********************************/

int orient3d_wrapped(const Point& p, const Point& q, const Point& r, const Point& query)
{
  return orient3d(
    p.x(), p.y(), p.z(),
    query.x(), query.y(), query.z(),
    q.x(), q.y(), q.z(),
    r.x(), r.y(), r.z()
  );
}

int orient3d_wrapped(const Triangle& tri, const Point& point)
{
  return orient3d(
    tri.ver0.x(), tri.ver0.y(), tri.ver0.z(),
    point.x(), point.y(), point.z(),
    tri.ver1.x(), tri.ver1.y(), tri.ver1.z(),
    tri.ver2.x(), tri.ver2.y(), tri.ver2.z()
  );
}

bool are_points_colinear(const Point& p0, const Point& p1, const Point& p2)
{
  return orient2d(p0.x(), p0.y(), p1.x(), p1.y(), p2.x(), p2.y()) == 0 &&
    orient2d(p0.y(), p0.z(), p1.y(), p1.z(), p2.y(), p2.z()) == 0 &&
    orient2d(p0.x(), p0.z(), p1.x(), p1.z(), p2.x(), p2.z()) == 0;
}

bool are_points_colinear(
  const Point& p0, const Point& p1, const Point& p2,
  const ExactPoint* ep0, const ExactPoint* ep1, const ExactPoint* ep2)
{
  if (ep0 || ep1 || ep2)
  {
    Point_3 pp0 = ep0 ? ep0->exact() : Point_3(p0.x(), p0.y(), p0.z());
    Point_3 pp1 = ep1 ? ep1->exact() : Point_3(p1.x(), p1.y(), p1.z());
    Point_3 pp2 = ep2 ? ep2->exact() : Point_3(p2.x(), p2.y(), p2.z());

    Triangle_3 tri(pp0, pp1, pp2);
    return K().is_degenerate_3_object()(tri);
  }
  else
  {
    return are_points_colinear(p0, p1, p2);
  }
}

int coplanar_orient2d(const Point& p0, const Point& p1, const Point& p2)
{
  int oxy_pqr = orient2d(p0.x(), p0.y(), p1.x(), p1.y(), p2.x(), p2.y());
  if (oxy_pqr != ZERO)
    return oxy_pqr;

  int oyz_pqr = orient2d(p0.y(), p0.z(), p1.y(), p1.z(), p2.y(), p2.z());
  if (oyz_pqr != ZERO)
    return oyz_pqr;

  return orient2d(p0.x(), p0.z(), p1.x(), p1.z(), p2.x(), p2.z());
}

/*********************************/
/********** Degenerate ***********/
/*********************************/

TriDeType triangle_degenerate_type(const Triangle& triangle)
{
  if (!are_points_colinear(triangle.ver0, triangle.ver1, triangle.ver2))
    return TriDeType::NOT_DE;
  else if (triangle[0] == triangle[1] && triangle[1] == triangle[2])
    return TriDeType::POINT_DE;
  else
    return TriDeType::SEG_DE;
}

/*********************************/
/*********** Accurate ************/
/*********************************/

Vec3d accurate_triangle_normal(const Triangle& triangle)
{
  double nvx, nvy, nvz;
  const double* t[3] = {
    triangle[0].data(),
    triangle[1].data(),
    triangle[2].data()
  };

  triangle_normal_exact(
    t[0][0], t[0][1], t[0][2],
    t[1][0], t[1][1], t[1][2],
    t[2][0], t[2][1], t[2][2],
    nvx, nvy, nvz
  );

  return Vec3d(nvx, nvy, nvz);
}

Vec3d accurate_cross_normalized(const Vec3d firstBegin, const Vec3d firstEnd, const Vec3d secondBegin, const Vec3d secondEnd)
{
  double nvx, nvy, nvz;
  const double* fb = firstBegin.data();
  const double* fe = firstEnd.data();
  const double* sb = secondBegin.data();
  const double* se = secondEnd.data();

  cross_product_normalized_exact(
    fb[0], fb[1], fb[2],
    fe[0], fe[1], fe[2],
    sb[0], sb[1], sb[2],
    se[0], se[1], se[2],
    nvx, nvy, nvz
  );

  return Vec3d(nvx, nvy, nvz);
}

}// namespace Geometry
}// namespace Cage