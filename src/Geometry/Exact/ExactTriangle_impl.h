#pragma

#include "ExactTriangle.h"
#include "IEEE754Union.h"

namespace Cage
{
namespace Geometry
{

template<typename BT>
ExactTriangle<BT>::ExactTriangle()
{
  ev0 = nullptr; ev1 = nullptr; ev2 = nullptr;
}

template<typename BT>
ExactTriangle<BT>::ExactTriangle(
  const Vec3d& fp0, const Vec3d& fp1, const Vec3d& fp2
)
{
  tri = BT(fp0, fp1, fp2);
  ev0 = nullptr; ev1 = nullptr; ev2 = nullptr;
}

template<>
inline ExactTriangle<Triangle>::ExactTriangle(
  unsigned idx, const Vec3d& fp0, const Vec3d& fp1, const Vec3d& fp2)
{
  tri = Triangle(fp0, fp1, fp2);
}

template<>
inline ExactTriangle<IndexedTriangle>::ExactTriangle(
  unsigned idx, const Vec3d& fp0, const Vec3d& fp1, const Vec3d& fp2)
{
  tri = IndexedTriangle(fp0, fp1, fp2, idx);
}

template<>
inline ExactTriangle<HeavyTriangle>::ExactTriangle(
  unsigned idx, const Vec3d& fp0, const Vec3d& fp1, const Vec3d& fp2)
{
  tri = HeavyTriangle(fp0, fp1, fp2, idx);
}

template<typename BT>
ExactTriangle<BT>::ExactTriangle(
  const ExactPoint* ep0, const ExactPoint* ep1, const ExactPoint* ep2,
  const Vec3d& fp0, const Vec3d& fp1, const Vec3d& fp2)
{
  tri = BT(fp0, fp1, fp2);
  ev0 = ep0; ev1 = ep1; ev2 = ep2;
}

template<>
inline ExactTriangle<Triangle>::ExactTriangle(
  unsigned idx,
  const ExactPoint* ep0, const ExactPoint* ep1, const ExactPoint* ep2,
  const Vec3d& fp0, const Vec3d& fp1, const Vec3d& fp2)
{
  tri = Triangle(fp0, fp1, fp2);
  ev0 = ep0; ev1 = ep1; ev2 = ep2;
}

template<>
inline ExactTriangle<IndexedTriangle>::ExactTriangle(
  unsigned idx,
  const ExactPoint* ep0, const ExactPoint* ep1, const ExactPoint* ep2,
  const Vec3d& fp0, const Vec3d& fp1, const Vec3d& fp2)
{
  tri = IndexedTriangle(fp0, fp1, fp2, idx);
  ev0 = ep0; ev1 = ep1; ev2 = ep2;
}

template<>
inline ExactTriangle<HeavyTriangle>::ExactTriangle(
  unsigned idx,
  const ExactPoint* ep0, const ExactPoint* ep1, const ExactPoint* ep2,
  const Vec3d& fp0, const Vec3d& fp1, const Vec3d& fp2)
{
  tri = HeavyTriangle(fp0, fp1, fp2, idx);
  switch (tri.obtuse_angle)
  {
  case -1:
  case 0:
    ev0 = ep0; ev1 = ep1; ev2 = ep2;
    break;
  case 1:
    ev0 = ep1; ev1 = ep2; ev2 = ep0;
    break;
  case 2:
    ev0 = ep2; ev1 = ep0; ev2 = ep1;
    break;
  }
}

template<typename BT>
bool ExactTriangle<BT>::allFloat()const
{
  return ev0 == nullptr && ev1 == nullptr && ev2 == nullptr;
}

/// @brief find min coords in exact points. 
/// NOTE: there are not always three exact points!
template<typename BT>
Point_3 ExactTriangle<BT>::minExactCoords()const
{
  ASSERT(!allFloat(), "No exact point.");
  Point_3 res;
  bool initialized = false;
  // process ev0
  if (ev0)
  {
    res = ev0->exact();
    initialized = true;
  }
  // process ev1
  if (ev1)
  {
    if (initialized)
    {
      minimize(res, ev1->exact());
    }
    else
    {
      res = ev1->exact();
      initialized = true;
    }
  }
  // process ev2
  if (ev2)
  {
    if (initialized)
    {
      minimize(res, ev2->exact());
    }
    else
    {
      res = ev2->exact();
      initialized = true;
    }
  }
  return res;
}

/// @brief find max coords in exact points. 
/// NOTE: there are not always three exact points!
template<typename BT>
Point_3 ExactTriangle<BT>::maxExactCoords()const
{
  ASSERT(!allFloat(), "No exact point.");
  Point_3 res;
  bool initialized = false;
  // process ev0
  if (ev0)
  {
    res = ev0->exact();
    initialized = true;
  }
  // process ev1
  if (ev1)
  {
    if (initialized)
    {
      maximize(res, ev1->exact());
    }
    else
    {
      res = ev1->exact();
      initialized = true;
    }
  }
  // process ev2
  if (ev2)
  {
    if (initialized)
    {
      maximize(res, ev2->exact());
    }
    else
    {
      res = ev2->exact();
      initialized = true;
    }
  }
  return res;
}

template<typename BT>
BoundingBox CalcBoxForExactTriangle::operator()(const ExactTriangle<BT>& exact_tri)const
{
  if (exact_tri.allFloat())
  {
    return BoundingBox(exact_tri.tri);
  }
  else
  {
    BoundingBox bbox;

    Point_3 exact_max = exact_tri.maxExactCoords();
    Point_3 exact_min = exact_tri.minExactCoords();

    bbox.max().x() = exact_max.approx().x().sup();
    bbox.max().y() = exact_max.approx().y().sup();
    bbox.max().z() = exact_max.approx().z().sup();

    bbox.min().x() = exact_min.approx().x().inf();
    bbox.min().y() = exact_min.approx().y().inf();
    bbox.min().z() = exact_min.approx().z().inf();

    if (exact_tri.ev0 == nullptr)
    {
      bbox.max().maximize(exact_tri.tri.ver0);
      bbox.min().minimize(exact_tri.tri.ver0);
    }
    if (exact_tri.ev1 == nullptr)
    {
      bbox.max().maximize(exact_tri.tri.ver1);
      bbox.min().minimize(exact_tri.tri.ver1);
    }
    if (exact_tri.ev2 == nullptr)
    {
      bbox.max().maximize(exact_tri.tri.ver2);
      bbox.min().minimize(exact_tri.tri.ver2);
    }

    return bbox;
  }
}

}// namespace Geometry
}// namespace Cage