#pragma once

#include <type_traits>
#include "Geometry/Basic/Triangle.h"
#include "Geometry/Basic/BoundingBox.h"
#include "Geometry/Exact/ExactPoint.h"

namespace Cage
{
namespace Geometry
{

template<typename BaseTriangle>
class ExactTriangle
{
public:
  typedef BaseTriangle BT;
  typedef ExactPoint EP;

  BT tri;
  const ExactPoint* ev0, * ev1, * ev2;

  ExactTriangle();

  ExactTriangle(const Vec3d& fp0, const Vec3d& fp1, const Vec3d& fp2);

  ExactTriangle(unsigned idx, const Vec3d& fp0, const Vec3d& fp1, const Vec3d& fp2);

  ExactTriangle(
    const ExactPoint* ep0, const ExactPoint* ep1, const ExactPoint* ep2,
    const Vec3d& fp0, const Vec3d& fp1, const Vec3d& fp2);

  ExactTriangle(
    unsigned idx,
    const ExactPoint* ep0, const ExactPoint* ep1, const ExactPoint* ep2,
    const Vec3d& fp0, const Vec3d& fp1, const Vec3d& fp2);

  bool allFloat()const;
  Point_3 minExactCoords()const;
  Point_3 maxExactCoords()const;
};

class ExactTriangleSplitPred
{
private:
  TriangleSplitPred split_pred;
public:
  ExactTriangleSplitPred() :split_pred(0) {}
  ExactTriangleSplitPred(size_t sd) :split_pred(sd) {}

  template<typename BaseTri>
  inline bool operator()(
    const ExactTriangle<BaseTri>& lhs, const ExactTriangle<BaseTri>& rhs)
  {
    return split_pred(lhs.tri, rhs.tri);
  }
};

class CalcBoxForExactTriangle
{
public:
  template<typename BaseTri>
  BoundingBox operator()(const ExactTriangle<BaseTri>& exact_tri)const;
};

}// namespace Geometry
}// namespace Cage

#include "ExactTriangle_impl.h"