#pragma once

#include "Geometry/Basic/Types.h"
#include "Geometry/Exact/CGALTypes.h"
#include "implicit_point.h"
#include "Utils/logger.hh"

namespace Cage
{
namespace Geometry
{
using namespace SimpleUtils;

class FloatPoint;
class RationalPoint;

class ExactPoint
{
private:
  Point_3 rp;     // rational point, exact predicates use this.
  bool same;   // true if rp and fp are same.
public:
  ExactPoint() {}
  // use default copy/operator=/move
  ~ExactPoint() {}

  ExactPoint(double x, double y, double z);
  ExactPoint(const Vec3d& src);
  ExactPoint(const Point_3& src);
  ExactPoint(const Vector_3& src);

  const Point_3& exact() const { return rp; }
  Point_3& exact() { return rp; }

  Vector_3 exactVec3();

  Vec3d approx() const;

  bool isSame() const { return same; }

  Point_3 round()const;
  void rounded(Vec3d& fp);

  void fromImplicitPoint(genericPoint* gp);
};

void minimize(Point_3& lhs, const Point_3& rhs);
void maximize(Point_3& lhs, const Point_3& rhs);

}// namespace Geometry
}// namespace Cage