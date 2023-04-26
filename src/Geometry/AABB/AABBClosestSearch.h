#pragma once

#include "AABBTraits.h"

namespace Cage
{
namespace Geometry
{
class AABBClosestSearchKernel
{
public:
  typedef HeavyTriangle Primitive;
  typedef TriangleSplitPred SplitPred;
  typedef DefaultCalcBox CalcBox;
};

class AABBClosestSearch : public AABBTree<AABBClosestSearchKernel>
{
public:
  typedef std::vector<HeavyTriangle> HeavyTriangles;
  typedef HeavyTriangles::iterator HeavyTriIter;
  typedef HeavyTriangle* HeavyTriPtr;
public:
  KdTree* m_p_search_tree = nullptr;

  void build();
  void build_kd_tree();

  void clear();
  void clear_search_tree();

  std::pair<Vec3d, HeavyTriPtr> best_hint(const Vec3d& query);
  std::pair<Vec3d, const HeavyTriPtr> closest_point(const Vec3d& query, const std::pair<Vec3d, HeavyTriPtr>& hint);
  std::pair<double, const HeavyTriPtr> closest_distance(const Vec3d& query, const std::pair<Vec3d, HeavyTriPtr>& hint);
};
}
// namespace Geometry
}// namespace Cage