#pragma once

#include "Geometry/Exact/Types.h"
#include "Geometry/AABB/AABBTraits.h"
#include "Mesh/SurfaceMesh/SurfaceMeshDefinition.h"

namespace Cage
{
namespace Geometry
{
using namespace SurfaceMesh;

class FaceTreeKernel
{
public:
  typedef HeavyTriangle Primitive;
  typedef TriangleSplitPred SplitPred;
  
  class CalcBox
  {
  public:
    template<typename Primitive>
    BoundingBox operator()(const Primitive& p){return BoundingBox(p);}
  };
};

/// Closest point search and intersection check.
class FaceTree : public AABBTree<FaceTreeKernel>
{
public:
  typedef std::vector<HeavyTriangle> HeavyTriangles;
  typedef HeavyTriangles::iterator HeavyTriIter;
  typedef HeavyTriangle* HeavyTriPtr;
public:
  std::unique_ptr<KdTree> m_p_search_tree;

  FaceTree() = default;
  FaceTree(HeavyTriangles&& primitives)
  {
    insert(std::move(primitives));
    build();
  }
  FaceTree(HeavyTriIter first, HeavyTriIter beyond)
  {
    insert(first, beyond);
    build();
  }

  void build();
  void build_kd_tree();

  void clear();
  void clear_search_tree();

  std::pair<Vec3d, HeavyTriPtr> best_hint(const Vec3d& query);
  std::pair<Vec3d, HeavyTriPtr> closest_point(const Vec3d& query, const std::pair<Vec3d, HeavyTriPtr>& hint);
  std::pair<double, HeavyTriPtr> closest_distance(const Vec3d& query, const std::pair<Vec3d, HeavyTriPtr>& hint);

  bool do_intersect(const Point& p, const Point& q, const Point& r);
public:
  // interfaces for openmesh
  bool pre_hint_set = false;
  std::pair<Vec3d, HeavyTriPtr> pre_hint;
  void set_hint(const SMeshT::Point& query)
  {
    pre_hint = best_hint(query);
    pre_hint_set = true;
  }
  void unset_hint() { pre_hint_set = false; }


  FaceTree(SMeshT& mesh);
  FaceTree(SMeshT& mesh, const std::vector<SMeshT::FaceHandle>& faces);

  std::pair<SMeshT::Point, SMeshT::FaceHandle>
    closest_point_and_face_handle(const SMeshT::Point& query);
  std::pair<double, SMeshT::FaceHandle>
    closest_distance_and_face_handle(const SMeshT::Point& query);
  SMeshT::Point closest_point(const SMeshT::Point& query);
  double  closest_distance(const SMeshT::Point& query);
};
}// namespace Geometry
}// namespace Cage