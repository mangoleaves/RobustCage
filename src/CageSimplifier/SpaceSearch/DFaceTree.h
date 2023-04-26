#pragma once

#include "Geometry/DynamicAABB/DAABBTree.h"
#include "Mesh/SurfaceMesh/SurfaceMeshDefinition.h"
#include "Geometry/Exact/Types.h"
#include <utility>
#include <stack>

namespace Cage
{
namespace Geometry
{
using namespace SurfaceMesh;

// foward declaration
class DFaceTree;

class DProjectionTrait
{
private:
  typedef HeavyTriangle* HeavyTriPtr;
  typedef DFaceTree* TreePtr;
private:
  Vec3d m_query;
  Vec3d m_closest_point;
  double m_square_distance;
  int m_closest_triangle;
  TreePtr m_tree;
public:
  DProjectionTrait(TreePtr tree, const Vec3d& query);
  DProjectionTrait(TreePtr tree, const Vec3d& query, const Vec3d& hint, int hint_triangle);

  bool intersection(int tri_idx);
  bool do_inter(int node_idx) const;

  inline const Vec3d& closest_point() const { return m_closest_point; }
  inline int primitive()const { return m_closest_triangle; }
  inline double square_distance() const { return m_square_distance; }
};

template<typename DTree>
class DTriInterTrait
{
private:
  typedef DTree* TreePtr;
private:
  TreePtr m_tree;
  ExactTriangle<Triangle> m_query;
  BoundingBox m_box_of_query;
  std::vector<int> m_ignore_primitives;
  bool m_intersected;
public:
  /// @brief ATTENTION: make sure "ignore" is ordered, binary search will be used.
  DTriInterTrait(
    TreePtr tree,
    const Point& p, const Point& q, const Point& r,
    const ExactPoint* ep, const ExactPoint* eq, const ExactPoint* er,
    std::vector<int>&& ignore)
    :m_tree(tree), m_query(ep, eq, er, p, q, r), m_intersected(false), m_ignore_primitives(std::move(ignore))
  {
    m_box_of_query = CalcBoxForExactTriangle()(m_query);
  }
  /// @brief ATTENTION: make sure "ignore" is ordered, binary search will be used.
  DTriInterTrait(
    TreePtr tree,
    const Point& p, const Point& q, const Point& r,
    const ExactPoint* ep, const ExactPoint* eq, const ExactPoint* er,
    std::vector<int>& ignore)
    :m_tree(tree), m_query(ep, eq, er, p, q, r), m_intersected(false), m_ignore_primitives(ignore)
  {
    m_box_of_query = CalcBoxForExactTriangle()(m_query);
  }

  bool intersection(int tri_idx);
  bool do_inter(int node_idx)const;
  bool result()const { return m_intersected; }
};

class DFaceTreeKernel
{
public:
  typedef ExactTriangle<HeavyTriangle> Primitive;
  typedef ExactTriangleSplitPred SplitPred;
  typedef CalcBoxForExactTriangle CalcBox;
};

class DFaceTree : public DAABBTree<DFaceTreeKernel>
{
  friend class DTriInterTrait<DFaceTree>;
  friend class DProjectionTrait;
public:
  std::pair<Vec3d, int> closest_point(const Vec3d& query);
  std::pair<double, int> closest_distance(const Vec3d& query);
  std::pair<Vec3d, int> closest_point(const Vec3d& query, std::pair<Vec3d, int> hint);
  std::pair<double, int> closest_distance(const Vec3d& query, std::pair<Vec3d, int> hint);

  bool do_intersect(
    const Point& p, const Point& q, const Point& r,
    const ExactPoint* ep, const ExactPoint* eq, const ExactPoint* er
  );

  template<typename IgnoreT>
  bool do_intersect(
    const Point& p, const Point& q, const Point& r,
    const ExactPoint* ep, const ExactPoint* eq, const ExactPoint* er,
    IgnoreT&& ignore_triangles)
  {
    DTriInterTrait<DFaceTree> trait(this, p, q, r, ep, eq, er, std::forward<IgnoreT>(ignore_triangles));
    traversal(trait);
    return trait.result();
  }

  inline virtual void insert(const Primitive& p) { DAABBTree::insert(p); }
  inline virtual void remove(int p_idx) { DAABBTree::remove(p_idx); }
  inline virtual void update(int p_idx, const Primitive& new_p) { DAABBTree::update(p_idx, new_p); }
public:
  // interfaces for OpenMesh

  // constructor for OpenMesh
  DFaceTree(SMeshT& mesh);
  DFaceTree(SMeshT& mesh, const std::vector<SMeshT::FaceHandle>& faces);

  // modification functions
  virtual void insert(SMeshT* mesh, SMeshT::FaceHandle fh);
  virtual void remove(SMeshT::FaceHandle fh);
  virtual void update(SMeshT* mesh, SMeshT::FaceHandle fh);

  // closest search hint
  bool pre_hint_set = false;
  std::pair<Vec3d, int> pre_hint;
  void set_hint(const SMeshT::Point& query)
  {
    pre_hint = closest_point(query);
    pre_hint_set = true;
  }
  void unset_hint() { pre_hint_set = false; }

  // closest search
  std::pair<SMeshT::Point, SMeshT::FaceHandle>
    closest_point_and_face_handle(const SMeshT::Point& query);
  std::pair<double, SMeshT::FaceHandle>
    closest_distance_and_face_handle(const SMeshT::Point& query);
};

class LightDFaceTreeKernel
{
public:
  typedef ExactTriangle<IndexedTriangle> Primitive;
  typedef ExactTriangleSplitPred SplitPred;
  typedef CalcBoxForExactTriangle CalcBox;
};

class LightDFaceTree : public DAABBTree<LightDFaceTreeKernel>
{
  friend class DTriInterTrait<LightDFaceTree>;
public:
  bool do_intersect(
    const Point& p, const Point& q, const Point& r,
    const ExactPoint* ep, const ExactPoint* eq, const ExactPoint* er
  );

  template<typename IgnoreT>
  bool do_intersect(
    const Point& p, const Point& q, const Point& r,
    const ExactPoint* ep, const ExactPoint* eq, const ExactPoint* er,
    IgnoreT&& ignore_triangles)
  {
    DTriInterTrait<LightDFaceTree> trait(this, p, q, r, ep, eq, er, std::forward<IgnoreT>(ignore_triangles));
    traversal(trait);
    return trait.result();
  }

  inline virtual void insert(const Primitive& p) { DAABBTree::insert(p); }
  inline virtual void remove(int p_idx) { DAABBTree::remove(p_idx); }
  inline virtual void update(int p_idx, const Primitive& new_p) { DAABBTree::update(p_idx, new_p); }
public:
  // interfaces for OpenMesh

  // constructor for OpenMesh
  LightDFaceTree(SMeshT& mesh);
  LightDFaceTree(SMeshT& mesh, const std::vector<SMeshT::FaceHandle>& faces);

  // modification functions
  virtual void insert(SMeshT* mesh, SMeshT::FaceHandle fh);
  virtual void remove(SMeshT::FaceHandle fh);
  virtual void update(SMeshT* mesh, SMeshT::FaceHandle fh);
};
}// namespace Geometry
}// namespace Cage