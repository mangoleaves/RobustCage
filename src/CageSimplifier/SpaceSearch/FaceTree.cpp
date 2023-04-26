#include "FaceTree.h"

namespace Cage
{
namespace Geometry
{
void FaceTree::build()
{
  AABBTree::build();
  if (m_primitives.size() > 1)
  {
    // build search tree
    build_kd_tree();
  }
}

void FaceTree::build_kd_tree()
{
  std::vector<Vec3d> points;
  std::vector<unsigned int> ids;
  points.reserve(m_primitives.size());
  ids.reserve(m_primitives.size());
  for (unsigned int i = 0;i < m_primitives.size();i++)
  {
    points.emplace_back(m_primitives[i].ver0);
    ids.emplace_back(i);
  }
  m_p_search_tree = std::make_unique<KdTree>(points, ids);
  if (m_p_search_tree == nullptr)
  {
    clear();
    throw std::bad_alloc();
  }
}

void FaceTree::clear()
{
  AABBTree::clear();
  clear_search_tree();
}

void FaceTree::clear_search_tree()
{
  if (m_p_search_tree)
  {
    m_p_search_tree = nullptr;
  }
}

std::pair<Vec3d, FaceTree::HeavyTriPtr>
FaceTree::best_hint(const Vec3d& query)
{
  auto result = m_p_search_tree->search_nearest_point(query);
  return std::make_pair(result.vec, &m_primitives[result.id]);
}

std::pair<Vec3d, FaceTree::HeavyTriPtr>
FaceTree::closest_point(const Vec3d& query, const std::pair<Vec3d, HeavyTriPtr>& hint)
{
  ProjectionTraits projection_traits(query, hint.first, hint.second);
  this->traversal(projection_traits);
  return std::pair<Vec3d, HeavyTriPtr>(projection_traits.closest_point(), projection_traits.primitive());
}

std::pair<double, FaceTree::HeavyTriPtr>
FaceTree::closest_distance(const Vec3d& query, const std::pair<Vec3d, HeavyTriPtr>& hint)
{
  ProjectionTraits projection_traits(query, hint.first, hint.second);
  this->traversal(projection_traits);
  return std::pair<double, HeavyTriPtr>(std::sqrt(projection_traits.square_distance()), projection_traits.primitive());
}

/// @brief do intersect with triangle pqr
bool FaceTree::do_intersect(const Point& p, const Point& q, const Point& r)
{
  TriInterTraits trait(p, q, r);
  traversal(trait);
  return trait.result();
}

/****** interfaces for OpenMesh *****/

FaceTree::FaceTree(SMeshT& mesh)
{
  std::vector<HeavyTriangle> triangles;
  triangles.reserve(mesh.n_faces());
  for (auto fh : mesh.faces())
  {
    auto fv_it = mesh.fv_begin(fh);
    auto& v0 = mesh.point(*fv_it); fv_it++;
    auto& v1 = mesh.point(*fv_it); fv_it++;
    auto& v2 = mesh.point(*fv_it);
    triangles.emplace_back(v0, v1, v2, fh.idx());
  }
  insert(std::move(triangles));
  build();
  pre_hint_set = false;
}

FaceTree::FaceTree(SMeshT& mesh, const std::vector<SMeshT::FaceHandle>& faces)
{
  std::vector<HeavyTriangle> triangles;
  triangles.reserve(faces.size());
  for (auto fh : faces)
  {
    auto fv_it = mesh.fv_begin(fh);
    auto& v0 = mesh.point(*fv_it); fv_it++;
    auto& v1 = mesh.point(*fv_it); fv_it++;
    auto& v2 = mesh.point(*fv_it);
    triangles.emplace_back(v0, v1, v2, fh.idx());
  }
  insert(std::move(triangles));
  build();
  pre_hint_set = false;
}


std::pair<SMeshT::Point, SMeshT::FaceHandle>
FaceTree::closest_point_and_face_handle(const SMeshT::Point& query)
{
  std::pair<Vec3d, HeavyTriPtr> result;
  if (pre_hint_set)
    result = closest_point(query, pre_hint);
  else
    result = closest_point(query, best_hint(query));

  return std::pair<SMeshT::Point, SMeshT::FaceHandle>(
    result.first,
    SMeshT::FaceHandle(result.second->index)
    );
}

std::pair<double, SMeshT::FaceHandle>
FaceTree::closest_distance_and_face_handle(const SMeshT::Point& query)
{
  std::pair<double, HeavyTriPtr> result;
  if (pre_hint_set)
    result = closest_distance(query, pre_hint);
  else
    result = closest_distance(query, best_hint(query));

  return std::pair<double, SMeshT::FaceHandle>(
    result.first,
    SMeshT::FaceHandle(result.second->index)
    );
}

SMeshT::Point FaceTree::closest_point(const SMeshT::Point& query)
{
  std::pair<Vec3d, HeavyTriPtr> result;
  if (pre_hint_set)
    result = closest_point(query, pre_hint);
  else
    result = closest_point(query, best_hint(query));

  return result.first;
}

double  FaceTree::closest_distance(const SMeshT::Point& query)
{
  std::pair<double, HeavyTriPtr> result;
  if (pre_hint_set)
    result = closest_distance(query, pre_hint);
  else
    result = closest_distance(query, best_hint(query));

  return result.first;
}
}// namespace Geometry
}// namespace Cage