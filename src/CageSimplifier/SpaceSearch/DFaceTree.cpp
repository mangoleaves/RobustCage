#include "DFaceTree.h"
#include "Geometry/Exact/TriTriIntersect.h"
#include "Utils/logger.hh"
#include <algorithm>

namespace Cage
{
namespace Geometry
{

/********************************/
/*       Projection Trait       */
/********************************/

DProjectionTrait::DProjectionTrait(TreePtr tree, const Vec3d& query)
  :m_tree(tree),
  m_query(query)
{
  if (!tree->m_primitives.empty())
  {
    m_closest_point = tree->m_primitives[0].tri.ver0;
    m_closest_triangle = 0;
    m_square_distance = (query - m_closest_point).squaredNorm();
  }
}

DProjectionTrait::DProjectionTrait(TreePtr tree, const Vec3d& query, const Vec3d& hint, int hint_triangle)
  :m_tree(tree),
  m_query(query),
  m_closest_point(hint),
  m_closest_triangle(hint_triangle),
  m_square_distance((query - hint).squaredNorm())
{}

bool DProjectionTrait::intersection(int tri_idx)
{
  Vec3d new_closest_point;
  auto& tri = m_tree->m_primitives[tri_idx];
  double new_square_distance = tri.tri.closest_point(m_query, new_closest_point);
  if (new_square_distance < m_square_distance)
  {
    m_closest_point = new_closest_point;
    m_square_distance = new_square_distance;
    m_closest_triangle = tri_idx;
  }
  return true;
}

bool DProjectionTrait::do_inter(int node_idx) const
{
  auto& bbox = m_tree->m_nodes[node_idx].m_bbox;
  return bbox.do_intersect(Sphere(m_query, m_square_distance));
}

/********************************/
/* Triangle Intersection Trait  */
/********************************/

template<typename DTree>
bool DTriInterTrait<DTree>::intersection(int tri_idx)
{
  if (!m_intersected &&
    !std::binary_search(m_ignore_primitives.begin(), m_ignore_primitives.end(), tri_idx))
  {
    auto& tri = m_tree->m_primitives[tri_idx];
    m_intersected = triangle_do_intersect(
      tri.tri.ver0, tri.tri.ver1, tri.tri.ver2,
      m_query.tri.ver0, m_query.tri.ver1, m_query.tri.ver2,
      tri.ev0, tri.ev1, tri.ev2,
      m_query.ev0, m_query.ev1, m_query.ev2);
  }
  return !m_intersected;
}

template<typename DTree>
bool DTriInterTrait<DTree>::do_inter(int node_idx)const
{
  auto& bbox = m_tree->m_nodes[node_idx].m_bbox;
  if (bbox.do_intersect(m_box_of_query))
    return bbox.do_intersect(m_query.tri);
  else
    return false;
}

/********************************/
/*       Dynamic Face Tree      */
/********************************/

std::pair<Vec3d, int> DFaceTree::closest_point(const Vec3d& query)
{
  DProjectionTrait projection_trait(this, query);
  this->traversal(projection_trait);
  return std::pair<Vec3d, int>(projection_trait.closest_point(), projection_trait.primitive());
}

std::pair<double, int> DFaceTree::closest_distance(const Vec3d& query)
{
  DProjectionTrait projection_trait(this, query);
  this->traversal(projection_trait);
  return std::pair<double, int>(std::sqrt(projection_trait.square_distance()), projection_trait.primitive());
}

std::pair<Vec3d, int> DFaceTree::closest_point(const Vec3d& query, std::pair<Vec3d, int> hint)
{
  DProjectionTrait projection_trait(this, query, hint.first, hint.second);
  this->traversal(projection_trait);
  return std::pair<Vec3d, int>(projection_trait.closest_point(), projection_trait.primitive());
}

std::pair<double, int> DFaceTree::closest_distance(const Vec3d& query, std::pair<Vec3d, int> hint)
{
  DProjectionTrait projection_trait(this, query, hint.first, hint.second);
  this->traversal(projection_trait);
  return std::pair<double, int>(std::sqrt(projection_trait.square_distance()), projection_trait.primitive());
}

bool DFaceTree::do_intersect(
  const Point& p, const Point& q, const Point& r,
  const ExactPoint* ep, const ExactPoint* eq, const ExactPoint* er)
{
  DTriInterTrait trait(this, p, q, r, ep, eq, er, std::vector<int>());
  traversal(trait);
  return trait.result();
}

DFaceTree::DFaceTree(SMeshT& mesh)
{
  Primitives triangles;
  triangles.reserve(mesh.n_faces());
  for (auto fh : mesh.faces())
  {
    auto fv_it = mesh.fv_begin(fh);
    auto& v0 = mesh.point(*fv_it);
    auto* ev0 = mesh.data(*fv_it).ep.get();
    fv_it++;
    auto& v1 = mesh.point(*fv_it);
    auto* ev1 = mesh.data(*fv_it).ep.get();
    fv_it++;
    auto& v2 = mesh.point(*fv_it);
    auto* ev2 = mesh.data(*fv_it).ep.get();
    triangles.emplace_back(fh.idx(), ev0, ev1, ev2, v0, v1, v2);
  }
  init_primitives(std::move(triangles));
  build();
  pre_hint_set = false;
}

DFaceTree::DFaceTree(SMeshT& mesh, const std::vector<SMeshT::FaceHandle>& faces)
{
  Primitives triangles;
  triangles.reserve(mesh.n_faces());
  for (auto fh : faces)
  {
    auto fv_it = mesh.fv_begin(fh);
    auto& v0 = mesh.point(*fv_it);
    auto* ev0 = mesh.data(*fv_it).ep.get();
    fv_it++;
    auto& v1 = mesh.point(*fv_it);
    auto* ev1 = mesh.data(*fv_it).ep.get();
    fv_it++;
    auto& v2 = mesh.point(*fv_it);
    auto* ev2 = mesh.data(*fv_it).ep.get();
    triangles.emplace_back(fh.idx(), ev0, ev1, ev2, v0, v1, v2);
  }
  init_primitives(std::move(triangles));
  build();
  pre_hint_set = false;
}

void DFaceTree::insert(SMeshT* mesh, SMeshT::FaceHandle fh)
{
  int fidx = fh.idx();
  auto fv_it = mesh->fv_begin(fh);
  auto& v0 = mesh->point(*fv_it);
  auto* ev0 = mesh->data(*fv_it).ep.get();
  fv_it++;
  auto& v1 = mesh->point(*fv_it);
  auto* ev1 = mesh->data(*fv_it).ep.get();
  fv_it++;
  auto& v2 = mesh->point(*fv_it);
  auto* ev2 = mesh->data(*fv_it).ep.get();
  Primitive tri(fh.idx(), ev0, ev1, ev2, v0, v1, v2);
  DAABBTree::insert(tri);
}

void DFaceTree::remove(SMeshT::FaceHandle fh)
{
  DAABBTree::remove(fh.idx());
}

void DFaceTree::update(SMeshT* mesh, SMeshT::FaceHandle fh)
{
  int fidx = fh.idx();
  auto fv_it = mesh->fv_begin(fh);
  auto& v0 = mesh->point(*fv_it);
  auto* ev0 = mesh->data(*fv_it).ep.get();
  fv_it++;
  auto& v1 = mesh->point(*fv_it);
  auto* ev1 = mesh->data(*fv_it).ep.get();
  fv_it++;
  auto& v2 = mesh->point(*fv_it);
  auto* ev2 = mesh->data(*fv_it).ep.get();
  Primitive tri(fh.idx(), ev0, ev1, ev2, v0, v1, v2);
  DAABBTree::update(fidx, tri);
}

std::pair<SMeshT::Point, SMeshT::FaceHandle>
DFaceTree::closest_point_and_face_handle(const SMeshT::Point& query)
{
  std::pair<Vec3d, int> result;
  if (pre_hint_set)
    result = closest_point(query, pre_hint);
  else
    result = closest_point(query);

  return std::pair<SMeshT::Point, SMeshT::FaceHandle>(
    result.first,
    SMeshT::FaceHandle(result.second));
}

std::pair<double, SMeshT::FaceHandle>
DFaceTree::closest_distance_and_face_handle(const SMeshT::Point& query)
{
  std::pair<double, int> result;
  if (pre_hint_set)
    result = closest_distance(query, pre_hint);
  else
    result = closest_distance(query);

  return std::pair<double, SMeshT::FaceHandle>(
    result.first,
    SMeshT::FaceHandle(result.second));
}

/********************************/
/*  Light  Dynamic Face Tree    */
/********************************/

bool LightDFaceTree::do_intersect(
  const Point& p, const Point& q, const Point& r,
  const ExactPoint* ep, const ExactPoint* eq, const ExactPoint* er)
{
  DTriInterTrait trait(this, p, q, r, ep, eq, er, std::vector<int>());
  traversal(trait);
  return trait.result();
}

LightDFaceTree::LightDFaceTree(SMeshT& mesh)
{
  Primitives triangles;
  triangles.reserve(mesh.n_faces());
  for (auto fh : mesh.faces())
  {
    auto fv_it = mesh.fv_begin(fh);
    auto& v0 = mesh.point(*fv_it);
    auto* ev0 = mesh.data(*fv_it).ep.get();
    fv_it++;
    auto& v1 = mesh.point(*fv_it);
    auto* ev1 = mesh.data(*fv_it).ep.get();
    fv_it++;
    auto& v2 = mesh.point(*fv_it);
    auto* ev2 = mesh.data(*fv_it).ep.get();
    triangles.emplace_back(fh.idx(), ev0, ev1, ev2, v0, v1, v2);
  }
  init_primitives(std::move(triangles));
  build();
}

LightDFaceTree::LightDFaceTree(SMeshT& mesh, const std::vector<SMeshT::FaceHandle>& faces)
{
  Primitives triangles;
  triangles.reserve(mesh.n_faces());
  for (auto fh : faces)
  {
    auto fv_it = mesh.fv_begin(fh);
    auto& v0 = mesh.point(*fv_it);
    auto* ev0 = mesh.data(*fv_it).ep.get();
    fv_it++;
    auto& v1 = mesh.point(*fv_it);
    auto* ev1 = mesh.data(*fv_it).ep.get();
    fv_it++;
    auto& v2 = mesh.point(*fv_it);
    auto* ev2 = mesh.data(*fv_it).ep.get();
    triangles.emplace_back(fh.idx(), ev0, ev1, ev2, v0, v1, v2);
  }
  init_primitives(std::move(triangles));
  build();
}

void LightDFaceTree::insert(SMeshT* mesh, SMeshT::FaceHandle fh)
{
  int fidx = fh.idx();
  auto fv_it = mesh->fv_begin(fh);
  auto& v0 = mesh->point(*fv_it);
  auto* ev0 = mesh->data(*fv_it).ep.get();
  fv_it++;
  auto& v1 = mesh->point(*fv_it);
  auto* ev1 = mesh->data(*fv_it).ep.get();
  fv_it++;
  auto& v2 = mesh->point(*fv_it);
  auto* ev2 = mesh->data(*fv_it).ep.get();
  Primitive tri(fh.idx(), ev0, ev1, ev2, v0, v1, v2);
  DAABBTree::insert(tri);
}

void LightDFaceTree::remove(SMeshT::FaceHandle fh)
{
  DAABBTree::remove(fh.idx());
}

void LightDFaceTree::update(SMeshT* mesh, SMeshT::FaceHandle fh)
{
  int fidx = fh.idx();
  auto fv_it = mesh->fv_begin(fh);
  auto& v0 = mesh->point(*fv_it);
  auto* ev0 = mesh->data(*fv_it).ep.get();
  fv_it++;
  auto& v1 = mesh->point(*fv_it);
  auto* ev1 = mesh->data(*fv_it).ep.get();
  fv_it++;
  auto& v2 = mesh->point(*fv_it);
  auto* ev2 = mesh->data(*fv_it).ep.get();
  Primitive tri(fh.idx(), ev0, ev1, ev2, v0, v1, v2);
  DAABBTree::update(fidx, tri);
}

}// namespace Geometry
}// namespace Cage