#include "VertexTree.h"

namespace Cage
{
namespace Geometry
{
VertexTree::VertexTree(SMeshT& mesh)
{
  std::vector<Vec3d> points;
  std::vector<unsigned int> ids;
  points.reserve(mesh.n_vertices());
  ids.reserve(mesh.n_vertices());

  for (auto vh : mesh.vertices())
  {
    points.emplace_back(mesh.point(vh));
    ids.emplace_back((unsigned)vh.idx());
  }
  KdTree::insert(points, ids);
  KdTree::build();
}

SMeshT::VertexHandle VertexTree::closest_vertex(const Vec3d& query)
{
  auto result = search_nearest_point(query);
  return SMeshT::VertexHandle(result.id);
}
}// namespace Geometry
}// namespace Cage