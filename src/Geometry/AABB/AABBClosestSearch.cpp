#include "AABBClosestSearch.h"

namespace Cage
{
namespace Geometry
{
void AABBClosestSearch::build()
{
  AABBTree::build();
  if (m_primitives.size() > 1)
  {
    // build search tree
    build_kd_tree();
  }
}

void AABBClosestSearch::build_kd_tree()
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
  m_p_search_tree = new KdTree(points, ids);
  if (m_p_search_tree == nullptr)
  {
    clear();
    throw std::bad_alloc();
  }
}

void AABBClosestSearch::clear()
{
  AABBTree::clear();
  clear_search_tree();
}

void AABBClosestSearch::clear_search_tree()
{
  if (m_p_search_tree)
  {
    delete m_p_search_tree;
    m_p_search_tree = nullptr;
  }
}

std::pair<Vec3d, AABBClosestSearch::HeavyTriPtr>
AABBClosestSearch::best_hint(const Vec3d& query)
{
  auto result = m_p_search_tree->search_nearest_point(query);
  return std::make_pair(result.vec, &m_primitives[result.id]);
}

std::pair<Vec3d, const AABBClosestSearch::HeavyTriPtr>
AABBClosestSearch::closest_point(const Vec3d& query, const std::pair<Vec3d, HeavyTriPtr>& hint)
{
  ProjectionTraits projection_traits(query, hint.first, hint.second);
  this->traversal(projection_traits);
  return std::pair<Vec3d, const HeavyTriPtr>(projection_traits.closest_point(), projection_traits.primitive());
}

std::pair<double, const AABBClosestSearch::HeavyTriPtr>
AABBClosestSearch::closest_distance(const Vec3d& query, const std::pair<Vec3d, HeavyTriPtr>& hint)
{
  ProjectionTraits projection_traits(query, hint.first, hint.second);
  this->traversal(projection_traits);
  return std::pair<double, const HeavyTriPtr>(std::sqrt(projection_traits.square_distance()), projection_traits.primitive());
}
}// namespace Geometry
}// namespace Cage