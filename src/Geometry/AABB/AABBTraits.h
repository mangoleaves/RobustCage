#pragma once

#include "Geometry/Basic/Types.h"
#include "Geometry/Exact/TriTriIntersect.h"
#include "Geometry/KdTree.h"
#include "AABBTree.h"

namespace Cage
{
namespace Geometry
{
class ProjectionTraits
{
private:
  typedef HeavyTriangle* HeavyTriPtr;
private:
  Vec3d m_query;
  Vec3d m_closest_point;
  double m_square_distance;
  HeavyTriPtr m_closest_triangle;
public:
  ProjectionTraits(const Vec3d& query, const Vec3d& hint, HeavyTriPtr hint_triangle)
    :m_query(query),
    m_closest_point(hint),
    m_closest_triangle(hint_triangle),
    m_square_distance((query - hint).squaredNorm())
  {}

  inline bool intersection(HeavyTriangle& tri)
  {
    Vec3d new_closest_point;
    double new_square_distance = tri.closest_point(m_query, new_closest_point);
    if (new_square_distance < m_square_distance)
    {
      m_closest_point = new_closest_point;
      m_square_distance = new_square_distance;
      m_closest_triangle = &tri;
    }
    return true;
  }

  inline bool do_inter(const BoundingBox& bbox) const
  {
    return bbox.do_intersect(Sphere(m_query, m_square_distance));
  }

  inline const Vec3d& closest_point() const { return m_closest_point; }
  inline HeavyTriPtr primitive()const { return m_closest_triangle; }
  inline double square_distance() const { return m_square_distance; }
};

template<typename Primitive>
class BoxInterTraits
{
private:
  typedef IndexedTriangle* TriPtr;
  typedef std::vector<int> Indices;
private:
  Primitive m_query;
  BoundingBox m_box_of_query;
  Indices m_indices;
public:
  BoxInterTraits(const Primitive& query)
    :m_query(query)
  {
    m_box_of_query = BoundingBox(m_query);
  }

  inline bool intersection(IndexedTriangle& tri) { m_indices.push_back(tri.index); return true; }

  inline bool do_inter(const BoundingBox& bbox) const
  {
    if (bbox.do_intersect(m_box_of_query))
      return bbox.do_intersect(m_query);
    else
      return false;
  }

  inline const Indices& result() const { return m_indices; }
};

class TriInterTraits
{
private:
  typedef IndexedTriangle* TriPtr;
private:
  Triangle m_query;
  BoundingBox m_box_of_query;
  bool m_intersected;
public:
  TriInterTraits(const Point& _p, const Point& _q, const Point& _r)
    :m_query(_p, _q, _r), m_intersected(false)
  {
    m_box_of_query = BoundingBox(m_query);
  }

  inline bool intersection(IndexedTriangle& tri)
  {
    if (!m_intersected)
    {
      m_intersected = triangle_do_intersect(
        tri.ver0, tri.ver1, tri.ver2,
        m_query.ver0, m_query.ver1, m_query.ver2,
        nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
    }
    // if havn't found an intersection, continue search.
    return !m_intersected;
  }

  inline bool do_inter(const BoundingBox& bbox) const
  {
    if (bbox.do_intersect(m_box_of_query))
      return bbox.do_intersect(m_query);
    else
      return false;
  }

  inline bool result() const { return m_intersected; }
};
}// namespace Geometry
}// namespace Cage