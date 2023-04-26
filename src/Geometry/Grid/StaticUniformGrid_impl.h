#include "StaticUniformGrid.h"
#include "omp.h"
#include <algorithm>
#include <map>

namespace Cage
{
namespace Geometry
{

template<typename Kernel>
StaticUniformGrid<Kernel>::StaticUniformGrid(PrimIter first, PrimIter beyond)
{
  insert(first, beyond);
  build();
}

template<typename Kernel>
void StaticUniformGrid<Kernel>::insert(PrimIter first, PrimIter beyond)
{
  clear();
  m_primitives.reserve(std::distance(first, beyond));
  while (first != beyond)
  {
    m_primitives.push_back(*first);
    first++;
  }
}

template<typename Kernel>
void StaticUniformGrid<Kernel>::insert(Primitives&& primitives)
{
  clear();
  m_primitives = std::move(primitives);
}

template<typename Kernel>
void StaticUniformGrid<Kernel>::build()
{
  // calculate box to bound all primitives
  m_prim_boxes.reserve(m_primitives.size());
  m_bbox = BoundingBox();
  auto calc_box = CalcBox();
  for (const Primitive& primitive : m_primitives)
  {
    m_prim_boxes.push_back(calc_box(primitive));
    m_bbox += m_prim_boxes.back();
  }
  // enlarge box a little
  double bdl = (m_bbox.max() - m_bbox.min()).length();
  m_bbox.enlarge(Vec3d(bdl, bdl, bdl) / (double)m_primitives.size());

  m_bbox_diag = m_bbox.max() - m_bbox.min();
  // partition box to grid
  best_grid_partition();

  // calculate voxel
  m_voxel.x() = m_bbox_diag.x() / m_grid_size.x();
  m_voxel.y() = m_bbox_diag.y() / m_grid_size.y();
  m_voxel.z() = m_bbox_diag.z() / m_grid_size.z();
  m_voxel_length = m_voxel.length();

  // allocate grid (one more tail grid)
  m_grid.resize(m_grid_size.x() * m_grid_size.y() * m_grid_size.z() + 1, nullptr);

  // insert all the primitives into the grid
  m_links.clear();
  m_links.reserve(m_primitives.size() + 1);
  for (int i = 0;i < m_primitives.size();i++)
  {
    const BoundingBox& bb = m_prim_boxes[i];
    // may need check whether bb intersects m_bbox
    BoundingBox3i bbi;
    Box2IBox(bb, bbi);
    int x, y, z;
    for (z = bbi.min().z();z <= bbi.max().z();++z)
    {
      int bz = z * m_grid_size.y();
      for (y = bbi.min().y();y <= bbi.max().y();++y)
      {
        int by = (y + bz) * m_grid_size.x();
        for (x = bbi.min().x();x <= bbi.max().x();++x)
          m_links.emplace_back(i, by + x);
      }
    }
  }

  m_links.emplace_back(-1, int(m_grid.size() - 1));

  std::sort(m_links.begin(), m_links.end());

  // build mapping relation
  auto link_iter = m_links.begin();
  for (int grid_idx = 0;grid_idx < m_grid.size();grid_idx++)
  {
    m_grid[grid_idx] = &(*link_iter);
    while (grid_idx == link_iter->gidx)
    {
      link_iter++;
      if (link_iter == m_links.end())
        break;
    }
  }

  // for better cache performance, sort primitives
  Primitives sorted_primitives; sorted_primitives.reserve(m_primitives.size());
  std::vector<BoundingBox> sorted_boxes; sorted_boxes.reserve(m_primitives.size());
  std::vector<int> new_idx(m_primitives.size(), -1);
  for (Link& link : m_links)
  {
    if (link.pidx >= 0)
    {
      if (new_idx[link.pidx] == -1)
      {
        sorted_primitives.push_back(m_primitives[link.pidx]);
        sorted_boxes.push_back(m_prim_boxes[link.pidx]);
        new_idx[link.pidx] = (int)sorted_primitives.size() - 1;
      }
      link.pidx = new_idx[link.pidx];
    }
  }
  m_primitives = std::move(sorted_primitives);
  m_prim_boxes = std::move(sorted_boxes);

  markers.resize(omp_get_max_threads(), Marker(m_primitives.size()));
}

template<typename Kernel>
void StaticUniformGrid<Kernel>::clear()
{
  m_primitives.clear();
  m_prim_boxes.clear();
}

template<typename Kernel>
void StaticUniformGrid<Kernel>::best_grid_partition()
{
  const size_t primitive_size = m_primitives.size();

  ASSERT(primitive_size > 0, "primitive size must be positive.");

  double box_diag_length = m_bbox_diag.length();
  double eps = box_diag_length * 1e-4;

  m_grid_size = Vec3i(1, 1, 1);
  if (m_bbox_diag.x() > eps)
  {
    if (m_bbox_diag.y() > eps)
    {
      if (m_bbox_diag.z() > eps)
      {
        double k = pow(primitive_size / (m_bbox_diag.x() * m_bbox_diag.y() * m_bbox_diag.z()), 1.0 / 3.0);
        m_grid_size.x() = int(m_bbox_diag.x() * k);
        m_grid_size.y() = int(m_bbox_diag.y() * k);
        m_grid_size.z() = int(m_bbox_diag.z() * k);
      }
      else
      {
        m_grid_size.x() = int(::sqrt(primitive_size * m_bbox_diag.x() / m_bbox_diag.y()));
        m_grid_size.y() = int(::sqrt(primitive_size * m_bbox_diag.y() / m_bbox_diag.x()));
      }
    }
    else
    {
      if (m_bbox_diag.z() > eps)
      {
        m_grid_size.x() = int(::sqrt(primitive_size * m_bbox_diag.x() / m_bbox_diag.z()));
        m_grid_size.z() = int(::sqrt(primitive_size * m_bbox_diag.z() / m_bbox_diag.x()));
      }
      else
        m_grid_size.x() = int(primitive_size);
    }
  }
  else
  {
    if (m_bbox_diag.y() > eps)
    {
      if (m_bbox_diag.z() > eps)
      {
        m_grid_size.y() = int(::sqrt(primitive_size * m_bbox_diag.y() / m_bbox_diag.z()));
        m_grid_size.z() = int(::sqrt(primitive_size * m_bbox_diag.z() / m_bbox_diag.y()));
      }
      else
        m_grid_size.y() = int(primitive_size);
    }
    else if (m_bbox_diag.z() > eps)
    {
      m_grid_size.z() = int(primitive_size);

    }
  }
  m_grid_size.x() = std::max(m_grid_size.x(), 1);
  m_grid_size.y() = std::max(m_grid_size.y(), 1);
  m_grid_size.z() = std::max(m_grid_size.z(), 1);
}

template<typename Kernel>
bool StaticUniformGrid<Kernel>::is_in_ex(const BoundingBox& box, const Vec3d& p)
{
  return box.min().x() <= p.x() && box.min().y() <= p.y() && box.min().z() <= p.z()
    && box.max().x() > p.x() && box.max().y() > p.y() && box.max().z() > p.z();
}

template<typename Kernel>
void StaticUniformGrid<Kernel>::find_grid(const Vec3i& ip, Cell& first, Cell& last)
{
  Cell* g = &m_grid[0] + (ip.x() + m_grid_size.x() * (ip.y() + m_grid_size.y() * ip.z()));
  first = *g;
  last = *(g + 1);
}

template<typename Kernel>
void StaticUniformGrid<Kernel>::P2IP(const Vec3d& p, Vec3i& pi)
{
  Vec3d t = p - m_bbox.min();
  pi.x() = int(t.x() / m_voxel.x());
  pi.y() = int(t.y() / m_voxel.y());
  pi.z() = int(t.z() / m_voxel.z());
}

template<typename Kernel>
void StaticUniformGrid<Kernel>::Box2IBox(const BoundingBox& b, BoundingBox3i& bi)
{
  P2IP(b.min(), bi.min());
  P2IP(b.max(), bi.max());
}

}// namespace Geometry
}// namespace Cage