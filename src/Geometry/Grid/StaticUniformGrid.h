#pragma once

#include "Geometry/Basic/BoundingBox.h"
#include "Geometry/Basic/BoundingBox3i.h"
#include "Utils/Marker.h"

namespace Cage
{
namespace Geometry
{
using namespace SimpleUtils;

/// @brief static uniform grid used for space searching.
/// @copyright refer to VCG library.  Visual Computing Laboratory;ISTI-CNR;Pisa.
template<typename Kernel>
class StaticUniformGrid
{
protected:
  typedef typename Kernel::Primitive Primitive;
  typedef typename Kernel::CalcBox CalcBox;

  typedef Primitive* PrimPtr;
  typedef std::vector<Primitive> Primitives;
  typedef typename Primitives::iterator PrimIter;

  class Link
  {
  public:
    Link() :pidx(-1), gidx(-1) {}
    Link(int _pidx, int _gidx) :pidx(_pidx), gidx(_gidx) {}
  public:
    inline bool operator<(const Link& rhs) const { return gidx < rhs.gidx; }
    inline bool operator<=(const Link& rhs) const { return gidx <= rhs.gidx; }
    inline bool operator>(const Link& rhs) const { return gidx > rhs.gidx; }
    inline bool operator>=(const Link& rhs) const { return gidx >= rhs.gidx; }
    inline bool operator==(const Link& rhs) const { return gidx == rhs.gidx; }
    inline bool operator!=(const Link& rhs) const { return gidx != rhs.gidx; }
  public:
    int pidx; // primitive index
    int gidx; // grid index
  };

  typedef Link* Cell;

public:
  StaticUniformGrid() = default;
  StaticUniformGrid(PrimIter first, PrimIter beyond);
public:
  void insert(PrimIter first, PrimIter beyond);
  void insert(Primitives&& primitives);
  void build();
  void clear();
protected:
  Primitives m_primitives;
  std::vector<BoundingBox> m_prim_boxes;

  BoundingBox m_bbox; // bound all primitives
  Vec3d m_bbox_diag;  // diagonal vector of m_bbox

  Vec3i m_grid_size;  // grid size on three dimension

  Vec3d m_voxel;      // the smallest unit of grid, m_voxel.x = m_bbox_diag.x / m_grid_size.x
  double m_voxel_length;

  std::vector<Link> m_links;
  std::vector<Cell> m_grid;

  std::vector<Marker> markers;  // for parallel searching, one marker for on thread.

  void best_grid_partition();

  bool is_in_ex(const BoundingBox& box, const Vec3d& p);

  void find_grid(const Vec3i& ip, Cell& first, Cell& last);

  void P2IP(const Vec3d& p, Vec3i& pi);
  void Box2IBox(const BoundingBox& b, BoundingBox3i& bi);
};

}// namespace Geometry
}// namespace Cage

#include "StaticUniformGrid_impl.h"