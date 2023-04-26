#include "FaceGrid.h"
#include "omp.h"

namespace Cage
{
namespace Geometry
{
std::pair<std::pair<Vec3d, double>, int> FaceGrid::closest_point(const Vec3d& query)
{
  double closest_distance_sqr = DBL_MAX;
  PrimPtr winner = nullptr;

  Marker& marker = markers[omp_get_thread_num()];
  marker.unmark_all();
  double new_radius_sqr = m_voxel_length * m_voxel_length;
  double radius_sqr;
  BoundingBox3i ibox_done, ibox_todo;
  Vec3d result;
  Cell first, last, l;

  if (is_in_ex(m_bbox, query))
  {
    Vec3i ip;
    P2IP(query, ip);
    find_grid(ip, first, last);
    for (l = first;l != last;++l) if (l->pidx >= 0)
    {
      Primitive& elem = m_primitives[l->pidx];
      Vec3d cp = elem.tri.closest_point(query);
      double distance_sqr = (cp - query).sqrnorm();
      if (distance_sqr < closest_distance_sqr)
      {
        winner = &elem;
        result = cp;
        closest_distance_sqr = distance_sqr;
        new_radius_sqr = closest_distance_sqr;
      }
      marker.mark(l->pidx);
    }
    ibox_done.min() = ip;
    ibox_done.max() = ip;
  }

  int ix, iy, iz;
  BoundingBox3i ibox(Vec3i(0, 0, 0), m_grid_size - Vec3i(1, 1, 1));
  do
  {
    radius_sqr = new_radius_sqr;
    BoundingBox box_todo(query);
    box_todo.enlarge(Vec3d(std::sqrt(radius_sqr), std::sqrt(radius_sqr), std::sqrt(radius_sqr)));
    Box2IBox(box_todo, ibox_todo);
    if (ibox_todo.do_intersect(ibox))
    {
      ibox_todo.intersect(ibox);
      for (ix = ibox_todo.min().x(); ix <= ibox_todo.max().x(); ix++)
      {
        for (iy = ibox_todo.min().y(); iy <= ibox_todo.max().y(); iy++)
        {
          for (iz = ibox_todo.min().z(); iz <= ibox_todo.max().z(); iz++)
          {
            // this test is to avoid to re-process already analyzed cells.
            if (ix<ibox_done.min().x() || ix>ibox_done.max().x() ||
              iy<ibox_done.min().y() || iy>ibox_done.max().y() ||
              iz<ibox_done.min().z() || iz>ibox_done.max().z())
            {
              find_grid(Vec3i(ix, iy, iz), first, last);
              for (l = first;l != last;++l) if (l->pidx >= 0)
              {
                Primitive& elem = m_primitives[l->pidx];
                if (!marker.is_marked(l->pidx))
                {
                  Vec3d cp = elem.tri.closest_point(query);
                  double distance_sqr = (cp - query).sqrnorm();
                  if (distance_sqr < closest_distance_sqr)
                  {
                    winner = &elem;
                    result = cp;
                    closest_distance_sqr = distance_sqr;
                  }
                  marker.mark(l->pidx);
                }
              }
            }
          }
        }
      }
    }
    if (!winner) new_radius_sqr = radius_sqr + m_voxel_length;
    else new_radius_sqr = closest_distance_sqr;
    ibox_done = ibox_todo;
  } while (closest_distance_sqr > radius_sqr);

  return std::make_pair(std::make_pair(result, std::sqrt(closest_distance_sqr)), winner->tri.index);
}

bool FaceGrid::do_intersect(const ExactTriangle<Triangle>& tri)
{
  Marker& marker = markers[omp_get_thread_num()];
  BoundingBox tri_box = CalcBox()(tri);

  if (!m_bbox.do_intersect(tri_box))
    return false;

  marker.unmark_all();

  BoundingBox3i ibox_todo;
  BoundingBox3i ibox(Vec3i(0, 0, 0), m_grid_size - Vec3i(1, 1, 1));
  Box2IBox(tri_box, ibox_todo);
  ibox_todo.enlarge(Vec3i(1, 1, 1));// enlarge ibox_todo for robustness.
  if (ibox_todo.do_intersect(ibox))
  {
    ibox_todo.intersect(ibox);
    int ix, iy, iz;
    Cell first, last, l;
    for (ix = ibox_todo.min().x(); ix <= ibox_todo.max().x(); ix++)
    {
      for (iy = ibox_todo.min().y(); iy <= ibox_todo.max().y(); iy++)
      {
        for (iz = ibox_todo.min().z(); iz <= ibox_todo.max().z(); iz++)
        {
          find_grid(Vec3i(ix, iy, iz), first, last);
          for (l = first;l != last;++l) if (l->pidx >= 0)
          {
            Primitive& elem = m_primitives[l->pidx];
            if (!marker.is_marked(l->pidx))
            {
              if (m_prim_boxes[l->pidx].do_intersect(tri_box))
              {
                const Vec3d& a = elem.tri.ver0, b = elem.tri.ver1, c = elem.tri.ver2, p = tri.tri.ver0, q = tri.tri.ver1, r = tri.tri.ver2;
                const ExactPoint* ea = elem.ev0, * eb = elem.ev1, * ec = elem.ev2, * ep = tri.ev0, * eq = tri.ev1, * er = tri.ev2;
                if (triangle_do_intersect(a, b, c, p, q, r, ea, eb, ec, ep, eq, er))
                {
                  return true;
                }
              }
              marker.mark(l->pidx);
            }
          }
        }
      }
    }
  }
  return false;
}

FaceGrid::FaceGrid(SMeshT& mesh)
{
  set_from_mesh(mesh);
  build();
}

FaceGrid::FaceGrid(SMeshT& mesh, const std::vector<FaceHandle>& faces)
{
  set_from_mesh(mesh, faces);
  build();
}

void FaceGrid::set_from_mesh(SMeshT& mesh)
{
  std::vector<ExactTriangle<HeavyTriangle>> triangles;triangles.reserve(mesh.n_faces());
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
  insert(std::move(triangles));
}

void FaceGrid::set_from_mesh(SMeshT& mesh, const std::vector<FaceHandle>& faces)
{
  std::vector<ExactTriangle<HeavyTriangle>> triangles;triangles.reserve(faces.size());
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
  insert(std::move(triangles));
}
}// namespace Geometry
}// namespace Cage