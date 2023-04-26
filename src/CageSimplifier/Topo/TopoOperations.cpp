#include "TopoOperations.h"

namespace Cage
{
namespace CageSimp
{

std::vector<FaceHandle>
faces_adjacent_to_edge(SMeshT* mesh, EdgeHandle e)
{
  std::vector<FaceHandle> one_ring_faces;
  HalfedgeHandle he = mesh->halfedge_handle(e, 0);
  if (mesh->is_boundary(e))
  {
    if (mesh->is_boundary(he))
      he = mesh->opposite_halfedge_handle(he);
    one_ring_faces = { mesh->face_handle(he) };
  }
  else
    one_ring_faces = { mesh->face_handle(he), mesh->opposite_face_handle(he) };
  return one_ring_faces;
}

void one_ring_faces_around_vertex(
  SMeshT* mesh,
  VertexHandle vh,
  std::set<FaceHandle>& faces)
{
  for (FaceHandle fh : mesh->vf_range(vh))
  {
    faces.insert(fh);
  }
}

void one_ring_faces_around_edge(
  SMeshT* mesh,
  HalfedgeHandle hh,
  std::set<FaceHandle>& faces)
{
  faces.clear();
  VertexHandle v_from = mesh->from_vertex_handle(hh);
  VertexHandle v_to = mesh->to_vertex_handle(hh);
  one_ring_faces_around_vertex(mesh, v_from, faces);
  one_ring_faces_around_vertex(mesh, v_to, faces);
}

void vertices_around_faces(
  SMeshT* mesh,
  const std::set<FaceHandle>& faces,
  std::set<VertexHandle>& vertices)
{
  for (FaceHandle fh : faces)
  {
    HalfedgeHandle hh = mesh->halfedge_handle(fh);
    vertices.insert(mesh->from_vertex_handle(hh));
    vertices.insert(mesh->to_vertex_handle(hh));
    vertices.insert(mesh->opposite_vh(hh));
  }
}

void extend_faces_by_one_ring(
  SMeshT* mesh,
  std::set<FaceHandle>& faces)
{
  std::set<VertexHandle> vertices;
  vertices_around_faces(mesh, faces, vertices);
  for (VertexHandle vh : vertices)
  {
    one_ring_faces_around_vertex(mesh, vh, faces);
  }
}

/// @brief given faces, collect n ring extended-faces.
/// @param [in] mesh
/// @param [in] faces input faces
/// @param [in] stencil_ring_size n ring, can be 1,2,3...
/// @param [out] extended_faces n ring faces, including input faces.
void extend_faces(
  SMeshT* mesh,
  const std::set<FaceHandle>& faces,
  int stencil_ring_size,
  std::set<FaceHandle>& extended_faces)
{
  extended_faces.clear();
  extended_faces = faces;
  for (int i = 0; i < stencil_ring_size; ++i)
  {
    extend_faces_by_one_ring(mesh, extended_faces);
  }
}

/// @brief find two ring vertices of given vertex handle.
/// then find one ring edges of two ring vertices.
std::vector<EdgeHandle> find_1rv_1re(SMeshT* mesh, VertexHandle v)
{
  std::set<EdgeHandle> surrounding_edges;
  for (VertexHandle vv : mesh->vv_range(v))
  {
    for (EdgeHandle ve : mesh->ve_range(vv))
      surrounding_edges.insert(ve);
  }
  return std::vector<EdgeHandle>(surrounding_edges.begin(), surrounding_edges.end());
}

void init_one_ring_faces(SMeshT* mesh, FaceHandle fh)
{
  std::set<FaceHandle>& one_ring_faces = mesh->data(fh).one_ring_faces;
  one_ring_faces.clear();
  for (VertexHandle fvh : mesh->fv_range(fh))
    for (FaceHandle vfh : mesh->vf_range(fvh))
      one_ring_faces.insert(vfh);
}

void init_one_ring_faces(SMeshT* mesh)
{
  for (FaceHandle fh : mesh->faces())
    init_one_ring_faces(mesh, fh);
}

/// @brief construct a local mesh depending on halfedges.
/// @param [in] mesh global mesh
/// @param [in] halfedges halfedges come fron collapse or relocate.
/// they are one ring halfedges of collapsed edge or relocating vertex.
/// @param [in] center_point position of center point.
/// @param [out] gf_map2_lf global face handle map to local face handle
/// @param [out] local_center_vh center vertex handle of local mesh
std::unique_ptr<SMeshT> construct_local_mesh(
  SMeshT* mesh,
  const std::vector<HalfedgeHandle>& halfedges,
  const Vec3d& center_point,
  VertexHandle& local_center_vh)
{
  auto local_mesh = std::make_unique<SMeshT>();

  local_center_vh = local_mesh->add_vertex(center_point);
  std::map<VertexHandle, VertexHandle> gv_map2_lv;
  for (HalfedgeHandle he : halfedges)
  {
    // add vertices to local mesh
    VertexHandle from_v = mesh->from_vertex_handle(he);
    VertexHandle to_v = mesh->to_vertex_handle(he);
    if (!gv_map2_lv.count(from_v))
      gv_map2_lv[from_v] = local_mesh->add_vertex(mesh->point(from_v));
    if (!gv_map2_lv.count(to_v))
      gv_map2_lv[to_v] = local_mesh->add_vertex(mesh->point(to_v));
    // add faces to local mesh
    local_mesh->add_face(gv_map2_lv[from_v], gv_map2_lv[to_v], local_center_vh);
  }
  pre_calculate_edge_length(local_mesh.get());
  pre_calculate_face_area(local_mesh.get());

  return local_mesh;
}

/// @brief a vector implementation of construct_local_mesh
std::vector<SMeshT> construct_local_meshes(
  SMeshT* mesh,
  const std::vector<HalfedgeHandle>& halfedges,
  const std::vector<Vec3d>& center_points,
  VertexHandle& local_center_vh)
{
  SMeshT local_mesh;
  Vec3d center_point(0., 0., 0.);
  local_center_vh = local_mesh.add_vertex(center_point);
  std::map<VertexHandle, VertexHandle> gv_map2_lv;
  for (HalfedgeHandle he : halfedges)
  {
    // add vertices to local mesh
    VertexHandle from_v = mesh->from_vertex_handle(he);
    VertexHandle to_v = mesh->to_vertex_handle(he);
    if (!gv_map2_lv.count(from_v))
      gv_map2_lv[from_v] = local_mesh.add_vertex(mesh->point(from_v));
    if (!gv_map2_lv.count(to_v))
      gv_map2_lv[to_v] = local_mesh.add_vertex(mesh->point(to_v));
    // add faces to local mesh
    local_mesh.add_face(gv_map2_lv[from_v], gv_map2_lv[to_v], local_center_vh);
  }

  std::vector<SMeshT> local_meshes(center_points.size(), local_mesh);
//#pragma omp parallel for
  for (int i = 0;i < (int)local_meshes.size();i++)
  {
    local_meshes[i].point(local_center_vh) = center_points[i];
    pre_calculate_edge_length(&local_meshes[i]);
    pre_calculate_face_area(&local_meshes[i]);
  }
  return local_meshes;
}
}// namespace CageSimp
}// namespace Cage