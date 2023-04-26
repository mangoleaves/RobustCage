#pragma once
#include "Mesh/SurfaceMesh/SurfaceMeshDefinition.h"

namespace Cage
{
namespace CageSimp
{
using namespace SurfaceMesh;
using OpenMesh::VertexHandle;
using OpenMesh::HalfedgeHandle;
using OpenMesh::EdgeHandle;
using OpenMesh::FaceHandle;

/********************************/
/*      neighbor search         */
/********************************/

std::vector<FaceHandle>
faces_adjacent_to_edge(SMeshT* mesh, EdgeHandle e);

void one_ring_faces_around_vertex(
  SMeshT* mesh,
  VertexHandle vh,
  std::set<FaceHandle>& faces);

void one_ring_faces_around_edge(
  SMeshT* mesh,
  HalfedgeHandle hh,
  std::set<FaceHandle>& faces);

void vertices_around_faces(
  SMeshT* mesh,
  const std::set<FaceHandle>& faces,
  std::set<VertexHandle>& vertices);

void extend_faces_by_one_ring(
  SMeshT* mesh,
  std::set<FaceHandle>& faces);

void extend_faces(
  SMeshT* mesh,
  const std::set<FaceHandle>& faces,
  int stencil_ring_size,
  std::set<FaceHandle>& extended_faces);

std::vector<EdgeHandle> find_1rv_1re(SMeshT* mesh, VertexHandle v);

void init_one_ring_faces(SMeshT* mesh, FaceHandle fh);
void init_one_ring_faces(SMeshT* mesh);

/**********************************/
/*       Local operations         */
/**********************************/

std::unique_ptr<SMeshT> construct_local_mesh(
  SMeshT* mesh,
  const std::vector<HalfedgeHandle>& halfedges,
  const Vec3d& center_point,
  VertexHandle& local_center_vh
);


std::vector<SMeshT> construct_local_meshes(
  SMeshT* mesh,
  const std::vector<HalfedgeHandle>& halfedges,
  const std::vector<Vec3d>& center_points,
  VertexHandle& local_center_vh);

}// namespace CageSimp
}// namespace Cage