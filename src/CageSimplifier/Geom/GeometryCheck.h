#pragma once
#include "Mesh/SurfaceMesh/SurfaceMeshDefinition.h"
#include "Geometry/GeometryMath.h"
#include "CageSimplifier/Topo/TopoOperations.h"
#include "CageSimplifier/SpaceSearch/FaceTree.h"
#include "CageSimplifier/SpaceSearch/DFaceTree.h"
#include "CageSimplifier/SpaceSearch/VertexTree.h"
#include "CageSimplifier/SpaceSearch/FaceGrid.h"
#include "Geometry/Exact/Predicates.h"
#include <set>

namespace Cage
{
namespace CageSimp
{
using namespace Geometry;
using namespace SurfaceMesh;
using Geometry::FaceTree;
using Geometry::DFaceTree;
using Geometry::VertexTree;

using OpenMesh::VertexHandle;
using OpenMesh::HalfedgeHandle;
using OpenMesh::EdgeHandle;
using OpenMesh::FaceHandle;

typedef const Vec3d& CR_Vec3d;
extern double cage_infinite_fp;

/*****************************/
/*   Sample Check Hausdorff  */
/*****************************/

/* samples points on faces and edges */

void get_nb_samples_on_faces(SMeshT* mesh, const std::vector<FaceHandle>& faces);
void get_nb_samples_on_edges(SMeshT* mesh, const std::vector<EdgeHandle>& edges);
void clear_nb_samples_on_faces(SMeshT* mesh, const std::vector<FaceHandle>& faces);
void clear_nb_samples_on_edges(SMeshT* mesh, const std::vector<EdgeHandle>& edges);
std::vector<Vec3d> generate_samples_on_face(SMeshT* mesh, FaceHandle f);
std::vector<Vec3d> generate_samples_on_edge(SMeshT* mesh, EdgeHandle e);

/* links between meshes */

void generate_face_out_links(SMeshT* out_mesh, SMeshT* in_mesh, FaceTree* in_tree, const std::vector<FaceHandle>& faces);
void generate_edge_out_links(SMeshT* out_mesh, SMeshT* in_mesh, FaceTree* in_tree, const std::vector<EdgeHandle>& edges);
void generate_vertex_out_links(SMeshT* out_mesh, SMeshT* in_mesh, FaceTree* in_tree, VertexHandle v);
void generate_out_links(SMeshT* out_mesh, SMeshT* in_mesh, FaceTree* in_tree);

double calc_face_in_error(SMeshT* in_mesh, const std::vector<FaceHandle>& faces);
#if USE_TREE_SEARCH
double calc_face_out_error(SMeshT* out_mesh, SMeshT* in_mesh, DFaceTree* in_tree, const std::vector<FaceHandle>& faces, double& threshold);
double calc_edge_out_error(SMeshT* out_mesh, SMeshT* in_mesh, DFaceTree* in_tree, const std::vector<EdgeHandle>& edges, double& threshold);
double calc_vertex_out_error(SMeshT* out_mesh, SMeshT* in_mesh, DFaceTree* in_tree, VertexHandle v);
double calc_out_error(SMeshT* out_mesh, SMeshT* in_mesh, DFaceTree* in_tree, double& threshold);
#else
double calc_face_out_error(SMeshT* out_mesh, SMeshT* in_mesh, FaceGrid* in_grid, const std::vector<FaceHandle>& faces, double& threshold);
double calc_edge_out_error(SMeshT* out_mesh, SMeshT* in_mesh, FaceGrid* in_grid, const std::vector<EdgeHandle>& edges, double& threshold);
double calc_vertex_out_error(SMeshT* out_mesh, SMeshT* in_mesh, FaceGrid* in_grid, VertexHandle v);
double calc_out_error(SMeshT* out_mesh, SMeshT* in_mesh, FaceGrid* in_grid, double& threshold);
#endif
void calc_out_surround_error(SMeshT* out_mesh, VertexHandle v);

LinkVector backup_local_in_links(SMeshT* mesh, const std::vector<FaceHandle>& faces);

size_t count_face_out_links(SMeshT* mesh);
size_t count_edge_out_links(SMeshT* mesh);
void clear_links(SMeshT* mesh);

/********************/
/*     Predicate    */
/********************/

int same_side_wrapped(CR_Vec3d a, CR_Vec3d b, CR_Vec3d c, CR_Vec3d p);

/*******************/
/*    Constraint   */
/*******************/

bool check_wrinkle(
  SMeshT* mesh,
  const std::vector<HalfedgeHandle>& halfedges,
  CR_Vec3d new_point
);

bool check_long_edge(
  SMeshT* mesh,
  const std::vector<HalfedgeHandle>& halfedges,
  CR_Vec3d new_point,
  double new_target_length
);

bool check_degenerate(
  SMeshT* mesh,
  const std::vector<HalfedgeHandle>& halfedges,
  CR_Vec3d new_point,
  const ExactPoint* new_ep
);

bool check_intersection(
  SMeshT* mesh,
  const std::vector<HalfedgeHandle>& halfedges,
  CR_Vec3d new_point,
  const ExactPoint* new_ep,
  DFaceTree* intersect_tree,
	FaceGrid* intersect_grid,
  LightDFaceTree* self_intersect_tree,
  const std::set<FaceHandle>& ignored_faces,
  bool check_selfinter = true,
  bool check_inter = true
  );

}// namespace CageSimp
}// namespace Cage