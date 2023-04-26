#pragma once

#include "Mesh/SurfaceMesh/SurfaceMeshDefinition.h"
#include "CageSimplifier/Topo/TopoOperations.h"
#include "CageSimplifier/Geom/GeometryCheck.h"

namespace Cage
{
namespace CageSimp
{
using namespace Geometry;
using namespace SurfaceMesh;
using OpenMesh::VertexHandle;
using OpenMesh::HalfedgeHandle;
using OpenMesh::EdgeHandle;
using OpenMesh::FaceHandle;

class EdgeFlipper
{
public:
  EdgeFlipper(
    SMeshT* _origin, SMeshT* _remeshing,
    DFaceTree* _origin_tree,
    LightDFaceTree* _remeshing_tree,
    FaceGrid* _origin_grid
  )
    :om(_origin), rm(_remeshing),
    ot(_origin_tree), lrt(_remeshing_tree),
    og(_origin_grid),
    initialized(false),
    f_update_links(false),
    f_update_target_length(false),
    f_update_normals(false)
  {}

  bool init(EdgeHandle e);
  void set_flags(
    bool _f_update_links, bool _f_update_target_length,
    bool _f_update_normals
  );

  bool try_flip_edge();
  bool try_flip_edge(EdgeHandle e);

  bool try_flip_edge_decrease_valence(EdgeHandle e);
  bool flip_will_decrease_valence()const;
  bool flip_will_cause_over_valence(size_t max_valence)const;

  bool flip_will_cause_small_large_angle()const;

  double local_Hausdorff_before_flipping()const;
  double local_Hausdorff_after_flipping()const;
private:
  bool initialized;
  // input
  SMeshT* om;
  SMeshT* rm;

  DFaceTree* ot;
  LightDFaceTree* lrt;
  FaceGrid* og;
  // flags
  bool f_update_links;
  bool f_update_target_length;
  bool f_update_normals;

  // ajdacent handles of flip edge
  // notations come from source code of OpenMesh
  HalfedgeHandle a0, b0, a1, a2, b1, b2;
  VertexHandle   va0, va1, vb0, vb1;
  FaceHandle     fa, fb;

  LinkVector faces_in_links;

  void clear();
  void init_handles(EdgeHandle e);

  bool flip_would_cause_intersection()const;

  void backup_in_links();
  void generate_in_links();

  void update_target_len();
  void update_length_and_area();
  void update_normals();
  void update_tree();
  void update_one_ring_faces();
};
}// namespace CageSimp
}// namespace Cage