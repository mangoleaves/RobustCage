#pragma once

#include "Config.hh"
#include "CageSimplifier/Topo/TopoOperations.h"
#include "CageSimplifier/Topo/VertexRelocater.h"
#include "CageSimplifier/Geom/GeometryCheck.h"

namespace Cage
{
namespace CageSimp
{

class RelocateStage
{
public:
  // input
  ParamRelocateStage* param;
  SMeshT* om;
  SMeshT* rm;
  VertexTree* vt;
  DFaceTree* ot;
  LightDFaceTree* lrt;
  FaceGrid* og;

  double original_diagonal_length;
  size_t candidate_points_size;
  bool allow_negtive;
  double max_distance_error;
public:
  RelocateStage(
    SMeshT* original, SMeshT* cage, ParamRelocateStage* p,
    VertexTree* vertex_tree, DFaceTree* original_tree, LightDFaceTree* remeshing_tree,
    FaceGrid* original_grid,
    double _original_diagonal_length);

  void update(size_t _candidate_points_size, bool _allow_negtive, double _max_distance_error);
  void do_relocate();
private:
  /* Relocate */
  std::vector<Vec3d> generate_candidate_points_for_relocate(VertexHandle vh, VertexRelocater& vertex_relocater);
  bool find_relocate_hausdorff_deviation(
    VertexRelocater& relocater, VertexHandle vh, double& local_hd_before, double& local_hd_after, Vec3d& new_point);
  void update_after_relocating(VertexHandle relocate_center);

  inline VertexRelocater new_vertex_relocater() { return VertexRelocater(om, rm, ot, lrt, vt, og); }
};

}// namespace CageSimp
}// namespace Cage