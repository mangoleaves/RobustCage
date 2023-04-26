#pragma once

#pragma once
#include <queue>
#include <set>
#include <vector>

#include "Mesh/SurfaceMesh/SurfaceMeshDefinition.h"
#include "Mesh/SurfaceMesh/MeshCurvature.h"
#include "Geometry/GeometryMath.h"
#include "Geom/GeometryCheck.h"
#include "Topo/EdgeCollapser.h"
#include "Topo/EdgeFlipper.h"
#include "Topo/VertexRelocater.h"
#include "Utils/logger.hh"
#include "Config.hh"

namespace Cage
{
namespace CageSimp
{
using namespace SimpleUtils;
using namespace Geometry;
using namespace SurfaceMesh;

class FastSimplifier
{
public:
  // input
  SMeshT* om;
  SMeshT* rm;
  // tree
  VertexTree* vt;
  DFaceTree* ot;
  LightDFaceTree* lrt;
  FaceGrid* og;
  // parameters
  ParamFastSimplifier* param;
  double original_diagonal_length;

  FastSimplifier(SMeshT* original, SMeshT* remeshing,
    VertexTree* vertex_tree, DFaceTree* original_tree, LightDFaceTree* remeshing_tree,
    FaceGrid* original_grid,
    ParamFastSimplifier* _param) :
    om(original), rm(remeshing), vt(vertex_tree), ot(original_tree), lrt(remeshing_tree), og(original_grid),
    param(_param)
  {}

  void simplify();
private:
  void calc_target_length_on_origin();
  void enlarge_target_length_on_origin();
  void find_target_length_on_remeshing();

  std::vector<size_t> update_state;
  struct EdgeLength
  {
    EdgeHandle eh;
    size_t state;
    double length;

    EdgeLength()=default;
    EdgeLength(EdgeHandle _eh, size_t _state, double _length) :
      eh(_eh), state(_state), length(_length)
    {}
    bool operator>(const EdgeLength& rhs)const { return length > rhs.length; }
  };
  typedef std::priority_queue<EdgeLength, std::vector<EdgeLength>, std::greater<EdgeLength>> EdgeQueue;
  EdgeQueue edges_to_collapse;

  std::queue<size_t> force_skip_vn;

  void initialize_edges_to_collapse();
  void update_after_collapsing(VertexHandle collapsed_center);
  size_t collapse_short_edges();

  size_t Laplacian_smooth();

  size_t equalize_valences();

  inline EdgeCollapser new_edge_collapser() { return EdgeCollapser(om, rm, ot, lrt, og); }
  inline EdgeFlipper new_edge_flipper() { return EdgeFlipper(om, rm, ot, lrt, og); }
  inline VertexRelocater new_vertex_relocater() { return VertexRelocater(om, rm, ot, lrt, vt, og); }
};
}// namespace CageSimp
}// namespace Cage