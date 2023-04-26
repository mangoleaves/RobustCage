#pragma once
#include <queue>
#include <set>
#include <vector>

#include "Mesh/SurfaceMesh/SurfaceMeshDefinition.h"
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
using namespace Geometry;
using namespace SimpleUtils;

class DegenerationRemover
{
public:
  // input
  SMeshT* om;
  SMeshT* rm;

  VertexTree* vt;
  DFaceTree* ot;
  LightDFaceTree* lrt;
  FaceGrid* og;
  // parameters
  double original_diagonal_length;
public:
  DegenerationRemover(
    SMeshT* om_, SMeshT* rm_,
    VertexTree* vt_, DFaceTree* ot_, LightDFaceTree* lrt_,
    FaceGrid* og_) :
    om(om_), rm(rm_),
    vt(vt_), ot(ot_), lrt(lrt_),
    og(og_)
  {}
  ~DegenerationRemover() {}

  void perform();
public:
  // almost degenerate
  bool is_face_almost_degenerate(FaceHandle fh);
  size_t eliminate_almost_degeneration();
  bool try_collapse_almost_degenerate_edge(EdgeCollapser& edge_collapser, EdgeHandle eh);

  inline EdgeCollapser new_edge_collapser() { return EdgeCollapser(om, rm, ot, lrt, og); }
  inline EdgeFlipper new_edge_flipper() { return EdgeFlipper(om, rm, ot, lrt, og); }
  inline VertexRelocater new_vertex_relocater() { return VertexRelocater(om, rm, ot, lrt, vt, og); }
};
}// namespace CageSimp
}// namespace Cage