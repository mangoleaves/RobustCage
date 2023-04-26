#pragma once

// volume mesh
#include "Mesh/VolumeMesh/VolumeMeshDefinition.hh"

namespace Cage
{
namespace CageInit
{
using namespace Geometry;
using namespace SimpleUtils;
namespace VM = VolumeMesh;

/// @brief an implementation of "QUALITY LOCAL REFINEMENT OF TETRAHEDRAL MESHES
/// BASED ON 8-SUBTETRAHEDRON SUBDIVISION"
class TetSubdivider
{
public:
  // input/output
  VM::VMeshT* VMesh;
  // if we only do global subdivide, it's better to store new elements in new mesh.
  VM::VMeshT sub_mesh;
public:
  TetSubdivider(VM::VMeshT* _VMesh);

  void subdivideAll();
  void subdivideNonManifoldVertex(VM::VertexHandle vh);
private:
  /***** global subdivide *****/
  std::vector<std::array<VM::VertexHandle, 4>> cell_vrts;

  inline VM::FaceHandle child(VM::FaceHandle fh, size_t idx) { return fh.idx() * 4 + idx; }
  inline VM::FaceHandle child(VM::FaceHandle fh, VM::VertexHandle vh)
  {
    auto& f = VMesh->face(fh);
    if (vh == VMesh->fromVertex(f.halfedges[0])) return child(fh, 0);
    else if (vh == VMesh->fromVertex(f.halfedges[1])) return child(fh, 1);
    else if (vh == VMesh->fromVertex(f.halfedges[2])) return child(fh, 2);
    else return child(fh, 4);
  }

  inline VM::EdgeHandle child(VM::EdgeHandle eh, size_t idx) { return eh.idx() * 2 + idx; }
  inline VM::EdgeHandle child(VM::EdgeHandle eh, VM::VertexHandle vh)
  {
    auto& e = VMesh->edge(eh);
    if (vh == e.vertices[0]) return child(eh, 0);
    else return child(eh, 1);
  }
  inline VM::VertexHandle& splitVertex(VM::EdgeHandle eh) { return sub_mesh.edge(child(eh, 0)).to(); }

  void generateVertices();

  void splitEdges();
  void updateProp(VM::EdgeHandle _parent, VM::EdgeHandle _e0, VM::EdgeHandle _e1);

  void subdivideFaces();
  void updateProp(VM::FaceHandle _parent, const std::vector<VM::FaceHandle>& _children);
  void updateProp(VM::FaceHandle _parent, VM::EdgeHandle _e0, VM::EdgeHandle _e1, VM::EdgeHandle _e2);

  void subdivideTets_1_12();

  void updateProp(VM::CellHandle _parent, const std::vector<VM::CellHandle>& _children);
};
}// namespace CageInit
}// namespace Cage