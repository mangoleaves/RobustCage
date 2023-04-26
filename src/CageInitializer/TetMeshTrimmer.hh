#pragma once
#include "TetRemover.hh"
#include "TetSubdivider.hh"
#include "VertexRounder.hh"

namespace Cage
{
namespace CageInit
{
namespace VM = VolumeMesh;

class TetMeshTrimmer
{
public:
  // input/output
  VM::VMeshT* VMesh;

  // sub components
  std::unique_ptr<TetRemover> tetRemover;
  std::unique_ptr<TetSubdivider> tetSubdivider;
  std::unique_ptr<VertexRounder> vertexRounder;
public:
  TetMeshTrimmer(VM::VMeshT* _VMesh);

  void trim();
private:
  size_t removeNonAdjacent2ConstraintTets();

  std::vector<VM::VertexHandle> findVrtsLeft();

#ifdef CHECK_MANIFOLD
  void detectAllNonManifold(std::vector<VM::EdgeHandle>& non_manifold_edges, std::vector<VM::VertexHandle>& non_manifold_vertices);
  bool isVertexManifold(VM::VertexHandle _vh, const std::vector<VM::FaceHandle>& adj_boundary_faces);
#endif
};
}// namespace CageInit
}// namespace Cage