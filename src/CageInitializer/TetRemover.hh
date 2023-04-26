#pragma once

// Simple mesh
#include "Mesh/VolumeMesh/VolumeMeshDefinition.hh"

namespace Cage
{
namespace CageInit
{
using namespace SimpleUtils;
namespace VM = VolumeMesh;

/// @brief Remove tetrahedron from a tetrahedral mesh.
/// We assume input tetrahedral mesh is 3-manifold.
/// TetRemover can keep mesh 3-manifold after removing tet.
/// See <<Removing Tetrahedra from manifold tetrahedralisation :
/// application to real-time surgical simulation>>
class TetRemover
{
public:
  // input/output
  VM::VMeshT* VMesh;
public:
  TetRemover(VM::VMeshT* _VMesh);

  void removeTet(VM::CellHandle _ch, bool _collectGarbage = false);
  void removeTet(const std::vector<VM::CellHandle>& _chs, bool _collectGarbage = true);
private:
  bool isTetRemovable(VM::CellHandle _ch);

  void updateBoundaryPropBeforeDeleting(VM::CellHandle _ch);
  void updateBoundaryPropBeforeDeleting(const std::vector<VM::CellHandle>& _chs);
  void updateBoundaryPropIfFaceIsDeletedOnce(VM::FaceHandle _fh);
};
}// namespace CageInit
}// namespace Cage