#pragma once
#include "Tetrahedralizer.hh"
#include "TetMeshTrimmer.hh"

namespace Cage
{
namespace CageInit
{
using namespace SimpleUtils;
namespace SM = SurfaceMesh;
namespace VM = VolumeMesh;

class CageInitializer
{
public:
  // input
  SM::SMeshT* SMesh;
  ParamCageInitializer* param;

  // sub components
  std::unique_ptr<Tetrahedralizer> tetrahedralizer;
  std::unique_ptr<TetRemover> tetRemover;
  std::unique_ptr<TetMeshTrimmer> tetMeshTrimmer;

  // middle results
  VM::VMeshT* outVMesh;

  // final results
  SM::SMeshT* outSMesh;
public:
  CageInitializer();
  CageInitializer(
    SM::SMeshT* _SMesh, ParamCageInitializer* _param,
    VM::VMeshT* _outVMesh, SM::SMeshT* _outSMesh
  );

  void generate();

private:
  BoundingBox bbox;
  double bbox_diag_length;
  void getBoundingBox();

  void tetrahedralizationPostProcess();
  void separateMeshToInOut();
  void retrieveBoundary(VM::VMeshT* vol_mesh, SM::SMeshT* boundary_mesh);
  void retrieveCage(VM::VMeshT* vol_mesh, SM::SMeshT* boundary_mesh);
};
}// namespace CageInit
}// namespace Cage 