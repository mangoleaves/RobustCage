#pragma once

#include "SurfaceMeshDefinition.h"

namespace Cage
{
namespace SurfaceMesh
{
using namespace OpenMesh;

class MeshConnectRegion
{
public:
  MeshConnectRegion() = default;
  MeshConnectRegion(SMeshT* _mesh) :mesh(_mesh) {}

  void find();
  inline const std::vector<std::vector<FaceHandle>>& get() { return regions; }
private:
  SMeshT* mesh;
  std::vector<std::vector<FaceHandle>> regions;
};

}// namespace SurfaceMesh
}// namespace Cage