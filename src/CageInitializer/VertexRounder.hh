#pragma once

#include "Mesh/VolumeMesh/VolumeMeshDefinition.hh"

namespace Cage
{
namespace CageInit
{
using namespace VolumeMesh;

class VertexRounder
{
public:
// Input
  VMeshT* VMesh;

  VertexRounder() = delete;
  VertexRounder(VMeshT* _VMesh):VMesh(_VMesh){}

  bool doRounding(VertexHandle vh);
  size_t doRounding(const std::vector<VertexHandle>& vrts);
  size_t doRounding();

  bool checkPositiveVolume(CellHandle ch, VertexHandle vh, const Point_3& approx_p);

  bool checkAllTets();

};
}// namespace CageInit
}// namespace Cage