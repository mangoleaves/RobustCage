#pragma once

#include "Geometry/AABB/AABBTraits.h"
#include "Mesh/SurfaceMesh/SurfaceMeshDefinition.h"

namespace Cage
{
namespace Geometry
{
using namespace SurfaceMesh;

class VertexTree : public KdTree
{
public:
  VertexTree(SMeshT& mesh);

  SMeshT::VertexHandle closest_vertex(const Vec3d& query);
};
}// namespace Geometry
}// namespace Cage