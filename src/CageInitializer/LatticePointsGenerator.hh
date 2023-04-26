#pragma once
#include <vector>
// SurfaceMesh
#include "Mesh/SurfaceMesh/SurfaceMeshDefinition.h"
// AABB Tree(closest point search)
#include "Geometry/AABB/AABBClosestSearch.h"
// mesh cluster
#include "Mesh/SurfaceMesh/MeshCluster.h"

#include "Utils/logger.hh"
#include "Config.hh"

namespace Cage
{
namespace CageInit
{
using namespace SimpleUtils;
using namespace Geometry;
namespace SM = SurfaceMesh;

class LatticePointsGenerator
{
public:
  // input
  SM::SMeshT* SMesh;
  ParamLatticePointsGenerator* param;
  // results
  std::vector<Vec3d> points;
public:
  LatticePointsGenerator(SM::SMeshT* _SMesh, ParamLatticePointsGenerator* _param);

  const std::vector<Vec3d>& generate();
private:
  std::unique_ptr<AABBClosestSearch> aabbTree;

  BoundingBox bbox;
  double length_threshold;
  double area_threshold;
  void getBoundingBox();

  void uniformlyDistributed();
  void offsetVertices();
  void offsetPoint(SM::FaceHandle fh, const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& n);

  void buildAABBTree();
  void releaseAABBTree();
};
}// namespace CageInit
}// namespace Cage