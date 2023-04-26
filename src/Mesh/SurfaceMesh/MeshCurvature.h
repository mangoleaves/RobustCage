#pragma once

#include "SurfaceMeshDefinition.h"

namespace Cage
{
namespace SurfaceMesh
{

class PrincipalCurvature
{
public:
  PrincipalCurvature() = default;
  PrincipalCurvature(SMeshT* _mesh);

  void compute_principal_curvature();
public:
  std::vector<double> principal_K1;
  std::vector<double> principal_K2;
  std::vector<Vec3d> principal_dir1;
  std::vector<Vec3d> principal_dir2;
private:
  SMeshT* mesh;
  std::vector<Vec3d> vertex_normal;

	std::vector<std::map<int, double>> cornerArea;
	std::vector<double> pointArea;
private:
  void compute_point_area();

  void rot_coord_sys(const Vec3d& old_u, const Vec3d& old_v,
    const Vec3d& new_norm,
    Vec3d& new_u, Vec3d& new_v);

  void proj_curv(const Vec3d& old_u, const Vec3d& old_v,
    double old_ku, double old_kuv, double old_kv,
    const Vec3d& new_u, const Vec3d& new_v,
    double& new_ku, double& new_kuv, double& new_kv);

  // Given a curvature tensor, find principal directions and curvatures
  // Makes sure that pdir1 and pdir2 are perpendicular to normal
  void diagonalize_curv(const Vec3d& old_u, const Vec3d& old_v,
    double ku, double kuv, double kv,
    const Vec3d& new_norm,
    Vec3d& pdir1, Vec3d& pdir2, double& vk1, double& vk2);
};

class GaussCurvature
{
public:
  GaussCurvature() = default;
  GaussCurvature(SMeshT* _mesh, bool _update);

  void update();
  std::vector<double> vertK();
  std::vector<double> edgeK();
  std::vector<double> faceK();
  void faceKMax(Vec3d* KMax);
private:
  SMeshT* mesh;
  std::vector<Vec3d> edge_normal;
  std::vector<std::array<double, 3>> f_alpha;
  std::vector<std::array<int, 3>> f_vId;

private:
  void calcFaceAlpha();
  void calcEdgeNormal();

  void calcMaxPrincpalCurvature(
    Vec3d dp1, Vec3d dp2,
    Vec3d dn1, Vec3d dn2, Vec3d& KMax);
};

}// namespace SurfaceMesh
}// namespace Cage