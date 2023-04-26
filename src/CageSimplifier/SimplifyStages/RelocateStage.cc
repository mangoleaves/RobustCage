#include "RelocateStage.hh"
#include <random>

namespace Cage
{
namespace CageSimp
{
RelocateStage::RelocateStage(
  SMeshT* original, SMeshT* cage, ParamRelocateStage* p,
  VertexTree* vertex_tree, DFaceTree* original_tree, LightDFaceTree* remeshing_tree,
  FaceGrid* original_grid,
  double _original_diagonal_length)
  :om(original), rm(cage), param(p),
  vt(vertex_tree), ot(original_tree), lrt(remeshing_tree), og(original_grid),
  original_diagonal_length(_original_diagonal_length)
{}

void RelocateStage::update(size_t _candidate_points_size, bool _allow_negtive, double _max_distance_error)
{
  candidate_points_size = _candidate_points_size;
  allow_negtive = _allow_negtive;
  max_distance_error = _max_distance_error;
}

std::vector<Vec3d> RelocateStage::generate_candidate_points_for_relocate(VertexHandle vh, VertexRelocater& vertex_relocater)
{
  // calculate radius, which is average length of adjacent edges.
  double radius = 0.0;
  for (EdgeHandle ve : rm->ve_range(vh))
    radius += rm->data(ve).edge_length;
  radius /= rm->valence(vh);
  radius *= 0.2;
  radius = std::min(radius, original_diagonal_length * 0.01);

  // tangential relaxation target
  Vec3d vertex_target;
  vertex_target = vertex_relocater.find_weighted_tangential_smooth_target();

  // local coordinate system on target point
  const Vec3d& vertex_normal = rm->normal(vh);
  Vec3d local_axis_x, local_axis_y;
  make_coordinate_system(vertex_normal, local_axis_x, local_axis_y);

  std::vector<Vec3d> points; points.reserve(candidate_points_size);
  // generate regular points
  points.push_back(vertex_target);

  // generate random points
  while (points.size() < candidate_points_size)
  {
    double height = (double)rand() / (double)RAND_MAX;
    double beta = (double)rand() / (double)RAND_MAX;
    double len = (double)rand() / (double)RAND_MAX;
    height = (height - 0.5) * original_diagonal_length * 0.005;
    beta = beta * 2 * M_PI;

    Vec3d new_point = vertex_target + (local_axis_x * cos(beta) + local_axis_y * sin(beta)) * len * radius + height * vertex_normal;
    points.push_back(new_point);
  }

  return points;
}

bool RelocateStage::find_relocate_hausdorff_deviation(
  VertexRelocater& vertex_relocater, VertexHandle vh, double& local_hd_before, double& local_hd_after, Vec3d& new_point)
{
  if (!vertex_relocater.init(vh))
    return false;

  local_hd_before = vertex_relocater.local_Hausdorff_before_relocating();

  // get candidate points
  std::vector<Vec3d> candidate_points = generate_candidate_points_for_relocate(vh, vertex_relocater);

  // construct local mesh
  VertexHandle local_center_v;
  std::vector<SMeshT> local_meshes = construct_local_meshes(rm, vertex_relocater.get_halfedges(), candidate_points, local_center_v);

  // find optimal point that minimize local hausdorff distance.
  size_t minimal_idx = 0;
  double minimal_local_hd = DBL_MAX;

#if USE_TREE_SEARCH
  ot->set_hint(candidate_points[0]);
#endif
#pragma omp parallel for schedule(dynamic)
  for (int i = 0;i < (int)candidate_points_size;i++)
  {
    double local_hd_i = vertex_relocater.local_Hausdorff_after_relocating(&local_meshes[i], candidate_points[i], minimal_local_hd);
  #pragma omp critical 
    if (local_hd_i < minimal_local_hd)
    {
      minimal_local_hd = local_hd_i;
      minimal_idx = i;
    }
  }
#if USE_TREE_SEARCH
  ot->unset_hint();
#endif
  // return result
  if (minimal_local_hd != DBL_MAX)
  {
    // found a available point
    local_hd_after = minimal_local_hd;
    new_point = candidate_points[minimal_idx];
    return true;
  }
  else return false;
}

void RelocateStage::update_after_relocating(VertexHandle relocate_center)
{
  // update in out error
  std::vector<FaceHandle> faces;
  std::vector<EdgeHandle> edges;
  for (FaceHandle vf : rm->vf_range(relocate_center)) faces.push_back(vf);
  for (EdgeHandle ve : rm->ve_range(relocate_center)) edges.push_back(ve);

  double infinite = DBL_MAX;
#ifdef USE_TREE_SEARCH
  calc_face_out_error(rm, om, ot, faces, infinite);
  calc_edge_out_error(rm, om, ot, edges, infinite);
  calc_vertex_out_error(rm, om, ot, relocate_center);
#else
  calc_face_out_error(rm, om, og, faces, infinite);
  calc_edge_out_error(rm, om, og, edges, infinite);
  calc_vertex_out_error(rm, om, og, relocate_center);
#endif
  calc_face_in_error(rm, faces);

  calc_out_surround_error(rm, relocate_center);
  for (VertexHandle vv : rm->vv_range(relocate_center))
    calc_out_surround_error(rm, vv);
}

void RelocateStage::do_relocate()
{
  auto vertex_relocater = new_vertex_relocater();
  vertex_relocater.set_flags(/*update_links*/true, /*update_target_length*/false, /*update_normals*/true, /*check_wrinkle*/false);
  size_t relocated_vertex_num = 0;
  for (size_t it = 0;it < param->smoothIter;it++)
  {
    for (VertexHandle vh : rm->vertices())
    {
      Vec3d new_point;
      double local_hd_before, local_hd_after = DBL_MAX;

      if (find_relocate_hausdorff_deviation(vertex_relocater, vh, local_hd_before, local_hd_after, new_point))
      {
        bool do_it = false;
        if (!allow_negtive && local_hd_before - local_hd_after >= 0.0)
        {
          do_it = true;
        }
        else if (allow_negtive&& local_hd_after < max_distance_error)
        {
          double local_hd_diff = (local_hd_after - local_hd_before) / local_hd_before;
          if (local_hd_diff <= 0.0)
            do_it = true;
          else if (local_hd_diff >= 1.0)
            do_it = false;
          else
            do_it = ((double)rand() / RAND_MAX) > local_hd_diff;
        }

        if (do_it)
        {
          vertex_relocater.relocate(new_point);
          relocated_vertex_num += 1;
          update_after_relocating(vh);
        }
      }
    }
  }
  Logger::user_logger->info("relocated {} vertices.", relocated_vertex_num);
}

}// namespace CageSimp
}// namespace Cage