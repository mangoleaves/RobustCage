#include "CollapseStage.hh"
#include <random>

#define Np 50
#define MUTATION_SCALE 0.9
#define CROSSOVER_RATE 0.9
#define CONVERGENCE_RATE 1e-4
#define MAX_CONSECUTIVE_ITER 5

namespace Cage
{
namespace CageSimp
{

CollapseStage::CollapseStage(
  SMeshT* original, SMeshT* cage, ParamCollapseStage* p,
  DFaceTree* original_tree, LightDFaceTree* remeshing_tree,
  FaceGrid* original_grid,
  double _original_diagonal_length)
  :om(original), rm(cage), param(p),
  ot(original_tree), lrt(remeshing_tree), og(original_grid),
  original_diagonal_length(_original_diagonal_length)
{}

void CollapseStage::update(size_t _candidate_points_size, bool _allow_negtive, double _max_distance_error)
{
  candidate_points_size = _candidate_points_size;
  allow_negtive = _allow_negtive;
  max_distance_error = _max_distance_error;
}

/// @brief generate random points around smooth target of edge.
std::vector<Vec3d> CollapseStage::generate_candidate_points_for_collapse(EdgeHandle e, EdgeCollapser& edge_collapser)
{
  // calculate radius, which is average length of adjacent edges.
  double radius = 0.0;
  HalfedgeHandle he = rm->halfedge_handle(e, 0);
  for (EdgeHandle ve : rm->ve_range(rm->to_vertex_handle(he)))
    radius += rm->data(ve).edge_length;
  for (EdgeHandle ve : rm->ve_range(rm->from_vertex_handle(he)))
    radius += rm->data(ve).edge_length;
  radius /= rm->valence(rm->to_vertex_handle(he)) + rm->valence(rm->from_vertex_handle(he));
  radius *= 0.2;
  radius = std::min(radius, original_diagonal_length * 0.01);

  // tangential relaxation target.
  auto& halfedges = edge_collapser.get_halfedges();

  Vec3d edge_midpoint = rm->calc_edge_midpoint(he);
  Vec3d new_vertex_normal, vertex_target;
  edge_collapser.predict_tangential_weighted_smooth_target(edge_midpoint, new_vertex_normal, vertex_target);

  // local coordinate system on target point.
  Vec3d local_axis_x, local_axis_y;
  make_coordinate_system(new_vertex_normal, local_axis_x, local_axis_y);

  std::vector<Vec3d> points; points.reserve(candidate_points_size);
  // generate regular points
  points.push_back(vertex_target);
  points.push_back(rm->calc_edge_midpoint(he));
  points.push_back(rm->point(rm->to_vertex_handle(he)));
  points.push_back(rm->point(rm->from_vertex_handle(he)));

  // generate random points
  while (points.size() < candidate_points_size)
  {
    double height = (double)rand() / (double)RAND_MAX;
    double beta = (double)rand() / (double)RAND_MAX;
    double len = (double)rand() / (double)RAND_MAX;
    height = (height - 0.5) * original_diagonal_length * 0.005;
    beta = beta * 2 * M_PI;

    Vec3d new_point = vertex_target + (local_axis_x * cos(beta) + local_axis_y * sin(beta)) * len * radius + height * new_vertex_normal;
    points.push_back(new_point);
  }
  points.resize(candidate_points_size);

  return points;
}

bool CollapseStage::find_collapse_hausdorff_deviation(
  EdgeHandle eh, double& local_hd_before, double& local_hd_after, Vec3d& new_point)
{
  auto edge_collapser = new_edge_collapser();

  if (!edge_collapser.init(eh))
    return false;

  // constrain valence
  HalfedgeHandle heh = rm->halfedge_handle(eh, 0);
  size_t valence_after_collapsing =
    (rm->valence(rm->to_vertex_handle(heh)) + rm->valence(rm->from_vertex_handle(heh))) - 3;
  if (valence_after_collapsing > param->maxValence)
    return false;

  local_hd_before = edge_collapser.local_Hausdorff_before_collapsing();

  // get candidate points
  std::vector<Vec3d> candidate_points = generate_candidate_points_for_collapse(eh, edge_collapser);

  // construct local mesh
  VertexHandle local_center_v;
  std::vector<SMeshT> local_meshes = construct_local_meshes(rm, edge_collapser.get_halfedges(), candidate_points, local_center_v);

  // find optimal point that minimize local hausdorff distance.
  size_t minimal_idx = 0;
  double minimal_local_hd = DBL_MAX;

#if USE_TREE_SEARCH
  ot->set_hint(candidate_points[0]);
#endif
#pragma omp parallel for schedule(dynamic)
  for (int i = 0;i < (int)candidate_points_size;i++)
  {
    double local_hd_i = edge_collapser.local_Hausdorff_after_collapsing(&local_meshes[i], candidate_points[i], minimal_local_hd);
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

void CollapseStage::initialize_collapse_edges_reward()
{
  // clear
  update_states.resize(rm->n_edges());
  std::fill(update_states.begin(), update_states.end(), 0);
  edges_to_collapse = CollapseEdgeRewardQueue();

  // calculate average edge length
  avg_edge_length = 0.0;
  for (EdgeHandle eh : rm->edges())
    avg_edge_length += rm->data(eh).edge_length;
  avg_edge_length /= rm->n_edges();

  // initialize for all edges
  for (EdgeHandle eh : rm->edges())
  {
    double local_hd_after, local_hd_before;
    Vec3d new_point;
    if (find_collapse_hausdorff_deviation(eh, local_hd_before, local_hd_after, new_point))
    {
      if (allow_negtive)
      {
        // we collapse closest first.
        if (local_hd_after < max_distance_error)
        {
          //edges_to_collapse.emplace(eh, 0, -local_hd_after, new_point);
          double error = local_hd_after + 0.1 * (rm->data(eh).edge_length - avg_edge_length);
          edges_to_collapse.emplace(eh, 0, -error, new_point);
        }
      }
      else
      {
        // we collapse largest gain of hd first.
        if (local_hd_before - local_hd_after >= 0.0)
        {
          //edges_to_collapse.emplace(eh, 0, local_hd_before - local_hd_after, new_point);
          double error = local_hd_before - local_hd_after + 0.1 * (avg_edge_length - rm->data(eh).edge_length);
          edges_to_collapse.emplace(eh, 0, error, new_point);
        }
      }
    }
  }
}

void CollapseStage::update_after_collapsing(VertexHandle collapsed_center)
{
  // update in out error
  std::vector<FaceHandle> faces;
  std::vector<EdgeHandle> edges;
  for (FaceHandle vf : rm->vf_range(collapsed_center)) faces.push_back(vf);
  for (EdgeHandle ve : rm->ve_range(collapsed_center)) edges.push_back(ve);

#ifdef USE_TREE_SEARCH
  calc_face_out_error(rm, om, ot, faces, cage_infinite_fp);
  calc_edge_out_error(rm, om, ot, edges, cage_infinite_fp);
  calc_vertex_out_error(rm, om, ot, collapsed_center);
#else
  calc_face_out_error(rm, om, og, faces, cage_infinite_fp);
  calc_edge_out_error(rm, om, og, edges, cage_infinite_fp);
  calc_vertex_out_error(rm, om, og, collapsed_center);
#endif
  calc_face_in_error(rm, faces);

  calc_out_surround_error(rm, collapsed_center);
  for (VertexHandle vv : rm->vv_range(collapsed_center))
    calc_out_surround_error(rm, vv);

  // update optimal vertex position
  std::vector<EdgeHandle> affected_edges = find_1rv_1re(rm, collapsed_center);
  for (EdgeHandle eh : affected_edges)
  {
    update_states[eh.idx()]++;

    double local_hd_after, local_hd_before;
    Vec3d new_point;
    if (find_collapse_hausdorff_deviation(eh, local_hd_before, local_hd_after, new_point))
    {
      if (allow_negtive)
      {
        // we collapse closest first.
        if (local_hd_after < max_distance_error)
        {
          //edges_to_collapse.emplace(eh, update_states[eh.idx()], -local_hd_after, new_point);
          double error = local_hd_after + 0.1 * (rm->data(eh).edge_length - avg_edge_length);
          edges_to_collapse.emplace(eh, update_states[eh.idx()], -error, new_point);
        }
      }
      else
      {
        // we collapse largest gain of hd first.
        if (local_hd_before - local_hd_after >= 0.0)
        {
          //edges_to_collapse.emplace(eh, update_states[eh.idx()], local_hd_before - local_hd_after, new_point);
          double error = local_hd_before - local_hd_after + 0.1 * (avg_edge_length - rm->data(eh).edge_length);
          edges_to_collapse.emplace(eh, update_states[eh.idx()], error, new_point);
        }
      }
    }
  }
}

void CollapseStage::do_collapse(size_t edge_num_to_collapse, size_t& total_collapsed_edge_num)
{
  auto edge_collapser = new_edge_collapser();
  edge_collapser.set_flags(
    /*update_links*/true, /*update_target_length*/false, /*update_normals*/true,
    /*check_wrinkle*/true, /*check_selfinter*/true, /*check_inter*/true);

  initialize_collapse_edges_reward();
  size_t n_vertices = rm->n_vertices();
  size_t collapsed_edge_num = 0;
  while (!edges_to_collapse.empty())
  {
    auto edge_reward = edges_to_collapse.top();
    edges_to_collapse.pop();

    // out of date
    if (edge_reward.state < update_states[edge_reward.eh.idx()])
      continue;
    if (rm->status(edge_reward.eh).deleted())
      continue;

    if (edge_collapser.try_collapse_edge(edge_reward.eh, edge_reward.new_point, nullptr))
    {
      VertexHandle center_v = edge_collapser.get_collapsed_center();
      update_after_collapsing(center_v);

      collapsed_edge_num++;
      total_collapsed_edge_num++;

      if (total_collapsed_edge_num == edge_num_to_collapse)
        break; // end collapsing

      //if (collapsed_edge_num % 1000 == 0)
      //  Logger::user_logger->info("collapsed {} edges.", collapsed_edge_num);
    }
  }
  Logger::user_logger->info("collapsed {} edges.", collapsed_edge_num);
  rm->garbage_collection();
  lrt->collect_garbage();
  init_one_ring_faces(rm);
}

}// namespace CageSimp
}// namespace Cage