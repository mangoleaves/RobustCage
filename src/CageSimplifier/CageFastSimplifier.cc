#include "CageFastSimplifier.hh"

namespace Cage
{
namespace CageSimp
{

void FastSimplifier::simplify()
{
  calc_target_length_on_origin();
  size_t force_skip_times = 3;
  for (size_t i = 0;i < force_skip_times;i++)
    force_skip_vn.push(param->targetVerticesNum * (force_skip_times - i));

  size_t it = 0;
  while (rm->n_vertices() > param->targetVerticesNum)
  {
    Logger::user_logger->info("fast simplify [iter {}].", ++it);
    find_target_length_on_remeshing();
    for (VertexHandle vh : rm->vertices())
    {
    #ifdef USE_TREE_SEARCH
      rm->data(vh).out_error = ot->closest_distance(rm->point(vh)).first;
    #else
      rm->data(vh).out_error = og->closest_point(rm->point(vh)).first.second;
    #endif
    }
    size_t n_collapse = collapse_short_edges();
    size_t n_equalize = equalize_valences();
    size_t n_smooth = Laplacian_smooth();
    enlarge_target_length_on_origin();
    if (n_collapse == 0)
      break;
  }
}

void FastSimplifier::calc_target_length_on_origin()
{
  // principal curvature
  PrincipalCurvature mesh_curvature(om);
  mesh_curvature.compute_principal_curvature();
  auto& K1 = mesh_curvature.principal_K1;
  auto& K2 = mesh_curvature.principal_K2;
  auto& dir1 = mesh_curvature.principal_dir1;
  auto& dir2 = mesh_curvature.principal_dir2;

  size_t V_N = om->n_vertices();
  std::vector<double> curvature; curvature.reserve(V_N);
  for (size_t i = 0; i < V_N; i++)
    curvature.push_back(std::max(abs(K1[i]), abs(K2[i])));

  // calculate target length on original mesh
  double epsilon = original_diagonal_length * 0.001;
  for (VertexHandle vh : om->vertices())
  {
    double length_temp = sqrt(6 * epsilon / curvature[vh.idx()] - 3 * epsilon * epsilon);
    if (length_temp < epsilon || isnan(length_temp))
      length_temp = epsilon;
    if (length_temp > 40 * epsilon)
      length_temp = 40 * epsilon;
    om->data(vh).target_length = length_temp;
  }
  for (EdgeHandle eh : om->edges())
  {
    HalfedgeHandle hh = om->halfedge_handle(eh, 0);
    VertexHandle v_from = om->from_vertex_handle(hh), v_to = om->to_vertex_handle(hh);
    om->data(eh).target_length = std::min(om->data(v_from).target_length,
      om->data(v_to).target_length);
  }

  // smooth target length on original mesh
  size_t max_smooth_iter = 10;
  for (size_t iter = 0;iter < max_smooth_iter;iter++)
  {
    for (VertexHandle vh : om->vertices())
    {
      double len_range = 0.0;
      for (VertexHandle vv : om->vv_range(vh))
        len_range += om->data(vv).target_length;
      om->data(vh).target_length = len_range / om->valence(vh);
    }
  }
}

void FastSimplifier::enlarge_target_length_on_origin()
{
  for (VertexHandle vh : om->vertices())
  {
    om->data(vh).target_length *= param->enlargeTargetLengthRatio;
  }
  for (EdgeHandle eh : om->edges())
  {
    HalfedgeHandle hh = om->halfedge_handle(eh, 0);
    VertexHandle v_from = om->from_vertex_handle(hh), v_to = om->to_vertex_handle(hh);
    om->data(eh).target_length = std::min(om->data(v_from).target_length,
      om->data(v_to).target_length);
  }
}

void FastSimplifier::find_target_length_on_remeshing()
{
  for (VertexHandle vh : rm->vertices())
  {
    VertexHandle cv = vt->closest_vertex(rm->point(vh));
    rm->data(vh).target_length = om->data(cv).target_length;
  }
  for (EdgeHandle eh : rm->edges())
  {
    HalfedgeHandle hh = rm->halfedge_handle(eh, 0);
    VertexHandle v_from = rm->from_vertex_handle(hh), v_to = rm->to_vertex_handle(hh);
    rm->data(eh).target_length = std::min(rm->data(v_from).target_length,
      rm->data(v_to).target_length);
  }
}

void FastSimplifier::initialize_edges_to_collapse()
{
  update_state.clear();
  update_state.resize(rm->n_edges(), 0);
  edges_to_collapse = EdgeQueue();

  for (EdgeHandle eh : rm->edges())
  {
    if (rm->data(eh).edge_length < (0.8 * rm->data(eh).target_length))
      edges_to_collapse.emplace(eh, 0, rm->data(eh).edge_length);
  }
}

void FastSimplifier::update_after_collapsing(VertexHandle collapsed_center)
{
  for (EdgeHandle veh : rm->ve_range(collapsed_center))
  {
    update_state[veh.idx()]++;
    if (rm->data(veh).edge_length < (0.8 * rm->data(veh).target_length))
      edges_to_collapse.emplace(veh, update_state[veh.idx()], rm->data(veh).edge_length);
  }
}

/// @brief collapse edges shorter than target length.
size_t FastSimplifier::collapse_short_edges()
{
  auto edge_collapser = new_edge_collapser();
  edge_collapser.set_flags(
    /*update_links*/false, /*update_target_length*/true,
    /*update_normals*/true, /*check_wrinkle*/false,
    /*check_selfinter*/true, /*check_inter*/true);
  initialize_edges_to_collapse();

  size_t n_vertices = rm->n_vertices();
  size_t collapsed_num = 0;
  while (!edges_to_collapse.empty())
  {
    auto edge = edges_to_collapse.top();
    edges_to_collapse.pop();

    // out of date
    if (edge.state < update_state[edge.eh.idx()])
      continue;
    if (rm->status(edge.eh).deleted())
      continue;

    if (edge_collapser.try_collapse_avoid_long_short_edge(edge.eh))
    {
      VertexHandle center_v = edge_collapser.get_collapsed_center();
      update_after_collapsing(center_v);
      collapsed_num++;
      n_vertices--;

      if (n_vertices <= force_skip_vn.front())
      {
        force_skip_vn.pop();
        break;
      }
    }
  }
  Logger::user_logger->info("collapsed [{}] edges.", collapsed_num);
  rm->garbage_collection();
  lrt->collect_garbage();
  init_one_ring_faces(rm);
  Logger::user_logger->info("[{}] vertices and [{}] faces remained.", rm->n_vertices(), rm->n_faces());
  return collapsed_num;
}

size_t FastSimplifier::Laplacian_smooth()
{
  auto vertex_relocater = new_vertex_relocater();
  vertex_relocater.set_flags(
    /*update_links*/false,
    /*update_target_length*/true,
    /*update_normals*/true,
    /*check_wrinkle*/false);

  size_t smooth_num = 0;
  for (size_t it = 0;it < param->smoothIter;it++)
  {
    size_t smooth_num_each_iter = 0;
    for (VertexHandle v : rm->vertices())
    {
      if (rm->status(v).deleted())
        continue;

      vertex_relocater.init(v);
      Vec3d new_point = vertex_relocater.find_smooth_target_with_target_length();
      if ((new_point - rm->point(v)).norm() < GeomThr::relocate_length_thr)
        continue; // too short relocate length

      if (vertex_relocater.try_relocate_vertex(new_point))
      {
        smooth_num_each_iter++;
      }
    }
    if (smooth_num_each_iter == 0)
      break;
    smooth_num += smooth_num_each_iter;
  }
  Logger::user_logger->info("smoothed [{}] vertices.", smooth_num);
  return smooth_num;
}

size_t FastSimplifier::equalize_valences()
{
  auto edge_flipper = new_edge_flipper();
  edge_flipper.set_flags( /*update_links*/false, /*update_target_length*/true, /*update_normals*/false);

  size_t flipped_num = 0;
  for (size_t it = 0;it < param->equalizeValenceIter;it++)
  {
    for (EdgeHandle e : rm->edges())
    {
      if (edge_flipper.try_flip_edge_decrease_valence(e))
      {
        flipped_num++;
      }
    }
  }

  Logger::user_logger->info("flipped [{}] edges.", flipped_num);
  return flipped_num;
}


}// namespace CageSimp
}// namespace Cage