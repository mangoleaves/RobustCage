#include "FlipStage.hh"

namespace Cage
{
namespace CageSimp
{

FlipStage::FlipStage(
  SMeshT* original, SMeshT* cage, ParamFlipStage* p,
  DFaceTree* original_tree, LightDFaceTree* remeshing_tree,
  FaceGrid* original_grid,
  double _original_diagonal_length)
  :om(original), rm(cage), param(p),
  ot(original_tree), lrt(remeshing_tree), og(original_grid),
  original_diagonal_length(_original_diagonal_length)
{ }

void FlipStage::update(bool _allow_negtive, double _max_distance_error)
{
  allow_negtive = _allow_negtive;
  max_distance_error = _max_distance_error;
}

void FlipStage::initialize_flip_edges_reward()
{
  // clear
  update_states.resize(rm->n_edges());
  std::fill(update_states.begin(), update_states.end(), 0);
  edges_to_flip = FlipEdgeRewardQueue();

  auto edge_flipper = new_edge_flipper();
  // initialize for all edges
  for (EdgeHandle eh : rm->edges())
  {
    if (!edge_flipper.init(eh))
      continue;
    if (edge_flipper.flip_will_cause_over_valence(param->maxValence))
      continue;
    if (!edge_flipper.flip_will_decrease_valence())
      continue;
    double local_hd_before = edge_flipper.local_Hausdorff_before_flipping();
    double local_hd_after = edge_flipper.local_Hausdorff_after_flipping();
    if (local_hd_after == DBL_MAX)
      continue;

    if (allow_negtive)
    {
      if (local_hd_after < max_distance_error)
        // decrease valence and won't cause out of distance error
        edges_to_flip.emplace(eh, 0, -local_hd_after);
    }
    else  // not allow_negtive
    {
      if (local_hd_before - local_hd_after >= 0.0)
        edges_to_flip.emplace(eh, 0, local_hd_before - local_hd_after);
    }
  }
}

void FlipStage::update_after_flipping(EdgeHandle flipped_edge)
{
  // two faces, five edges and four vertices need to be updated.
  HalfedgeHandle flip_he = rm->halfedge_handle(flipped_edge, 0);
  HalfedgeHandle flip_he_opp = rm->halfedge_handle(flipped_edge, 1);
  // update in out error
  std::vector<FaceHandle> faces = { rm->face_handle(flip_he), rm->face_handle(flip_he_opp) };
  std::vector<EdgeHandle> edges = { flipped_edge };
  std::vector<VertexHandle> vertices = {
    rm->to_vertex_handle(flip_he), rm->to_vertex_handle(flip_he_opp) ,
    rm->opposite_vh(flip_he), rm->opposite_vh(flip_he_opp) };

#ifdef USE_TREE_SEARCH
  calc_face_out_error(rm, om, ot, faces, cage_infinite_fp);
  calc_edge_out_error(rm, om, ot, edges, cage_infinite_fp);
#else
  calc_face_out_error(rm, om, og, faces, cage_infinite_fp);
  calc_edge_out_error(rm, om, og, edges, cage_infinite_fp);
#endif
  calc_face_in_error(rm, faces);
  for (VertexHandle v : vertices)
    calc_out_surround_error(rm, v);

  // update affected edges 
  std::vector<EdgeHandle> affected_edges = {
    rm->edge_handle(rm->next_halfedge_handle(flip_he)),
    rm->edge_handle(rm->prev_halfedge_handle(flip_he)),
    rm->edge_handle(rm->next_halfedge_handle(flip_he_opp)),
    rm->edge_handle(rm->prev_halfedge_handle(flip_he_opp)),
  };

  auto edge_flipper = new_edge_flipper();
  for (EdgeHandle eh : affected_edges)
  {
    update_states[eh.idx()]++;

    if (!edge_flipper.init(eh))
      continue;
    if (edge_flipper.flip_will_cause_over_valence(param->maxValence))
      continue;
    if (!edge_flipper.flip_will_decrease_valence())
      continue;
    double local_hd_before = edge_flipper.local_Hausdorff_before_flipping();
    double local_hd_after = edge_flipper.local_Hausdorff_after_flipping();
    if (local_hd_after == DBL_MAX)
      continue;

    if (allow_negtive)
    {
      if (local_hd_after < max_distance_error)
        // decrease valence and won't cause out of distance error
        edges_to_flip.emplace(eh, update_states[eh.idx()], -local_hd_after);
    }
    else  // not allow_negtive
    {
      if (local_hd_before - local_hd_after >= 0.0)
        edges_to_flip.emplace(eh, update_states[eh.idx()], local_hd_before - local_hd_after);
    }
  }
}

void FlipStage::do_flip()
{
  auto edge_flipper = new_edge_flipper();
  edge_flipper.set_flags(/*update_links*/true, /*update_target_length*/false, /*update_normals*/true);

  initialize_flip_edges_reward();
  size_t flipped_edge_num = 0;
  while (!edges_to_flip.empty())
  {
    auto edge_reward = edges_to_flip.top();
    edges_to_flip.pop();

    // out of date
    if (edge_reward.state < update_states[edge_reward.eh.idx()])
      continue;

    if (edge_flipper.try_flip_edge(edge_reward.eh))
    {
      update_after_flipping(edge_reward.eh);
      flipped_edge_num++;
    }
  }
  Logger::user_logger->info("flipped {} edges.", flipped_edge_num);
}
}// namespace CageSimp
}// namespace Cage