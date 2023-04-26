#include "DegenerationRemover.h"

namespace Cage
{
namespace CageSimp
{

void DegenerationRemover::perform()
{
  eliminate_almost_degeneration();

  rm->garbage_collection();
  lrt->collect_garbage();
  init_one_ring_faces(rm);
}

bool DegenerationRemover::is_face_almost_degenerate(FaceHandle fh)
{
#define HE(p) p->first
#define COS(p) p->second

  auto sectors_angle = face_sector_angle(rm, fh);
  // 1. find min and max angle
  // angle_a > angle_b === cos(angle_a) < cos(angle_b) when angle in [0, pi]
  auto max_angle = sectors_angle.begin(), min_angle = sectors_angle.begin();
  for (auto sa = sectors_angle.begin() + 1;sa != sectors_angle.end();sa++)
  {
    if (COS(sa) < COS(max_angle))max_angle = sa;
    else if (COS(sa) > COS(min_angle))min_angle = sa;
  }
  HalfedgeHandle he_opp_max_angle = rm->prev_halfedge_handle(HE(max_angle));
  EdgeHandle e_opp_max_angle = rm->edge_handle(he_opp_max_angle);
  HalfedgeHandle he_opp_min_angle = rm->prev_halfedge_handle(HE(min_angle));
  EdgeHandle e_opp_min_angle = rm->edge_handle(he_opp_min_angle);
  // 2. for different cases, use different operations.
  if (COS(min_angle) > GeomThr::cos_de_thr &&
    rm->data(e_opp_min_angle).edge_length <
    rm->data(e_opp_max_angle).edge_length * GeomThr::edge_length_ratio_de_thr
    )
  {
    // case (1): Has a small angle and a short edge
    return true;
  }
  else if (COS(max_angle) < -GeomThr::cos_de_thr)
  {
    // case (2): Has a large angle but no short edge
    return true;
  }
  return false;
#undef HE
#undef COS
}

/// @brief try to eliminate small angles by collapsing opposite edge.
/// try to eliminate large angles by flip opposite edge.
size_t DegenerationRemover::eliminate_almost_degeneration()
{
  auto edge_flipper = new_edge_flipper();
  edge_flipper.set_flags(
    /*update_links*/false,
    /*update_target_length*/false,
    /*update_normals*/false
  );
  auto edge_collapser = new_edge_collapser();
  edge_collapser.set_flags(
    /*update_links*/false,
    /*update_target_length*/false,
    /*update_normals*/false,
    /*check_wrinkle*/false,
    /*check_selfinter*/true,
    /*check_inter*/true
  );

#define HE(p) p->first
#define COS(p) p->second
  // record number of eliminated degeneration each iteration
  size_t eliminate_case = 1;
  size_t total_eliminate_case = 0;
  while (eliminate_case > 0)
  {
    eliminate_case = 0;
    for (FaceHandle fh : rm->faces())
    {
      if (rm->status(fh).deleted())
        continue;

      auto sectors_angle = face_sector_angle(rm, fh);
      // 1. find min and max angle
      // angle_a > angle_b === cos(angle_a) < cos(angle_b) when angle in [0, pi]
      auto max_angle = sectors_angle.begin(), min_angle = sectors_angle.begin();
      for (auto sa = sectors_angle.begin() + 1;sa != sectors_angle.end();sa++)
      {
        if (COS(sa) < COS(max_angle))max_angle = sa;
        else if (COS(sa) > COS(min_angle))min_angle = sa;
      }
      HalfedgeHandle he_opp_max_angle = rm->prev_halfedge_handle(HE(max_angle));
      EdgeHandle e_opp_max_angle = rm->edge_handle(he_opp_max_angle);
      HalfedgeHandle he_opp_min_angle = rm->prev_halfedge_handle(HE(min_angle));
      EdgeHandle e_opp_min_angle = rm->edge_handle(he_opp_min_angle);
      // 2. for different cases, use different operations.
      if (COS(min_angle) > GeomThr::cos_de_thr &&
        rm->data(e_opp_min_angle).edge_length <
        rm->data(e_opp_max_angle).edge_length * GeomThr::edge_length_ratio_de_thr
        )
      {
        // case (1): Has a small angle and a short edge, try to collapse short edge.
        if (try_collapse_almost_degenerate_edge(edge_collapser, e_opp_min_angle))
          eliminate_case++;
      }
      else if (COS(max_angle) < -GeomThr::cos_de_thr)
      {
        // case (2): Has a large angle but no short edge, try to flip loggest edge.
        edge_flipper.init(e_opp_max_angle);
        if (!edge_flipper.flip_will_cause_small_large_angle())
        {
          if (edge_flipper.try_flip_edge(e_opp_max_angle))
            eliminate_case++;
        }
      }
    }
    total_eliminate_case += eliminate_case;
  }
#undef HE
#undef COS

  Logger::user_logger->info("eliminate [{}] almost degenerate cases.", total_eliminate_case);
  return total_eliminate_case;
}

bool DegenerationRemover::try_collapse_almost_degenerate_edge(EdgeCollapser& edge_collapser, EdgeHandle eh)
{
  bool collapsed = edge_collapser.try_collapse_almost_degenerate_edge(eh);
  if (!collapsed)
  {
    // if fail to collapse because a adjacent vertex whose valence is 3,
    // we try to eliminate this vertex.
    HalfedgeHandle he = rm->halfedge_handle(eh, 0);
    VertexHandle vh = rm->opposite_vh(he);
    VertexHandle opp_vh = rm->opposite_he_opposite_vh(he);

    std::vector<VertexHandle> valence_3_vertices;
    if (rm->valence(opp_vh) == 3)
      valence_3_vertices.push_back(opp_vh);
    if (rm->valence(vh) == 3)
      valence_3_vertices.push_back(vh);

    bool adj_collapsed = !valence_3_vertices.empty();
    for (VertexHandle vh : valence_3_vertices)
    {
      HalfedgeHandle he_out = *rm->voh_begin(vh);
      adj_collapsed = adj_collapsed && edge_collapser.try_collapse_edge(
        rm->edge_handle(he_out),
        rm->point(rm->to_vertex_handle(he_out)),
        rm->data(rm->to_vertex_handle(he_out)).ep.get()
      );
    }
    if (adj_collapsed)
    {
      // we try to collapse edge.
      collapsed = edge_collapser.try_collapse_almost_degenerate_edge(eh);
    }
  }
  return collapsed;
}
}// namespace CageSimp
}// namespace Cage