#include "EdgeFlipper.h"

namespace Cage
{
namespace CageSimp
{

bool EdgeFlipper::init(EdgeHandle e)
{
  clear();

  if (!rm->is_flip_ok(e))
    return false;

  init_handles(e);

  initialized = true;
  return true;
}

void EdgeFlipper::set_flags(
  bool _f_update_links, bool _f_update_target_length,
  bool _f_update_normals
)
{
  f_update_links = _f_update_links;
  f_update_target_length = _f_update_target_length;
  f_update_normals = _f_update_normals;
}

/// @brief A general purpose function. try flip edge. it will
/// check intersection constriant.
bool EdgeFlipper::try_flip_edge()
{
  ASSERT(initialized, "edge flipper not initialized.");

  // check and backup before flipping
  if (flip_would_cause_intersection())
    return false;
  if (f_update_links)
    backup_in_links();

  // real flip
  rm->flip(rm->edge_handle(a0));

  // update after flipping
  update_one_ring_faces();
  if (f_update_target_length)
    update_target_len();
  update_length_and_area();
  if (f_update_normals)
    update_normals();
  update_tree();
  if (f_update_links)
    generate_in_links();
  return true;
}

bool EdgeFlipper::try_flip_edge(EdgeHandle e)
{
  if (!init(e))
    return false;
  return try_flip_edge();
}

bool EdgeFlipper::try_flip_edge_decrease_valence(EdgeHandle e)
{
  if (!init(e))
    return false;

  if (!flip_will_decrease_valence())
    return false;

  return try_flip_edge();
}

bool EdgeFlipper::flip_will_cause_small_large_angle()const
{
  double la0 = (rm->point(va1) - rm->point(vb1)).length();
  double lb0 = la0;
  double la1 = rm->data(rm->edge_handle(a1)).edge_length;
  double lb1 = rm->data(rm->edge_handle(b1)).edge_length;
  double la2 = rm->data(rm->edge_handle(a2)).edge_length;
  double lb2 = rm->data(rm->edge_handle(b2)).edge_length;

  double la0_sqr = la0 * la0;
  double lb0_sqr = la0_sqr;
  double la1_sqr = la1 * la1;
  double lb1_sqr = lb1 * lb1;
  double la2_sqr = la2 * la2;
  double lb2_sqr = lb2 * lb2;

  double sa0 = (la0_sqr + la2_sqr - lb1_sqr) / (2. * la0 * la2);
  if (sa0 < -GeomThr::cos_de_thr || sa0 > GeomThr::cos_de_thr) return true;
  double sa1 = (la1_sqr + lb0_sqr - lb2_sqr) / (2. * la1 * lb0);
  if (sa1 < -GeomThr::cos_de_thr || sa0 > GeomThr::cos_de_thr) return true;
  double sa2 = (la2_sqr + lb1_sqr - la0_sqr) / (2. * la2 * lb1);
  if (sa2 < -GeomThr::cos_de_thr || sa0 > GeomThr::cos_de_thr) return true;
  double sb0 = (lb0_sqr + lb2_sqr - la1_sqr) / (2. * lb0 * lb2);
  if (sb0 < -GeomThr::cos_de_thr || sa0 > GeomThr::cos_de_thr) return true;
  double sb1 = (lb1_sqr + la0_sqr - la2_sqr) / (2. * lb1 * la0);
  if (sb1 < -GeomThr::cos_de_thr || sa0 > GeomThr::cos_de_thr) return true;
  double sb2 = (lb2_sqr + la1_sqr - lb0_sqr) / (2. * lb2 * la1);
  if (sb2 < -GeomThr::cos_de_thr || sa0 > GeomThr::cos_de_thr) return true;

  return false;
}

bool EdgeFlipper::flip_will_decrease_valence()const
{
  int target_val_a0 = rm->is_boundary(va0) ? 4 : 6;
  int target_val_a1 = rm->is_boundary(va1) ? 4 : 6;
  int target_val_b0 = rm->is_boundary(vb0) ? 4 : 6;
  int target_val_b1 = rm->is_boundary(vb1) ? 4 : 6;

  int deviation_pre =
    std::abs((int)(rm->valence(va0)) - target_val_a0) +
    std::abs((int)(rm->valence(va1)) - target_val_a1) +
    std::abs((int)(rm->valence(vb0)) - target_val_b0) +
    std::abs((int)(rm->valence(vb1)) - target_val_b1);
  int deviation_post =
    std::abs((int)(rm->valence(va0)) - 1 - target_val_a0) +
    std::abs((int)(rm->valence(va1)) + 1 - target_val_a1) +
    std::abs((int)(rm->valence(vb0)) - 1 - target_val_b0) +
    std::abs((int)(rm->valence(vb1)) + 1 - target_val_b1);

  return deviation_post < deviation_pre;
}

bool EdgeFlipper::flip_will_cause_over_valence(size_t max_valence)const
{
  return rm->valence(va1) + 1 > max_valence || rm->valence(vb1) + 1 > max_valence;
}

double EdgeFlipper::local_Hausdorff_before_flipping()const
{
  ASSERT(initialized, "flipper not initialized.");

  double hd = 0.0;
  // two faces in/out error
  hd = std::max(rm->data(fa).out_error, rm->data(fa).in_error);
  if (rm->data(fb).in_error > hd) hd = rm->data(fb).in_error;
  if (rm->data(fb).out_error > hd) hd = rm->data(fb).out_error;
  if (rm->data(rm->edge_handle(a0)).out_error > hd) hd = rm->data(rm->edge_handle(a0)).out_error;
  // four edges out error
  if (rm->data(rm->edge_handle(a1)).out_error > hd) hd = rm->data((rm->edge_handle(a1))).out_error;
  if (rm->data(rm->edge_handle(a2)).out_error > hd) hd = rm->data(rm->edge_handle(a2)).out_error;
  if (rm->data(rm->edge_handle(b1)).out_error > hd) hd = rm->data(rm->edge_handle(b1)).out_error;
  if (rm->data(rm->edge_handle(b2)).out_error > hd) hd = rm->data(rm->edge_handle(b2)).out_error;
  // four vertices out error
  if (rm->data(va0).out_error > hd) hd = rm->data(va0).out_error;
  if (rm->data(va1).out_error > hd) hd = rm->data(va1).out_error;
  if (rm->data(vb0).out_error > hd) hd = rm->data(vb0).out_error;
  if (rm->data(vb1).out_error > hd) hd = rm->data(va1).out_error;
  return hd;
}

/// @brief calculate local hausdorff distance after flipping.
/// @return DBL_MAX if flipping will violate constraints.
double EdgeFlipper::local_Hausdorff_after_flipping()const
{
  ASSERT(initialized, "flipper not initialized.");

  if (flip_would_cause_intersection())
    return DBL_MAX;

  double hd = 0.0;

  // hausdorff of "in"
  auto& pa0 = rm->point(va0);
  auto& pa1 = rm->point(va1);
  auto& pb0 = rm->point(vb0);
  auto& pb1 = rm->point(vb1);
  Geometry::HeavyTriangle tri_a(pa1, pb0, pb1);
  Geometry::HeavyTriangle tri_b(pa0, pa1, pb1);

#define hd_for_links(links)\
  for (Link& link : links) {\
    double dis_a = (tri_a.closest_point(link.first) - link.first).squaredNorm();\
    double dis_b = (tri_b.closest_point(link.first) - link.first).squaredNorm();\
    hd = std::max(hd, std::sqrt(std::min(dis_a, dis_b)));}

  hd_for_links(rm->data(fa).face_in_links);
  hd_for_links(rm->data(fb).face_in_links);
#undef hd_for_links

  // hausdorff of  "out"
  SMeshT local_rm;
  VertexHandle la0 = local_rm.add_vertex(pa0), la1 = local_rm.add_vertex(pa1),
    lb0 = local_rm.add_vertex(pb0), lb1 = local_rm.add_vertex(pb1);
  local_rm.add_face(la0, la1, lb1);
  local_rm.add_face(lb0, lb1, la1);

#ifdef USE_TREE_SEARCH
  hd = std::max(hd, calc_out_error(&local_rm, om, ot, cage_infinite_fp));
#else
  hd = std::max(hd, calc_out_error(&local_rm, om, og, cage_infinite_fp));
#endif
  return hd;
}

void EdgeFlipper::clear()
{
  initialized = false;
  faces_in_links.clear();
}

void EdgeFlipper::init_handles(EdgeHandle e)
{
  a0 = rm->halfedge_handle(e, 0);
  b0 = rm->halfedge_handle(e, 1);

  a1 = rm->next_halfedge_handle(a0);
  a2 = rm->next_halfedge_handle(a1);

  b1 = rm->next_halfedge_handle(b0);
  b2 = rm->next_halfedge_handle(b1);

  va0 = rm->to_vertex_handle(a0);
  va1 = rm->to_vertex_handle(a1);

  vb0 = rm->to_vertex_handle(b0);
  vb1 = rm->to_vertex_handle(b1);

  fa = rm->face_handle(a0);
  fb = rm->face_handle(b0);
}

bool EdgeFlipper::flip_would_cause_intersection()const
{
  auto& pa0 = rm->point(va0);
  auto& pa1 = rm->point(va1);
  auto& pb0 = rm->point(vb0);
  auto& pb1 = rm->point(vb1);
  const auto* epa0 = rm->data(va0).ep.get();
  const auto* epa1 = rm->data(va1).ep.get();
  const auto* epb0 = rm->data(vb0).ep.get();
  const auto* epb1 = rm->data(vb1).ep.get();
  // 0. check degenerate
  {
    // check triangle va1-vb0-vb1
    if (are_points_colinear(pa1, pb0, pb1, epa1, epb0, epb1))
      return true;
    // check triangle va0-va1-vb1
    if (are_points_colinear(pa0, pa1, pb1, epa0, epa1, epb1))
      return true;
  }
  // 1. check intersect with original mesh.
  {
    // check triangle va1-vb0-vb1
    if (ot->do_intersect(pa1, pb0, pb1, epa1, epb0, epb1))
      return true;
    // check triangle va0-va1-vb1
    if (ot->do_intersect(pa0, pa1, pb1, epa0, epa1, epb1))
      return true;
  }
  // 2. check overlap triangles with a common edge.
  {
  #define wrapped_check(edge, p0, p1, p2, ep0, ep1, ep2)\
  if (!rm->is_boundary(rm->opposite_halfedge_handle(edge))) {\
    auto opp_vh = rm->opposite_vh(rm->opposite_halfedge_handle(edge));\
    auto& p_oppo_edge = rm->point(opp_vh);\
    auto* ep_opp_edge = rm->data(opp_vh).ep.get();\
    if (Geometry::triangle_do_overlap(p0, p1, p2, p_oppo_edge, ep0, ep1, ep2, ep_opp_edge))\
      return true;\
  }

    wrapped_check(a1, pb1, pa0, pa1, epb1, epa0, epa1);
    wrapped_check(a2, pb1, pa1, pb0, epb1, epa1, epb0);
    wrapped_check(b1, pa1, pb0, pb1, epa1, epb0, epb1);
    wrapped_check(b2, pa1, pb1, pa0, epa1, epb1, epa0);
  #undef wrapped_check
  }
  // 3. check intersect triangles with a common vertex.
  {
    { // (1) check faces have common vertex va1
      FaceHandle fab_ignore = rm->face_handle(a0);
      FaceHandle fb_ignore = rm->opposite_face_handle(a1);
      FaceHandle fa_ignore = rm->opposite_face_handle(a2);

      for (FaceHandle vf : rm->vf_range(va1))
      {
        if (vf == fab_ignore) continue;

        VertexHandle to_v, from_v;
        for (auto fhe : rm->fh_range(vf))
        {
          to_v = rm->to_vertex_handle(fhe);
          from_v = rm->from_vertex_handle(fhe);
          if (to_v != va1 && from_v != va1) break;
        }
        Vec3d& to_p = rm->point(to_v);
        Vec3d& from_p = rm->point(from_v);
        ExactPoint* to_ep = rm->data(to_v).ep.get();
        ExactPoint* from_ep = rm->data(from_v).ep.get();

        if (vf != fa_ignore &&
          Geometry::triangle_do_intersect(pb0, pb1, pa1, from_p, to_p, epb0, epb1, epa1, from_ep, to_ep))
          return true;

        if (vf != fb_ignore &&
          Geometry::triangle_do_intersect(pb1, pa0, pa1, from_p, to_p, epb1, epa0, epa1, from_ep, to_ep))
          return true;
      }
    }
    {// (2) check faces have common vertex vb1
      FaceHandle fab_ignore = rm->face_handle(b0);
      FaceHandle fb_ignore = rm->opposite_face_handle(b2);
      FaceHandle fa_ignore = rm->opposite_face_handle(b1);

      for (FaceHandle vf : rm->vf_range(vb1))
      {
        if (vf == fab_ignore) continue;

        VertexHandle to_v, from_v;
        for (auto fhe : rm->fh_range(vf))
        {
          to_v = rm->to_vertex_handle(fhe);
          from_v = rm->from_vertex_handle(fhe);
          if (to_v != vb1 && from_v != vb1) break;
        }
        Vec3d& to_p = rm->point(to_v);
        Vec3d& from_p = rm->point(from_v);
        ExactPoint* to_ep = rm->data(to_v).ep.get();
        ExactPoint* from_ep = rm->data(from_v).ep.get();

        if (vf != fa_ignore &&
          Geometry::triangle_do_intersect(pa1, pb0, pb1, from_p, to_p, epa1, epb0, epb1, from_ep, to_ep))
          return true;

        if (vf != fb_ignore &&
          Geometry::triangle_do_intersect(pa0, pa1, pb1, from_p, to_p, epa0, epa1, epb1, from_ep, to_ep))
          return true;
      }
    }
    {// (3) check faces have common vertex vb0
      std::set<FaceHandle> one_ring_faces;
      for (FaceHandle vf : rm->vf_range(vb0))one_ring_faces.insert(vf);
      std::set<FaceHandle> ignore_faces(
        { fa,fb,rm->opposite_face_handle(a2),rm->opposite_face_handle(b1) });
      std::vector<FaceHandle> check_faces(one_ring_faces.size());
      auto it = std::set_difference(
        one_ring_faces.begin(), one_ring_faces.end(),
        ignore_faces.begin(), ignore_faces.end(),
        check_faces.begin());
      check_faces.resize(it - check_faces.begin());
      for (FaceHandle f : check_faces)
      {
        VertexHandle to_v, from_v;
        for (auto fhe : rm->fh_range(f))
        {
          to_v = rm->to_vertex_handle(fhe);
          from_v = rm->from_vertex_handle(fhe);
          if (to_v != vb0 && from_v != vb0) break;
        }
        Vec3d& to_p = rm->point(to_v);
        Vec3d& from_p = rm->point(from_v);
        ExactPoint* to_ep = rm->data(to_v).ep.get();
        ExactPoint* from_ep = rm->data(from_v).ep.get();
        if (Geometry::triangle_do_intersect(pb1, pa1, pb0, from_p, to_p, epb1, epa1, epb0, from_ep, to_ep))
          return true;
      }
    }
    {// (4) check faces have common vertex va0
      std::set<FaceHandle> one_ring_faces;
      for (FaceHandle vf : rm->vf_range(va0))one_ring_faces.insert(vf);
      std::set<FaceHandle> ignore_faces(
        { fa,fb,rm->opposite_face_handle(a1),rm->opposite_face_handle(b2) });
      std::vector<FaceHandle> check_faces(one_ring_faces.size());
      auto it = std::set_difference(
        one_ring_faces.begin(), one_ring_faces.end(),
        ignore_faces.begin(), ignore_faces.end(),
        check_faces.begin());
      check_faces.resize(it - check_faces.begin());
      for (FaceHandle f : check_faces)
      {
        VertexHandle to_v, from_v;
        for (auto fhe : rm->fh_range(f))
        {
          to_v = rm->to_vertex_handle(fhe);
          from_v = rm->from_vertex_handle(fhe);
          if (to_v != va0 && from_v != va0) break;
        }
        Vec3d& to_p = rm->point(to_v);
        Vec3d& from_p = rm->point(from_v);
        ExactPoint* to_ep = rm->data(to_v).ep.get();
        ExactPoint* from_ep = rm->data(from_v).ep.get();
        if (Geometry::triangle_do_intersect(pa1, pb1, pa0, from_p, to_p, epa1, epb1, epa0, from_ep, to_ep))
          return true;
      }
    }
  }
  // 4. check intersect non-adjacent triangles
  {
    //std::set<FaceHandle> one_ring_faces;
    //one_ring_faces.insert(fa);
    //one_ring_faces.insert(fb);
    //extend_faces_by_one_ring(rm, one_ring_faces);
    std::set<FaceHandle> one_ring_faces = rm->data(fa).one_ring_faces;
    one_ring_faces.insert(rm->data(fb).one_ring_faces.begin(), rm->data(fb).one_ring_faces.end());

    std::vector<int> _one_ring_faces(one_ring_faces.size());
    std::transform(one_ring_faces.begin(), one_ring_faces.end(),
      _one_ring_faces.begin(), [&](FaceHandle fh) {return fh.idx();});

    // check triangle va1-vb0-vb1
    if (lrt->do_intersect(pa1, pb0, pb1, epa1, epb0, epb1, _one_ring_faces))
      return true;
    // check triangle va0-va1-vb1
    if (lrt->do_intersect(pa0, pa1, pb1, epa0, epa1, epb1, std::move(_one_ring_faces)))
      return true;
  }
  return false;
}

void EdgeFlipper::backup_in_links()
{
  faces_in_links = backup_local_in_links(rm, { fa, fb });
}

void EdgeFlipper::generate_in_links()
{
  auto& pa0 = rm->point(va0);
  auto& pa1 = rm->point(va1);
  auto& pb0 = rm->point(vb0);
  auto& pb1 = rm->point(vb1);

  Geometry::HeavyTriangle tri_a(pa1, pb0, pb1);
  Geometry::HeavyTriangle tri_b(pa0, pa1, pb1);

  for (Link& link : faces_in_links)
  {
    Vec3d& first = link.first;
    auto cp_a = tri_a.closest_point(first);
    auto cp_b = tri_b.closest_point(first);
    if ((cp_a - first).squaredNorm() < (cp_b - first).squaredNorm())
    {
      link.second = cp_a;
      rm->data(fa).face_in_links.push_back(link);
    }
    else
    {
      link.second = cp_b;
      rm->data(fb).face_in_links.push_back(link);
    }
  }
}

void EdgeFlipper::update_target_len()
{
  // update edge target length
  rm->data(rm->edge_handle(a0)).target_length = std::min(
    rm->data(va1).target_length,
    rm->data(vb1).target_length
  );
}

void EdgeFlipper::update_length_and_area()
{
  EdgeHandle flipped_edge = rm->edge_handle(a0);
  rm->data(flipped_edge).edge_length = rm->calc_edge_length(flipped_edge);

  rm->data(fa).face_area = calc_face_area(rm, fa);
  rm->data(fb).face_area = calc_face_area(rm, fb);
}

void EdgeFlipper::update_normals()
{
  rm->update_normal(fa);
  rm->update_normal(fb);
  rm->update_normal(va0);
  rm->update_normal(va1);
  rm->update_normal(vb0);
  rm->update_normal(vb1);
}

void EdgeFlipper::update_tree()
{
  lrt->update(rm, fa);
  lrt->update(rm, fb);
}

void EdgeFlipper::update_one_ring_faces()
{
  for (FaceHandle vfh : rm->vf_range(va0))
    rm->data(vfh).one_ring_faces.erase(fa);
  for (FaceHandle vfh : rm->vf_range(vb0))
    rm->data(vfh).one_ring_faces.erase(fb);
  for (FaceHandle vfh : rm->vf_range(va1))
  {
    rm->data(vfh).one_ring_faces.insert(fa);
    rm->data(vfh).one_ring_faces.insert(fb);
  }
  for (FaceHandle vfh : rm->vf_range(vb1))
  {
    rm->data(vfh).one_ring_faces.insert(fa);
    rm->data(vfh).one_ring_faces.insert(fb);
  }
  init_one_ring_faces(rm, fa);
  init_one_ring_faces(rm, fb);
}

}// namespace CageSimp
}// namespace Cage