#include "KdTree.h"

namespace Cage
{
namespace Geometry
{
// region KdBox
void KdBox::update_from_point_pointers(KdVec3dPtrIter begin, KdVec3dPtrIter end)
{
  if (begin == end)
    return;

  minBound = (**begin).vec;
  maxBound = (**begin).vec;
  for (begin++; begin != end; begin++)
  {
    (*this) += ((**begin).vec);
  }
}
// end region KdBox

// region PointContainer 
void PointContainer::split(PointContainer& c_low, Separator& sep, bool sliding)
{
  c_low.m_bbox = m_bbox;
  m_build_dim = sep.cut_dim;
  c_low.m_build_dim = sep.cut_dim;

  KdVec3dPtrIter it = std::partition(m_begin, m_end,
    [&](KdVec3dPtr v) {return v->vec[sep.cut_dim] < sep.cut_val; });
  // now [begin,it) are lower and [it,end) are upper
  if (sliding)    // avoid empty lists
  {
    if (it == m_begin)
    {
      KdVec3dPtrIter min_elt = std::min_element(m_begin, m_end,
        [&](KdVec3dPtr lhs, KdVec3dPtr rhs) {return lhs->vec[sep.cut_dim] < rhs->vec[sep.cut_dim]; });
      if (min_elt != it)
      {
        std::iter_swap(min_elt, it);
      }
      sep.cut_val = (**it)[sep.cut_dim];
      it++;
    }
    if (it == m_end)
    {
      KdVec3dPtrIter max_elt = std::max_element(m_begin, m_end,
        [&](KdVec3dPtr lhs, KdVec3dPtr rhs) {return lhs->vec[sep.cut_dim] < rhs->vec[sep.cut_dim]; });
      it--;
      if (max_elt != it)
      {
        std::iter_swap(max_elt, it);
      }
      sep.cut_val = (**it)[sep.cut_dim];
      it++;
    }
  }

  c_low.m_begin = m_begin;
  c_low.m_end = it;
  m_begin = it;
  // adjusting boxes
  m_bbox.min()[sep.cut_dim] =  sep.cut_val;
  m_tbox.update_from_point_pointers(m_begin, m_end);
  c_low.m_bbox.max()[sep.cut_dim] = sep.cut_val;
  c_low.m_tbox.update_from_point_pointers(c_low.m_begin, c_low.m_end);
}
// endregion PointContainer

// region KdTree
void KdTree::insert(const std::vector<Vec3d>& points, const std::vector<unsigned int>& ids)
{
  pts.reserve(points.size());
  for (size_t i = 0; i < points.size(); i++)
    pts.emplace_back(points[i], ids[i]);
}
void KdTree::build()
{
  data.reserve(pts.size());
  for (size_t i = 0; i < pts.size(); i++)
  {
    data.emplace_back(&pts[i]);
  }

  PointContainer c(data.begin(), data.end());
  bbox = c.bbox();
  if (c.size() <= bucket_size)
  {
    tree_root = create_leaf_node(c);
  }
  else
  {
    tree_root = new_internal_node();
    create_internal_node(tree_root, c);
  }

  KdVec3ds pts_tmp;
  pts_tmp.reserve(pts.size());
  for (size_t i = 0; i < pts.size(); i++)
  {
    pts_tmp.emplace_back(*data[i]);
  }
  for (size_t i = 0; i < leaf_nodes.size(); i++)
  {
    std::ptrdiff_t tmp = leaf_nodes[i].data - pts.begin();
    leaf_nodes[i].data = pts_tmp.begin() + tmp;
  }
  pts.swap(pts_tmp);

  data.clear();
}

const KdVec3d& KdTree::search_nearest_point(const Vec3d& query)
{
  OrthogonalNearestSeach ons(this, KdVec3d(query));
  KdVec3dIter cp = ons.closest_point_iter();
  return *cp;
}

KdTree::NodePtr KdTree::create_leaf_node(PointContainer& c)
{
  LeafNode node(c.size());
  std::ptrdiff_t tmp = c.begin() - data.begin();
  node.data = pts.begin() + tmp;

  leaf_nodes.emplace_back(node);
  return &(leaf_nodes.back());
}

KdTree::NodePtr KdTree::new_internal_node()
{
  internal_nodes.emplace_back();
  return &(internal_nodes.back());
}

void KdTree::create_internal_node(NodePtr n, PointContainer& c)
{
  InternalPtr nh = static_cast<InternalPtr>(n);

  Separator sep;
  PointContainer c_low;

  split(sep, c, c_low);
  nh->set_separator(sep);

  handle_extended_node(nh, c, c_low);

  if (c_low.size() > bucket_size)
  {
    nh->m_lower_ch = new_internal_node();
    create_internal_node(nh->m_lower_ch, c_low);
  }
  else
  {
    nh->m_lower_ch = create_leaf_node(c_low);
  }

  if (c.size() > bucket_size)
  {
    nh->m_upper_ch = new_internal_node();
    create_internal_node(nh->m_upper_ch, c);
  }
  else
  {
    nh->m_upper_ch = create_leaf_node(c);
  }
}

void KdTree::split(Separator& sep, PointContainer& c_origin, PointContainer& c_low)
{
  size_t cutdim = c_origin.m_bbox.longest_axis();

  // Fix degenerated cases
  if (c_origin.m_tbox.min()[cutdim] != c_origin.m_tbox.max()[cutdim])
  {
    sep = Separator(cutdim,
      (c_origin.m_bbox.max()[cutdim] + c_origin.m_bbox.min()[cutdim]) / 2.0);
  }
  else
  {
    cutdim = c_origin.m_tbox.longest_axis();
    sep = Separator(cutdim,
      (c_origin.m_tbox.max()[cutdim] + c_origin.m_tbox.min()[cutdim]) / 2.0);
  }

  double max_span_upper = c_origin.m_tbox.max()[cutdim];
  double max_span_lower = c_origin.m_tbox.min()[cutdim];
  if (max_span_upper <= sep.cut_val)
  {
    sep.cut_val = max_span_upper;
  }
  if (max_span_lower >= sep.cut_val)
  {
    sep.cut_val = max_span_lower;
  }

  c_origin.split(c_low, sep, true);
}

void KdTree::handle_extended_node(InternalPtr nh, PointContainer& c, PointContainer& c_low)
{
  size_t cut_dim = nh->m_cut_dim;
  if (c_low.size() > 0)
  {
    nh->m_lower_low_val = c_low.m_tbox.min()[cut_dim];
    nh->m_lower_high_val = c_low.m_tbox.max()[cut_dim];
  }
  else
  {
    nh->m_lower_low_val = nh->m_cut_val;
    nh->m_lower_high_val = nh->m_cut_val;
  }
  if (c.size() > 0)
  {
    nh->m_upper_low_val = c.m_tbox.min()[cut_dim];
    nh->m_upper_high_val = c.m_tbox.max()[cut_dim];
  }
  else
  {
    nh->m_upper_low_val = nh->m_cut_val;
    nh->m_upper_high_val = nh->m_cut_val;
  }
}

KdTree::OrthogonalNearestSeach::OrthogonalNearestSeach(KdTree* tree, const KdVec3d& query)
  :m_tree(tree), m_query(query), m_dists(0, 0, 0), m_square_distance(DBL_MAX)
{
  if (m_tree->empty()) return;

  double distance_to_root = m_tree->bbox.distance_to_point(m_query.vec, m_dists.vec);
  compute_nearest_neightbors_orthogonally(m_tree->tree_root, distance_to_root);
}

void KdTree::OrthogonalNearestSeach::compute_nearest_neightbors_orthogonally(NodePtr N, double rd)
{
  if (N->is_leaf)
  {
    LeafPtr node = static_cast<LeafPtr>(N);
    if (node->n > 0)
    {
      search_nearest_in_leaf(node);
    }
  }
  else
  {
    InternalPtr node = static_cast<InternalPtr>(N);
    size_t cut_dim = node->m_cut_dim;

    NodePtr best_ch, other_ch;

    double offset;
    double val = m_query[cut_dim];
    double diff1 = val - node->m_upper_low_val;
    double diff2 = val - node->m_lower_high_val;

    if (diff1 + diff2 < 0)
    {
      offset = diff1;
      best_ch = node->m_lower_ch;
      other_ch = node->m_upper_ch;
    }
    else
    {
      offset = diff2;
      best_ch = node->m_upper_ch;
      other_ch = node->m_lower_ch;
    }

    compute_nearest_neightbors_orthogonally(best_ch, rd);
    double dst = m_dists[cut_dim];
    double new_rd = new_square_distance(rd, dst, offset);
    m_dists[cut_dim] = offset;
    if (new_rd < m_square_distance)
    {
      compute_nearest_neightbors_orthogonally(other_ch, new_rd);
    }
    m_dists[cut_dim] = dst;
  }
}

void KdTree::OrthogonalNearestSeach::search_nearest_in_leaf(LeafPtr node)
{
  for (KdVec3dIter begin = node->data, end = node->data + node->n; begin != end; begin++)
  {
    double distance = ((*begin).vec - m_query.vec).squaredNorm();
    if (distance < m_square_distance)
    {
      m_result = begin;
      m_square_distance = distance;
    }
  }
}

double KdTree::OrthogonalNearestSeach::new_square_distance(double dist, double old_off, double new_off)
{
  return dist + (new_off * new_off - old_off * old_off);
}
// endregion KdTree
}// namespace Geometry
}// namespace Cage