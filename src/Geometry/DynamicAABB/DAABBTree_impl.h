#pragma once
#include "DAABBTree.h"
#include <numeric>
#include <cassert>

namespace Cage
{
namespace Geometry
{
using namespace SimpleUtils;

template<typename Kernel>
void DAABBTree<Kernel>::init_primitives(Primitives&& primitives)
{
  clear();
  m_primitives = std::move(primitives);
  m_primitive_map2_node.resize(m_primitives.size());
  m_primitive_deleted.resize(m_primitives.size(), false);
}

template<typename Kernel>
void DAABBTree<Kernel>::init_primitives(PrimitiveIter first, PrimitiveIter beyond)
{
  clear();
  m_primitives.reserve(std::distance(first, beyond));
  while (first != beyond)
  {
    m_primitives.push_back(*first);
    first++;
  }
  m_primitive_map2_node.resize(m_primitives.size());
  m_primitive_deleted.resize(m_primitives.size(), false);
}

template<typename Kernel>
void DAABBTree<Kernel>::build()
{
  if (m_primitives.size() > 1)
  {
    std::vector<int> ids(m_primitives.size());
    std::iota(ids.begin(), ids.end(), 0);
    m_nodes.reserve(m_primitives.size() * 2);
    m_node_deleted.reserve(m_primitives.size() * 2);
    expand(ids.begin(), ids.end(), /*parent_idx*/-1);
  }
}

template<typename Kernel>
void DAABBTree<Kernel>::rebuild()
{
  m_nodes.clear();
  m_node_deleted.clear();
  build();
}

template<typename Kernel>
void DAABBTree<Kernel>::clear()
{
  m_primitives.clear();
  m_nodes.clear();
  m_primitive_deleted.clear();
  m_node_deleted.clear();
  m_primitive_map2_node.clear();
}

template<typename Kernel>
void DAABBTree<Kernel>::insert(const Primitive& p)
{
  m_primitives.push_back(p);
  m_primitive_deleted.push_back(false);
  m_primitive_map2_node.emplace_back();
  BoundingBox bbox = CalcBox()(p);
  int p_idx = (int)m_primitives.size() - 1;

  switch (m_primitives.size())
  {
  case 1: // only one just added primitive, no node exists.
  {
    m_nodes.emplace_back();
    m_node_deleted.push_back(false);
    auto& n = m_nodes.back();
    n.m_bbox = bbox;
    n.m_nb_primitives = 1;
    n.m_parent = -1;
    n.m_p_left_child = p_idx;
    n.m_p_right_child = -1;
  }
  break;
  default:
  {
    insert_to_node(p_idx, bbox, 0);
  }
  break;
  }
}

template<typename Kernel>
void DAABBTree<Kernel>::remove(int p_idx)
{
  m_primitive_deleted[p_idx] = true;
  // remove it from m_primitives when collecting garbage.
  // now remove it from tree(leaf node).
  int node_idx = mapto(p_idx);
  remove_from_node(p_idx, node_idx);
}

/// @brief After change geometry of primitive, update its bounding box
/// and propogate change to ancestors.
template<typename Kernel>
void DAABBTree<Kernel>::update(int p_idx, const Primitive& new_p)
{
  auto& p = primitive(p_idx);
  auto& n = node(mapto(p_idx));
  ASSERT(n.m_p_left_child == p_idx, "mapto error.");
  p = new_p;
  n.m_bbox = CalcBox()(p);
  update_node(n.m_parent);
}

/// @brief really delete marked primitives and nodes.
/// @attention If primitives are corresponding to a mesh,
/// garbage collection should be called after the mesh collecting garbage.
/// Only in this way can we keep order of primitives right!
template<typename Kernel>
void DAABBTree<Kernel>::collect_garbage()
{
  {// remove deleted faces
    int front = 0;
    int back = (int)m_primitive_deleted.size() - 1;

    while (true)
    {
      // find 1st deleted and last un-deleted
      while (!m_primitive_deleted[front] && front < back)++front;
      while (m_primitive_deleted[back] && front < back)--back;
      if (front >= back)break;

      // swap front and back
      std::swap(m_primitive_deleted[front], m_primitive_deleted[back]);
      primitive(front) = primitive(back);
      mapto(front) = mapto(back);
      auto& pn = node(mapto(back));
      pn.m_p_left_child = front;
    }
    m_primitives.resize(m_primitive_deleted[front] ? front : front + 1);
    m_primitive_map2_node.resize(m_primitives.size());
    m_primitive_deleted.resize(m_primitives.size());
  }
  rebuild();
  /*
  {// remove deleted nodes
    int front = 0;
    int back = (int)m_node_deleted.size() - 1;

    while (true)
    {
      // find 1st deleted and last un-deleted
      while (!m_node_deleted[front] && front < back)++front;
      while (m_node_deleted[back] && front < back)--back;
      if (front >= back)break;

      // swap front and back
      std::swap(m_node_deleted[front], m_node_deleted[back]);
      swap_node(front, back);
    }
    m_nodes.resize(m_node_deleted[front] ? front : front + 1);
    m_node_deleted.resize(m_nodes.size());
  }
  */
}

template<typename Kernel>
template<typename Trait>
void DAABBTree<Kernel>::traversal(Trait& traits)
{
  switch (m_primitives.size())
  {
  case 0:
    break;
  case 1:
    traits.intersection(0);
    break;
  default:
    traversal_node(traits, 0); // root node's index is always 0.
  }
}

template<typename Kernel>
template<typename Trait>
bool DAABBTree<Kernel>::traversal_node(Trait& trait, int node_idx)
{
  // trait.intersection will check intersection with primitive.
  // trait.do_inter will check intersection with box.
  bool go_next = true;
  Node& n = node(node_idx);
  switch (n.m_nb_primitives)
  {
  case 1:
    return trait.intersection(n.m_p_left_child);
  default:
    if (trait.do_inter(n.m_p_left_child))
    {
      go_next = traversal_node(trait, n.m_p_left_child);
      if (go_next && trait.do_inter(n.m_p_right_child))
        go_next = traversal_node(trait, n.m_p_right_child);
      return go_next;
    }
    else if (trait.do_inter(n.m_p_right_child))
    {
      return traversal_node(trait, n.m_p_right_child);
    }
    return true;
  }
}


template<typename Kernel>
BoundingBox DAABBTree<Kernel>::calc_bbox(IndexIter first, IndexIter beyond)
{
  BoundingBox bbox = CalcBox()(primitive(*first));
  for (++first; first != beyond; ++first)
  {
    bbox += CalcBox()(primitive(*first));
  }
  return bbox;
}

template<typename Kernel>
int DAABBTree<Kernel>::split_primitives(IndexIter first, IndexIter beyond, const BoundingBox& box)
{
  IndexIter middle = first + (beyond - first) / 2;
  size_t split_dim = box.longest_axis();

  SplitPred pred(split_dim);
  auto wrapped_pred = [&](int id0, int id1)
  {
    return pred(primitive(id0), primitive(id1));
  };

  std::nth_element(first, middle, beyond, wrapped_pred);
  return (int)split_dim;
}

template<typename Kernel>
int DAABBTree<Kernel>::expand(IndexIter first, IndexIter beyond, int parent_idx)
{
  m_nodes.emplace_back();
  m_node_deleted.push_back(false);
  DAABBNode& n = m_nodes.back();
  int node_idx = (int)m_nodes.size() - 1;

  n.m_parent = parent_idx;
  n.m_nb_primitives = std::distance(first, beyond);

  n.m_bbox = calc_bbox(first, beyond);
  n.m_split_dim = split_primitives(first, beyond, n.m_bbox);

  switch (n.m_nb_primitives)
  {
  case 1:
    n.m_p_left_child = *first;
    n.m_p_right_child = -1;
    mapto(*first) = node_idx;
    break;
  default:
    IndexIter middle = first + n.m_nb_primitives / 2;
    node(node_idx).m_p_left_child = expand(first, middle, node_idx);
    node(node_idx).m_p_right_child = expand(middle, beyond, node_idx);
    break;
  }
  return node_idx;
}

template<typename Kernel>
void DAABBTree<Kernel>::insert_to_node(int p_idx, BoundingBox& pbox, int node_idx)
{
  auto& n = node(node_idx);
  n.m_nb_primitives += 1;
  switch (n.m_nb_primitives)
  {
  case 2:
  {
    // re-split three primitives
    std::vector<int> indices = { n.m_p_left_child, p_idx };
    auto first = indices.begin(), beyond = indices.end();
    n.m_bbox = calc_bbox(first, beyond);
    n.m_split_dim = split_primitives(first, beyond, n.m_bbox);

    IndexIter middle = first + n.m_nb_primitives / 2;
    node(node_idx).m_p_left_child = expand(first, middle, node_idx);
    node(node_idx).m_p_right_child = expand(middle, beyond, node_idx);
  }
  break;
  default:
  {
    auto& lnode = node(n.m_p_left_child);
    auto& rnode = node(n.m_p_right_child);

    int sd = n.m_split_dim;
    double left_max = lnode.m_bbox.max_coord(sd);
    double right_min = rnode.m_bbox.min_coord(sd);
    double insert_min = pbox.min_coord(sd);
    double insert_max = pbox.max_coord(sd);
    n.m_bbox += pbox;

    if (insert_max <= right_min)
      insert_to_node(p_idx, pbox, n.m_p_left_child);
    else if (insert_min >= left_max)
      insert_to_node(p_idx, pbox, n.m_p_right_child);
    else
    {
      // choose a branch that makes tree more balance.
      if (lnode.m_nb_primitives > rnode.m_nb_primitives)
        insert_to_node(p_idx, pbox, n.m_p_right_child);
      else
        insert_to_node(p_idx, pbox, n.m_p_left_child);
    }
  }
  break;
  }
}

template<typename Kernel>
void DAABBTree<Kernel>::remove_from_node(int p_idx, int node_idx)
{
  auto& n = node(node_idx);
  n.m_p_left_child = -1;
  n.m_nb_primitives = 0;

  // decrease nb_primitive in ancestors
  int parent_idx = n.m_parent;
  while (parent_idx != -1)
  {
    node(parent_idx).m_nb_primitives -= 1;
    parent_idx = node(parent_idx).m_parent;
  }

  // nb_primitives == 0
  remove_node(node_idx);
}

/// @brief remove a node. removed node must have no primitives.
template<typename Kernel>
void DAABBTree<Kernel>::remove_node(int node_idx)
{
  ASSERT(node(node_idx).m_nb_primitives == 0, "remove a non-empty node.");

  m_node_deleted[node_idx] = true;
  // remove node from its parent.
  // move its sibling up.
  int parent_idx = node(node_idx).m_parent;
  node(node_idx).clear();
  if (parent_idx != -1)
  {
    auto& p_node = node(parent_idx);
    if (p_node.m_p_left_child == node_idx)
    {
      p_node.m_p_left_child = -1;
      moveup_node(p_node.m_p_right_child);
    }
    else
    {
      p_node.m_p_right_child = -1;
      moveup_node(p_node.m_p_left_child);
    }
  }
}

/// @brief parent of this node only has one child, which is this node.
/// so replace parent with this node.
template<typename Kernel>
void DAABBTree<Kernel>::moveup_node(int node_idx)
{
  assert(node(node_idx).m_parent != -1);

  auto& n = node(node_idx);
  int parent_idx = n.m_parent;
  auto& p_node = node(parent_idx);
  int grandpa_idx = p_node.m_parent;

  p_node = n;
  p_node.m_parent = grandpa_idx;
  n.clear();

  if (p_node.m_nb_primitives >= 2)
  {
    // has child nodes, update parent in child nodes.
    node(p_node.m_p_left_child).m_parent = parent_idx;
    node(p_node.m_p_right_child).m_parent = parent_idx;
  }
  else if (p_node.m_nb_primitives == 1)
  {
    mapto(p_node.m_p_left_child) = parent_idx;
  }
  update_node(grandpa_idx);
  remove_node(node_idx);
}

template<typename Kernel>
void DAABBTree<Kernel>::swap_node(int node0, int node1)
{
  auto& n0 = node(node0);
  auto& n1 = node(node1);

  // swap parent's child index
  if (n0.m_parent != -1)
  {
    if (node(n0.m_parent).m_p_left_child == node0)
      node(n0.m_parent).m_p_left_child = node1;
    else
      node(n0.m_parent).m_p_right_child = node1;
  }
  if (n1.m_parent != -1)
  {
    if (node(n1.m_parent).m_p_left_child == node1)
      node(n1.m_parent).m_p_left_child = node0;
    else
      node(n1.m_parent).m_p_right_child = node0;
  }
  // swap children's parent index
  if (n0.m_nb_primitives >= 2)
  {
    node(n0.m_p_left_child).m_parent = node1;
    node(n0.m_p_right_child).m_parent = node1;
  }
  else if (n0.m_nb_primitives == 1)
  {
    mapto(n0.m_p_left_child) = node1;
  }

  if (n1.m_nb_primitives >= 2)
  {
    node(n1.m_p_left_child).m_parent = node0;
    node(n1.m_p_right_child).m_parent = node0;
  }
  else if (n1.m_nb_primitives == 1)
  {
    mapto(n1.m_p_left_child) = node0;
  }
  // swap two nodes
  std::swap(n0, n1);
}

template<typename Kernel>
void DAABBTree<Kernel>::update_node(int node_idx)
{
  if (node_idx != -1)
  {
    auto& cur_node = node(node_idx);
    BoundingBox new_box = 
      node(cur_node.m_p_left_child).m_bbox +
      node(cur_node.m_p_right_child).m_bbox;
    if (!(cur_node.m_bbox == new_box))
    {
      cur_node.m_bbox = new_box;
      update_node(cur_node.m_parent);
    }
  }
}


}// namespace Geometry
}// namespace Cage