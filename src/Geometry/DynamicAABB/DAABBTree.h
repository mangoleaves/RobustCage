#pragma once

#include "DAABBNode.h"
#include <vector>

namespace Cage
{
namespace Geometry
{
template<typename Kernel>
class DAABBTree
{
protected:
  typedef Kernel K;

  typedef typename K::Primitive Primitive;
  typedef std::vector<Primitive> Primitives;
  typedef typename Primitives::iterator PrimitiveIter;
  typedef Primitive* PrimitivePtr;

  typedef typename K::SplitPred SplitPred;

  typedef typename K::CalcBox CalcBox;

  typedef std::vector<int> Indices;
  typedef typename Indices::iterator IndexIter;

  typedef DAABBNode Node;
  typedef std::vector<Node> Nodes;
public:
  Primitives m_primitives;
  Nodes m_nodes;
  std::vector<int> m_primitive_map2_node;

  std::vector<uint8_t> m_primitive_deleted;
  std::vector<uint8_t> m_node_deleted;
public:
  DAABBTree() {}

  DAABBTree(PrimitiveIter first, PrimitiveIter second)
  {
    init_primitives(first, second);
    build();
  }

  DAABBTree(Primitives&& primitives)
  {
    init_primitives(primitives);
    build();
  }

  void init_primitives(Primitives&& primitives);
  void init_primitives(PrimitiveIter first, PrimitiveIter second);
  void build();
  void rebuild();

  void clear();

  inline size_t size() const { return m_primitives.size(); }
  inline bool empty() const { return m_primitives.empty(); }

  virtual void insert(const Primitive& p);
  virtual void remove(int p_idx);
  virtual void update(int p_idx, const Primitive& new_p);
  void collect_garbage();

  template<typename Trait>
  void traversal(Trait& trait);
protected:
  inline Node& node(int node_idx) { return m_nodes[node_idx]; }
  inline Primitive& primitive(int p_idx) { return m_primitives[p_idx]; }
  inline int& mapto(int p_idx) { return m_primitive_map2_node[p_idx]; }

  // build functions
  BoundingBox calc_bbox(IndexIter first, IndexIter beyond);
  int split_primitives(IndexIter first, IndexIter beyond, const BoundingBox& box);
  int expand(IndexIter first, IndexIter beyond, int parent_idx);

  // insert functions
  void insert_to_node(int p_idx, BoundingBox& pbox, int node_idx);

  // delete functions
  void remove_from_node(int p_idx, int node_idx);
  void remove_node(int node_idx);
  void moveup_node(int node_idx);

  // update functions
  void swap_node(int n0, int n1);
  void update_node(int node_idx);

  // traversal functions
  template<typename Trait>
  bool traversal_node(Trait& trait, int node_idx);
};
}// namespace Geometry
}// namespace Cage
#include "DAABBTree_impl.h"