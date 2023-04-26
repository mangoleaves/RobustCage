#pragma once

#include "Basic/Types.h"
#include "Basic/BoundingBox.h"
#include <algorithm>
#include <cfloat>
#include <deque>
#include <vector>

namespace Cage
{
namespace Geometry
{
class KdVec3d
{
public:
  Vec3d vec;
  unsigned int id;
public:
  KdVec3d()
    :vec()
  {}
  KdVec3d(const double& v0, const double& v1, const double& v2)
    :vec(v0, v1, v2)
  {}
  KdVec3d(const double* v)
    :vec(v)
  {}
  KdVec3d(const Vec3d& v)
    :vec(v)
  {}
  KdVec3d(const Vec3d& v, unsigned int id)
    :vec(v), id(id)
  {}
  KdVec3d& operator=(const KdVec3d& v) = default;
  inline double& operator[](size_t dim) { return vec[dim]; }
  inline const double& operator[](size_t dim)const { return vec[dim]; }
};

typedef std::vector<KdVec3d>       KdVec3ds;
typedef KdVec3ds::iterator         KdVec3dIter;
typedef KdVec3d* KdVec3dPtr;
typedef std::vector<KdVec3dPtr>    KdVec3dPtrs;
typedef KdVec3dPtrs::iterator      KdVec3dPtrIter;

class KdBox : public BoundingBox
{
public:
  KdBox() = default;
  KdBox(const KdVec3d& minB, const KdVec3d& maxB)
    :BoundingBox(minB.vec, maxB.vec)
  {}

  void update_from_point_pointers(KdVec3dPtrIter begin, KdVec3dPtrIter end);
};

class Separator
{
public:
  size_t cut_dim;
  double cut_val;

  Separator() = default;
  Separator(size_t cut_dim, double cut_value)
    :cut_dim(cut_dim), cut_val(cut_value)
  {}
};

class PointContainer
{
public:
  KdVec3dPtrIter m_begin, m_end;
  size_t m_build_dim;
  KdBox m_bbox, m_tbox;
  PointContainer() {}

  PointContainer(KdVec3dPtrIter begin, KdVec3dPtrIter end)
    :m_begin(begin), m_end(end)
  {
    if (begin != end)
    {
      m_bbox = KdBox(**begin, **begin);
      m_bbox.update_from_point_pointers(begin + 1, end);
    }
    m_tbox = m_bbox;
    m_build_dim = m_bbox.longest_axis();
  }

  inline size_t size()const { return std::distance(m_begin, m_end); }
  inline KdVec3dPtrIter begin()const { return m_begin; }
  inline KdVec3dPtrIter end()const { return m_end; }
  inline const KdBox& bbox()const { return m_bbox; }
  inline const KdBox& tbox()const { return m_tbox; }

  void split(PointContainer& c_low, Separator& sep, bool sliding = false);
};

class KdTree_Node
{
public:
  bool is_leaf;
  KdTree_Node(bool leaf)
    :is_leaf(leaf)
  {}
};

class KdTree_InternalNode : public KdTree_Node
{
public:
  typedef KdTree_Node     Node;
  typedef Node* NodePtr;

  size_t m_cut_dim;
  double m_cut_val;
  double m_lower_low_val, m_lower_high_val, m_upper_low_val, m_upper_high_val;
  NodePtr m_lower_ch, m_upper_ch;

public:
  KdTree_InternalNode()
    :KdTree_Node(false)
  {}

  inline void set_separator(const Separator& sep)
  {
    m_cut_dim = sep.cut_dim;
    m_cut_val = sep.cut_val;
  }
};

class KdTree_LeafNode : public KdTree_Node
{
public:
  size_t n;
  KdVec3dIter data;
public:
  KdTree_LeafNode()
    :KdTree_Node(true)
  {}

  KdTree_LeafNode(size_t n)
    :KdTree_Node(true), n(n)
  {}
};

class KdTree
{
public:
  typedef KdTree_Node             Node;
  typedef Node* NodePtr;
  typedef KdTree_InternalNode     InternalNode;
  typedef InternalNode* InternalPtr;
  typedef KdTree_LeafNode         LeafNode;
  typedef LeafNode* LeafPtr;
public:
  std::deque<InternalNode>    internal_nodes;
  std::deque<LeafNode>        leaf_nodes;

  NodePtr       tree_root;
  KdBox         bbox;
  KdVec3ds      pts;
  KdVec3dPtrs   data;

  size_t bucket_size = 10;
public:
  KdTree() = default;
  KdTree(const std::vector<Vec3d>& points, const std::vector<unsigned int>& ids)
  {
    insert(points, ids);
    build();
  }

  void insert(const std::vector<Vec3d>& points, const std::vector<unsigned int>& ids);
  void build();
  inline bool empty() { return pts.empty(); }
  const KdVec3d& search_nearest_point(const Vec3d& query);
  NodePtr create_leaf_node(PointContainer& c);
  NodePtr new_internal_node();
  void create_internal_node(NodePtr n, PointContainer& c);
  void split(Separator& sep, PointContainer& c_origin, PointContainer& c_low);
  void handle_extended_node(InternalPtr nh, PointContainer& c, PointContainer& c_low);
  class OrthogonalNearestSeach
  {
  public:
    KdVec3d m_dists;
    KdVec3d m_query;
    KdTree* m_tree;

    KdVec3dIter m_result;
    double m_square_distance;
    OrthogonalNearestSeach(KdTree* tree, const KdVec3d& query);
    inline KdVec3dIter closest_point_iter() { return m_result; }

    void compute_nearest_neightbors_orthogonally(NodePtr N, double rd);
    void search_nearest_in_leaf(LeafPtr node);
    double new_square_distance(double dist, double old_off, double new_off);
  };
};
}// namespace Geometry
}// namespace Cage