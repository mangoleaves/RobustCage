#pragma once

#include <queue>

#include "Config.hh"
#include "CageSimplifier/Topo/TopoOperations.h"
#include "CageSimplifier/Topo/EdgeCollapser.h"
#include "CageSimplifier/Geom/GeometryCheck.h"

namespace Cage
{
namespace CageSimp
{

class CollapseStage
{
public:
  // input
  ParamCollapseStage* param;
  SMeshT* om;
  SMeshT* rm;
  DFaceTree* ot;
  LightDFaceTree* lrt;
  FaceGrid* og;

  double original_diagonal_length;
  size_t candidate_points_size;
  bool allow_negtive;
  double max_distance_error;
public:
  CollapseStage(
    SMeshT* original, SMeshT* cage, ParamCollapseStage* p,
    DFaceTree* original_tree, LightDFaceTree* remeshing_tree,
    FaceGrid* original_grid,
    double _original_diagonal_length);

  void do_collapse(size_t edge_num_to_collapse, size_t& total_collapsed_edge_num);
  void update(size_t _candidate_points_size, bool _allow_negtive, double _max_distance_error);
private:
  double avg_edge_length;
  std::vector<size_t> update_states;

  struct CollapseEdgeReward
  {
    EdgeHandle eh;
    size_t state;
    double reward;
    Vec3d new_point;

    CollapseEdgeReward() = default;
    CollapseEdgeReward(EdgeHandle _eh, size_t _state, double _reward, const Vec3d& _new_point) :
      eh(_eh), state(_state), reward(_reward), new_point(_new_point)
    {}

    bool operator<(const CollapseEdgeReward& rhs)const { return reward < rhs.reward; }
  };
  typedef std::priority_queue<CollapseEdgeReward> CollapseEdgeRewardQueue;
  CollapseEdgeRewardQueue edges_to_collapse;

  std::vector<Vec3d> generate_candidate_points_for_collapse(EdgeHandle e, EdgeCollapser& edge_collapser);
  bool find_collapse_hausdorff_deviation(
    EdgeHandle eh, double& local_hd_before, double& local_hd_after, Vec3d& new_point);
  void initialize_collapse_edges_reward();
  void update_after_collapsing(VertexHandle collapsed_center);


  inline EdgeCollapser new_edge_collapser() { return EdgeCollapser(om, rm, ot, lrt, og); }
};

}// namespace CageSimp
}// namespace Cage