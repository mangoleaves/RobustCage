#pragma once

#include <queue>

#include "Config.hh"
#include "CageSimplifier/Topo/TopoOperations.h"
#include "CageSimplifier/Topo/EdgeFlipper.h"
#include "CageSimplifier/Geom/GeometryCheck.h"

namespace Cage
{
namespace CageSimp
{
class FlipStage
{
public:
  ParamFlipStage* param;
  SMeshT* om;
  SMeshT* rm;
  DFaceTree* ot;
  LightDFaceTree* lrt;
  FaceGrid* og;

  double original_diagonal_length;
  bool allow_negtive;
  double max_distance_error;
public:
  FlipStage(
    SMeshT* original, SMeshT* cage, ParamFlipStage* p,
    DFaceTree* orginal_tree, LightDFaceTree* remeshing_tree,
    FaceGrid* original_grid,
    double _original_diagonal_length);

  void do_flip();
  void update(bool _allow_negtive, double _max_distance_error);
private:
  std::vector<size_t> update_states;

  struct FlipEdgeReward
  {
    EdgeHandle eh;
    size_t state;
    double reward;

    FlipEdgeReward() = default;
    FlipEdgeReward(EdgeHandle _eh, size_t _state, double _reward) :
      eh(_eh), state(_state), reward(_reward)
    {}

    bool operator<(const FlipEdgeReward& rhs)const { return reward < rhs.reward; }
  };
  typedef std::priority_queue<FlipEdgeReward> FlipEdgeRewardQueue;
  FlipEdgeRewardQueue edges_to_flip;

  void initialize_flip_edges_reward();
  void update_after_flipping(EdgeHandle flipped_edge);

  inline EdgeFlipper new_edge_flipper() { return EdgeFlipper(om, rm, ot, lrt, og); }
};
}// namespace CageSimp
}// namespace Cage