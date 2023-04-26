#pragma once

#include "DegenerationRemover.h"
#include "CageFastSimplifier.hh"
#include "SimplifyStages/CollapseStage.hh"
#include "SimplifyStages/RelocateStage.hh"
#include "SimplifyStages/FlipStage.hh"

namespace Cage
{
namespace CageSimp
{

class CageSimplifier
{
public:
  // input
  ParamCageSimplifier* param;
  SMeshT* om;
  SMeshT* rm;

  // sub-components
  std::unique_ptr<DegenerationRemover> degeneration_remover;
  std::unique_ptr<FastSimplifier> fast_simplifier;

  std::unique_ptr<CollapseStage> collapse_stage;
  std::unique_ptr<RelocateStage> relocate_stage;
  std::unique_ptr<FlipStage> flip_stage;
public:
  CageSimplifier(SMeshT* original, SMeshT* cage, ParamCageSimplifier* p);

  void simplify();
private:
  /* Bounding error */

  double original_diagonal_length;
  void calc_diagonal_length();

  /* tree */

  std::unique_ptr<VertexTree> vt;
  std::unique_ptr<DFaceTree> ot;         // static tree for input mesh.
  std::unique_ptr<FaceTree> rt;          // static tree, only used for initialization links.
  std::unique_ptr<LightDFaceTree> lrt;   // light dynamic tree

  /* grid */
  std::unique_ptr<FaceGrid> og;

  /* simplification */

  // allow (after HD - before HD) be negtive
  bool allow_negtive;
  // Hausdorff distance
  double max_distance_error;

  void update_strategy();

  void simplify_to_target_num();
};
}// namespace CageSimp
}// namespace Cage