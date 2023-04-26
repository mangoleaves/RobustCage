#include "CageSimplifier.hh"

namespace Cage
{
namespace CageSimp
{
CageSimplifier::CageSimplifier(SMeshT* original, SMeshT* cage, ParamCageSimplifier* p)
  :om(original), rm(cage), param(p)
{}

void CageSimplifier::simplify()
{
  Logger::user_logger->info("initializing simplifier.");

  pre_calculate_edge_length(om);
  pre_calculate_edge_length(rm);
  pre_calculate_face_area(om);
  pre_calculate_face_area(rm);
  init_one_ring_faces(rm);

  vt = std::make_unique<VertexTree>(*om);
  ot = std::make_unique<DFaceTree>(*om);
  og = std::make_unique<FaceGrid>(*om);
  lrt = std::make_unique<LightDFaceTree>(*rm);

  if (!rm->has_face_normals())
    rm->request_face_normals();
  if (!rm->has_vertex_normals())
    rm->request_vertex_normals();
  rm->update_normals();

  degeneration_remover = std::make_unique<DegenerationRemover>(
    om, rm, vt.get(), ot.get(), lrt.get(), og.get());
  fast_simplifier = std::make_unique<FastSimplifier>(
    om, rm, vt.get(), ot.get(), lrt.get(), og.get(), &param->paramFastSimplifier);
  calc_diagonal_length();
  if (fast_simplifier->param->targetVerticesNum == 0)
    fast_simplifier->param->targetVerticesNum = (size_t)(om->n_vertices() * 3);

  Logger::user_logger->info("begin simplifying.");

  degeneration_remover->perform();
  fast_simplifier->simplify();

  #ifdef OUTPUT_MIDDLE_RESULT
    OpenMesh::IO::write_mesh(*rm, param->fileOutPath + std::to_string(param->cageLabel) + "_fast_simplify.obj",
      OpenMesh::IO::Options::Default, 15);
  #endif
  degeneration_remover = nullptr;
  fast_simplifier = nullptr;
  vt = nullptr;

  // initialize hausdorff distance.
  rt = std::make_unique<FaceTree>(*rm);
  generate_out_links(om, rm, rt.get());
  rt = nullptr;
  calc_face_in_error(rm, std::vector<FaceHandle>(rm->faces_begin(), rm->faces_end()));
  #ifdef USE_TREE_SEARCH
  calc_out_error(rm, om, ot.get(), cage_infinite_fp);
  #else
  calc_out_error(rm, om, og.get(), cage_infinite_fp);
  #endif

  collapse_stage = std::make_unique<CollapseStage>(
    om, rm, &param->paramCollapse,
    ot.get(), lrt.get(), og.get(), original_diagonal_length);

  relocate_stage = std::make_unique<RelocateStage>(
    om, rm, &param->paramRelocate,
    vt.get(), ot.get(), lrt.get(), og.get(), original_diagonal_length);

  flip_stage = std::make_unique<FlipStage>(
    om, rm, &param->paramFlip,
    ot.get(), lrt.get(), og.get(), original_diagonal_length);

  simplify_to_target_num();
  Logger::user_logger->info("simplification done.");
}

void CageSimplifier::calc_diagonal_length()
{
  Vec3d ptMin(DBL_MAX, DBL_MAX, DBL_MAX);
  Vec3d ptMax(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  for (const auto& vh : om->vertices())
  {
    const auto& point = om->point(vh);
    ptMin.minimize(point);
    ptMax.maximize(point);
  }
  original_diagonal_length = (ptMax - ptMin).norm();
  degeneration_remover->original_diagonal_length = original_diagonal_length;
  fast_simplifier->original_diagonal_length = original_diagonal_length;
}

void CageSimplifier::update_strategy()
{
  if (!allow_negtive)
  {
    collapse_stage->update(6, allow_negtive, max_distance_error);
    relocate_stage->update(6, allow_negtive, max_distance_error);
  }
  else
  {
    collapse_stage->update(10, allow_negtive, max_distance_error);
    relocate_stage->update(10, allow_negtive, max_distance_error);
  }
  flip_stage->update(allow_negtive, max_distance_error);
}

void CageSimplifier::simplify_to_target_num()
{
  if (rm->n_vertices() <= param->targetVerticesNum)
    return;
  Logger::user_logger->info("collapsing to target vertices number {}.", param->targetVerticesNum);

  // goal
  const size_t edge_num_to_collapse = rm->n_vertices() - param->targetVerticesNum;

  // strategy: force collapsing procedure to jump out and do smoothing.
  std::queue<size_t> force_skip_en;
  size_t force_skip_times = 3;
  for (size_t i = 0;i < force_skip_times - 1;i++)
    force_skip_en.push(edge_num_to_collapse / force_skip_times);
  force_skip_en.push(edge_num_to_collapse - (edge_num_to_collapse / force_skip_times) * (force_skip_times - 1));

  // initialize parameters during simplification
  allow_negtive = false;
  max_distance_error = original_diagonal_length * param->initError;

  // record statistics during simplification
  size_t total_cen = 0;               // cen: Collapsed Edges Number
  size_t last_iter_cen = 0;

  size_t error_relax_iter = 0;
  size_t iter = 0;
  size_t no_collapsed_iter = 0;
  while (true)  // when collapsed_edge_num == edge_num_to_collapse, end loop immediately.
  {
    update_strategy();
    collapse_stage->do_collapse(force_skip_en.front(), total_cen);
    flip_stage->do_flip();
    relocate_stage->do_relocate();
    size_t collapsed_this_iter = total_cen - last_iter_cen;
    last_iter_cen = total_cen;
    // forced to jump out,
    // continue if still have edges to collapse |OR| break if reach the target vertices number.
    if (force_skip_en.front() == total_cen)
    {
      // don't relax distance error, do next iteration or stop.
      last_iter_cen = 0;
      total_cen = 0;
      force_skip_en.pop();
      if (force_skip_en.empty())
        break;
      else
        continue;
    }
    iter++;
    // if reach max iter, break.
    if (iter == param->maxIter)
      break;
    // little edges are collapsed, perhaps need relax valence constraint.
    if (collapsed_this_iter <= 10 && error_relax_iter >= param->maxErrorRelaxIter)
      param->paramCollapse.maxValence += 1;

    // no edge is collapsed.
    if (collapsed_this_iter == 0)
    {
      no_collapsed_iter++;
      if (no_collapsed_iter == 5)
        break;
    }
    else no_collapsed_iter = 0;   // reset

    // no edge is able to collapse, perhaps need relax distance error.
    if (iter % param->relaxErrorIterStep == 0)
    {
      if (!allow_negtive)
        allow_negtive = true;
      else if (error_relax_iter < param->maxErrorRelaxIter)
        max_distance_error += original_diagonal_length * param->errorStep;

      error_relax_iter++;
      Logger::user_logger->info("[{}]current maximal distance error {}", error_relax_iter, max_distance_error);
      Logger::user_logger->info("{} vertices and {} faces remained.", rm->n_vertices(), rm->n_faces());

    #ifdef OUTPUT_MIDDLE_RESULT
      OpenMesh::IO::write_mesh(*rm, param->fileOutPath + std::to_string(param->cageLabel) + "_simplify_iter_" + std::to_string(iter) + ".obj",
        OpenMesh::IO::Options::Default, 15);
    #endif
    }
  }
  rm->garbage_collection();
  lrt->collect_garbage();
}

}// namespace CageSimp
}// namespace Cage