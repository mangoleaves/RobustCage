#include "VertexRounder.hh"

namespace Cage
{
namespace CageInit
{

bool VertexRounder::doRounding(VertexHandle vh)
{
  // constrained point, skip.
  if (VMesh->vertex(vh).prop.is_constraint)
    return false;
  // exact point and float point are same, skip.
  if (VMesh->exact_point(vh).isSame())
    return false;
  // get rounded point
  Point_3 round_ep = VMesh->exact_point(vh).round();
  // for each connected cell, check volume after rounding
  bool all_positive = true;
  for (CellHandle ch : VMesh->getConnCells(vh))
  {
    if (!checkPositiveVolume(ch, vh, round_ep))
    {
      all_positive = false;
      break;
    }
  }
  if (all_positive)
  {
    // do real rounding
    VMesh->exact_point(vh).rounded(VMesh->point(vh));
    return true;
  }
  else return false;
}

size_t VertexRounder::doRounding(const std::vector<VertexHandle>& vrts)
{
  Logger::user_logger->info("rounding vertex.");
  size_t rounded_vertex = 0;
  for (VertexHandle vh : vrts)
  {
    rounded_vertex += doRounding(vh);
  }
  Logger::user_logger->info("rounded {} vertices.", rounded_vertex);
  return rounded_vertex;
}

size_t VertexRounder::doRounding()
{
  Logger::user_logger->info("rounding vertex.");
  size_t rounded_vertex = 0;
  for (size_t vidx = 0;vidx < VMesh->nVertices();vidx++)
  {
    rounded_vertex += doRounding(VertexHandle(vidx));
  }
  Logger::user_logger->info("rounded {} vertices.", rounded_vertex);
  return rounded_vertex;
}

bool VertexRounder::checkPositiveVolume(CellHandle ch, VertexHandle vh, const Point_3& approx_p)
{
  HalfFaceHandle opp_hfh = VMesh->halffaceOppositeVertex(ch, vh);
  std::array<VertexHandle, 3> fvs = VMesh->findHFV(opp_hfh);
  CGAL::Sign ori = CGAL::orientation(
    VMesh->exact_point(fvs[0]).exact(), VMesh->exact_point(fvs[1]).exact(),
    VMesh->exact_point(fvs[2]).exact(), approx_p);
  return ori == CGAL::ON_POSITIVE_SIDE;
}

bool VertexRounder::checkAllTets()
{
  // first we count rounded points.
  size_t rounded_vertices = 0;
  size_t unconstrained_vertices = 0;
  for (size_t vidx = 0;vidx < VMesh->vertices.size();vidx++)
  {
    VertexHandle vh(vidx);
    if (!VMesh->vertex(vh).prop.is_constraint)
    {
      if (VMesh->exact_point(vh).isSame())
        rounded_vertices++;
      unconstrained_vertices++;
    }
  }
  Logger::user_logger->info("rounded_vertices ratio = {}", (double)rounded_vertices / (double)unconstrained_vertices);
  // second we check tets' volume.
  for (size_t cidx = 0;cidx < VMesh->cells.size();cidx++)
  {
    CellHandle ch(cidx);
    for (HalfFaceHandle cf : VMesh->cell(ch).halffaces)
    {
      std::array<VertexHandle, 3> fvs = VMesh->findHFV(cf);
      VertexHandle opp_vh = VMesh->vertexOppositeFace(ch, cf);
      CGAL::Sign ori = CGAL::orientation(
        VMesh->exact_point(fvs[0]).exact(), VMesh->exact_point(fvs[1]).exact(),
        VMesh->exact_point(fvs[2]).exact(), VMesh->exact_point(opp_vh).exact());
      if (ori != CGAL::ON_POSITIVE_SIDE)
      {
        Logger::user_logger->error("detected negtive volume tet.");
        return false;
      }
    }
  }
  return true;
}

}// namespace CageInit
}// namespace Cage