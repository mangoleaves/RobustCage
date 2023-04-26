#include "TetRemover.hh"
#include <map>
#include <set>
#include <unordered_map>

namespace Cage
{
namespace CageInit
{

/// @param [inout] _VMesh volume tetrahedral mesh.
/// @param [in] _param parameters control behaviors.
TetRemover::TetRemover(VM::VMeshT* _VMesh)
{
  VMesh = _VMesh;
}

void TetRemover::removeTet(VM::CellHandle _ch, bool _collectGarbage)
{
  updateBoundaryPropBeforeDeleting(_ch);
  VMesh->quickDeleteCell(_ch);
  if (_collectGarbage)
    VMesh->collectGarbage();
}

void TetRemover::removeTet(const std::vector<VM::CellHandle>& _chs, bool _collectGarbage)
{
  if (_chs.empty())return;

  updateBoundaryPropBeforeDeleting(_chs);
  for (auto& ch : _chs)
    VMesh->quickDeleteCell(ch);
  if (_collectGarbage)
    VMesh->collectGarbage();
}

/// @brief detect whether a tet is removable in topological sense.
bool TetRemover::isTetRemovable(VM::CellHandle _ch)
{
  std::vector<VM::HalfFaceHandle> facesOnBoundary;
  for (auto hfh : VMesh->cell(_ch).halffaces)
  {
    if (VMesh->face(hfh).prop.is_boundary)
      facesOnBoundary.emplace_back(hfh);
  }

  switch (facesOnBoundary.size())
  {
  case 1:
  {
    // if T has exactly one triangle on the surface of the tetrahedralisation, it
    // is safely removable ifand only if the vertex opposite to that triangle is
    // not on the surface.
    auto vhOpposite = VMesh->vertexOppositeFace(_ch, facesOnBoundary[0]);
    return !VMesh->vertex(vhOpposite).prop.is_on_boundary;
  }
  break;
  case 2:
  {
    // if T has exactly two triangles on the surface(T makes an angle in the
    // surface), it is safely removable if and only if the edge opposite to those
    // two triangles does not belong to the surface. Otherwise the removal of
    // the tetrahedron creates a singularity located on that edge.
    auto ehOpposite = VMesh->edgeOppositeTwoFace(_ch, facesOnBoundary[0], facesOnBoundary[1]);
    return !VMesh->edge(ehOpposite).prop.is_on_boundary;
  }
  break;
  case 0: // inside mesh
  case 3: // pyramid lying on the mesh
  case 4: // single floating tetrahedron
    return true;
    break;
  default:
    return false;
    break;
  }
}

/// @brief Update whether faces, edges and vertices of a cell are on boundary
/// after deleting the cell.
/// Update is executed before real deleting. but we can foresee that 
/// number of conn_cells of each face of the cell will decrease by exact one.
void TetRemover::updateBoundaryPropBeforeDeleting(VM::CellHandle _ch)
{
  for (auto& hfh : VMesh->cell(_ch).halffaces)
    updateBoundaryPropIfFaceIsDeletedOnce(faceHdl(hfh));
}

/// @brief Update whether faces, edges and vertices of a cell are on boundary
/// after deleting the cells.
/// Update is executed before real deleting. but we can foresee that 
/// number of conn_cells of each face.
void TetRemover::updateBoundaryPropBeforeDeleting(const std::vector<VM::CellHandle>& _chs)
{
  // 1. find all affected faces.
  // affected faces are faces adjacent to the cells.
  // their boundary property may changed after deleting the cells.
  std::vector<VM::HalfFaceHandle> affectedFaces;
  affectedFaces.reserve(_chs.size() * 4);

  for (auto& ch : _chs)
  {
    auto& c = VMesh->cell(ch);
    affectedFaces.insert(affectedFaces.end(),
      c.halffaces.begin(), c.halffaces.end());
  }

  std::sort(affectedFaces.begin(), affectedFaces.end());
  // 2. for each affected face, determine whether
  //    it will become boundary after deleting.
  size_t i = 0;
  while (i < affectedFaces.size() - 1)
  {
    auto& cur_fh = affectedFaces[i];
    auto& next_fh = affectedFaces[i + 1];
    if (faceHdl(cur_fh) == faceHdl(next_fh))
    {
      // the face has two conn_cells, and two conn_cells are deleted.
      // so this face will be deleted later.
      i += 2;
    }
    else
    {
      // the face is deleted only once.
      updateBoundaryPropIfFaceIsDeletedOnce(faceHdl(cur_fh));
      i += 1;
    }
  }
  // last one face
  if (i < affectedFaces.size())
    updateBoundaryPropIfFaceIsDeletedOnce(faceHdl(affectedFaces[i]));
}

/// @brief update boundary property of face, assuming one connected cell
/// of face is deleted.
void TetRemover::updateBoundaryPropIfFaceIsDeletedOnce(VM::FaceHandle _fh)
{
  auto& face = VMesh->face(_fh);
  switch (face.nConnCells())
  {
  case 2:// the face will become boundary later.
  {
    face.prop.is_boundary = true;
    // all edges and vertices of the face is on boundary.
    for (auto& eh : face.halfedges)
    {
      auto& edge = VMesh->edge(eh);
      edge.prop.is_on_boundary = true;
      VMesh->vertex(edge.from()).prop.is_on_boundary = true;
      VMesh->vertex(edge.to()).prop.is_on_boundary = true;
    }
  }
  break;
  case 1:// the face will be removed later.
  default:// impossible case in manifold mesh.
    break;
  }
}
}// namespace CageInit
}// namespace Cage