#pragma once
#include <map>
#include <set>
#include <string>
#include "Utils/logger.hh"
#include "MeshComponents.hh"

namespace Cage
{
namespace VolumeMesh
{
using namespace SimpleUtils;
using namespace Geometry;

template<
  typename VertexProp_,
  typename EdgeProp_,
  typename FaceProp_,
  typename CellProp_
>
class TetMeshKernel : public MeshKernel<
  VertexProp_, Vertex,
  EdgeProp_, Edge,
  FaceProp_, Face,
  CellProp_, Cell
>
{};

template<typename Kernel>
class TetMesh
{
public:
  typedef Kernel K;
  typedef typename K::VertexProp VP;
  typedef typename K::EdgeProp EP;
  typedef typename K::FaceProp FP;
  typedef typename K::CellProp CP;

  typedef typename K::Vertex VT;
  typedef typename K::Edge ET;
  typedef typename K::Face FT;
  typedef typename K::Cell CT;
public:
  std::vector<VT> vertices;
  std::vector<ET> edges;
  std::vector<FT> faces;
  std::vector<CT> cells;

  /***** constructor *****/
  TetMesh() = default;
  TetMesh(TetMesh& source);
  TetMesh(TetMesh&& source);
  TetMesh& operator=(TetMesh& source);
  TetMesh& operator=(TetMesh&& source);

  /***** query existance *****/

  inline bool exist(VertexHandle vh) { return vh.isValid() && vh.idx() < vertices.size(); }
  inline bool exist(EdgeHandle eh) { return eh.isValid() && eh.idx() < edges.size(); }
  inline bool exist(HalfEdgeHandle heh) { return heh.isValid() && edgeHdl(heh).idx() < edges.size(); }
  inline bool exist(FaceHandle fh) { return fh.isValid() && fh.idx() < faces.size(); }
  inline bool exist(HalfFaceHandle hfh) { return hfh.isValid() && faceHdl(hfh).idx() < faces.size(); }
  inline bool exist(CellHandle ch) { return ch.isValid() && ch.idx() < cells.size(); }

  inline bool deleted(CellHandle ch) { return is_cell_deleted.at(ch.idx()); }

  /***** query size *****/
  inline size_t nVertices() { return vertices.size(); }
  inline size_t nEdges() { return edges.size(); }
  inline size_t nFaces() { return faces.size(); }
  inline size_t nCells() { return cells.size(); }

  /***** query instance by handle *****/
  // attention: won't check whether handle exists!

  inline ExactPoint& exact_point(VertexHandle vh) { return vertices[vh.idx()].ep; }
  inline Vec3d& point(VertexHandle vh) { return vertices[vh.idx()].fp; }
  inline VT& vertex(VertexHandle vh) { return vertices[vh.idx()]; }
  inline ET& edge(EdgeHandle eh) { return edges[eh.idx()]; }
  inline ET& edge(HalfEdgeHandle heh) { return edges[edgeHdl(heh).idx()]; }
  inline FT& face(FaceHandle fh) { return faces[fh.idx()]; }
  inline FT& face(HalfFaceHandle hfh) { return faces[faceHdl(hfh).idx()]; }
  inline CT& cell(CellHandle ch) { return cells[ch.idx()]; }

  /***** query adjacent *****/

  bool isAdjacent(EdgeHandle eh, VertexHandle vh);
  bool isAdjacent(EdgeHandle e0, EdgeHandle e1);
  bool isAdjacent(HalfEdgeHandle heh, VertexHandle vh);
  bool isAdjacent(HalfEdgeHandle he0, HalfEdgeHandle he1);
  bool isAdjacent(FaceHandle fh, VertexHandle vh);
  bool isAdjacent(FaceHandle fh, EdgeHandle eh);
  bool isAdjacent(FaceHandle f0, FaceHandle f1);

  VertexHandle toVertex(HalfEdgeHandle heh);
  VertexHandle fromVertex(HalfEdgeHandle heh);

  std::array<VertexHandle, 3> findHFV(HalfFaceHandle hfh);
  SubHalfEdges findHFHE(HalfFaceHandle hfh);

  std::array<VertexHandle, 3> findFV(FaceHandle fh);

  std::array<VertexHandle, 4> findCV(CellHandle ch);
  std::array<EdgeHandle, 6> findCE(CellHandle ch);

  /***** query common *****/
  EdgeHandle commonEdge(FaceHandle f0, FaceHandle f1);

  /***** reorder elements *****/
  void sortHalfEdges(SubHalfEdges& hehs);

  /***** add elements *****/

  VertexHandle addVertex(const Vec3d& point);
  VertexHandle addVertex(const ExactPoint& point);

  EdgeHandle addEdge(VertexHandle v0, VertexHandle v1);

  FaceHandle addFace( SubHalfEdges& hehs, bool sortEdges = true);

  CellHandle addCell(SubHalfFaces& hfhs);

  /***** delete elements *****/
  bool quickDeleteCell(CellHandle ch);

  /***** fixed size conn_cell of vertices for less memory *****/
  /***** it won't update with local changes, only be re-build after global changes *****/
  std::vector<uint32_t> conn_cells_begin;
  std::vector<CellHandle> conn_cells;
  bool valid_conn_cell;

  void updateVertexConnCells(std::vector<std::array<VertexHandle, 4>>& cell_vrts);
  void updateFaceConnCells();
  std::vector<CellHandle> getConnCells(VertexHandle vh);

  /***** garbage collection *****/
  void collectGarbage();

  std::vector<uint8_t> is_cell_deleted;

  void clear();
};
}// namespace VolumeMesh
}// namespace Cage

#include "TetMesh_impl.hh"

namespace Cage
{
namespace VolumeMesh
{
typedef TetMeshKernel<VoidProp, VoidProp, VoidProp, VoidProp> TetMeshKernelDT;
typedef TetMesh<TetMeshKernelDT> TetMeshDT;  // default type
}// namespace VolumeMesh
}// namespace Cage