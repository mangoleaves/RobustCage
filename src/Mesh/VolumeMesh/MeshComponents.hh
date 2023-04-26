#pragma once
#include <array>
#include <vector>
#include "MeshHandle.hh"
#include "Geometry/Exact/ExactPoint.h"

namespace Cage
{
namespace VolumeMesh
{
using namespace Geometry;
using namespace SimpleUtils;

struct VoidProp {};

template<typename VertexProp>
class Vertex
{
public:
  // geometry implementation
  Vec3d fp;     // float point, error tolerant algorithms use this.
  ExactPoint ep;

  VertexProp prop;

  // constructor
  Vertex() = default;
  Vertex(const Vec3d& p);
  Vertex(const ExactPoint& p);
  // copy
  Vertex(const Vertex& v);
  Vertex& operator=(const Vertex& v);
  // move
  Vertex(Vertex&& v);
  Vertex& operator=(Vertex&& v);
};

typedef std::array<VertexHandle, 2> SubVerts;

template<typename EdgeProp>
class Edge
{
public:
  SubVerts vertices;

  EdgeProp prop;

  // constructor
  Edge();
  Edge(const SubVerts& sub_verts);
  Edge(VertexHandle v0, VertexHandle v1);
  // copy
  Edge(const Edge& e);
  Edge& operator=(const Edge& e);
  // move
  Edge(Edge&& e);
  Edge& operator=(Edge&& e);

  inline VertexHandle& from() { return vertices[0]; }
  inline VertexHandle from() const { return vertices[0]; }
  inline VertexHandle& to() { return vertices[1]; }
  inline VertexHandle to() const { return vertices[1]; }
};


typedef std::array<EdgeHandle, 3> SubEdges;
typedef std::array<HalfEdgeHandle, 3> SubHalfEdges;

template<typename FaceProp>
class Face
{
public:
  SubHalfEdges halfedges;
  // Halfface Sequence->[0], Halfface Reverse->[1]
  std::array<CellHandle, 2> conn_cells;

  FaceProp prop;

  Face();
  Face(const SubHalfEdges& hehs);
  // copy
  Face(const Face& f);
  Face& operator=(const Face& f);
  // move
  Face(Face&& f);
  Face& operator=(Face&& f);

  SubEdges edges();
  SubHalfEdges reverseHalfEdges();

  void addConnCell(CellHandle ch, HalfOrientation ori);
  void removeConnCell(CellHandle ch);
  void updateConnCell(const std::vector<size_t>& old_mapto_new);
  const std::vector<CellHandle> connCells()const;
  size_t nConnCells()const;
};

typedef std::array<FaceHandle, 4> SubFaces;
typedef std::array<HalfFaceHandle, 4> SubHalfFaces;

template<typename CellProp>
class Cell
{
public:
  SubHalfFaces halffaces;

  CellProp prop;

  // constructor
  Cell();
  Cell(const SubHalfFaces& hfhs);
  // copy
  Cell(const Cell& c);
  Cell& operator=(const Cell& c);
  // move
  Cell(Cell&& c);
  Cell& operator=(Cell&& c);

  SubFaces faces();
};

template<
  typename VertexProp_,
  template<typename VP_> typename Vertex_,
  typename EdgeProp_,
  template<typename EP_> typename Edge_,
  typename FaceProp_,
  template<typename FP_> typename Face_,
  typename CellProp_,
  template<typename CP_> typename Cell_
>
class MeshKernel
{
public:
  typedef VertexProp_ VertexProp;
  typedef Vertex_<VertexProp> Vertex;
  typedef EdgeProp_ EdgeProp;
  typedef Edge_<EdgeProp> Edge;
  typedef FaceProp_ FaceProp;
  typedef Face_<FaceProp> Face;
  typedef CellProp_ CellProp;
  typedef Cell_<CellProp> Cell;
};
}// namespace VolumeMesh
}// namespace Cage

#include "MeshComponents_impl.hh"