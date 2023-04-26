#pragma once
#include "MeshComponents.hh"

namespace Cage
{
namespace VolumeMesh
{

template<typename VertexProp>
Vertex<VertexProp>::Vertex(const Vec3d& p)
{
  fp = p;
  ep = ExactPoint(fp);
}

template<typename VertexProp>
Vertex<VertexProp>::Vertex(const ExactPoint& p)
{
  ep = p;
  fp = ep.approx();
}

template<typename VertexProp>
Vertex<VertexProp>::Vertex(const Vertex& v)
{
  fp = v.fp;
  ep = v.ep;

  prop = v.prop;
}

template<typename VertexProp>
Vertex<VertexProp>& Vertex<VertexProp>::operator=(const Vertex& v)
{
  fp = v.fp;
  ep = v.ep;

  prop = v.prop;
  return *this;
}

template<typename VertexProp>
Vertex<VertexProp>::Vertex(Vertex&& v)
{
  fp = std::move(v.fp);
  ep = std::move(v.ep);

  prop = std::move(v.prop);
}

template<typename VertexProp>
Vertex<VertexProp>& Vertex<VertexProp>::operator=(Vertex&& v)
{
  fp = std::move(v.fp);
  ep = std::move(v.ep);

  prop = std::move(v.prop);
  return *this;
}

template<typename EdgeProp>
Edge<EdgeProp>::Edge()
{
  vertices[0].setInvalid();
  vertices[1].setInvalid();
}


template<typename EdgeProp>
Edge<EdgeProp>::Edge(const SubVerts& sub_verts)
{
  vertices = sub_verts;
}

template<typename EdgeProp>
Edge<EdgeProp>::Edge(VertexHandle v0, VertexHandle v1)
{
  vertices[0] = v0;
  vertices[1] = v1;
}

template<typename EdgeProp>
Edge<EdgeProp>::Edge(const Edge& e)
{
  vertices = e.vertices;
  prop = e.prop;
}

template<typename EdgeProp>
Edge<EdgeProp>& Edge<EdgeProp>::operator=(const Edge& e)
{
  vertices = e.vertices;
  prop = e.prop;
  return *this;
}

template<typename EdgeProp>
Edge<EdgeProp>::Edge(Edge&& e)
{
  vertices = std::move(e.vertices);
  prop = std::move(e.prop);
}

template<typename EdgeProp>
Edge<EdgeProp>& Edge<EdgeProp>::operator=(Edge&& e)
{
  vertices = std::move(e.vertices);
  prop = std::move(e.prop);
  return *this;
}

template<typename FaceProp>
Face<FaceProp>::Face()
{
  conn_cells[0].setInvalid();
  conn_cells[1].setInvalid();
}

template<typename FaceProp>
Face<FaceProp>::Face(const SubHalfEdges& hehs)
{
  halfedges = hehs;
  conn_cells[0].setInvalid();
  conn_cells[1].setInvalid();
}

template<typename FaceProp>
Face<FaceProp>::Face(const Face& f)
{
  halfedges = f.halfedges;
  conn_cells = f.conn_cells;

  prop = f.prop;
}

template<typename FaceProp>
Face<FaceProp>& Face<FaceProp>::operator=(const Face& f)
{
  halfedges = f.halfedges;
  conn_cells = f.conn_cells;

  prop = f.prop;
  return *this;
}

template<typename FaceProp>
Face<FaceProp>::Face(Face&& f)
{
  halfedges = std::move(f.halfedges);
  conn_cells = std::move(f.conn_cells);

  prop = std::move(f.prop);
}

template<typename FaceProp>
Face<FaceProp>& Face<FaceProp>::operator=(Face&& f)
{
  halfedges = std::move(f.halfedges);
  conn_cells = std::move(f.conn_cells);

  prop = std::move(f.prop);
  return *this;
}

template<typename FaceProp>
SubEdges Face<FaceProp>::edges()
{
  return SubEdges({ edgeHdl(halfedges[0]), edgeHdl(halfedges[1]), edgeHdl(halfedges[2]) });
}

template<typename FaceProp>
SubHalfEdges Face<FaceProp>::reverseHalfEdges()
{
  return SubHalfEdges({ reverse(halfedges[2]), reverse(halfedges[1]), reverse(halfedges[0]) });
}

template<typename FaceProp>
void Face<FaceProp>::addConnCell(CellHandle ch, HalfOrientation ori)
{
  ASSERT(!conn_cells[(size_t)ori].isValid(), "Add duplicate conn cell to face.");
  conn_cells[(size_t)ori] = ch;
}

template<typename FaceProp>
void Face<FaceProp>::removeConnCell(CellHandle ch)
{
  if (conn_cells[0] == ch)conn_cells[0].setInvalid();
  else if (conn_cells[1] == ch)conn_cells[1].setInvalid();
}

template<typename FaceProp>
void Face<FaceProp>::updateConnCell(const std::vector<size_t>& old_mapto_new)
{
  if (conn_cells[0].isValid())
    conn_cells[0] = CellHandle(old_mapto_new[conn_cells[0].idx()]);
  if (conn_cells[1].isValid())
    conn_cells[1] = CellHandle(old_mapto_new[conn_cells[1].idx()]);
}

template<typename FaceProp>
const std::vector<CellHandle> Face<FaceProp>::connCells()const
{
  std::vector<CellHandle> result;
  if (conn_cells[0].isValid()) result.push_back(conn_cells[0]);
  if (conn_cells[1].isValid()) result.push_back(conn_cells[1]);
  return result;
}

template<typename FaceProp>
size_t Face<FaceProp>::nConnCells()const
{
  size_t result = 0;
  if (conn_cells[0].isValid()) result++;
  if (conn_cells[1].isValid()) result++;
  return result;
}

template<typename CellProp>
Cell<CellProp>::Cell()
{}

template<typename CellProp>
Cell<CellProp>::Cell(const SubHalfFaces& hfhs)
{
  halffaces = hfhs;
}

template<typename CellProp>
Cell<CellProp>::Cell(const Cell& c)
{
  halffaces = c.halffaces;
  prop = c.prop;
}

template<typename CellProp>
Cell<CellProp>& Cell<CellProp>::operator=(const Cell& c)
{
  halffaces = c.halffaces;
  prop = c.prop;
  return *this;
}

template<typename CellProp>
Cell<CellProp>::Cell(Cell&& c)
{
  halffaces = std::move(c.halffaces);
  prop = std::move(c.prop);
}

template<typename CellProp>
Cell<CellProp>& Cell<CellProp>::operator=(Cell&& c)
{
  halffaces = std::move(c.halffaces);
  prop = std::move(c.prop);
  return *this;
}

template<typename CellProp>
SubFaces Cell<CellProp>::faces()
{
  return SubFaces({
    faceHdl(halffaces[0]), faceHdl(halffaces[1]),
    faceHdl(halffaces[2]), faceHdl(halffaces[3])
    });
}
}// namespace VolumeMesh
}// namespace Cage