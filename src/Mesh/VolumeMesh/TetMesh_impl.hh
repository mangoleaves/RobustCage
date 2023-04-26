#pragma once
#include "TetMesh_impl.hh"
#include "Geometry/Exact/Predicates.h"
#include <algorithm>

namespace Cage
{
namespace VolumeMesh
{

template<typename Kernel>
TetMesh<Kernel>::TetMesh(TetMesh& source)
{
  vertices = source.vertices;
  edges = source.edges;
  faces = source.faces;
  cells = source.cells;

  is_cell_deleted = source.is_cell_deleted;
  conn_cells = source.conn_cells;
  conn_cells_begin = source.conn_cells_begin;
  valid_conn_cell = source.valid_conn_cell;
}

template<typename Kernel>
TetMesh<Kernel>::TetMesh(TetMesh&& source)
{
  vertices = std::move(source.vertices);
  edges = std::move(source.edges);
  faces = std::move(source.faces);
  cells = std::move(source.cells);

  is_cell_deleted = std::move(source.is_cell_deleted);
  conn_cells = std::move(source.conn_cells);
  conn_cells_begin = std::move(source.conn_cells_begin);
  valid_conn_cell = source.valid_conn_cell;
}

template<typename Kernel>
TetMesh<Kernel>& TetMesh<Kernel>::operator=(TetMesh& source)
{
  vertices = source.vertices;
  edges = source.edges;
  faces = source.faces;
  cells = source.cells;

  is_cell_deleted = source.is_cell_deleted;
  conn_cells = source.conn_cells;
  conn_cells_begin = source.conn_cells_begin;
  valid_conn_cell = source.valid_conn_cell;
  return *this;
}

template<typename Kernel>
TetMesh<Kernel>& TetMesh<Kernel>::operator=(TetMesh&& source)
{
  vertices = std::move(source.vertices);
  edges = std::move(source.edges);
  faces = std::move(source.faces);
  cells = std::move(source.cells);

  is_cell_deleted = std::move(source.is_cell_deleted);
  conn_cells = std::move(source.conn_cells);
  conn_cells_begin = std::move(source.conn_cells_begin);
  valid_conn_cell = source.valid_conn_cell;
  return *this;
}

/// @brief query whether edge and vertex are adjacent.
template<typename Kernel>
bool TetMesh<Kernel>::isAdjacent(EdgeHandle eh, VertexHandle vh)
{
  ASSERT(exist(eh) && exist(vh), "one of elements doesn't exist.");

  const ET& e = edge(eh);
  return e.from() == vh || e.to() == vh;
}

/// @brief query whether two edges are adjacent.
/// if two edges have common vertex, they are adjacent.
template<typename Kernel>
bool TetMesh<Kernel>::isAdjacent(EdgeHandle e0, EdgeHandle e1)
{
  ASSERT(exist(e0) && exist(e1), "one of edges doesn't exist.");
  ASSERT(e0 != e1, "two edges are same.");
  VertexHandle e0v0 = edge(e0).from();
  VertexHandle e0v1 = edge(e0).to();
  VertexHandle e1v0 = edge(e1).from();
  VertexHandle e1v1 = edge(e1).to();

  return e0v0 == e1v0 || e0v0 == e1v1 || e0v1 == e1v0 || e0v1 == e1v1;
}

/// @brief query whether two halfedges are adjacent.
template<typename Kernel>
bool TetMesh<Kernel>::isAdjacent(HalfEdgeHandle heh, VertexHandle vh)
{
  ASSERT(exist(heh) && exist(vh), "one of elements doesn't exist.");
  return toVertex(heh) == vh || fromVertex(heh) == vh;
}

/// @brief query whether two halfedges are adjacent.
/// there is a strong requirement that
/// he0's from-vertex is connected to he1's to-vertex
/// or he0's to-vertex is connected to he1's from-vertex.
template<typename Kernel>
bool TetMesh<Kernel>::isAdjacent(HalfEdgeHandle he0, HalfEdgeHandle he1)
{
  ASSERT(exist(he0) && exist(he1), "one of edges doesn't exist.");
  ASSERT(edgeHdl(he0) != edgeHdl(he1), "two halfedge are same edge.");
  return toVertex(he0) == fromVertex(he1) || fromVertex(he0) == toVertex(he1);
}

/// @brief query whether face and vertex are adjacent.
template<typename Kernel>
bool TetMesh<Kernel>::isAdjacent(FaceHandle fh, VertexHandle vh)
{
  ASSERT(exist(fh) && exist(vh), "one of elements doesn't exist.");

  for (EdgeHandle eh : face(fh).edges())
  {
    ET& e = edge(eh);
    if (e.from() == vh || e.to() == vh)
      return true;
  }
  return false;
}

/// @brief query whether face and edge are adjacent.
template<typename Kernel>
bool TetMesh<Kernel>::isAdjacent(FaceHandle fh, EdgeHandle eh)
{
  ASSERT(exist(fh) && exist(eh), "one of elements doesn't exist.");
  for (HalfEdgeHandle feh : face(fh).halfedges)
  {
    if (edgeHdl(feh) == eh)
      return true;
  }
  return false;
}

/// @brief query whether two faces are adjacent.
template<typename Kernel>
bool TetMesh<Kernel>::isAdjacent(FaceHandle f0, FaceHandle f1)
{
  ASSERT(exist(f0) && exist(f1), "one of elements doesn't exist.");
  SubEdges fe0 = face(f0).edges();
  SubEdges fe1 = face(f1).edges();
  for (EdgeHandle e0 : fe0)
    for (EdgeHandle e1 : fe1)
      if (e0 == e1)
        return true;
  return false;
}

/// @brief query to-vertex by halfedge.
template<typename Kernel>
VertexHandle TetMesh<Kernel>::toVertex(HalfEdgeHandle heh)
{
  ASSERT(exist(heh), "halfedge handle doesn't exist.");
  if (orientation(heh) == Sequence)
    return edges[edgeHdl(heh).idx()].to();
  else  // Reverse
    return edges[edgeHdl(heh).idx()].from();
}

/// @brief query from-vertex by halfedge.
template<typename Kernel>
VertexHandle TetMesh<Kernel>::fromVertex(HalfEdgeHandle heh)
{
  ASSERT(exist(heh), "halfedge handle doesn't exist.");
  if (orientation(heh) == Sequence)
    return edges[edgeHdl(heh).idx()].from();
  else  // Reverse
    return edges[edgeHdl(heh).idx()].to();
}

/// @brief query adjacent halfedges of a halfface.
template<typename Kernel>
std::array<VertexHandle, 3> TetMesh<Kernel>::findHFV(HalfFaceHandle hfh)
{
  ASSERT(exist(hfh), "face doesn't exist.");

  std::array<VertexHandle, 3> vhs;
  size_t idx = 0;
  if (orientation(hfh) == Sequence)
  {
    for (HalfEdgeHandle he : face(hfh).halfedges)
      vhs[idx++] = fromVertex(he);
  }
  else
  {
    for (HalfEdgeHandle he : face(hfh).reverseHalfEdges())
      vhs[idx++] = fromVertex(he);
  }
  return vhs;
}

/// @brief query adjacent halfedges of a halfface.
template<typename Kernel>
SubHalfEdges TetMesh<Kernel>::findHFHE(HalfFaceHandle hfh)
{
  ASSERT(exist(hfh), "face doesn't exist.");

  if (orientation(hfh) == Sequence)
    return face(hfh).halfedges;
  else
    return face(hfh).reverseHalfEdges();
}

/// @brief query adjacent vertices of a face.
template<typename Kernel>
std::array<VertexHandle, 3> TetMesh<Kernel>::findFV(FaceHandle fh)
{
  ASSERT(exist(fh), "face doesn't exist.");
  std::array<VertexHandle, 3> vrts; size_t i = 0;
  for (HalfEdgeHandle heh : face(fh).halfedges)
    vrts[i++] = fromVertex(heh);
  return vrts;
}

/// @brief query adjacent vertices of a cell.
template<typename Kernel>
std::array<VertexHandle, 4> TetMesh<Kernel>::findCV(CellHandle ch)
{
  ASSERT(exist(ch), "cell doesn't exist.");
  std::set<VertexHandle> vrts;

  for (HalfFaceHandle hfh : cell(ch).halffaces)
  {
    for (HalfEdgeHandle heh : face(hfh).halfedges)
    {
      vrts.insert(fromVertex(heh));
    }
  }
  ASSERT(vrts.size() == 4, "cell must have 4 vertices.");
  std::array<VertexHandle, 4> result; size_t i = 0;
  for (VertexHandle vh : vrts) result[i++] = vh;
  return result;
}

/// @brief query adjacent edges of a cell.
template<typename Kernel>
std::array<EdgeHandle, 6> TetMesh<Kernel>::findCE(CellHandle ch)
{
  ASSERT(exist(ch), "cell doesn't exist.");
  std::set<EdgeHandle> es;

  for (HalfFaceHandle hfh : cell(ch).halffaces)
  {
    for (HalfEdgeHandle heh : face(hfh).halfedges)
    {
      es.insert(edgeHdl(heh));
    }
  }
  ASSERT(es.size() == 6, "cell must have 6 edges.");
  std::array<EdgeHandle, 6> result; size_t i = 0;
  for (EdgeHandle eh : es) result[i++] = eh;
  return result;
}

/// @brief find common edge between two faces.
template<typename Kernel>
EdgeHandle TetMesh<Kernel>::commonEdge(FaceHandle f0, FaceHandle f1)
{
  ASSERT(exist(f0) && exist(f1) && f0 != f1, "faces doesn't exist.");
  SubEdges edges0 = face(f0).edges();
  SubEdges edges1 = face(f1).edges();
  std::sort(edges0.begin(), edges0.end());
  std::sort(edges1.begin(), edges1.end());
  std::vector<EdgeHandle> result(std::min(edges0.size(), edges1.size()));
  std::set_intersection(edges0.begin(), edges0.end(), edges1.begin(), edges1.end(), result.begin());
  ASSERT(result.size() <= 1, "find more than one common edges.");
  return result.front();
}

/// @brief sort halfedges. in vector, adjacent halfedges are connected by common vertex.
template<typename Kernel>
void TetMesh<Kernel>::sortHalfEdges(SubHalfEdges& hehs)
{
  // get all from-vertices and to-vertices
  std::vector<VertexHandle> toVs, fromVs;
  toVs.reserve(hehs.size());    fromVs.reserve(hehs.size());
  for (HalfEdgeHandle heh : hehs)
  {
    toVs.push_back(toVertex(heh));
    fromVs.push_back(fromVertex(heh));
  }
  // adjust halfedges' orientation to form a closed and manifold loop
  for (size_t i = 0;i < hehs.size() - 1;i++)
  {
    if (toVs[i] != fromVs[i + 1])
    {
      if (fromVs[i] == fromVs[i + 1])   // reverse halfedge[i]
      {
        hehs[i] = reverse(hehs[i]);
        std::swap(toVs[i], fromVs[i]);
      }
      else if (toVs[i] == toVs[i + 1])// reverse halfedge[i+1]
      {
        hehs[i + 1] = reverse(hehs[i + 1]);
        std::swap(toVs[i + 1], fromVs[i + 1]);
      }
      else  // reverse both
      {
        hehs[i] = reverse(hehs[i]);
        std::swap(toVs[i], fromVs[i]);
        hehs[i + 1] = reverse(hehs[i + 1]);
        std::swap(toVs[i + 1], fromVs[i + 1]);
      }
    }
  }
}

template<typename Kernel>
VertexHandle TetMesh<Kernel>::addVertex(const Vec3d& point)
{
  vertices.emplace_back(point);
  return VertexHandle(vertices.size() - 1);
}

/// @brief add a vertex into mesh
/// @return success or failure.
template<typename Kernel>
VertexHandle TetMesh<Kernel>::addVertex(const ExactPoint& point)
{
  vertices.emplace_back(point);
  return VertexHandle(vertices.size() - 1);
}

/// @brief add an edge into mesh. two vertices of edge must exist.
/// @return success or failure.
template<typename Kernel>
EdgeHandle TetMesh<Kernel>::addEdge(
  VertexHandle v0, VertexHandle v1)
{
  ASSERT(exist(v0) && exist(v1), "Add Edge: one of vertices doesn't exist.");

  edges.emplace_back(v0, v1);
  EdgeHandle eh(edges.size() - 1);
  return eh;
}

/// @brief add a face into mesh. edges of face must exist.
/// @return success or failure.
template<typename Kernel>
FaceHandle TetMesh<Kernel>::addFace(SubHalfEdges& hehs, bool sortEdges)
{
  {// check existance.
    bool all_exist = std::find_if_not(hehs.begin(), hehs.end(),
      [&](HalfEdgeHandle heh) {return exist(heh);}) == hehs.end();
    ASSERT(all_exist, "one of edges doesn't exist.");
  }

  if (sortEdges)
  {
    sortHalfEdges(hehs);
  }

  faces.emplace_back(hehs);
  FaceHandle fh(faces.size() - 1);
  return fh;
}

/// @brief add a cell into mesh. faces of cell must exist.
/// @return success or failure.
template<typename Kernel>
CellHandle TetMesh<Kernel>::addCell(SubHalfFaces& hfhs)
{
  {// check existance.
    bool all_exist = std::find_if_not(hfhs.begin(), hfhs.end(),
      [&](HalfFaceHandle hfh) {return exist(hfh);}) == hfhs.end();
    ASSERT(all_exist, "one of faces doesn't exist.");
  }

  cells.emplace_back(hfhs);
  is_cell_deleted.push_back(false);

  CellHandle ch(cells.size() - 1);
  std::array<VertexHandle, 4> cell_vertices = findCV(ch);
  for (HalfFaceHandle hfh : hfhs)
  {
    face(hfh).addConnCell(ch, orientation(hfh));
  }
  valid_conn_cell = false;
  return ch;
}

/// @brief quickly delete cell.
/// 1. mark cell deleted, then wait garbage collection.
/// 2. remove cell from it's vertices' conn_cells. 
template<typename Kernel>
bool TetMesh<Kernel>::quickDeleteCell(CellHandle ch)
{
  ASSERT(exist(ch), "the cell to be deleted doesn't exist.");
  ASSERT(!is_cell_deleted[ch.idx()], "repeat delete cell.");

  is_cell_deleted[ch.idx()] = true;
  // remove ch from it's faces conn_cells.
  for (HalfFaceHandle hfh : cell(ch).halffaces)
  {
    face(hfh).removeConnCell(ch);
  }
  // remove ch from it's vertices conn_cells.
  valid_conn_cell = false;
  return true;
}

template<typename Kernel>
void TetMesh<Kernel>::updateVertexConnCells(std::vector<std::array<VertexHandle, 4>>& cell_vrts)
{
  // we add conn_cell to faces and vertices
  size_t n_cells = nCells();
  if (cell_vrts.empty())
  {
    cell_vrts.reserve(n_cells);
    for (int cidx = 0;cidx < n_cells;cidx++)
    {
      cell_vrts.push_back(findCV(CellHandle(cidx)));
    }
  }
  // count
  std::vector<uint32_t> vc_cnt(nVertices(), 0);
  for (size_t cidx = 0;cidx < n_cells;cidx++)
  {
    const std::array<VertexHandle, 4>& cell_vertices = cell_vrts[cidx];
    for (VertexHandle vh : cell_vertices)
      vc_cnt[vh.idx()]++;
  }
  // alloc space for conn_cells 
  conn_cells_begin.clear();
  conn_cells_begin.resize(nVertices() + 1, 0);
  uint32_t begin = 0;
  for (size_t vidx = 0;vidx < vc_cnt.size();vidx++)
  {
    conn_cells_begin[vidx] = begin;
    begin += vc_cnt[vidx];
  }
  conn_cells_begin.back() = begin;
  conn_cells.clear();
  conn_cells.resize(conn_cells_begin.back());
  vc_cnt.clear(); vc_cnt.resize(nVertices(), 0);
  // add conn_cell
  for (size_t cidx = 0;cidx < n_cells;cidx++)
  {
    CellHandle ch(cidx);
    const std::array<VertexHandle, 4>& cell_vertices = cell_vrts[ch.idx()];

    for (VertexHandle vh : cell_vertices)
    {
      conn_cells[conn_cells_begin[vh.idx()] + vc_cnt[vh.idx()]] = ch;
      vc_cnt[vh.idx()]++;
    }
  }
  valid_conn_cell = true;
}

template<typename Kernel>
void TetMesh<Kernel>::updateFaceConnCells()
{
  // clear
  for (size_t fidx = 0;fidx < nFaces();fidx++)
  {
    faces[fidx].conn_cells[0].setInvalid();
    faces[fidx].conn_cells[1].setInvalid();
  }
  // rebuild
  for (size_t cidx = 0;cidx < nCells();cidx++)
  {
    CellHandle ch(cidx);
    for (HalfFaceHandle hfh : cell(ch).halffaces)
      face(hfh).addConnCell(ch, orientation(hfh));
  }
}

template<typename Kernel>
std::vector<CellHandle> TetMesh<Kernel>::getConnCells(VertexHandle vh)
{
  ASSERT(valid_conn_cell, "invalid conn cell.");

  std::vector<CellHandle> results;
  size_t cc_begin = conn_cells_begin[vh.idx()];
  size_t cc_end = conn_cells_begin[vh.idx() + 1];
  results.reserve(cc_end - cc_begin);
  for (size_t cc = cc_begin;cc < cc_end;cc++)
    results.push_back(conn_cells[cc]);
  return results;
}

template<typename Kernel>
void TetMesh<Kernel>::collectGarbage()
{
  // vertices, edges and faces reference count.
  const size_t
    n_old_vertices = vertices.size(),
    n_old_edges = edges.size(),
    n_old_faces = faces.size(),
    n_old_cells = cells.size();
  std::vector<size_t>
    v_tag(n_old_vertices, 0),
    e_tag(n_old_edges, 0),
    f_tag(n_old_faces, 0),
    c_tag(n_old_cells, 0);

  // add reference count to sub-elements.
  for (size_t cidx = 0;cidx < cells.size();cidx++)
  {
    if (is_cell_deleted[cidx])
      continue;

    c_tag[cidx] = 1;
    CT& c = cells[cidx];
    for (HalfFaceHandle hfh : c.halffaces)
    {
      FaceHandle fh = faceHdl(hfh);
      f_tag[fh.idx()]++;
      FT& f = faces[fh.idx()];
      for (HalfEdgeHandle heh : f.halfedges)
      {
        EdgeHandle eh = edgeHdl(heh);
        e_tag[eh.idx()]++;
        ET& e = edges[eh.idx()];

        v_tag[e.vertices[0].idx()]++;
        v_tag[e.vertices[1].idx()]++;
      }
    }
  }
  // those who has 0 reference count will be deleted.
  // now, we see x_rag as a map, mapping old index to new index.
#define update_tag(tag) \
  for (size_t idx = 0, new_idx = 0;idx < tag.size();idx++) { \
    if (tag[idx] > 0) {\
      tag[idx] = new_idx;\
      new_idx++;\
    } else { \
      tag[idx] = tag.size(); }\
  }

  update_tag(v_tag);
  update_tag(e_tag);
  update_tag(f_tag);
  update_tag(c_tag);
#undef update_tag

  // really delete element that has 0 reference count .
  size_t n_remain_vertices = 0;
  for (size_t i = 0;i < n_old_vertices;i++)
  {
    if (v_tag[i] < n_old_vertices)
    {
      vertices[v_tag[i]] = std::move(vertices[i]);
      n_remain_vertices++;
    }
  }
  vertices.resize(n_remain_vertices);
  size_t n_remain_edges = 0;
  for (size_t i = 0;i < n_old_edges;i++)
  {
    if (e_tag[i] < n_old_edges)
    {
      edges[e_tag[i]] = std::move(edges[i]);
      n_remain_edges++;
    }
  }
  edges.resize(n_remain_edges);
  size_t n_remain_faces = 0;
  for (size_t i = 0;i < n_old_faces;i++)
  {
    if (f_tag[i] < n_old_faces)
    {
      faces[f_tag[i]] = std::move(faces[i]);
      n_remain_faces++;
    }
  }
  faces.resize(n_remain_faces);
  size_t n_remain_cells = 0;
  for (size_t i = 0;i < n_old_cells;i++)
  {
    if (c_tag[i] < n_old_cells)
    {
      cells[n_remain_cells] = cells[i];
      n_remain_cells++;
    }
  }
  cells.resize(n_remain_cells);

  // correct index in all elements
  for (ET& e : edges)
  {
    e.vertices[0] = VertexHandle(v_tag[e.vertices[0].idx()]);
    e.vertices[1] = VertexHandle(v_tag[e.vertices[1].idx()]);
  }
  for (FT& f : faces)
  {
    for (HalfEdgeHandle& heh : f.halfedges)
    {
      EdgeHandle eh = edgeHdl(heh);
      HalfOrientation ori = orientation(heh);
      eh = EdgeHandle(e_tag[eh.idx()]);
      heh = halfEdgeHdl(eh, ori);
    }
  }
  for (CT& c : cells)
  {
    for (HalfFaceHandle& hfh : c.halffaces)
    {
      FaceHandle fh = faceHdl(hfh);
      HalfOrientation ori = orientation(hfh);
      fh = FaceHandle(f_tag[fh.idx()]);
      hfh = halfFaceHdl(fh, ori);
    }
  }
  updateFaceConnCells();
  is_cell_deleted.resize(cells.size());
  std::fill(is_cell_deleted.begin(), is_cell_deleted.end(), false);
}

template<typename Kernel>
void TetMesh<Kernel>::clear()
{
  vertices.clear(); vertices.shrink_to_fit();
  edges.clear();    edges.shrink_to_fit();
  faces.clear();    faces.shrink_to_fit();
  cells.clear();    cells.shrink_to_fit();
  is_cell_deleted.clear();  is_cell_deleted.shrink_to_fit();
  conn_cells.clear();       conn_cells.shrink_to_fit();
  conn_cells_begin.clear(); conn_cells_begin.shrink_to_fit();
}
}// namespace VolumeMesh
}// namespace Cage