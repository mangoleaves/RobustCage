#include "TetMeshTrimmer.hh"
#include "Graph/UndiGraph.hh"
#include "Geometry/GeometryMath.h"
#include <queue>

namespace Cage
{
namespace CageInit
{
using namespace Geometry;

TetMeshTrimmer::TetMeshTrimmer(VM::VMeshT* _VMesh)
{
  VMesh = _VMesh;
}

void TetMeshTrimmer::trim()
{
  tetRemover = std::make_unique<TetRemover>(VMesh);
  tetSubdivider = std::make_unique<TetSubdivider>(VMesh);
  vertexRounder = std::make_unique<VertexRounder>(VMesh);
  std::vector<VM::VertexHandle> vrts_left;

  std::vector<std::array<VM::VertexHandle, 4>> cell_vrts;
  VMesh->updateVertexConnCells(cell_vrts);
  vertexRounder->doRounding();
  removeNonAdjacent2ConstraintTets();

  Logger::user_logger->info("subdividing all tetrahedrons and remove non-adjacent tetrahedrons.(round 1)");
  tetSubdivider->subdivideAll();
  vrts_left = findVrtsLeft();
  vertexRounder->doRounding(vrts_left);
  removeNonAdjacent2ConstraintTets();

  Logger::user_logger->info("subdividing all tetrahedrons and remove non-adjacent tetrahedrons.(round 2)");
  tetSubdivider->subdivideAll();
  tetSubdivider.reset();
  vrts_left = findVrtsLeft();
  vertexRounder->doRounding(vrts_left);
  removeNonAdjacent2ConstraintTets();

  Logger::user_logger->info("final tetrahedral mesh has {} cells, {} faces, {} edges, {} vertices.",
    VMesh->nCells(), VMesh->nFaces(), VMesh->nEdges(), VMesh->nVertices());

#ifdef CHECK_MANIFOLD
  std::vector<VM::EdgeHandle> non_manifold_edges;
  std::vector<VM::VertexHandle> non_manifold_vertices;

  detectAllNonManifold(non_manifold_edges, non_manifold_vertices);
  Logger::user_logger->info("detected {} non-manifold edges, {} non-manifold vertices outside.",
    non_manifold_edges.size(), non_manifold_vertices.size());
  if (!non_manifold_vertices.empty())
    throw std::logic_error("error in manifold.");

  vertexRounder->checkAllTets();
#endif
}

/// @brief remove tets that aren't adjacent to constraint.
/// @return number of deleted tetrahedron
size_t TetMeshTrimmer::removeNonAdjacent2ConstraintTets()
{
  REP_FUNC;
  std::vector<VM::CellHandle> cells_to_delete;
  cells_to_delete.reserve(VMesh->cells.size());
  for (size_t cidx = 0; cidx < VMesh->cells.size();cidx++)
  {
    VM::CellHandle ch(cidx);
    if (!VMesh->deleted(ch) && !VMesh->cell(ch).prop.is_adjacent_constraint)
      cells_to_delete.push_back(ch);
  }
  tetRemover->removeTet(cells_to_delete);
  return cells_to_delete.size();
}

/// @brief tets adjacent to constraint are left.
/// vertices belonging to these tets are left 
/// if they are not adjacent to constraint.
std::vector<VM::VertexHandle> TetMeshTrimmer::findVrtsLeft()
{
  std::vector<VM::VertexHandle> vrts_left;
  for (size_t vidx = 0;vidx < VMesh->vertices.size();vidx++)
  {
    VM::VertexHandle vh(vidx);
    if (VMesh->vertex(vh).prop.is_constraint)
      continue;
    // find tets adjacent to constraint in conn_cells
    for (VM::CellHandle ch : VMesh->getConnCells(vh))
    {
      if (VMesh->cell(ch).prop.is_adjacent_constraint)
      {
        vrts_left.push_back(vh);
        break;
      }
    }
  }
  return vrts_left;
}

#ifdef CHECK_MANIFOLD
/// @brief find all non-manifold vertices and edges, mark them.
void TetMeshTrimmer::detectAllNonManifold(
  std::vector<VM::EdgeHandle>& non_manifold_edges,
  std::vector<VM::VertexHandle>& non_manifold_vertices)
{
  REP_FUNC;
  non_manifold_edges.clear();
  non_manifold_vertices.clear();

  // 1. set all elements' flag to manifold.
  for (auto& vertex : VMesh->vertices)
    vertex.prop.is_manifold = true;
  for (auto& edge : VMesh->edges)
    edge.prop.is_manifold = true;

  // 2. find non-manifold edge and non-manifold face.
  std::vector<size_t> edge_adj_boundary_faces(VMesh->nEdges(), 0);
  std::vector<std::vector<VM::FaceHandle>> vertex_adj_boundary_faces(VMesh->nVertices());

  for (size_t fidx = 0;fidx < VMesh->faces.size();fidx++)
  {
    VM::FaceHandle fh(fidx);
    auto& face = VMesh->face(fh);
    if (face.prop.is_boundary && !face.prop.is_constraint)
    {
      for (VM::EdgeHandle eh : face.edges())
        edge_adj_boundary_faces[eh.idx()]++;
      for (VM::VertexHandle vh : VMesh->findFV(fh))
        vertex_adj_boundary_faces[vh.idx()].push_back(fh);
    }
  }

  for (size_t eidx = 0;eidx < VMesh->edges.size();eidx++)
  {
    VM::EdgeHandle eh(eidx);
    VMesh->edges[eidx].prop.is_manifold = edge_adj_boundary_faces[eidx] <= 2;
    if (!VMesh->edges[eidx].prop.is_manifold)
      non_manifold_edges.push_back(eh);
  }

  // 3. find non-manifold vertex.
  for (size_t vidx = 0;vidx < VMesh->vertices.size();vidx++)
  {
    VM::VertexHandle vh(vidx);
    VMesh->vertices[vidx].prop.is_manifold = isVertexManifold(vh, vertex_adj_boundary_faces[vidx]);
    if (!VMesh->vertices[vidx].prop.is_manifold)
      non_manifold_vertices.push_back(vh);
  }
}

/// @brief estimate whether an edge is manifold in tet mesh.
/// @details @par An embedded tetrahedralisation is manifold if and only if the
/// neighborhood, within the surface, of every surface vertex is connected and if
/// each surface edge is adjacent to exactly two surface triangles.
/// @par If an edge on the mesh surface is adjacent to more than two surface
/// triangles, then the tetrahedralisation is said to have a singularity located on
/// that edge.

//bool TetMeshTrimmer::isEdgeManifold(VM::EdgeHandle _eh);

/// @brief estimate whether a vertex is manifold in tet mesh.
/// @details @par An embedded tetrahedralisation is manifold if and only if the
/// neighborhood, within the surface, of every surface vertex is connected and if
/// each surface edge is adjacent to exactly two surface triangles.
/// @par If the neighborhood of a vertex is not connected (i.e. if its neighborhood
/// has several connected components), then the tetrahedralisation is said to have
/// a topological singularity located on that vertex.
bool TetMeshTrimmer::isVertexManifold(VM::VertexHandle _vh, const std::vector<VM::FaceHandle>& adj_boundary_faces)
{
  if (!VMesh->vertex(_vh).prop.is_on_boundary)
    return true;
  if (adj_boundary_faces.empty())
    return true;
  Graph::UndiGraph local_graph;
  auto& adj_mat = local_graph.adjacentMat;
  adj_mat.resize(adj_boundary_faces.size());

  // set neighbor relation in local graph
  for (size_t i = 0;i < adj_boundary_faces.size() - 1;i++)
  {
    for (size_t j = i + 1;j < adj_boundary_faces.size();j++)
    {
      if (VMesh->isAdjacent(adj_boundary_faces[i], adj_boundary_faces[j]))
      {
        adj_mat[i].push_back((int)j);
        adj_mat[j].push_back((int)i);
      }
    }
  }
  // detect whether local graph is strongly connected by using DFS.
  return local_graph.isConnected();
}
#endif
}// namespace CageInit
}// namespace Cage