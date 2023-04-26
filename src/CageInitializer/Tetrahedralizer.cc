#include <map>
#include <unordered_map>
#include <unordered_set>
#include "Tetrahedralizer.hh"
#include "Geometry/GeometryMath.h"

namespace Cage
{
namespace CageInit
{

/// @param [in] _SMesh surface triangle mesh.
/// @param [in] _param parameters control behaviors of tetrahedralizer.
/// @param [out] _VMesh volume tetrahedral mesh.
Tetrahedralizer::Tetrahedralizer(SM::SMeshT* _SMesh, ParamTetrahedralizer* _param, VM::VMeshT* _VMesh)
{
  SMesh = _SMesh;
  param = _param;
  VMesh = _VMesh;
}

/// @brief generate tet mesh from tri mesh depend on "VolumeMesher".
void Tetrahedralizer::tetrahedralize()
{
  REP_FUNC;
  {Logger::user_logger->info("tetrahedralizing...");}
  // step 1: generate lattice points arond triangle surface mesh.
  //   because we want more regular tets.
  generateLatticePoints();
  // step 2: mainly call makePolyhedralMesh of "VolumeMesher"
  //   for getting a bsp complex.
  bspComplex = std::unique_ptr<BSPcomplex>(bspDivide());
  // step 3: generate tet mesh from bsp complex.
  decomposeBSP();
  {Logger::user_logger->info("tetrahedralizing done!");}
}

/// @brief generate lattice points to form more regular tets around surface mesh.
void Tetrahedralizer::generateLatticePoints()
{
  REP_FUNC;
  // Generate lattice points
  latticePointsGenerator = std::make_unique<LatticePointsGenerator>
    (SMesh, &param->paramLatticePointsGenerator);

  latticePoints = latticePointsGenerator->generate();

  // delete pointer
  latticePointsGenerator.reset();

  Logger::user_logger->info("total number of points after adding lattice points is {}",
    SMesh->n_vertices() + latticePoints.size());
}

/// @brief given triangle mesh, generate a BSP dividing result.
/// @return bspComplex: result of BSP dividing.
BSPcomplex* Tetrahedralizer::bspDivide()
{
  REP_FUNC;

  // construct input for makePolyhedralMesh.
  uint32_t nPnts = (uint32_t)(SMesh->n_vertices() + latticePoints.size());
  double* coords = (double*)malloc(nPnts * 3 * sizeof(double));
  uint32_t nTri = (uint32_t)SMesh->n_faces();
  uint32_t* triIdx = (uint32_t*)malloc(nTri * 3 * sizeof(uint32_t));

  size_t v_offset = 0;
  for (SM::VertexHandle vh : SMesh->vertices())
  {
    memcpy(coords + v_offset, SMesh->point(vh).data(), 3 * sizeof(double));
    v_offset += 3;
  }
  for (Vec3d& p : latticePoints)
  {
    memcpy(coords + v_offset, p.data(), 3 * sizeof(double));
    v_offset += 3;
  }

  for (SM::FaceHandle fh : SMesh->faces())
  {
    auto vh_iter = SMesh->fv_begin(fh);
    triIdx[fh.idx() * 3] = vh_iter->idx(); vh_iter++;
    triIdx[fh.idx() * 3 + 1] = vh_iter->idx(); vh_iter++;
    triIdx[fh.idx() * 3 + 2] = vh_iter->idx(); vh_iter++;
  }

  double* coords_B = nullptr;
  uint32_t nCoords_B = 0;
  uint32_t* triIdx_B = nullptr;
  uint32_t nTriIdx_B = 0;

  {Logger::user_logger->info("making polyhedral mesh...");}

  bool verbose = Logger::user_logger->should_log(spdlog::level::info);

  auto* bspComplex = makePolyhedralMesh(
    "input triangle mesh", coords, nPnts, triIdx, nTri,
    /*fileB_name=*/nullptr, coords_B, nCoords_B, triIdx_B, nTriIdx_B, /*bool_opcode=*/'0',    // unused parameters
    /*free_mem=*/true, /*verbose=*/verbose, /*logging=*/false
  );
  {Logger::user_logger->info("total number of points after making polyhedral is {}", bspComplex->vertices.size());}
  {Logger::user_logger->info("making polyhedral done!");}
  return bspComplex;
}

/// @brief decompose bsp complex to tet mesh.
/// @param [in] bspComplex bsp complex fron "VolumeMesher"
/// @result tet mesh will be saved in member "polyMesh"
void Tetrahedralizer::decomposeBSP()
{
  REP_FUNC;
  // step 1. triangualte all bsp faces.
  {Logger::user_logger->info("triangulating BSP faces...");}
  bspComplex->triangulateAllBspFaces();
  {Logger::user_logger->info("triangulating done!");}

  // step 2. decompose all bsp cells to tets.
  {Logger::user_logger->info("decomposing BSP cells...");}
  std::vector<DecomposeType> decomposeType;
  std::vector<uint32_t> decomposeVrt;
  std::map<size_t, std::vector<uint32_t>> barycenterVrts;
  findDecomposeType(decomposeType, decomposeVrt, barycenterVrts);
  convertBSP2TetMesh(decomposeType, decomposeVrt, barycenterVrts);
  {Logger::user_logger->info("decomposing BSP cells done!");}
}

/// @brief find each cell's decomposition type:
///   1. no decompose
///   2. from one vertex of cell
///   3. from barycenter of cell
///  @sa BSPcomplex::makeTetrahedra()
void Tetrahedralizer::findDecomposeType(
  std::vector<DecomposeType>& decomposeType,
  std::vector<uint32_t>& decomposeVrt,
  std::map<size_t, std::vector<uint32_t>>& barycenterVrts)
{
  REP_FUNC;
  auto& cells = bspComplex->cells;
  decomposeType.resize(cells.size(), NoDecompose);
  decomposeVrt.resize(cells.size(), UINT32_MAX);
  size_t nNoDecompose = 0, nFromVrt = 0, nFromBarycenter = 0; // only used for logging

  for (size_t cellIdx = 0;cellIdx < cells.size();cellIdx++)
  {
    auto& cell = cells[cellIdx];

    if (cell.faces.size() > 4)
    {
      uint64_t numCellEdges = bspComplex->count_cellEdges(cell);
      uint64_t numCellVrts = bspComplex->count_cellVertices(cell, &numCellEdges);
      std::vector<uint32_t> cellVrts(numCellVrts, UINT32_MAX);
      bspComplex->list_cellVertices(cell, numCellEdges, cellVrts);

      // check whether cell can be decomposed from a vertex.
      bool needBarycenter = true;
      for (uint32_t vrtIdx : cellVrts)
      {
        if (bspComplex->cell_is_tetrahedrizable_from_v(cell, vrtIdx))
        {
          decomposeType[cellIdx] = FromVertex;
          decomposeVrt[cellIdx] = vrtIdx;
          needBarycenter = false;
          nFromVrt += 1;
          break;
        }
      }

      if (needBarycenter)
      {
        decomposeType[cellIdx] = FromBarycenter;
        decomposeVrt[cellIdx] = (uint32_t)(bspComplex->vertices.size() - 1);
        barycenterVrts[cellIdx] = cellVrts;
        nFromBarycenter += 1;
      }
    }
    else nNoDecompose += 1;
  }
  {// log
    Logger::user_logger->info(
      "decompose {} cells from vertex, {} cells from barycenter.\n"\
      "the remaining {} cells are not decomposed.",
      nFromVrt, nFromBarycenter, nNoDecompose);
  }
}

/// @brief convert bsp complex to tet mesh, do decomposition in this procedure.
/// @param [in] bspComplex bsp complex
/// @result tet, saved as class member
/// @attention bsp complex will be deleted!
void Tetrahedralizer::convertBSP2TetMesh(
  const std::vector<DecomposeType>& decomposeType,
  const std::vector<uint32_t>& decomposeVrt,
  const std::map<size_t, std::vector<uint32_t>>& barycenterVrts)
{
  REP_FUNC;

  // pre calculate number of new edges, faces and cells
  size_t nCellToDecompose = bspComplex->cells.size();
  size_t nNewCell = 0, nNewFace = 0, nNewEdge = 0;
  for (size_t cellIdx = 0;cellIdx < nCellToDecompose;cellIdx++)
  {
    if (decomposeType[cellIdx] != DecomposeType::NoDecompose)
    {
      // count number of cell's face, edge and vertex
      size_t nLocalFace = bspComplex->cells[cellIdx].faces.size();
      size_t nLocalEdge = 0;
      for (size_t& faceIdx : bspComplex->cells[cellIdx].faces)
        nLocalEdge += bspComplex->faces[faceIdx].edges.size();
      nLocalEdge /= 2;
      size_t nLocalVertex = nLocalEdge - nLocalFace + 2;
      // accumulate
      nNewCell += nLocalFace;
      nNewFace += nLocalEdge;
      nNewEdge += nLocalVertex;
    }
  }
  // pre reserve memory
  VMesh->vertices.reserve(bspComplex->vertices.size());
  VMesh->edges.reserve(bspComplex->edges.size() + nNewEdge);
  VMesh->faces.reserve(bspComplex->faces.size() + nNewFace);
  VMesh->cells.reserve(bspComplex->cells.size() + nNewCell);

  for (auto* vertex : bspComplex->vertices)
  {
    ExactPoint point;
    point.fromImplicitPoint(vertex);
    VMesh->addVertex(point);
  }

  for (auto& edge : bspComplex->edges)
  {
    VMesh->addEdge(edge.vertices[0], edge.vertices[1]);
  }

  for (auto& face : bspComplex->faces)
  {
    VM::SubEdges ehs; size_t i = 0;
    // collect all halfedges
    for (auto& eIdx : face.edges)
      ehs[i++] = VM::EdgeHandle(eIdx);
    // add
    auto fh = adjustEdgeOriAddFace(ehs);
    VMesh->face(fh).prop.is_constraint = face.colour == BLACK_A;  // see BSPcomplex::fill_face_colour
  }
  faceEntered.clear(); faceEntered.resize(bspComplex->faces.size(), false);
  adjFaces.clear(); adjFaces.resize(bspComplex->faces.size());

  // intialize procudure of adding cell
  cellEntered.clear(); cellEntered.resize(nCellToDecompose, false);

  // choose one tet, determine it's faces orientation.
  size_t initCellIdx = 0;
  for (size_t i = 0;i < bspComplex->cells.size();i++)
  {
    if (decomposeType[i] == DecomposeType::NoDecompose)
    {
      initCellIdx = i;break;
    }
  }
  auto& initFaces = bspComplex->cells[initCellIdx].faces;
  VM::SubFaces fhs = {
    VM::FaceHandle(initFaces[0]), VM::FaceHandle(initFaces[1]),
    VM::FaceHandle(initFaces[2]), VM::FaceHandle(initFaces[3]) };
  auto ch = adjustFaceOriAddCell(fhs);
  if (ch.isValid())
    VMesh->cell(ch).prop.is_inside = bspComplex->cells[initCellIdx].place == INTERNAL_A; // see BSPcomplex::makeTetrahedra
  cellEntered[initCellIdx] = true;

  // Add four incident cells of init cell to queue
  for (VM::HalfFaceHandle hfh : VMesh->cell(ch).halffaces)
  {
    VM::FaceHandle fh = faceHdl(hfh);
    BSPface& bspFace = bspComplex->faces[fh.idx()];
    for (size_t i = 0;i < 2;i++)
    {
      if (bspFace.conn_cells[i] != UINT64_MAX && !cellEntered[bspFace.conn_cells[i]])
      {
        cellsToAdd.emplace(bspFace.conn_cells[i], reverse(hfh));
        cellEntered[bspFace.conn_cells[i]] = true;
      }
    }
  }
  // wide first search incident cells, and add them.
  size_t processed_n = 1;
  while (!cellsToAdd.empty())
  {
    CellToAdd cellToAdd = cellsToAdd.front();
    cellsToAdd.pop();
    processed_n++;
    switch (decomposeType[cellToAdd.cellIdx])
    {
    case DecomposeType::FromVertex:
      decomposeCellFromVertex(decomposeVrt, cellToAdd);
      break;
    case DecomposeType::FromBarycenter:
      decomposeCellFromBarycenter(decomposeVrt, barycenterVrts, cellToAdd);
      break;
    case DecomposeType::NoDecompose:
      addCellDirectly(cellToAdd);
      break;
    }
  }
  ASSERT(processed_n == bspComplex->cells.size(), "remain cells not processed.");
}

void Tetrahedralizer::addCellDirectly(const CellToAdd& cellToAdd)
{
  size_t cellIdx = cellToAdd.cellIdx;
  VM::HalfFaceHandle enterFace = cellToAdd.enterFace;

  auto& faces = bspComplex->cells[cellIdx].faces;
  VM::SubFaces fhs = {
    VM::FaceHandle(faces[0]), VM::FaceHandle(faces[1]),
    VM::FaceHandle(faces[2]), VM::FaceHandle(faces[3])
  };
  VM::SubHalfFaces hfhs = {
    halfFaceHdl(fhs[0], VM::Sequence), halfFaceHdl(fhs[1], VM::Sequence),
    halfFaceHdl(fhs[2], VM::Sequence), halfFaceHdl(fhs[3], VM::Sequence)
  };
  // adjust enter face's orientation
  auto enterFaceIt = std::find(fhs.begin(), fhs.end(), faceHdl(enterFace));
  auto enterHalfFaceIt = hfhs.begin() + (enterFaceIt - fhs.begin());
  *enterHalfFaceIt = enterFace;
  // adjust remain faces' orientation
  VM::SubHalfEdges fhes = VMesh->findHFHE(enterFace);
  for (VM::HalfEdgeHandle& he : fhes)
    he = reverse(he);
  for (auto halfFaceIt = hfhs.begin();halfFaceIt != hfhs.end();halfFaceIt++)
  {
    if (halfFaceIt == enterHalfFaceIt)
      continue;

    VM::SubHalfEdges adj_fhes = VMesh->findHFHE(*halfFaceIt);
    if (std::find_if(adj_fhes.begin(), adj_fhes.end(),
      [&](VM::HalfEdgeHandle he) { return std::find(fhes.begin(), fhes.end(), he) != fhes.end(); })
      == adj_fhes.end())
    {
      *halfFaceIt = reverse(*halfFaceIt);
    }
  }

  // add cell
  VM::CellHandle newCell = VMesh->addCell(hfhs);
  if (newCell.isValid())
    VMesh->cell(newCell).prop.is_inside = bspComplex->cells[cellIdx].place == INTERNAL_A; // see BSPcomplex::makeTetrahedra

  // Add adjacent cells to queue
  std::vector<VM::HalfFaceHandle> halffaces(
    VMesh->cell(newCell).halffaces.begin(),
    VMesh->cell(newCell).halffaces.end());
  walkToAdjacentCells(cellToAdd, halffaces);
}

/// @brief determine cell's faces' orientation by enter face and enter orientation
std::vector<VM::HalfFaceHandle>
Tetrahedralizer::findCellHalffaces(const CellToAdd& cellToAdd)
{
  size_t cellIdx = cellToAdd.cellIdx;
  VM::HalfFaceHandle enterFace = cellToAdd.enterFace;

  BSPcell& cell = bspComplex->cells[cellIdx];

  std::vector<VM::HalfFaceHandle> cellHalffaces;  cellHalffaces.reserve(cell.faces.size());
  // initialize
  for (size_t faceIdx : cell.faces) faceEntered[faceIdx] = false;

  // find adjacent relation
  for (size_t i = 0;i < cell.faces.size() - 1;i++)
  {
    for (size_t j = i + 1;j < cell.faces.size();j++)
    {
      auto& face_i = VMesh->face(VM::FaceHandle(cell.faces[i]));
      auto& face_j = VMesh->face(VM::FaceHandle(cell.faces[j]));
      for (size_t ei = 0;ei < 3;ei++)
      {
        for (size_t ej = 0;ej < 3;ej++)
        {
          if (edgeHdl(face_i.halfedges[ei]) == edgeHdl(face_j.halfedges[ej]))
          {
            adjFaces[cell.faces[i]][ei] = cell.faces[j];
            adjFaces[cell.faces[j]][ej] = cell.faces[i];
            goto FOUND_COMMON_EDGE;
          }
        }
      }
    FOUND_COMMON_EDGE:;
    }
  }
  // determine enter face's orientation
  cellHalffaces.push_back(enterFace);
  faceEntered[faceHdl(enterFace).idx()] = true;
  // add adjacent faces to queue
  facesToCheck = std::queue<FaceToCheck>();
  VM::HalfOrientation enterFaceOri = orientation(enterFace);
  for (size_t ei = 0;ei < 3;ei++)
  {
    VM::FaceHandle adjFace = adjFaces[faceHdl(enterFace).idx()][ei];
    VM::HalfEdgeHandle enterHalfedge;
    if (enterFaceOri == VM::Sequence)
      enterHalfedge = reverse(VMesh->face(enterFace).halfedges[ei]);
    else
      enterHalfedge = VMesh->face(enterFace).halfedges[ei];
    facesToCheck.emplace(adjFace, enterHalfedge);
    faceEntered[adjFace.idx()] = true;
  }

  // wide first search to determine remain faces' orientation
  while (!facesToCheck.empty())
  {
    FaceToCheck faceToCheck = facesToCheck.front();
    facesToCheck.pop();

    // determine orientation
    VM::FaceHandle fh = faceToCheck.fh;
    VM::HalfEdgeHandle enterEdge = faceToCheck.enterEdge;
    VM::HalfFaceHandle hfh = halfFaceHdl(fh, VM::Reverse);

    for (VM::HalfEdgeHandle he : VMesh->face(fh).halfedges)
    {
      if (he == enterEdge) // Sequence
      {
        hfh = halfFaceHdl(fh, VM::Sequence);
        break;
      }
    }
    cellHalffaces.push_back(hfh);
    // add adjacent faces to queue
    VM::HalfOrientation ori = orientation(hfh);
    for (size_t ei = 0;ei < 3;ei++)
    {
      VM::FaceHandle adjFace = adjFaces[fh.idx()][ei];
      if (faceEntered[adjFace.idx()])
        continue;
      VM::HalfEdgeHandle enterHalfedge;
      if (ori == VM::Sequence)
        enterHalfedge = reverse(VMesh->face(fh).halfedges[ei]);
      else
        enterHalfedge = VMesh->face(fh).halfedges[ei];
      facesToCheck.emplace(adjFace, enterHalfedge);
      faceEntered[adjFace.idx()] = true;
    }
  }

  return cellHalffaces;
}

void Tetrahedralizer::decomposeCellFromVertex(
  const std::vector<uint32_t>& decomposeVrt,
  const CellToAdd& cellToAdd
)
{
  size_t cellIdx = cellToAdd.cellIdx;
  VM::HalfFaceHandle enterFace = cellToAdd.enterFace;

  BSPcell& cell = bspComplex->cells[cellIdx];
  VM::VertexHandle start_v(decomposeVrt[cellIdx]);

  // second EdgeHandle is edge which connects first VertexHandle and start_v.
  std::map<VM::VertexHandle, VM::EdgeHandle> vrts2Edge;
  // second FaceHandle is face which contains first EdgeHandle and start_v.
  std::map<VM::EdgeHandle, VM::FaceHandle> edges2Face;

  std::vector<VM::HalfFaceHandle> cellHalffaces = findCellHalffaces(cellToAdd);

  // we are going to find all faces that are not adjacent to start_v,
  // and connect them to start_v to form new tets.

  // find faces adjacent or not to start_v
  std::vector<VM::HalfFaceHandle> nonAdjacentFaces;
  for (VM::HalfFaceHandle hfh : cellHalffaces)
  {
    VM::FaceHandle fh = faceHdl(hfh);
    std::array<VM::VertexHandle, 3> face_vrts = VMesh->findFV(fh);
    if (std::find(face_vrts.begin(), face_vrts.end(), start_v) != face_vrts.end())
    {
      // This face is adjacent to start_v. So,
      // a. find an edge of this face that is not adjacent to start_v,
      // and add the pair to edges2Face.
      // b. find two edges of this face that are adjacent to start_v,
      // and add the pairs to vrts2Edge.
      for (VM::EdgeHandle eh : VMesh->face(fh).edges())
      {
        auto& edge = VMesh->edge(eh);
        if (edge.vertices[0] == start_v || edge.vertices[1] == start_v)
        {
          // adjacent edge
          VM::VertexHandle end_vh = edge.vertices[0] == start_v ? edge.vertices[1] : edge.vertices[0];
          vrts2Edge[end_vh.idx()] = eh.idx();
        }
        else
        {
          // non-adjacent edge
          edges2Face[eh.idx()] = fh.idx();
        }
      }
    }
    else
    {
      nonAdjacentFaces.push_back(hfh);
    }
  }

  // for each non-adjacent faces, add edge and face.
  for (VM::HalfFaceHandle hfh : nonAdjacentFaces)
  {
    // all edges of face are not adjacent to start_v.
    // NOTE: vertices of non-adjacent edge may be connected to start_v.
    VM::FaceHandle fh = faceHdl(hfh);
    auto& face = VMesh->face(fh);
    // we are going to connect edges of fh to start_v, to form a new face.
    for (VM::EdgeHandle eh : face.edges())
    {
      if (edges2Face.count(eh.idx()))
        continue;
      auto& edge = VMesh->edge(eh);

      // we are going to connect vertices of eh to start_v.
      // if already connected, find the linking edge.
      for (VM::VertexHandle vh : edge.vertices)
      {
        if (vrts2Edge.count(vh.idx()))
          continue;

        vrts2Edge[vh.idx()] = VMesh->addEdge(vh, start_v);
      }
      // add new face if no face contains eh and start_v.
      if (!edges2Face.count(eh.idx()))
      {
        // connect non-adjacent edges with start_v
        VM::SubEdges ehs = { eh, vrts2Edge[edge.vertices[0]], vrts2Edge[edge.vertices[1]] };
        edges2Face[eh.idx()] = adjustEdgeOriAddFace(ehs);
      }
    }
  }

  addSubCells(cellToAdd, vrts2Edge, edges2Face, nonAdjacentFaces);
  walkToAdjacentCells(cellToAdd, cellHalffaces);
}

void Tetrahedralizer::decomposeCellFromBarycenter(
  const std::vector<uint32_t>& decomposeVrt,
  const std::map<size_t, std::vector<uint32_t>>& barycenterVrts,
  const CellToAdd& cellToAdd)
{
  size_t cellIdx = cellToAdd.cellIdx;
  VM::HalfFaceHandle enterFace = cellToAdd.enterFace;
  BSPcell& cell = bspComplex->cells[cellIdx];
  // barycenter
  const auto& cellVrts = barycenterVrts.at(cellIdx);
  Vector_3 barycenter(0., 0., 0.);
  for (uint32_t vidx : cellVrts)
    barycenter += VMesh->exact_point(VM::VertexHandle(vidx)).exactVec3();
  barycenter /= (double)cellVrts.size();
  VM::VertexHandle center_v = VMesh->addVertex(ExactPoint(barycenter));
  // second EdgeHandle is edge which connects first VertexHandle and start_v.
  std::map<VM::VertexHandle, VM::EdgeHandle> vrts2Edge;
  // second FaceHandle is face which contains first EdgeHandle and start_v.
  std::map<VM::EdgeHandle, VM::FaceHandle> edges2Face;

  std::vector<VM::HalfFaceHandle> cellHalffaces = findCellHalffaces(cellToAdd);

  // prepare vrts2Edge and edges2Face
  for (VM::HalfFaceHandle hfh : cellHalffaces)
  {
    VM::FaceHandle fh = faceHdl(hfh);
    for (VM::EdgeHandle eh : VMesh->face(fh).edges())
    {
      if (edges2Face.count(eh))
        continue;

      // connect vertices of eh to start_v.
      VM::VMeshT::ET& edge = VMesh->edge(eh);
      for (VM::VertexHandle& vh : edge.vertices)
      {
        if (vrts2Edge.count(vh))
          continue;
        vrts2Edge[vh] = VMesh->addEdge(vh, center_v);
      }
      // connect eh to start_v.
      VM::SubEdges ehs = { eh, vrts2Edge[edge.from()], vrts2Edge[edge.to()] };
      edges2Face[eh] = adjustEdgeOriAddFace(ehs);
    } // after adding edges and faces.
  }

  addSubCells(cellToAdd, vrts2Edge, edges2Face, cellHalffaces);
  walkToAdjacentCells(cellToAdd, cellHalffaces);
}

void Tetrahedralizer::addSubCells(
  const CellToAdd& cellToAdd,
  const std::map<VM::VertexHandle, VM::EdgeHandle>& vrts2Edge,
  const std::map<VM::EdgeHandle, VM::FaceHandle>& edges2Face,
  const std::vector<VM::HalfFaceHandle>& halffaces)
{
  size_t cellIdx = cellToAdd.cellIdx;
  bool inside = bspComplex->cells[cellIdx].place == INTERNAL_A;

  for (VM::HalfFaceHandle hfh : halffaces)
  {
    VM::FaceHandle fh = faceHdl(hfh);
    auto& face = VMesh->face(fh);
    VM::SubFaces fhs = { fh,
      edges2Face.at(edgeHdl(face.halfedges[0])),
      edges2Face.at(edgeHdl(face.halfedges[1])),
      edges2Face.at(edgeHdl(face.halfedges[2])) };

    // adjust face orientation and add cell.

    // get halffaces
    VM::SubHalfFaces hfhs =
    { hfh, halfFaceHdl(fhs[1], VM::Sequence),
      halfFaceHdl(fhs[2], VM::Sequence), halfFaceHdl(fhs[3], VM::Sequence), };

    VM::SubHalfEdges reversed_fhes;
    if (orientation(hfh) == VM::Sequence)
      reversed_fhes = { reverse(face.halfedges[0]), reverse(face.halfedges[1]), reverse(face.halfedges[2]) };
    else
      reversed_fhes = { face.halfedges[0], face.halfedges[1], face.halfedges[2] };

    for (size_t i = 1;i < 4;i++)
    {
      VM::SubHalfEdges adj_fhes = VMesh->findHFHE(hfhs[i]);
      if (std::find(adj_fhes.begin(), adj_fhes.end(), reversed_fhes[i - 1]) == adj_fhes.end())
      {
        hfhs[i] = reverse(hfhs[i]);
      }

      ASSERT(!VMesh->face(hfhs[i]).conn_cells[orientation(hfhs[i])].isValid(), "duplicate conn cell.");
    }
    VM::CellHandle newCell = VMesh->addCell(hfhs);
    // set tet's property.
    if (newCell.isValid())
      VMesh->cell(newCell).prop.is_inside = inside;
  }
}

void Tetrahedralizer::walkToAdjacentCells(
  const CellToAdd& cellToAdd,
  const std::vector<VM::HalfFaceHandle>& halffaces)
{
  size_t cellIdx = cellToAdd.cellIdx;

  for (VM::HalfFaceHandle hfh : halffaces)
  {
    VM::FaceHandle fh = faceHdl(hfh);
    BSPface& bspFace = bspComplex->faces[fh.idx()];
    for (size_t i = 0;i < 2;i++)
    {
      if (bspFace.conn_cells[i] != UINT64_MAX
        && bspFace.conn_cells[i] != cellIdx
        && !cellEntered[bspFace.conn_cells[i]])
      {
        cellsToAdd.emplace(bspFace.conn_cells[i], reverse(hfh));
        cellEntered[bspFace.conn_cells[i]] = true;
      }
    }
  }
}

/// @brief adjust edges' orientation to form a closed and manifold loop, then add face
VM::FaceHandle Tetrahedralizer::adjustEdgeOriAddFace(VM::SubEdges& ehs)
{
  // get all halfedges
  VM::SubHalfEdges hehs; size_t i = 0;
  for (VM::EdgeHandle& eh : ehs)
    hehs[i++] = halfEdgeHdl(eh, VM::Sequence);
  VMesh->sortHalfEdges(hehs);
  return VMesh->addFace(hehs, /*sortEdge*/false);
}

/// @brief adjust faces' orientation to form a closed and manifold cell, then add cell 
VM::CellHandle Tetrahedralizer::adjustFaceOriAddCell(VM::SubFaces& fhs)
{
  const std::vector<genericPoint*>& exact_points = bspComplex->vertices;
  // find all vertices
  std::set<VM::VertexHandle> vrts;
  for (VM::FaceHandle& fh : fhs)
  {
    for (VM::HalfEdgeHandle& heh : VMesh->face(fh).halfedges)
    {
      for (VM::VertexHandle& vh : VMesh->edge(heh).vertices)
      {
        vrts.insert(vh);
        if (vrts.size() == 4)
          goto STAGE_1_END;
      }
    }
  }
STAGE_1_END:;

  // get halffaces
  VM::SubHalfFaces hfhs =
  { halfFaceHdl(fhs[0], VM::Sequence), halfFaceHdl(fhs[1], VM::Sequence),
    halfFaceHdl(fhs[2], VM::Sequence), halfFaceHdl(fhs[3], VM::Sequence), };
  // adjust first face's normal to point to center.
  std::array<VM::VertexHandle, 3> fvs = VMesh->findHFV(hfhs[0]);
  // find opposite vertex
  VM::VertexHandle opp_v =
    *std::find_if(vrts.begin(), vrts.end(), [&](VM::VertexHandle v) {
    return std::find(fvs.begin(), fvs.end(), v) == fvs.end(); });

  int exact_ori = genericPoint::orient3D(
    *exact_points[fvs[0].idx()], *exact_points[opp_v.idx()],
    *exact_points[fvs[1].idx()], *exact_points[fvs[2].idx()]
  );
  if (exact_ori == -1)
    hfhs[0] = reverse(hfhs[0]);

  // adjust remain faces' orientation
  VM::SubHalfEdges fhes = VMesh->findHFHE(hfhs[0]);
  for (VM::HalfEdgeHandle& he : fhes)
    he = reverse(he);
  for (size_t i = 1;i < 4;i++)
  {
    VM::SubHalfEdges adj_fhes = VMesh->findHFHE(hfhs[i]);
    if (std::find_if(adj_fhes.begin(), adj_fhes.end(),
      [&](VM::HalfEdgeHandle he) { return std::find(fhes.begin(), fhes.end(), he) != fhes.end(); })
      == adj_fhes.end())
    {
      hfhs[i] = reverse(hfhs[i]);
    }
  }

  // add cell
  return VMesh->addCell(hfhs);
}

}// namespace CageInit
}// namespace Cage