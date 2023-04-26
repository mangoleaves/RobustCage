#pragma once
#include <queue>
// VolumeMesher
#include "BSP.h"
// Simple mesh
#include "Mesh/SurfaceMesh/SurfaceMeshDefinition.h"
#include "Mesh/VolumeMesh/VolumeMeshDefinition.hh"

#include "LatticePointsGenerator.hh"

namespace Cage
{
namespace CageInit
{
using namespace SimpleUtils;
namespace VM = VolumeMesh;
namespace SM = SurfaceMesh;

class Tetrahedralizer
{
public:
  /***** used for decomposition *****/
  enum DecomposeType : uint8_t
  {
    NoDecompose = 0,
    FromVertex = 1,
    FromBarycenter = 2
  };
public:
  // input
  SM::SMeshT* SMesh;
  ParamTetrahedralizer* param;

  // sub components
  std::unique_ptr<LatticePointsGenerator> latticePointsGenerator;

  // middle results
  std::vector<Vec3d> latticePoints;
  std::unique_ptr<BSPcomplex> bspComplex;

  // results 
  VM::VMeshT* VMesh;
public:
  Tetrahedralizer(SM::SMeshT* _SMesh, ParamTetrahedralizer* _param, VM::VMeshT* _VMesh);

  void tetrahedralize();
private:
  void generateLatticePoints();

  BSPcomplex* bspDivide();

  /// @name decompose
  /// @{
  void decomposeBSP();
  void findDecomposeType(
    std::vector<DecomposeType>& decomposeType,
    std::vector<uint32_t>& decomposeVrt,
    std::map<size_t, std::vector<uint32_t>>& barycenterVrts
  );
  void convertBSP2TetMesh(
    const std::vector<DecomposeType>& decomposeType,
    const std::vector<uint32_t>& decomposeVrt,
    const std::map<size_t, std::vector<uint32_t>>& barycenterVrts
  );

  // data assist orientation determination
  std::vector<uint8_t> faceEntered;
  struct FaceToCheck
  {
    VM::FaceHandle fh;
    VM::HalfEdgeHandle enterEdge;
    FaceToCheck() = default;
    FaceToCheck(VM::FaceHandle f, VM::HalfEdgeHandle ee) :
      fh(f), enterEdge(ee)
    {}
  };
  std::queue<FaceToCheck> facesToCheck;
  std::vector<std::array<VM::FaceHandle, 3>> adjFaces;

  std::vector<uint8_t> cellEntered;
  struct CellToAdd
  {
    size_t cellIdx;
    VM::HalfFaceHandle enterFace;
    CellToAdd() = default;
    CellToAdd(size_t c, VM::HalfFaceHandle ef) :
      cellIdx(c), enterFace(ef)
    {}
  };
  std::queue<CellToAdd> cellsToAdd;

  void addCellDirectly(const CellToAdd& cellToAdd);
  std::vector<VM::HalfFaceHandle>
    findCellHalffaces(const CellToAdd& cellToAdd);
  void decomposeCellFromVertex(
    const std::vector<uint32_t>& decomposeVrt,
    const CellToAdd& cellToAdd
  );
  void decomposeCellFromBarycenter(
    const std::vector<uint32_t>& decomposeVrt,
    const std::map<size_t, std::vector<uint32_t>>& barycenterVrts,
    const CellToAdd& cellToAdd
  );
  void addSubCells(
    const CellToAdd& cellToAdd,
    const std::map<VM::VertexHandle, VM::EdgeHandle>& vrts2Edge,
    const std::map<VM::EdgeHandle, VM::FaceHandle>& edges2Face,
    const std::vector<VM::HalfFaceHandle>& halffaces
  );
  void walkToAdjacentCells(
    const CellToAdd& cellToAdd,
    const std::vector<VM::HalfFaceHandle>& halffaces
  );

  // add elements with extra process
  VM::FaceHandle adjustEdgeOriAddFace(VM::SubEdges& ehs);
  VM::CellHandle adjustFaceOriAddCell(VM::SubFaces& fhs);
  /// @}
};
}// namespace CageInit
}// namespace Cage