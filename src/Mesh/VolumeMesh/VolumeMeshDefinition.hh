#pragma once

// Simple mesh
#include "TetMesh.hh"

namespace Cage
{
namespace VolumeMesh
{
struct VVertexProp
{
  bool is_constraint;
  bool is_on_boundary;
#ifdef CHECK_MANIFOLD
  bool is_manifold;
#endif
  VVertexProp();
};

struct VEdgeProp
{
  bool is_constraint;
  bool is_on_boundary;
#ifdef CHECK_MANIFOLD
  bool is_manifold;
#endif
  VEdgeProp();
};

struct VFaceProp
{
  bool is_constraint;
  bool is_boundary;
  VFaceProp();
};

struct VCellProp
{
  bool is_inside;
  bool is_adjacent_constraint;
  VCellProp();
};

typedef TetMeshKernel<VVertexProp, VEdgeProp, VFaceProp, VCellProp> VMeshKernel;

class VMeshT : public TetMesh<VMeshKernel>
{
public:
  /***** only used in tetrahedral mesh *****/

  /*****    query opposite elements    *****/
  VertexHandle vertexOppositeEdge(FaceHandle _fh, EdgeHandle _eh);
  VertexHandle vertexOppositeFace(CellHandle _ch, HalfFaceHandle _fh);

  HalfEdgeHandle edgeOppositeVertex(HalfFaceHandle _fh, VertexHandle _vh);
  EdgeHandle edgeOppositeTwoFace(CellHandle _ch, HalfFaceHandle _f0, HalfFaceHandle _f1);

  FaceHandle faceOppositeVertex(CellHandle _ch, VertexHandle _vh);
  HalfFaceHandle halffaceOppositeVertex(CellHandle _ch, VertexHandle _vh);
};
#ifdef OUTPUT_VMESH
void find_cut_cells(VMeshT& mesh, const BoundingBox& bbox, VMeshT& cut_mesh, std::vector<VertexHandle>& vertices, std::vector<FaceHandle>& faces);
void write_cut_cells_to_file(VMeshT& mesh, std::vector<VertexHandle>& vertices, std::vector<FaceHandle>& faces, std::string file_name);
void write_vmesh_boundary(VMeshT& mesh, std::string file_name);
#endif
}// namespace VolumeMesh
}// namespace Cage