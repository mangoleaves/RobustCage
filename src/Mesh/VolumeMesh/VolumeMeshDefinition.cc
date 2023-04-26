#include "VolumeMeshDefinition.hh"
#include "CageInitializer/TetRemover.hh"

namespace Cage
{
namespace VolumeMesh
{
VVertexProp::VVertexProp()
{
  is_constraint = false;
  is_on_boundary = false;
#ifdef CHECK_MANIFOLD
  is_manifold = true;
#endif
}

VEdgeProp::VEdgeProp()
{
  is_constraint = false;
  is_on_boundary = false;
#ifdef CHECK_MANIFOLD
  is_manifold = true;
#endif
}

VFaceProp::VFaceProp()
{
  is_constraint = false;
  is_boundary = false;
}

VCellProp::VCellProp()
{
  is_inside = false;
  is_adjacent_constraint = false;
}

VertexHandle VMeshT::vertexOppositeEdge(FaceHandle _fh, EdgeHandle _eh)
{
  ET& e = edge(_eh);
  VertexHandle ev0 = e.vertices[0];
  VertexHandle ev1 = e.vertices[1];
  std::array<VertexHandle, 3> fvs = findFV(_fh);
  for (VertexHandle fv : fvs)
  {
    if (fv != ev0 && fv != ev1)
      return fv;
  }
  return VertexHandle::invalidHandle();
}

/// @brief find the vertex opposite to the triangle in a tetrahedron.
VertexHandle VMeshT::vertexOppositeFace(CellHandle _ch, HalfFaceHandle _fh)
{
  CT& c = cell(_ch);
  // find another face
  HalfFaceHandle adj_fh = *std::find_if(c.halffaces.begin(), c.halffaces.end(),
    [&](HalfFaceHandle& fh) {return fh != _fh;});
  // find vertex belongs to adj_fh but not to _fh.
  std::array<VertexHandle, 3> fh_vertices = findFV(faceHdl(_fh));
  for (HalfEdgeHandle& heh : face(adj_fh).halfedges)
  {
    VertexHandle adj_vh = fromVertex(heh);
    if (std::find_if(fh_vertices.begin(), fh_vertices.end(),
      [&](VertexHandle& vh) {return vh == adj_vh;}) == fh_vertices.end())
    {
      return adj_vh;  // it's a non-adjacent vertex.
    }
  }
  return VertexHandle::invalidHandle();
}

/// @brief find the halfedge opposite to the vertex in a halfface.
HalfEdgeHandle VMeshT::edgeOppositeVertex(HalfFaceHandle _fh, VertexHandle _vh)
{
  SubHalfEdges fhehs = findHFHE(_fh);
  for (HalfEdgeHandle fheh : fhehs)
  {
    if (!isAdjacent(edgeHdl(fheh), _vh))
      return fheh;
  }
  return HalfEdgeHandle::invalidHandle();
}

/// @brief find the edge opposite to two triangles in a tetrahedron.
EdgeHandle VMeshT::edgeOppositeTwoFace(CellHandle _ch, HalfFaceHandle _f0, HalfFaceHandle _f1)
{
  CT& c = cell(_ch);
  // find the other two faces
  auto hf_iter = std::find_if(c.halffaces.begin(), c.halffaces.end(),
    [&](HalfFaceHandle& fh) {return fh != _f0 && fh != _f1;});
  HalfFaceHandle f2 = *(hf_iter++);
  HalfFaceHandle f3 = *std::find_if(hf_iter, c.halffaces.end(),
    [&](HalfFaceHandle& fh) {return fh != _f0 && fh != _f1;});
  // find common edge between f2 and f3
  return commonEdge(faceHdl(f2), faceHdl(f3));
}

/// @brief find the face opposite to a vertex in a tetrahedron.
FaceHandle VMeshT::faceOppositeVertex(CellHandle _ch, VertexHandle _vh)
{
  return faceHdl(halffaceOppositeVertex(_ch, _vh));
}

/// @brief find the halfface opposite to a vertex in a tetrahedron.
HalfFaceHandle VMeshT::halffaceOppositeVertex(CellHandle _ch, VertexHandle _vh)
{
  CT& c = cell(_ch);
  for (auto hfh : c.halffaces)
  {
    if (!isAdjacent(faceHdl(hfh), _vh))
      return hfh;
  }
  return HalfFaceHandle::invalidHandle();
}

#ifdef OUTPUT_VMESH
void find_cut_cells(VMeshT& mesh, const BoundingBox& bbox, VMeshT& cut_mesh, std::vector<VertexHandle>& cut_vertices, std::vector<FaceHandle>& cut_faces)
{
  // 1. find cut vertices
  for (size_t vidx = 0;vidx < mesh.vertices.size();vidx++)
  {
    VertexHandle vh(vidx);
    auto& p = mesh.point(vh);
    if (bbox.do_intersect(p))
    {
      cut_vertices.push_back(vh);
    }
  }
  // 2. find cut cells
  std::vector<CellHandle> uncut_cells;
  for (size_t cidx = 0;cidx < mesh.cells.size();cidx++)
  {
    CellHandle ch(cidx);
    auto cv = mesh.findCV(ch);
    size_t cut_vrt_cnt = 0;
    cut_vrt_cnt += std::binary_search(cut_vertices.begin(), cut_vertices.end(), cv[0]);
    cut_vrt_cnt += std::binary_search(cut_vertices.begin(), cut_vertices.end(), cv[1]);
    cut_vrt_cnt += std::binary_search(cut_vertices.begin(), cut_vertices.end(), cv[2]);
    cut_vrt_cnt += std::binary_search(cut_vertices.begin(), cut_vertices.end(), cv[3]);
    if (cut_vrt_cnt < 1)
    {
      uncut_cells.push_back(ch);
    }
  }
  // 3. remove uncut cells
  cut_mesh = mesh;
  CageInit::TetRemover tet_remover(&cut_mesh);
  tet_remover.removeTet(uncut_cells);
  // 4. retrieve boundary faces
  cut_vertices.clear();
  for (size_t vidx = 0;vidx < cut_mesh.vertices.size();vidx++)
  {
    if (cut_mesh.vertices[vidx].prop.is_on_boundary)
      cut_vertices.push_back(VertexHandle(vidx));
  }
  for (size_t fidx = 0;fidx < cut_mesh.faces.size();fidx++)
  {
    if (cut_mesh.faces[fidx].prop.is_boundary)
      cut_faces.push_back(FaceHandle(fidx));
  }
}

void write_cut_cells_to_file(VMeshT& mesh, std::vector<VertexHandle>& vertices, std::vector<FaceHandle>& faces, std::string file_name)
{
  std::fstream fout;
  fout.open(file_name, std::fstream::out | std::fstream::binary);
  ASSERT(fout.is_open(), "fail to open file");

  // map vertices to sequence
  std::unordered_map<int, int> v2v;
  for (VertexHandle vh : vertices)
  {
    v2v[vh.idx()] = (int)v2v.size();
  }

  // save vertices to file
  int vn = (int)vertices.size();
  fout.write((char*)(&vn), sizeof(int));
  for (VertexHandle vh : vertices)
  {
    auto& p = mesh.point(vh);
    fout.write((char*)(&p.x()), sizeof(double));
    fout.write((char*)(&p.y()), sizeof(double));
    fout.write((char*)(&p.z()), sizeof(double));
  }

  // save faces to file
  int fn = (int)faces.size();
  fout.write((char*)(&fn), sizeof(int));
  for (FaceHandle fh : faces)
  {
    auto fv = mesh.findFV(fh);
    fout.write((char*)&(v2v[fv[0].m_idx]), sizeof(int));
    fout.write((char*)&(v2v[fv[1].m_idx]), sizeof(int));
    fout.write((char*)&(v2v[fv[2].m_idx]), sizeof(int));
  }

  fout.close();
}

void write_vmesh_boundary(VMeshT& mesh, std::string file_name)
{
  std::fstream fout;
  fout.open(file_name, std::fstream::out | std::fstream::binary);
  ASSERT(fout.is_open(), "fail to open file");

  // write vertices
  //int vn = (int)mesh.vertices.size();
  int vn = 0;
  for (size_t vidx = 0;vidx < mesh.vertices.size();vidx++)
    vn += mesh.vertices[vidx].prop.is_on_boundary && !mesh.vertices[vidx].prop.is_constraint;
  fout.write((char*)(&vn), sizeof(int));

  std::vector<int> v_map2_b; v_map2_b.resize(mesh.vertices.size());
  int output_vn = 0;

  for (size_t vidx = 0;vidx < mesh.vertices.size();vidx++)
  {
    if (mesh.vertices[vidx].prop.is_on_boundary && !mesh.vertices[vidx].prop.is_constraint)
    {
      v_map2_b[vidx] = output_vn++;
      VertexHandle vh(vidx);
      auto& p = mesh.point(vh);
      fout.write((char*)(&p.x()), sizeof(double));
      fout.write((char*)(&p.y()), sizeof(double));
      fout.write((char*)(&p.z()), sizeof(double));
    }
  }

  // write faces
  //int fn = (int)mesh.faces.size();
  int fn = 0;
  for (size_t fidx = 0;fidx < mesh.faces.size();fidx++)
    fn += mesh.faces[fidx].prop.is_boundary && !mesh.faces[fidx].prop.is_constraint;
  fout.write((char*)(&fn), sizeof(int));

  for (size_t fidx = 0;fidx < mesh.faces.size();fidx++)
  {
    if (mesh.faces[fidx].prop.is_boundary && !mesh.faces[fidx].prop.is_constraint)
    {
      FaceHandle fh(fidx);
      auto fv = mesh.findFV(fh);
      //fout.write((char*)&(fv[0].m_idx), sizeof(int));
      //fout.write((char*)&(fv[1].m_idx), sizeof(int));
      //fout.write((char*)&(fv[2].m_idx), sizeof(int));
      fout.write((char*)&(v_map2_b[fv[0].m_idx]), sizeof(int));
      fout.write((char*)&(v_map2_b[fv[1].m_idx]), sizeof(int));
      fout.write((char*)&(v_map2_b[fv[2].m_idx]), sizeof(int));
    }
  }

  // write tets
  /*
  int cn = (int)mesh.cells.size();
  fout.write((char*)(&cn), sizeof(int));
  for (size_t cidx = 0;cidx < mesh.cells.size();cidx++)
  {
    CellHandle ch(cidx);
    auto cf = mesh.cells[cidx].faces();
    fout.write((char*)&(cf[0].m_idx), sizeof(int));
    fout.write((char*)&(cf[1].m_idx), sizeof(int));
    fout.write((char*)&(cf[2].m_idx), sizeof(int));
    fout.write((char*)&(cf[3].m_idx), sizeof(int));
  }
  */
  fout.close();

#ifdef CHECK_MANIFOLD
  // write_singular_vertices
  fout.open(file_name + "sv.txt", std::fstream::out);
  ASSERT(fout.is_open(), "fail to open file");

  // write vertices
  for (size_t vidx = 0;vidx < mesh.vertices.size();vidx++)
  {
    VertexHandle vh(vidx);
    if (!mesh.vertex(vh).prop.is_manifold && mesh.vertex(vh).prop.is_on_boundary && !mesh.vertex(vh).prop.is_constraint)
    {
      fout << v_map2_b[vidx] << std::endl;
    }
  }

  fout.close();

  fout.open(file_name + "se.txt", std::fstream::out);
  ASSERT(fout.is_open(), "fail to open file");
  // write edges
  for (size_t eidx = 0;eidx < mesh.edges.size();eidx++)
  {
    EdgeHandle eh(eidx);
    auto& edge = mesh.edge(eh);
    if (!edge.prop.is_manifold && edge.prop.is_on_boundary && !edge.prop.is_constraint)
    {
      fout << v_map2_b[edge.from().idx()] << " " << v_map2_b[edge.to().idx()] << std::endl;
    }
  }
#endif
}
#endif
}// namespace VolumeMesh
}// namespace Cage