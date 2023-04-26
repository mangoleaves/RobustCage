#include "TetSubdivider.hh"

namespace Cage
{
namespace CageInit
{

TetSubdivider::TetSubdivider(VM::VMeshT* _VMesh)
{
  VMesh = _VMesh;
}

/// @brief Subdivide each tetrahedron by SUB8.
void TetSubdivider::subdivideAll()
{
  REP_FUNC;
  sub_mesh.clear();
  {Logger::user_logger->info("memory used: {} MB", getMegabytesUsed());}
  {Logger::user_logger->info("generating vertices.");}
  generateVertices();
  {Logger::user_logger->info("memory used: {} MB", getMegabytesUsed());}
  {Logger::user_logger->info("splitting all edges.");}
  splitEdges();
  {Logger::user_logger->info("memory used: {} MB", getMegabytesUsed());}
  {Logger::user_logger->info("subdividing all faces.");}
  subdivideFaces();
  {Logger::user_logger->info("memory used: {} MB", getMegabytesUsed());}
  {Logger::user_logger->info("subdividing all tetrahedrons.");}
  subdivideTets_1_12();
  //subdivideTets_1_8();
  //subdivideTets_1_4();
  {Logger::user_logger->info("memory used: {} MB", getMegabytesUsed());}

  cell_vrts.clear();
  *VMesh = std::move(sub_mesh);
}

void TetSubdivider::generateVertices()
{
  size_t vertex_idx_offset = VMesh->nVertices();
  size_t n_edges_to_split = VMesh->nEdges();
  // vertices from original mesh
  // keep order and property
  sub_mesh.vertices.resize(VMesh->nVertices());
#pragma omp parallel for
  for (int vidx = 0;vidx < vertex_idx_offset;vidx++)
  {
    VM::VMeshT::VT& dst = sub_mesh.vertices[vidx];
    VM::VMeshT::VT& src = VMesh->vertices[vidx];
    dst.fp = src.fp;
    dst.ep = src.ep;
    dst.prop = src.prop;
  }
  // vertices from midpoint of all edges
  sub_mesh.vertices.resize(vertex_idx_offset + n_edges_to_split);
# pragma omp parallel for
  for (int eidx = 0;eidx < n_edges_to_split;eidx++)
  {
    VM::EdgeHandle eh(eidx);
    // add mid point
    ExactPoint mid_point((VMesh->exact_point(VMesh->edge(eh).vertices[0]).exactVec3() +
      VMesh->exact_point(VMesh->edge(eh).vertices[1]).exactVec3()) * 0.5);
    VM::VMeshT::VT& mid_v = sub_mesh.vertices[vertex_idx_offset + eidx];
    mid_v = VM::VMeshT::VT(mid_point);
    // update property
    VM::VMeshT::ET& parent_e = VMesh->edge(eh);
    mid_v.prop.is_constraint = parent_e.prop.is_constraint;
    mid_v.prop.is_on_boundary = parent_e.prop.is_on_boundary;
  }
}

/// @brief split marked edges to two sub edges.
/// edge connecting mid point to vertices[0] is stored in sub_children[0].
/// edge connecting mid point to vertices[1] is stored in sub_children[1].
void TetSubdivider::splitEdges()
{
  size_t vertex_idx_offset = VMesh->nVertices();
  size_t n_edges_to_split = VMesh->nEdges();
  //sub_edges.resize(n_edges_to_split);
  sub_mesh.edges.resize(n_edges_to_split * 2);
# pragma omp parallel for
  for (int eidx = 0;eidx < n_edges_to_split;eidx++)
  {
    VM::EdgeHandle eh(eidx);
    VM::VertexHandle mid_vh(vertex_idx_offset + eidx);
    // add two splitted edges
    sub_mesh.edges[eidx * 2] = VM::VMeshT::ET(VMesh->edge(eh).vertices[0], mid_vh);
    sub_mesh.edges[eidx * 2 + 1] = VM::VMeshT::ET(VMesh->edge(eh).vertices[1], mid_vh);

    updateProp(eh, VM::EdgeHandle(eidx * 2), VM::EdgeHandle(eidx * 2 + 1));
  }
}

void TetSubdivider::updateProp(VM::EdgeHandle _parent, VM::EdgeHandle _e0, VM::EdgeHandle _e1)
{
  VM::VMeshT::ET& parent_e = VMesh->edge(_parent);
  VM::VMeshT::ET& child_e0 = sub_mesh.edge(_e0);
  VM::VMeshT::ET& child_e1 = sub_mesh.edge(_e1);

  // set parent child relation
  //child(_parent, 0) = _e0;
  //child(_parent, 1) = _e1;

  // set inherit properties
  child_e0.prop.is_constraint = parent_e.prop.is_constraint;
  child_e0.prop.is_on_boundary = parent_e.prop.is_on_boundary;
  child_e1.prop.is_constraint = parent_e.prop.is_constraint;
  child_e1.prop.is_on_boundary = parent_e.prop.is_on_boundary;
#ifdef CHECK_MANIFOLD
  child_e0.prop.is_manifold = parent_e.prop.is_manifold;
  child_e1.prop.is_manifold = parent_e.prop.is_manifold;
#endif
}

/// @brief subdivide all affected regular faces to sub faces
void TetSubdivider::subdivideFaces()
{
  size_t n_faces_to_subdivide = VMesh->nFaces();
  size_t edge_idx_offset = sub_mesh.nEdges();
  //sub_faces.resize(n_faces_to_subdivide);
  sub_mesh.edges.resize(edge_idx_offset + n_faces_to_subdivide * 3);
  sub_mesh.faces.resize(n_faces_to_subdivide * 4);
#pragma omp parallel for
  for (int fidx = 0;fidx < n_faces_to_subdivide;fidx++)
  {
    VM::FaceHandle fh(fidx);
    // Use notations in fig 2 in ref-paper.
    // Suppose we are subdividing face t0t1t2.
    VM::EdgeHandle t0_t1 = edgeHdl(VMesh->face(fh).halfedges[0]);
    VM::EdgeHandle t1_t2 = edgeHdl(VMesh->face(fh).halfedges[1]);
    VM::EdgeHandle t2_t0 = edgeHdl(VMesh->face(fh).halfedges[2]);

    VM::VertexHandle t0 = VMesh->fromVertex(VMesh->face(fh).halfedges[0]);
    VM::VertexHandle t1 = VMesh->fromVertex(VMesh->face(fh).halfedges[1]);
    VM::VertexHandle t2 = VMesh->fromVertex(VMesh->face(fh).halfedges[2]);

    VM::VertexHandle t01 = splitVertex(t0_t1);
    VM::VertexHandle t12 = splitVertex(t1_t2);
    VM::VertexHandle t02 = splitVertex(t2_t0);

    // add three interior edges
    sub_mesh.edges[edge_idx_offset + fidx * 3] = VM::VMeshT::ET(t01, t02);
    sub_mesh.edges[edge_idx_offset + fidx * 3 + 1] = VM::VMeshT::ET(t12, t01);
    sub_mesh.edges[edge_idx_offset + fidx * 3 + 2] = VM::VMeshT::ET(t02, t12);
    VM::EdgeHandle t01_t02(edge_idx_offset + fidx * 3);
    VM::EdgeHandle t12_t01(edge_idx_offset + fidx * 3 + 1);
    VM::EdgeHandle t02_t12(edge_idx_offset + fidx * 3 + 2);

    std::vector<VM::FaceHandle> child_fhs;  child_fhs.resize(4);
    // face 0: t0 -> t01 -> t02
    VM::SubHalfEdges hehs;
    hehs = {
      halfEdgeHdl(child(t0_t1, t0), VM::Sequence),
      halfEdgeHdl(t01_t02, VM::Sequence),
      halfEdgeHdl(child(t2_t0, t0), VM::Reverse) };
    sub_mesh.faces[fidx * 4] = VM::VMeshT::FT(hehs);
    child_fhs[0] = VM::FaceHandle(fidx * 4);

    // face 1: t1 -> t12 -> t01 
    hehs = {
      halfEdgeHdl(child(t1_t2, t1), VM::Sequence),
      halfEdgeHdl(t12_t01, VM::Sequence),
      halfEdgeHdl(child(t0_t1, t1), VM::Reverse) };
    sub_mesh.faces[fidx * 4 + 1] = VM::VMeshT::FT(hehs);
    child_fhs[1] = VM::FaceHandle(fidx * 4 + 1);

    // face 2: t2 -> t02 -> t12
    hehs = {
      halfEdgeHdl(child(t2_t0, t2), VM::Sequence),
      halfEdgeHdl(t02_t12, VM::Sequence),
      halfEdgeHdl(child(t1_t2, t2), VM::Reverse) };
    sub_mesh.faces[fidx * 4 + 2] = VM::VMeshT::FT(hehs);
    child_fhs[2] = VM::FaceHandle(fidx * 4 + 2);

    // face 3: t01 -> t12 -> t02
    hehs = {
      halfEdgeHdl(t12_t01, VM::Reverse),
      halfEdgeHdl(t02_t12, VM::Reverse),
      halfEdgeHdl(t01_t02, VM::Reverse) };
    sub_mesh.faces[fidx * 4 + 3] = VM::VMeshT::FT(hehs);
    child_fhs[3] = VM::FaceHandle(fidx * 4 + 3);

    updateProp(fh, child_fhs);
    updateProp(fh, t01_t02, t12_t01, t02_t12);
  }
}

/// @brief after subdividing faces, update properties.
void TetSubdivider::updateProp(VM::FaceHandle _parent, const std::vector<VM::FaceHandle>& _children)
{
  // set child and properties
  VM::VMeshT::FT& parent_f = VMesh->face(_parent);
  for (size_t i = 0;i < _children.size();i++)
  {
    VM::VMeshT::FT& child_f = sub_mesh.face(_children[i]);
    //child(_parent, i) = _children[i];
    child_f.prop.is_boundary = parent_f.prop.is_boundary;
    child_f.prop.is_constraint = parent_f.prop.is_constraint;
  }
}

void TetSubdivider::updateProp(VM::FaceHandle _parent, VM::EdgeHandle _e0, VM::EdgeHandle _e1, VM::EdgeHandle _e2)
{
  VM::VMeshT::FT& parent_f = VMesh->face(_parent);
  sub_mesh.edge(_e0).prop.is_constraint = parent_f.prop.is_constraint;
  sub_mesh.edge(_e0).prop.is_on_boundary = parent_f.prop.is_boundary;
  sub_mesh.edge(_e1).prop.is_constraint = parent_f.prop.is_constraint;
  sub_mesh.edge(_e1).prop.is_on_boundary = parent_f.prop.is_boundary;
  sub_mesh.edge(_e2).prop.is_constraint = parent_f.prop.is_constraint;
  sub_mesh.edge(_e2).prop.is_on_boundary = parent_f.prop.is_boundary;
}

#define ADD_EDGE(edge_name, v0, v1) \
  VM::EdgeHandle edge_name(edge_idx_offset + cidx * 6 + local_edge_idx++);\
  sub_mesh.edge(edge_name) = VM::VMeshT::ET(v0, v1);
#define ADD_FACE(face_name) \
  VM::FaceHandle face_name(face_idx_offset + cidx * 16 + local_face_idx++);\
  sub_mesh.face(face_name) = VM::VMeshT::FT(hehs);
#define ADD_CELL() \
  just_add_cell = VM::CellHandle(cidx * 12 + local_cell_idx++);\
  sub_mesh.cell(just_add_cell) = VM::VMeshT::CT(hfhs);\
  children.push_back(just_add_cell);

void TetSubdivider::subdivideTets_1_12()
{
  size_t vertex_idx_offset = sub_mesh.nVertices();
  size_t edge_idx_offset = sub_mesh.nEdges();
  size_t face_idx_offset = sub_mesh.nFaces();
  size_t n_tets_to_subdivide = VMesh->nCells();
  // we will add extra
  // n_tets_to_subdivide vertices,
  // n_tets_to_subdivide * 6 edges,
  // n_tets_to_subdivide * 16 faces
  // and new
  // n_tets_to_subdivide * 12 tets.
  sub_mesh.vertices.resize(vertex_idx_offset + n_tets_to_subdivide);
  sub_mesh.edges.resize(edge_idx_offset + n_tets_to_subdivide * 6);
  sub_mesh.faces.resize(face_idx_offset + n_tets_to_subdivide * 16);
  sub_mesh.cells.resize(n_tets_to_subdivide * 12);
  sub_mesh.is_cell_deleted.resize(n_tets_to_subdivide * 12, false);
  cell_vrts.resize(n_tets_to_subdivide * 12);

#pragma omp parallel for
  for (int cidx = 0;cidx < n_tets_to_subdivide;cidx++)
  {
    VM::CellHandle ch(cidx);

    auto& halffaces = VMesh->cell(ch).halffaces;
    // all notations come from fig.2 in ref-paper
    // in 12-subtetrahedron subdivision, we add a new point into tetrahedron,
    // which is berycenter of tet, called t0123.
    VM::HalfFaceHandle f_t0_t1_t2, f_t0_t1_t3, f_t1_t2_t3, f_t0_t2_t3;

    VM::EdgeHandle e_t0_t1, e_t0_t2, e_t0_t3, e_t1_t2, e_t1_t3, e_t2_t3;

    VM::VertexHandle v_t0, v_t1, v_t2, v_t3;
    VM::VertexHandle v_t01, v_t02, v_t03, v_t12, v_t13, v_t23;
    VM::VertexHandle v_center;

    // we need to fix order of faces, edges and vertices.
    f_t0_t1_t2 = halffaces[0];
    std::array<VM::HalfEdgeHandle, 3> hehs_012 = VMesh->findHFHE(f_t0_t1_t2);
    e_t0_t1 = edgeHdl(hehs_012[0]);
    e_t1_t2 = edgeHdl(hehs_012[1]);
    e_t0_t2 = edgeHdl(hehs_012[2]);
    v_t0 = VMesh->fromVertex(hehs_012[0]);
    v_t1 = VMesh->fromVertex(hehs_012[1]);
    v_t2 = VMesh->fromVertex(hehs_012[2]);

    for (auto hfh_iter = halffaces.begin() + 1;hfh_iter != halffaces.end();hfh_iter++)
    {
      if (VMesh->isAdjacent(faceHdl(*hfh_iter), e_t0_t1))
      {
        f_t0_t1_t3 = *hfh_iter;
        for (VM::HalfEdgeHandle heh : VMesh->face(f_t0_t1_t3).halfedges)
        {
          VM::EdgeHandle eh = edgeHdl(heh);
          if (eh != e_t0_t1)
          {
            if (VMesh->isAdjacent(eh, v_t0))
            {
              e_t0_t3 = eh;
              v_t3 = VMesh->edge(e_t0_t3).from() == v_t0 ?
                VMesh->edge(e_t0_t3).to() : VMesh->edge(e_t0_t3).from();
            }
            else e_t1_t3 = eh;
          }
        }
      }
      else if (VMesh->isAdjacent(faceHdl(*hfh_iter), e_t1_t2))
      {
        f_t1_t2_t3 = *hfh_iter;
        for (VM::HalfEdgeHandle heh : VMesh->face(f_t1_t2_t3).halfedges)
        {
          VM::EdgeHandle eh = edgeHdl(heh);
          if (eh != e_t1_t2)
          {
            if (VMesh->isAdjacent(eh, v_t2))
            {
              e_t2_t3 = eh;
              break;
            }
          }
        }
      }
      else
      {
        f_t0_t2_t3 = *hfh_iter;
      }
    }
    v_t01 = splitVertex(e_t0_t1);
    v_t02 = splitVertex(e_t0_t2);
    v_t03 = splitVertex(e_t0_t3);
    v_t12 = splitVertex(e_t1_t2);
    v_t13 = splitVertex(e_t1_t3);
    v_t23 = splitVertex(e_t2_t3);

    std::vector<VM::CellHandle> children; children.reserve(12);

    size_t local_edge_idx = 0;
    size_t local_face_idx = 0;
    size_t local_cell_idx = 0;
    VM::CellHandle just_add_cell;

    // Add four tets at corner and corresponding four faces.
    VM::HalfOrientation ori_012, ori_013, ori_123, ori_023;
    ori_012 = VM::orientation(f_t0_t1_t2);
    ori_013 = VM::orientation(f_t0_t1_t3);
    ori_123 = VM::orientation(f_t1_t2_t3);
    ori_023 = VM::orientation(f_t0_t2_t3);

    VM::SubHalfEdges hehs;
    VM::SubHalfFaces hfhs;
    // Add face t01_t02_t03 and tet t0_t01_t02_t03
    VM::HalfFaceHandle f_t0_t01_t02 = VM::halfFaceHdl(child(faceHdl(f_t0_t1_t2), v_t0), ori_012);
    VM::HalfFaceHandle f_t0_t01_t03 = halfFaceHdl(child(faceHdl(f_t0_t1_t3), v_t0), ori_013);
    VM::HalfFaceHandle f_t0_t02_t03 = halfFaceHdl(child(faceHdl(f_t0_t2_t3), v_t0), ori_023);
    VM::HalfEdgeHandle e_t01_t02 = sub_mesh.edgeOppositeVertex(f_t0_t01_t02, v_t0);
    VM::HalfEdgeHandle e_t01_t03 = sub_mesh.edgeOppositeVertex(f_t0_t01_t03, v_t0);
    VM::HalfEdgeHandle e_t02_t03 = sub_mesh.edgeOppositeVertex(f_t0_t02_t03, v_t0);
    hehs = { e_t01_t02 , e_t02_t03 , e_t01_t03 };
    ADD_FACE(f_t01_t02_t03);
    hfhs = { f_t0_t01_t02 , f_t0_t01_t03 , f_t0_t02_t03 , halfFaceHdl(f_t01_t02_t03, VM::Reverse) };
    ADD_CELL();
    cell_vrts[just_add_cell.idx()] = { v_t0, v_t01, v_t02, v_t03 };

    // Add face t01_t12_t13 and tet t1_t01_t12_t13
    VM::HalfFaceHandle f_t1_t01_t12 = halfFaceHdl(child(faceHdl(f_t0_t1_t2), v_t1), ori_012);
    VM::HalfFaceHandle f_t1_t01_t13 = halfFaceHdl(child(faceHdl(f_t0_t1_t3), v_t1), ori_013);
    VM::HalfFaceHandle f_t1_t12_t13 = halfFaceHdl(child(faceHdl(f_t1_t2_t3), v_t1), ori_123);
    VM::HalfEdgeHandle e_t01_t12 = sub_mesh.edgeOppositeVertex(f_t1_t01_t12, v_t1);
    VM::HalfEdgeHandle e_t01_t13 = sub_mesh.edgeOppositeVertex(f_t1_t01_t13, v_t1);
    VM::HalfEdgeHandle e_t12_t13 = sub_mesh.edgeOppositeVertex(f_t1_t12_t13, v_t1);
    hehs = { e_t01_t12 , e_t01_t13 , e_t12_t13 };
    ADD_FACE(f_t01_t12_t13);
    hfhs = { f_t1_t01_t12 , f_t1_t01_t13 , f_t1_t12_t13 , halfFaceHdl(f_t01_t12_t13, VM::Reverse) };
    ADD_CELL();
    cell_vrts[just_add_cell.idx()] = { v_t1, v_t01, v_t12, v_t13 };

    // Add face t02_t12_t23 and tet t2_t02_t12_t23 
    VM::HalfFaceHandle f_t2_t02_t12 = halfFaceHdl(child(faceHdl(f_t0_t1_t2), v_t2), ori_012);
    VM::HalfFaceHandle f_t2_t02_t23 = halfFaceHdl(child(faceHdl(f_t0_t2_t3), v_t2), ori_023);
    VM::HalfFaceHandle f_t2_t12_t23 = halfFaceHdl(child(faceHdl(f_t1_t2_t3), v_t2), ori_123);
    VM::HalfEdgeHandle e_t02_t12 = sub_mesh.edgeOppositeVertex(f_t2_t02_t12, v_t2);
    VM::HalfEdgeHandle e_t02_t23 = sub_mesh.edgeOppositeVertex(f_t2_t02_t23, v_t2);
    VM::HalfEdgeHandle e_t12_t23 = sub_mesh.edgeOppositeVertex(f_t2_t12_t23, v_t2);
    hehs = { e_t02_t12 , e_t12_t23 , e_t02_t23 };
    ADD_FACE(f_t02_t12_t23);
    hfhs = { f_t2_t02_t12 , f_t2_t02_t23 , f_t2_t12_t23 , halfFaceHdl(f_t02_t12_t23 , VM::Reverse) };
    ADD_CELL();
    cell_vrts[just_add_cell.idx()] = { v_t2, v_t02, v_t12, v_t23 };

    // Add face t03_t13_t23 and tet t3_t03_t13_t23
    VM::HalfFaceHandle f_t3_t03_t13 = halfFaceHdl(child(faceHdl(f_t0_t1_t3), v_t3), ori_013);
    VM::HalfFaceHandle f_t3_t03_t23 = halfFaceHdl(child(faceHdl(f_t0_t2_t3), v_t3), ori_023);
    VM::HalfFaceHandle f_t3_t13_t23 = halfFaceHdl(child(faceHdl(f_t1_t2_t3), v_t3), ori_123);
    VM::HalfEdgeHandle e_t03_t13 = sub_mesh.edgeOppositeVertex(f_t3_t03_t13, v_t3);
    VM::HalfEdgeHandle e_t03_t23 = sub_mesh.edgeOppositeVertex(f_t3_t03_t23, v_t3);
    VM::HalfEdgeHandle e_t13_t23 = sub_mesh.edgeOppositeVertex(f_t3_t13_t23, v_t3);
    hehs = { e_t03_t13 , e_t03_t23 , e_t13_t23 };
    ADD_FACE(f_t03_t13_t23);
    hfhs = { f_t3_t03_t13 , f_t3_t03_t23 , f_t3_t13_t23 , halfFaceHdl(f_t03_t13_t23 , VM::Reverse) };
    ADD_CELL();
    cell_vrts[just_add_cell.idx()] = { v_t3, v_t03, v_t13, v_t23 };

    // Add barycenter point, add eight center tets and their faces.

    // connect barycenter and edges' center.
    v_center = VM::VertexHandle(vertex_idx_offset + cidx);
    sub_mesh.vertex(v_center) = VM::VMeshT::VT(
      ExactPoint((sub_mesh.exact_point(v_t0).exactVec3() + sub_mesh.exact_point(v_t1).exactVec3()
        + sub_mesh.exact_point(v_t2).exactVec3() + sub_mesh.exact_point(v_t3).exactVec3()) * 0.25));
    sub_mesh.vertex(v_center).prop.is_constraint = false;
    sub_mesh.vertex(v_center).prop.is_on_boundary = false;
    ADD_EDGE(e_center_t01, v_center, v_t01);
    ADD_EDGE(e_center_t02, v_center, v_t02);
    ADD_EDGE(e_center_t03, v_center, v_t03);
    ADD_EDGE(e_center_t12, v_center, v_t12);
    ADD_EDGE(e_center_t13, v_center, v_t13);
    ADD_EDGE(e_center_t23, v_center, v_t23);

    // Add faces adjacent to barycenter
    hehs = { halfEdgeHdl(e_center_t01, VM::Sequence), e_t01_t02, halfEdgeHdl(e_center_t02, VM::Reverse) };
    ADD_FACE(f_center_t01_t02);
    hehs = { halfEdgeHdl(e_center_t01, VM::Sequence), reverse(e_t01_t03), halfEdgeHdl(e_center_t03, VM::Reverse) };
    ADD_FACE(f_center_t01_t03);
    hehs = { halfEdgeHdl(e_center_t01, VM::Sequence), reverse(e_t01_t12), halfEdgeHdl(e_center_t12, VM::Reverse) };
    ADD_FACE(f_center_t01_t12);
    hehs = { halfEdgeHdl(e_center_t01, VM::Sequence), e_t01_t13, halfEdgeHdl(e_center_t13, VM::Reverse) };
    ADD_FACE(f_center_t01_t13);
    hehs = { halfEdgeHdl(e_center_t02, VM::Sequence), e_t02_t03, halfEdgeHdl(e_center_t03, VM::Reverse) };
    ADD_FACE(f_center_t02_t03);
    hehs = { halfEdgeHdl(e_center_t02, VM::Sequence), e_t02_t12, halfEdgeHdl(e_center_t12, VM::Reverse) };
    ADD_FACE(f_center_t02_t12);
    hehs = { halfEdgeHdl(e_center_t02, VM::Sequence), reverse(e_t02_t23), halfEdgeHdl(e_center_t23, VM::Reverse) };
    ADD_FACE(f_center_t02_t23);
    hehs = { halfEdgeHdl(e_center_t03, VM::Sequence), reverse(e_t03_t13), halfEdgeHdl(e_center_t13, VM::Reverse) };
    ADD_FACE(f_center_t03_t13);
    hehs = { halfEdgeHdl(e_center_t03, VM::Sequence), e_t03_t23, halfEdgeHdl(e_center_t23, VM::Reverse) };
    ADD_FACE(f_center_t03_t23);
    hehs = { halfEdgeHdl(e_center_t12, VM::Sequence), reverse(e_t12_t13), halfEdgeHdl(e_center_t13, VM::Reverse) };
    ADD_FACE(f_center_t12_t13);
    hehs = { halfEdgeHdl(e_center_t12, VM::Sequence), e_t12_t23, halfEdgeHdl(e_center_t23, VM::Reverse) };
    ADD_FACE(f_center_t12_t23);
    hehs = { halfEdgeHdl(e_center_t13, VM::Sequence), reverse(e_t13_t23), halfEdgeHdl(e_center_t23, VM::Reverse) };
    ADD_FACE(f_center_t13_t23);

    // Add center cells
    std::array<VM::CellHandle, 8> center_conn_cells;
    // tet center-t01-t02-t03
    hfhs = {
      halfFaceHdl(f_center_t01_t02, VM::Reverse), halfFaceHdl(f_center_t01_t03, VM::Sequence),
      halfFaceHdl(f_center_t02_t03, VM::Reverse), halfFaceHdl(f_t01_t02_t03, VM::Sequence)
    };
    ADD_CELL();
    cell_vrts[just_add_cell.idx()] = { v_center, v_t01, v_t02, v_t03 };
    center_conn_cells[0] = just_add_cell;
    // tet center-t01-t02-t12
    hfhs = {
      halfFaceHdl(f_center_t01_t02, VM::Sequence), halfFaceHdl(f_center_t01_t12, VM::Reverse),
      halfFaceHdl(f_center_t02_t12, VM::Sequence), halfFaceHdl(child(faceHdl(f_t0_t1_t2),3), ori_012)
    };
    ADD_CELL();
    cell_vrts[just_add_cell.idx()] = { v_center, v_t01, v_t02, v_t12 };
    center_conn_cells[1] = just_add_cell;
    // tet center-t01-t03-t13
    hfhs = {
      halfFaceHdl(f_center_t01_t03, VM::Reverse), halfFaceHdl(f_center_t01_t13, VM::Sequence),
      halfFaceHdl(f_center_t03_t13, VM::Reverse), halfFaceHdl(child(faceHdl(f_t0_t1_t3),3), ori_013)
    };
    ADD_CELL();
    cell_vrts[just_add_cell.idx()] = { v_center, v_t01, v_t03, v_t13 };
    center_conn_cells[2] = just_add_cell;
    // tet center-t01-t12-t13
    hfhs = {
      halfFaceHdl(f_center_t01_t12, VM::Sequence), halfFaceHdl(f_center_t01_t13, VM::Reverse),
      halfFaceHdl(f_center_t12_t13, VM::Sequence), halfFaceHdl(f_t01_t12_t13, VM::Sequence)
    };
    ADD_CELL();
    cell_vrts[just_add_cell.idx()] = { v_center, v_t01, v_t12, v_t13 };
    center_conn_cells[3] = just_add_cell;
    // tet center-t02-t03-t23
    hfhs = {
      halfFaceHdl(f_center_t02_t03, VM::Sequence), halfFaceHdl(f_center_t02_t23, VM::Reverse),
      halfFaceHdl(f_center_t03_t23, VM::Sequence), halfFaceHdl(child(faceHdl(f_t0_t2_t3),3), ori_023)
    };
    ADD_CELL();
    cell_vrts[just_add_cell.idx()] = { v_center, v_t02, v_t03, v_t23 };
    center_conn_cells[4] = just_add_cell;
    // tet center-t02-t12-t23
    hfhs = {
      halfFaceHdl(f_center_t02_t12, VM::Reverse), halfFaceHdl(f_center_t02_t23, VM::Sequence),
      halfFaceHdl(f_center_t12_t23, VM::Reverse), halfFaceHdl(f_t02_t12_t23, VM::Sequence)
    };
    ADD_CELL();
    cell_vrts[just_add_cell.idx()] = { v_center, v_t02, v_t12, v_t23 };
    center_conn_cells[5] = just_add_cell;
    // tet center-t03-t13-t23
    hfhs = {
      halfFaceHdl(f_center_t03_t13, VM::Sequence), halfFaceHdl(f_center_t03_t23, VM::Reverse),
      halfFaceHdl(f_center_t13_t23, VM::Sequence), halfFaceHdl(f_t03_t13_t23, VM::Sequence)
    };
    ADD_CELL();
    cell_vrts[just_add_cell.idx()] = { v_center, v_t03, v_t13, v_t23 };
    center_conn_cells[6] = just_add_cell;
    // tet center-t12-t13-t23
    hfhs = {
      halfFaceHdl(f_center_t12_t13, VM::Reverse), halfFaceHdl(f_center_t12_t23, VM::Sequence),
      halfFaceHdl(f_center_t13_t23, VM::Reverse), halfFaceHdl(child(faceHdl(f_t1_t2_t3),3), ori_123)
    };
    ADD_CELL();
    cell_vrts[just_add_cell.idx()] = { v_center, v_t12, v_t13, v_t23 };
    center_conn_cells[7] = just_add_cell;

    updateProp(ch, children);

    auto checkPositiveVolume = [&](VM::CellHandle ch, VM::VertexHandle vh, const Point_3& approx_p)
    {
      VM::HalfFaceHandle opp_hfh = sub_mesh.halffaceOppositeVertex(ch, vh);
      std::array<VM::VertexHandle, 3> fvs = sub_mesh.findHFV(opp_hfh);
      CGAL::Sign ori = CGAL::orientation(
        sub_mesh.exact_point(fvs[0]).exact(), sub_mesh.exact_point(fvs[1]).exact(),
        sub_mesh.exact_point(fvs[2]).exact(), approx_p);
      return ori == CGAL::ON_POSITIVE_SIDE;
    };

    // exact point and float point are same, skip.
    if (!sub_mesh.exact_point(v_center).isSame())
    {
      // get rounded point
      Point_3 round_ep = sub_mesh.exact_point(v_center).round();
      // for each connected cell, check volume after rounding
      bool all_positive = true;
      for (VM::CellHandle ch : center_conn_cells)
      {
        if (!checkPositiveVolume(ch, v_center, round_ep))
        {
          all_positive = false;
          break;
        }
      }
      if (all_positive)
      {
        // do real rounding
        sub_mesh.exact_point(v_center).rounded(sub_mesh.point(v_center));
      }
    }
  }
  VMesh->clear();

  sub_mesh.updateVertexConnCells(cell_vrts);
  sub_mesh.updateFaceConnCells();
}
#undef ADD_EDGE
#undef ADD_FACE
#undef ADD_CELL

void TetSubdivider::updateProp(VM::CellHandle _parent, const std::vector<VM::CellHandle>& _children)
{
  const VM::VMeshT::CT& parent = VMesh->cell(_parent);
  for (size_t i = 0;i < _children.size();i++)
  {
    VM::VMeshT::CT& child = sub_mesh.cell(_children[i]);
    child.prop.is_inside = parent.prop.is_inside;
    if (!parent.prop.is_adjacent_constraint)
    {
      child.prop.is_adjacent_constraint = false;
    }
    else
    {
      bool is_constraint = false;
      for (VM::HalfFaceHandle hfh : child.halffaces)
      {
        if (sub_mesh.face(hfh).prop.is_constraint)
        {
          is_constraint = true;
          goto CONSTRAINT_END;
        }
        for (VM::HalfEdgeHandle heh : sub_mesh.face(hfh).halfedges)
        {
          if (sub_mesh.edge(heh).prop.is_constraint)
          {
            is_constraint = true;
            goto CONSTRAINT_END;
          }
          for (VM::VertexHandle vh : sub_mesh.edge(heh).vertices)
          {
            if (sub_mesh.vertex(vh).prop.is_constraint)
            {
              is_constraint = true;
              goto CONSTRAINT_END;
            }
          }
        }
      }
    CONSTRAINT_END:
      child.prop.is_adjacent_constraint = is_constraint;
    }
  }
}
}// namespace CageInit
}// namespace Cage