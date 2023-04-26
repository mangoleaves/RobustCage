#pragma once
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include "Geometry/Basic/Types.h"
#include "Geometry/Exact/ExactPoint.h"
#include "Utils/logger.hh"
#include "Utils/Marker.h"

namespace Cage
{
namespace SurfaceMesh
{
typedef Geometry::Vec3d Vec3d;
typedef Geometry::ExactPoint ExactPoint;

struct MeshTraits : public OpenMesh::DefaultTraits
{
  typedef Geometry::Vec3d Point;
  typedef Geometry::Vec3d Normal;

  typedef typename std::pair<Point, Point> Link;
  typedef typename std::vector<Link> LinkVector;

  VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
  FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
  EdgeAttributes(OpenMesh::Attributes::Status);
  HalfedgeAttributes(OpenMesh::Attributes::Status);

  FaceTraits
  {
  public:
    std::set<OpenMesh::FaceHandle> one_ring_faces;

    size_t samples_num = 0;
    LinkVector face_in_links;

    double face_area = 0.0;

    double in_error = 0.0;
    double out_error = 0.0;
  };

  EdgeTraits
  {
  public:
    size_t samples_num = 0;

    double edge_length = 0.0;

    double target_length = 0.0;

    double out_error = 0.0;
  };

  HalfedgeTraits
  {
  };

  VertexTraits
  {
  public:
    std::unique_ptr<ExactPoint> ep;

    double target_length;

    double out_error;
    double out_surround_error;

    VertexT() : ep(nullptr), target_length(0.0), out_error(0.0), out_surround_error(0.0)
    {}
    // copy constructor
    VertexT(const VertexT& rhs)
    {
      if (rhs.ep) ep = std::make_unique<ExactPoint>(*rhs.ep);
      target_length = rhs.target_length;
      out_error = rhs.out_error;
      out_surround_error = rhs.out_surround_error;
    }
    VertexT& operator=(const VertexT& rhs)
    {
      if (rhs.ep) ep = std::make_unique<ExactPoint>(*rhs.ep);
      target_length = rhs.target_length;
      out_error = rhs.out_error;
      out_surround_error = rhs.out_surround_error;
      return *this;
    }
    // move constructor
    VertexT(VertexT&& rhs)
    {
      if (rhs.ep)  ep = std::move(rhs.ep);
      target_length = rhs.target_length;
      out_error = rhs.out_error;
      out_surround_error = rhs.out_surround_error;
    }
    VertexT& operator=(VertexT&& rhs)
    {
      if (rhs.ep)  ep = std::move(rhs.ep);
      target_length = rhs.target_length;
      out_error = rhs.out_error;
      out_surround_error = rhs.out_surround_error;
      return *this;
    }
  };
};

typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> SMeshT;
typedef typename MeshTraits::Link        Link;
typedef typename MeshTraits::LinkVector  LinkVector;

using OpenMesh::VertexHandle;
using OpenMesh::EdgeHandle;
using OpenMesh::HalfedgeHandle;
using OpenMesh::FaceHandle;

void pre_calculate_edge_length(SMeshT* mesh);
void pre_calculate_face_area(SMeshT* mesh);

std::array<std::pair<HalfedgeHandle, double>, 3>
face_sector_angle(SMeshT* mesh, FaceHandle fh);

double calc_face_area(SMeshT* mesh, FaceHandle fh);
double calc_mesh_area(SMeshT* mesh);

/// graph vertex color
void graph_color_vertex(const SMeshT& mesh, std::vector<std::vector<int>>& phase);
}// namespace SurfaceMesh
}// namespace Cage