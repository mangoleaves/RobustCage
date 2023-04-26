#pragma once

#include "SurfaceMeshDefinition.h"
#include "MeshConnectRegion.h"
#include <queue>

namespace Cage
{
namespace SurfaceMesh
{

using namespace OpenMesh;

class LloydCluster
{
public:
  LloydCluster() = default;
  LloydCluster(SMeshT* _mesh) :mesh(_mesh) {}

  void run(size_t cluster_size, size_t max_iteration);
  inline size_t getClusterNum() { return K; }
  inline const std::vector<int>& getPartition() { return partition; }
  inline const std::vector<Vec3d>& getCentroid() { return face_centroid; }
  std::vector<FaceHandle> centerFaces();
private:
  typedef std::vector<FaceHandle> SeedTriangles;

  struct Proxy
  {
    size_t idx;
    size_t size;
    Vec3d point;
    Vec3d normal;

    Proxy() = default;
    Proxy(size_t i, size_t s, const Vec3d& p, const Vec3d& n) : idx(i), size(s), point(p), normal(n) {}
  };

  typedef std::vector<Proxy> Proxys;

  struct QueueElem
  {
    FaceHandle fh;
    size_t label;
    size_t state;
    double priority;

    QueueElem() = default;
    QueueElem(FaceHandle f, size_t l, size_t s, double p) :fh(f), label(l), state(s), priority(p) {}

    bool operator<(const QueueElem& rhs)const { return priority < rhs.priority; }
    bool operator>(const QueueElem& rhs)const { return priority > rhs.priority; }
  };

  typedef std::priority_queue<QueueElem, std::vector<QueueElem>, greater<QueueElem>> TriangleQueue;

  // Mesh

  SMeshT* mesh;
  std::vector<Vec3d> face_centroid;

  // Lloyd iteration

  size_t K;                     // cluster number
  size_t S;                     // cluster size
  Proxys proxys;                // proxy planes
  SeedTriangles seed_tris;      // seed triangles
  std::vector<int> partition;   // great than or equal to 0: proxy id. equal to -1: unclassified.

  double CalcError(OpenMesh::SmartFaceHandle& fh, Proxy& proxy);
  double CalcErrorWithBarrier(OpenMesh::SmartFaceHandle& fh, Proxy& proxy);
  void InitialSeeding();
  void GeometryPartitioning();
  bool ProxyFitting();
  void GenerateSeeds();
};

}// namespace SurfaceMesh
}// namespace Cage