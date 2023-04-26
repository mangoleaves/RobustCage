#include "MeshCluster.h"
#include "MeshConnectRegion.h"

namespace Cage
{
namespace SurfaceMesh
{
void LloydCluster::run(size_t cluster_size, size_t max_iteration)
{
  S = cluster_size;
  InitialSeeding();
  GeometryPartitioning();
  ProxyFitting();
  for (int i = 0; i < max_iteration; i++)
  {
    GenerateSeeds();
    GeometryPartitioning();
    if (!ProxyFitting())
      break;
  }
}

std::vector<FaceHandle> LloydCluster::centerFaces()
{
  // for each proxy, find triangles that are closest to proxy as center.
  std::vector<double> min_distance(K, DBL_MAX);
  std::vector<FaceHandle> cluster_centers(K);
  for (auto fh : mesh->faces())
  {
    int proxy_idx = partition[fh.idx()];
    if (proxy_idx == -1)
      continue;
    double distance = (face_centroid[fh.idx()] - proxys[proxy_idx].point).sqrnorm();
    if (distance < min_distance[proxy_idx])
    {
      min_distance[proxy_idx] = distance;
      cluster_centers[proxy_idx] = fh;
    }
  }
  // add unclassfied faces to result
  for(auto fh : mesh->faces())
  {
    int proxy_idx = partition[fh.idx()];
    if (proxy_idx == -1)
      cluster_centers.push_back(fh);
  }
  return cluster_centers;
}

void LloydCluster::InitialSeeding()
{
  // compute normal, centroid and area.
  mesh->request_face_normals();
  mesh->update_face_normals();

  face_centroid.clear(); face_centroid.reserve(mesh->n_faces());
  for (auto fh : mesh->faces())
  {
    face_centroid.push_back(mesh->calc_face_centroid(fh));
  }

  partition.clear();
  partition.resize(mesh->n_faces(), -1);

  // find connecting regions
  MeshConnectRegion MCR(mesh);
  MCR.find();
  const std::vector<std::vector<FaceHandle>>& connect_regions = MCR.get();

  K = 0;
  for (const std::vector<FaceHandle>& region : connect_regions)
  {
    // compute cluster number depending on cluster size
    size_t num_face = region.size();
    size_t region_K = num_face / S;
    if (region_K <= 3)
      continue; // if a connect region is too small, we skip it.
    K += region_K;

    // initialize seed triangles uniformly
    for (size_t i = 0, k = 0; k < region_K; i += S, k++)
    {
      auto fh = region[i];
      seed_tris.push_back(fh);
      size_t proxy_idx = proxys.size();
      proxys.push_back(Proxy(proxy_idx, 1, face_centroid[fh.idx()], mesh->normal(fh)));
      partition[fh.idx()] = (int)proxy_idx;
    }
  }
}

double LloydCluster::CalcError(OpenMesh::SmartFaceHandle& fh, Proxy& proxy)
{
  return mesh->data(fh).face_area * (mesh->normal(fh) - proxy.normal).sqrnorm();
}

double LloydCluster::CalcErrorWithBarrier(OpenMesh::SmartFaceHandle& fh, Proxy& proxy)
{
  double geo_error = mesh->data(fh).face_area * (mesh->normal(fh) - proxy.normal).sqrnorm();
  if (proxy.size <= S)
    return geo_error;
  else
    return geo_error + (proxy.size - S);  // barrier
}

void LloydCluster::GeometryPartitioning()
{
  TriangleQueue triq;
  std::vector<size_t> state(mesh->n_faces(), 0);
  // push neighbor face of seed triangles into queue.
  // label is the partition idx of seed triangle.
  for (auto fh : seed_tris)
  {
    for (auto ff : mesh->ff_range(fh))
      triq.emplace(ff, partition[fh.idx()], (size_t)0, CalcError(ff, proxys[partition[fh.idx()]]));
  }
  while (!triq.empty())
  {
    QueueElem minElem = triq.top();
    triq.pop();
    if (minElem.state < state[minElem.fh.idx()]) // outdated
      continue;
    else if (partition[minElem.fh.idx()] >= 0)
      // element is already put into cluster, skip.
      continue;
    else
    {
      // put element with minimal error to corresponding cluster.
      // then update queue.
      partition[minElem.fh.idx()] = (int)minElem.label;
      proxys[minElem.label].size += 1;
      for (auto ff : mesh->ff_range(minElem.fh)) if (partition[ff.idx()] == -1)
      {
        state[ff.idx()] += 1;
        triq.emplace(ff, partition[minElem.fh.idx()], state[ff.idx()], CalcErrorWithBarrier(ff, proxys[partition[minElem.fh.idx()]]));
      }
    }
  }
}

bool LloydCluster::ProxyFitting()
{
  std::vector<double> region_area(K, 0.0);
  std::vector<Vec3d> region_center(K, Vec3d(0, 0, 0));
  std::vector<Vec3d> region_normal(K, Vec3d(0, 0, 0));

  // compute area, weighted center and weighted normal of each region.
  for (auto fh : mesh->faces())
  {
    int proxy_idx = partition[fh.idx()];
    if (proxy_idx == -1)
      continue;

    double area = mesh->data(fh).face_area;
    region_area[proxy_idx] += area;
    region_center[proxy_idx] += area * face_centroid[fh.idx()];
    region_normal[proxy_idx] += area * mesh->normal(fh);
  }
  for (size_t i = 0;i < K;i++)
  {
    region_center[i] /= region_area[i];
    region_normal[i] /= region_area[i];
    region_normal[i].normalize();
  }
  // update proxy
  bool is_updated = false;
  for (int i = 0; i < K; i++)
  {
    if ((proxys[i].point - region_center[i]).length() > 1e-5)
    {
      proxys[i].point = region_center[i];
      is_updated = true;
    }
    if ((proxys[i].normal - region_normal[i]).length() > 1e-5)
    {
      proxys[i].normal = region_normal[i];
      is_updated = true;
    }
  }
  return is_updated;
}

void LloydCluster::GenerateSeeds()
{
  // for each proxy, find triangles that are closest to proxy as new seeds.
  std::vector<double> min_error(K, DBL_MAX);
  for (auto fh : mesh->faces())
  {
    int proxy_idx = partition[fh.idx()];
    if (proxy_idx == -1)
      continue;
    double error = CalcError(fh, proxys[proxy_idx]);
    if (error < min_error[proxy_idx])
    {
      min_error[proxy_idx] = error;
      seed_tris[proxy_idx] = fh;
    }
  }
  // reset proxys and partition
  for (auto fh : mesh->faces())
  {
    if (partition[fh.idx()] != -1 &&
      fh.idx() != seed_tris[partition[fh.idx()]].idx())
    {
      partition[fh.idx()] = -1;
    }
  }
  for (auto& proxy : proxys)
  {
    proxy.size = 1;
  }
}
}// namespace SurfaceMesh
}// namespace Cage