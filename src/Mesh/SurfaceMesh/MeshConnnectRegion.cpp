#include "MeshConnectRegion.h"
#include <queue>

namespace Cage
{
namespace SurfaceMesh
{

void MeshConnectRegion::find()
{
  std::vector<uint8_t> visited(mesh->n_faces(), false);

  auto FindUnvisited = [&](FaceHandle& fh) {
    for (size_t i = 0;i < visited.size();i++)
    {
      if (!visited[i])
      {
        fh = FaceHandle((int)i);
        return true;
      }
    }
    return false;
  };

  FaceHandle walk_fh;
  while (FindUnvisited(walk_fh))
  {
    regions.emplace_back();

    std::queue<FaceHandle> adj_fhs;
    adj_fhs.push(walk_fh);
    visited[walk_fh.idx()] = true;

    while (!adj_fhs.empty())
    {
      walk_fh = adj_fhs.front();
      adj_fhs.pop();

      regions.back().push_back(walk_fh);
      for (auto ff : mesh->ff_range(walk_fh))
      {
        if (!visited[ff.idx()])
        {
          adj_fhs.push(ff);
          visited[ff.idx()] = true;
        }
      }
    }
  }
}

}// namespace SurfaceMesh
}// namespace Cage