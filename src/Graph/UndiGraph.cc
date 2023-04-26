#include "UndiGraph.hh"

namespace Cage
{
namespace Graph
{

bool UndiGraph::isConnected()
{
  isVisited.clear(); isVisited.resize(adjacentMat.size(), false);
  DFSConnection(0);
  return std::find(isVisited.begin(), isVisited.end(), false) == isVisited.end();
}

/// @brief by deep first search, find all vertices connected to v[idx].
void UndiGraph::DFSConnection(int idx)
{
  isVisited[idx] = true;
  for (int adj_idx : adjacentMat[idx])
  {
    if (!isVisited[adj_idx])
      DFSConnection(adj_idx);
  }
}

}// namespace Graph
}// namespace Cage