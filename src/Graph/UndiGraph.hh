#pragma once
#include <vector>
#include <cstdint>

namespace Cage
{
namespace Graph
{
/// @brief Undirected Graph
class UndiGraph
{
public:
  std::vector<std::vector<int>> adjacentMat;

  bool isConnected();
private:
  std::vector<uint8_t> isVisited;

  void DFSConnection(int idx);
};
}// namespace Graph
}// namespace Cage