#pragma once

#include <vector>
#include <cstdint>
#include <cstddef>

namespace Cage
{
namespace SimpleUtils
{

class Marker
{
public:
  Marker() = default;
  Marker(size_t primitive_size) { m_mark.resize(primitive_size, false); }

  void mark(int idx);
  
  inline bool is_marked(int idx) { return m_mark[idx]; }

  void unmark_all();
private:
  std::vector<uint8_t> m_mark;
  std::vector<int> m_marked;
};

}// namespace SimpleUtils
}// namespace Cage