#include "Marker.h"

namespace Cage
{
namespace SimpleUtils
{

void Marker::mark(int idx)
{
  if (!is_marked(idx))
  {
    m_mark[idx] = true;
    m_marked.push_back(idx);
  }
}

void Marker::unmark_all()
{
  for (int idx : m_marked)
    m_mark[idx] = false;
  m_marked.clear();
}

}// namespace SimpleUtils
}// namespace Cage