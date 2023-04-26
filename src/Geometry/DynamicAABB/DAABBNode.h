#pragma once

#include "Geometry/Basic/BoundingBox.h"

namespace Cage
{
namespace Geometry
{
class DAABBNode
{
public:
  BoundingBox m_bbox;
  size_t m_nb_primitives;
  int m_parent, m_p_left_child, m_p_right_child;
  int m_split_dim;
public:
  DAABBNode() :
    m_bbox(),
    m_nb_primitives(0),
    m_parent(-1),
    m_p_left_child(-1),
    m_p_right_child(-1),
    m_split_dim(-1)
  {}

  ~DAABBNode() {}

  void clear()
  {
    m_nb_primitives = 0;
    m_parent = -1;
    m_p_left_child = -1;
    m_p_right_child = -1;
  }
};
}// namespace Geometry
}// namespace Cage