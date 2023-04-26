#include "Algorithms.h"

namespace Cage
{
namespace Geometry
{

int find_obtuse_angle(const Triangle& triangle)
{
  for (int i = 0; i < 3; i++)
  {
    if (dot_sign(
      triangle[(i + 1) % 3] - triangle[i],
      triangle[(i + 2) % 3] - triangle[i]) == -1)
    {
      return i;
    }
  }
  return -1;
}

Vec3d calc_vertical_vec2vec(const Vec3d& vec)
{
  if (vec[0] != 0.0)
  {
    return Vec3d(-(vec[1] + vec[2]) / vec[0], 1.0, 1.0).normalized();
  }
  else	// vec must be in y-z-plane or z-axis
  {
    return Vec3d(1.0, 0, 0);
  }
}

/// @brief find which octant the vec lies in.
int find_octant(const Vec3d& vec)
{
  static int octant[2][2][2] = { {{6, 2}, {5, 1}}, {{7, 3},{4, 0}} };
  int xsign = vec[0] >= 0 ? 1 : 0;
  int ysign = vec[1] >= 0 ? 1 : 0;
  int zsign = vec[2] >= 0 ? 1 : 0;
  return octant[xsign][ysign][zsign];
}

}// namespace Geometry
}// namespace Cage