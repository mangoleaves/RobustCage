#pragma once

#include "Basic/Types.h"
#include "Basic/Triangle.h"
#include "Exact/Predicates.h"

namespace Cage
{
namespace Geometry
{
int find_obtuse_angle(const Triangle& triangle);

Vec3d calc_vertical_vec2vec(const Vec3d& vec);

int find_octant(const Vec3d& vec);
}// namespace Geometry
}// namespace Cage