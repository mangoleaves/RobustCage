#pragma once

#include "Geometry/Basic/Types.h"
#include "Geometry/Exact/Predicates.h"

namespace Cage
{
namespace Geometry
{
typedef const Point& CR_FP;       // Const Reference to Float Point

bool do_intersect(CR_FP a, CR_FP b, CR_FP c, CR_FP p, CR_FP q);

}// namespace Geometry
}// namespace Cage