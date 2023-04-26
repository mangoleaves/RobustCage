#pragma once

#include "Geometry/Basic/Triangle.h"
#include "Geometry/Exact/Predicates.h"

namespace Cage
{
namespace Geometry
{
typedef const Point& CR_FP;       // Const Reference to Float Point
typedef const ExactPoint* CP_EP;  // Const Pointer to Exact Point
typedef const Point_3& CR_P3;     // Const Reference to Point_3

bool triangle_do_intersect(
  CR_FP p, CR_FP q, CR_FP r,
  CR_FP a, CR_FP b, CR_FP c,
  CP_EP ep, CP_EP eq, CP_EP er,
  CP_EP ea, CP_EP eb, CP_EP ec
);

bool triangle_do_overlap(
  CR_FP p, CR_FP common_a, CR_FP common_b, CR_FP q,
  CP_EP ep, CP_EP common_ea, CP_EP common_eb, CP_EP eq
);

bool triangle_do_intersect(
  CR_FP p, CR_FP q, CR_FP common_v, CR_FP a, CR_FP b,
  CP_EP ep, CP_EP eq, CP_EP common_ev, CP_EP ea, CP_EP eb
);
}// namespace Geometry
}// namespace Cage