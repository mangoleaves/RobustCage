#pragma once

#include "CGAL/Exact_predicates_exact_constructions_kernel.h"

namespace Cage
{
namespace Geometry
{

typedef CGAL::Exact_predicates_exact_constructions_kernel K;

typedef K::Epeck::FT ET;

typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Triangle_2 Triangle_2;
typedef K::Intersect_2 Intersect_2;

typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef K::Segment_3 Segment_3;
typedef K::Plane_3 Plane_3;
typedef K::Triangle_3 Triangle_3;
typedef K::Iso_cuboid_3 Iso_cuboid_3;
typedef K::Intersect_3 Intersect_3;
typedef K::Tetrahedron_3 Tetrahedron_3;

}// namespace Geometry
}// namespace Cage