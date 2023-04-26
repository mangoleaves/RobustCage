#include "TriSegIntersect.h"

namespace Cage
{
namespace Geometry
{

// forward declaration
bool do_intersect_coplanar(CR_FP A, CR_FP B, CR_FP C, CR_FP p, CR_FP q);

/// @brief test whether triangle(abc) and segment(pq) intersect.
/// ref to CGAL/Segment_3_Triangle_3_do_intersect.h
bool do_intersect(CR_FP a, CR_FP b, CR_FP c, CR_FP p, CR_FP q)
{
  int orip = orient3d_wrapped(a, b, c, p);
  int oriq = orient3d_wrapped(a, b, c, q);

  switch (orip)
  {
  case POSITIVE:
    switch (oriq)
    {
    case POSITIVE:
      // the segment lies in the positive open halfspaces defined by the
      // triangle's supporting plane
      return false;
    case NEGATIVE:
      // p sees the triangle in counterclockwise order
      return orient3d_wrapped(p, q, a, b) != POSITIVE
        && orient3d_wrapped(p, q, b, c) != POSITIVE
        && orient3d_wrapped(p, q, c, a) != POSITIVE;
    case ZERO:
      // q belongs to the triangle's supporting plane
      // p sees the triangle in counterclockwise order
      return orient3d_wrapped(p, q, a, b) != POSITIVE
        && orient3d_wrapped(p, q, b, c) != POSITIVE
        && orient3d_wrapped(p, q, c, a) != POSITIVE;
    default: // should not happen.
      CGAL_kernel_assertion(false);
      return false;
    }
  case NEGATIVE:
    switch (oriq)
    {
    case POSITIVE:
      // q sees the triangle in counterclockwise order
      return orient3d_wrapped(q, p, a, b) != POSITIVE
        && orient3d_wrapped(q, p, b, c) != POSITIVE
        && orient3d_wrapped(q, p, c, a) != POSITIVE;
    case NEGATIVE:
      // the segment lies in the negative open halfspaces defined by the
      // triangle's supporting plane
      return false;
    case ZERO:
      // q belongs to the triangle's supporting plane
      // p sees the triangle in clockwise order
      return orient3d_wrapped(q, p, a, b) != POSITIVE
        && orient3d_wrapped(q, p, b, c) != POSITIVE
        && orient3d_wrapped(q, p, c, a) != POSITIVE;

    default: // should not happen.
      CGAL_kernel_assertion(false);
      return false;
    }
  case ZERO: // p belongs to the triangle's supporting plane
    switch (oriq)
    {
    case POSITIVE:
      // q sees the triangle in counterclockwise order
      return orient3d_wrapped(q, p, a, b) != POSITIVE
        && orient3d_wrapped(q, p, b, c) != POSITIVE
        && orient3d_wrapped(q, p, c, a) != POSITIVE;
    case NEGATIVE:
      // q sees the triangle in clockwise order
      return orient3d_wrapped(p, q, a, b) != POSITIVE
        && orient3d_wrapped(p, q, b, c) != POSITIVE
        && orient3d_wrapped(p, q, c, a) != POSITIVE;
    case ZERO:
      // the segment is coplanar with the triangle's supporting plane
      // we test whether the segment intersects the triangle in the common
      // supporting plane
      return do_intersect_coplanar(a, b, c, p, q);

    default: // should not happen.
      ASSERT(false, "");
      return false;
    }
  default: // should not happen.
    ASSERT(false, "");
    return false;
  }
}

bool do_intersect_coplanar(CR_FP A, CR_FP B, CR_FP C, CR_FP p, CR_FP q)
{
  const Vec3d* a = &A;
  const Vec3d* b = &B;
  const Vec3d* c = &C;

  // Determine the orientation of the triangle in the common plane

  if (coplanar_orient2d(A, B, C) != POSITIVE)
  {
    // The triangle is not counterclockwise oriented
    // swap two vertices.
    b = &C;
    c = &B;
  }

  // Test whether the segment's supporting line intersects the
  // triangle in the common plane

  int pqa = coplanar_orient2d(p, q, *a);
  int pqb = coplanar_orient2d(p, q, *b);
  int pqc = coplanar_orient2d(p, q, *c);

  switch (pqa)
  {
  case POSITIVE:
    switch (pqb)
    {
    case POSITIVE:
      if (pqc == POSITIVE)
        return false;

      // the triangle lies in the positive halfspace
      // defined by the segment's supporting line.
    // c is isolated on the negative side
      return coplanar_orient2d(*b, *c, q) != NEGATIVE
        && coplanar_orient2d(*c, *a, p) != NEGATIVE;
    case NEGATIVE:
      if (pqc == POSITIVE) // b is isolated on the negative side
        return coplanar_orient2d(*a, *b, q) != NEGATIVE
        && coplanar_orient2d(*b, *c, p) != NEGATIVE;
      // a is isolated on the positive side
      return coplanar_orient2d(*a, *b, q) != NEGATIVE
        && coplanar_orient2d(*c, *a, p) != NEGATIVE;
    case ZERO:
      if (pqc == POSITIVE) // b is isolated on the negative side
        return coplanar_orient2d(*a, *b, q) != NEGATIVE
        && coplanar_orient2d(*b, *c, p) != NEGATIVE;
      // a is isolated on the positive side
      return coplanar_orient2d(*a, *b, q) != NEGATIVE
        && coplanar_orient2d(*c, *a, p) != NEGATIVE;
    default:// should not happen.
      CGAL_kernel_assertion(false);
      return false;
    }
  case NEGATIVE:
    switch (pqb)
    {
    case POSITIVE:
      if (pqc == POSITIVE) // a is isolated on the negative side
        return coplanar_orient2d(*a, *b, p) != NEGATIVE
        && coplanar_orient2d(*c, *a, q) != NEGATIVE;
      // b is isolated on the positive side
      return coplanar_orient2d(*a, *b, p) != NEGATIVE
        && coplanar_orient2d(*b, *c, q) != NEGATIVE;
    case NEGATIVE:
      if (pqc == NEGATIVE)
        return false;

      // the triangle lies in the negative halfspace
      // defined by the segment's supporting line.
      // c is isolated on the positive side
      return coplanar_orient2d(*b, *c, p) != NEGATIVE
        && coplanar_orient2d(*c, *a, q) != NEGATIVE;
    case ZERO:
      if (pqc == NEGATIVE) // b is isolated on the positive side
        return coplanar_orient2d(*a, *b, p) != NEGATIVE
        && coplanar_orient2d(*b, *c, q) != NEGATIVE;
      // a is isolated on the negative side
      return coplanar_orient2d(*a, *b, p) != NEGATIVE
        && coplanar_orient2d(*c, *a, q) != NEGATIVE;

    default:// should not happen.
      CGAL_kernel_assertion(false);
      return false;
    }
  case ZERO:
    switch (pqb)
    {
    case POSITIVE:
      if (pqc == POSITIVE) // a is isolated on the negative side
        return coplanar_orient2d(*a, *b, p) != NEGATIVE
        && coplanar_orient2d(*c, *a, q) != NEGATIVE;
      // b is isolated on the positive side
      return coplanar_orient2d(*a, *b, p) != NEGATIVE
        && coplanar_orient2d(*b, *c, q) != NEGATIVE;
    case NEGATIVE:
      if (pqc == NEGATIVE) // a is isolated on the positive side
        return coplanar_orient2d(*a, *b, q) != NEGATIVE
        && coplanar_orient2d(*c, *a, p) != NEGATIVE;
      // b is isolated on the negative side
      return coplanar_orient2d(*a, *b, q) != NEGATIVE
        && coplanar_orient2d(*b, *c, p) != NEGATIVE;
    case ZERO:
      if (pqc == POSITIVE) // c is isolated on the positive side
        return coplanar_orient2d(*b, *c, p) != NEGATIVE
        && coplanar_orient2d(*c, *a, q) != NEGATIVE;
      // c is isolated on the negative side
      return coplanar_orient2d(*b, *c, q) != NEGATIVE
        && coplanar_orient2d(*c, *a, p) != NEGATIVE;
      // case pqc == COLLINEAR is impossible since the triangle is
      // assumed to be non flat

    default:// should not happen.
      CGAL_kernel_assertion(false);
      return false;

    }
  default:// should not happen.
    CGAL_kernel_assertion(false);
    return false;
  }
}

}// namespace Geometry
}// namespace Cage