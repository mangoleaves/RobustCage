#include "TriTriIntersect.h"
#include "TriSegIntersect.h"

namespace Cage
{
namespace Geometry
{

// fowrad declaration.
bool coplanar_triangle_do_intersect(
  CR_FP p, CR_FP q, CR_FP r,
  CR_FP a, CR_FP b, CR_FP c
);

bool _intersection_test_edge(
  CR_FP p, CR_FP q, CR_FP r,
  CR_FP a, CR_FP b, CR_FP c
);

bool _intersection_test_vertex(
  CR_FP p, CR_FP q, CR_FP r,
  CR_FP a, CR_FP b, CR_FP c
);

/// @brief check whether two triangles intersect.
/// p, q, r define first triangle. a, b, c define second triangle.
/// order of vertices defines triangle's orientation.
/// @see ref paper: "Faster Triangle-Triangle Intersection Tests".
/// ref code: CGAL/Intersections_3/internal/Triangle_3_Triangle_3_do_intersect.h
bool triangle_do_intersect(
  CR_FP p, CR_FP q, CR_FP r,
  CR_FP a, CR_FP b, CR_FP c
)
{
  const Point* s_min1;
  const Point* t_min1;
  const Point* s_max1;
  const Point* t_max1;
  // check of triA respect to triB
  const int dp = orient3d_wrapped(a, b, c, p);
  const int dq = orient3d_wrapped(a, b, c, q);
  const int dr = orient3d_wrapped(a, b, c, r);

  switch (dp)
  {
  case POSITIVE:
    if (dq == POSITIVE)
    {
      if (dr == POSITIVE)
        // pqr lies in the open positive halfspace induced by
        // the plane of triangle(a,b,c)
        return false;
      // r is isolated on the negative side of the plane
      s_min1 = &q; t_min1 = &r; s_max1 = &r; t_max1 = &p;
    }
    else
    {
      if (dr == POSITIVE)
      {
        // q is isolated on the negative side of the plane
        s_min1 = &p; t_min1 = &q; s_max1 = &q; t_max1 = &r;
      }
      else
      {
        // p is isolated on the positive side of the plane
        s_min1 = &p; t_min1 = &q; s_max1 = &r; t_max1 = &p;
      }
    }
    break;
  case NEGATIVE:
    if (dq == NEGATIVE)
    {
      if (dr == NEGATIVE)
        // pqr lies in the open negative halfspace induced by
        // the plane of triangle(a,b,c)
        return false;
      // r is isolated on the positive side of the plane
      s_min1 = &r; t_min1 = &p; s_max1 = &q; t_max1 = &r;
    }
    else
    {
      if (dr == NEGATIVE)
      {
        // q is isolated on the positive side of the plane
        s_min1 = &q; t_min1 = &r; s_max1 = &p; t_max1 = &q;
      }
      else
      {
        // p is isolated on the negative side of the plane
        s_min1 = &r; t_min1 = &p; s_max1 = &p; t_max1 = &q;
      }
    }
    break;
  case ZERO:
    switch (dq)
    {
    case POSITIVE:
      if (dr == POSITIVE)
      {
        // p is isolated on the negative side of the plane
        s_min1 = &r; t_min1 = &p; s_max1 = &p; t_max1 = &q;
      }
      else
      {
        // q is isolated on the positive side of the plane
        s_min1 = &q; t_min1 = &r; s_max1 = &p; t_max1 = &q;
      }
      break;
    case NEGATIVE:
      if (dr == NEGATIVE)
      {
        // p is isolated on the positive side of the plane
        s_min1 = &p; t_min1 = &q; s_max1 = &r; t_max1 = &p;
      }
      else
      {
        // q is isolated on the negative side of the plane
        s_min1 = &p; t_min1 = &q; s_max1 = &q; t_max1 = &r;
      }
      break;
    case ZERO:
      switch (dr)
      {
      case POSITIVE:
        // r is isolated on the positive side of the plane
        s_min1 = &r; t_min1 = &p; s_max1 = &q; t_max1 = &r;
        break;
      case NEGATIVE:
        // r is isolated on the negative side of the plane
        s_min1 = &q; t_min1 = &r; s_max1 = &r; t_max1 = &p;
        break;
      case ZERO:
        return coplanar_triangle_do_intersect(p, q, r, a, b, c);
      }
      break;
    }
    break;
  }

  const Point* s_min2;
  const Point* t_min2;
  const Point* s_max2;
  const Point* t_max2;

  // Compute distance signs  of a, b and c to the plane of triangle(p,q,r)
  const int da = orient3d_wrapped(p, q, r, a);
  const int db = orient3d_wrapped(p, q, r, b);
  const int dc = orient3d_wrapped(p, q, r, c);


  switch (da)
  {
  case POSITIVE:
    if (db == POSITIVE)
    {
      if (dc == POSITIVE)
        // abc lies in the open positive halfspace induced by
        // the plane of triangle(p,q,r)
        return false;
      // c is isolated on the negative side of the plane
      s_min2 = &b; t_min2 = &c; s_max2 = &c; t_max2 = &a;
    }
    else
    {
      if (dc == POSITIVE)
      {
        // b is isolated on the negative side of the plane
        s_min2 = &a; t_min2 = &b; s_max2 = &b; t_max2 = &c;
      }
      else
      {
        // a is isolated on the positive side of the plane
        s_min2 = &a; t_min2 = &b; s_max2 = &c; t_max2 = &a;
      }
    }
    break;
  case NEGATIVE:
    if (db == NEGATIVE)
    {
      if (dc == NEGATIVE)
        // abc lies in the open negative halfspace induced by
        // the plane of triangle(p,q,r)
        return false;
      // c is isolated on the positive side of the plane
      s_min2 = &c; t_min2 = &a; s_max2 = &b; t_max2 = &c;
    }
    else
    {
      if (dc == NEGATIVE)
      {
        // b is isolated on the positive side of the plane
        s_min2 = &b; t_min2 = &c; s_max2 = &a; t_max2 = &b;
      }
      else
      {
        // a is isolated on the negative side of the plane
        s_min2 = &c; t_min2 = &a; s_max2 = &a; t_max2 = &b;
      }
    }
    break;
  case ZERO:
    switch (db)
    {
    case POSITIVE:
      if (dc == POSITIVE)
      {
        // a is isolated on the negative side of the plane
        s_min2 = &c; t_min2 = &a; s_max2 = &a; t_max2 = &b;
      }
      else
      {
        // b is isolated on the positive side of the plane
        s_min2 = &b; t_min2 = &c; s_max2 = &a; t_max2 = &b;
      }
      break;
    case NEGATIVE:
      if (dc == NEGATIVE)
      {
        // a is isolated on the positive side of the plane
        s_min2 = &a; t_min2 = &b; s_max2 = &c; t_max2 = &a;
      }
      else
      {
        // b is isolated on the negative side of the plane
        s_min2 = &a; t_min2 = &b; s_max2 = &b; t_max2 = &c;
      }
      break;
    case ZERO:
      switch (dc)
      {
      case POSITIVE:
        // c is isolated on the positive side of the plane
        s_min2 = &c; t_min2 = &a; s_max2 = &b; t_max2 = &c;
        break;
      case NEGATIVE:
        // c is isolated on the negative side of the plane
        s_min2 = &b; t_min2 = &c; s_max2 = &c; t_max2 = &a;
        break;
      case ZERO:
        // Supposed unreachable code
        // since the triangles are assumed to be non-flat
        return coplanar_triangle_do_intersect(p, q, r, a, b, c);
      }
      break;
    }
    break;
  }

  return  orient3d_wrapped(*s_min1, *t_min1, *s_min2, *t_min2) != POSITIVE
    && orient3d_wrapped(*s_max1, *t_max1, *t_max2, *s_max2) != POSITIVE;
}

/// @brief triangle pqr and triangle abc are coplanar. check whether
/// they are intersect.
/// assume that pqr and abc are counter-clock-wise.
bool coplanar_triangle_do_intersect(
  CR_FP p, CR_FP q, CR_FP r,
  CR_FP a, CR_FP b, CR_FP c
)
{
  if (coplanar_orient2d(a, b, p) != NEGATIVE)
  {
    if (coplanar_orient2d(b, c, p) != NEGATIVE)
    {
      if (coplanar_orient2d(c, a, p) != NEGATIVE)
        // p is inside triangle abc
        return true;
      // p sees ac
      return _intersection_test_edge(p, q, r, a, b, c);
    }
    if (coplanar_orient2d(c, a, p) != NEGATIVE)//p sees bc
      return _intersection_test_edge(p, q, r, c, a, b);
    // p sees c
    return _intersection_test_vertex(p, q, r, a, b, c);

  }
  if (coplanar_orient2d(b, c, p) != NEGATIVE)
  {
    if (coplanar_orient2d(c, a, p) != NEGATIVE) //p sees ab
      return _intersection_test_edge(p, q, r, b, c, a);
    // p sees a
    return _intersection_test_vertex(p, q, r, b, c, a);
  }
  // p sees b
  return _intersection_test_vertex(p, q, r, c, a, b);
}

bool _intersection_test_edge(
  CR_FP p, CR_FP q, CR_FP r,
  CR_FP a, CR_FP b, CR_FP c
)
{
  if (coplanar_orient2d(c, a, q) != NEGATIVE)
  {  //pq straddles (ac)
    if (coplanar_orient2d(p, a, q) != NEGATIVE)
      return coplanar_orient2d(p, q, c) != NEGATIVE;

    return coplanar_orient2d(q, r, a) != NEGATIVE
      && coplanar_orient2d(r, p, a) != NEGATIVE;
  }

  if (coplanar_orient2d(c, a, r) != NEGATIVE)
  {
    // pr and qr straddle line (pr)
    return coplanar_orient2d(p, a, r) != NEGATIVE
      && (coplanar_orient2d(p, r, c) != NEGATIVE
        || coplanar_orient2d(q, r, c) != NEGATIVE);
  }

  return false; //ac separes
}

bool _intersection_test_vertex(
  CR_FP p, CR_FP q, CR_FP r,
  CR_FP a, CR_FP b, CR_FP c
)
{
  if (coplanar_orient2d(c, a, q) != NEGATIVE)
  {
    if (coplanar_orient2d(c, b, q) != POSITIVE)
    {
      if (coplanar_orient2d(p, a, q) == POSITIVE)
        return  coplanar_orient2d(p, b, q) != POSITIVE;

      return coplanar_orient2d(p, a, r) != NEGATIVE
        && coplanar_orient2d(q, r, a) != NEGATIVE;
    }

    if (coplanar_orient2d(p, b, q) != POSITIVE)
      return coplanar_orient2d(c, b, r) != POSITIVE
      && coplanar_orient2d(q, r, b) != NEGATIVE;
    return false;

  }

  if (coplanar_orient2d(c, a, r) != NEGATIVE)
  { //qr straddles (ac)
    if (coplanar_orient2d(q, r, c) != NEGATIVE)
      return (coplanar_orient2d(p, a, r) != NEGATIVE);

    return coplanar_orient2d(q, r, b) != NEGATIVE
      && coplanar_orient2d(c, r, b) != NEGATIVE;
  }
  return false; // ca separes
}

/// @brief test whether two triangle overlap.
/// Two triangle are coplanar and adjacent. They have a common edge.
/// @note Filter Version
bool triangle_do_overlap(
  CR_FP p,
  CR_FP common_a, CR_FP common_b,
  CR_FP q)
{
  if (orient3d_wrapped(p, common_a, common_b, q) == ZERO)
  {
    // q is on plane of p-a-b
    int ori_p = coplanar_orient2d(common_a, common_b, p);
    int ori_q = coplanar_orient2d(common_a, common_b, q);
    return ori_p * ori_q >= 0;
    // include cases:
    // ori_p == 0 or ori_q == 0, degenerate.
    // ori_p == ori_q, same side.
    // ori_p != ori_q, different side.
  }
  else
  {
    // q is not coplanar with pab, won't overlap.
    return false;
  }
}

/// @brief test whether two triangle overlap.
/// Two triangle are coplanar and adjacent. They have a common edge.
/// @note CGAL Version
bool triangle_do_overlap(
  CR_P3 p,
  CR_P3 common_a, CR_P3 common_b,
  CR_P3 q)
{
  if (CGAL::orientation(p, common_a, common_b, q) == CGAL::ON_ORIENTED_BOUNDARY)
  {
    // q is on plane of p-a-b
    int ori_p = CGAL::coplanar_orientation(common_a, common_b, p);
    int ori_q = CGAL::coplanar_orientation(common_a, common_b, q);
    return ori_p * ori_q >= 0;
    // include cases:
    // ori_p == 0 or ori_q == 0, degenerate.
    // ori_p == ori_q, same side.
    // ori_p != ori_q, different side.
  }
  else
  {
    // q is not coplanar with pab, won't overlap.
    return false;
  }
}

/// @brief test whether two triangle intersect. they have no common point.
/// this is a wrapper to CGAL/Fast filter.
bool triangle_do_intersect(
  CR_FP p, CR_FP q, CR_FP r,
  CR_FP a, CR_FP b, CR_FP c,
  CP_EP ep, CP_EP eq, CP_EP er,
  CP_EP ea, CP_EP eb, CP_EP ec)
{
  if (ep || eq || er || ea || eb || ec)
  {
    // if any one point is exact point, we gonna use CGAL's predicates.
    Point_3 pp, pq, pr, pa, pb, pc;
    pp = ep ? ep->exact() : Point_3(p.x(), p.y(), p.z());
    pq = eq ? eq->exact() : Point_3(q.x(), q.y(), q.z());
    pr = er ? er->exact() : Point_3(r.x(), r.y(), r.z());
    pa = ea ? ea->exact() : Point_3(a.x(), a.y(), a.z());
    pb = eb ? eb->exact() : Point_3(b.x(), b.y(), b.z());
    pc = ec ? ec->exact() : Point_3(c.x(), c.y(), c.z());
    Triangle_3 tri0(pp, pq, pr), tri1(pa, pb, pc);

    return CGAL::do_intersect(tri0, tri1);
  }
  else
  {
    // all points are float point, we gonna use fast filter.
    return triangle_do_intersect(p, q, r, a, b, c);
  }
}

/// @brief test whether two triangle intersect. they have a common edge.
/// this is a wrapper to CGAL/Fast filter.
bool triangle_do_overlap(
  CR_FP p, CR_FP common_a, CR_FP common_b, CR_FP q,
  CP_EP ep, CP_EP common_ea, CP_EP common_eb, CP_EP eq)
{
  if (ep || eq || common_ea || common_eb)
  {
    // if any one point is exact point, we gonna use CGAL's predicates.
    Point_3 pp, common_pa, common_pb, pq;
    pp = ep ? ep->exact() : Point_3(p.x(), p.y(), p.z());
    pq = eq ? eq->exact() : Point_3(q.x(), q.y(), q.z());
    common_pa = common_ea ? common_ea->exact() : Point_3(common_a.x(), common_a.y(), common_a.z());
    common_pb = common_eb ? common_eb->exact() : Point_3(common_b.x(), common_b.y(), common_b.z());

    return triangle_do_overlap(pp, common_pa, common_pb, pq);
  }
  else
  {
    // all points are float point, we gonna use fast filter.
    return triangle_do_overlap(p, common_a, common_b, q);
  }
}

/// @brief test whether two triangle intersect. they have a common vertex.
/// this is a wrapper to CGAL/Fast filter.
bool triangle_do_intersect(
  CR_FP p, CR_FP q, CR_FP common_v, CR_FP a, CR_FP b,
  CP_EP ep, CP_EP eq, CP_EP common_ev, CP_EP ea, CP_EP eb)
{
  if (ep || eq || common_ev || ea || eb)
  {
    // if any one point is exact point, we gonna use CGAL's predicates.
    Point_3 pp, pq, common_pv, pa, pb;
    pp = ep ? ep->exact() : Point_3(p.x(), p.y(), p.z());
    pq = eq ? eq->exact() : Point_3(q.x(), q.y(), q.z());
    pa = ea ? ea->exact() : Point_3(a.x(), a.y(), a.z());
    pb = eb ? eb->exact() : Point_3(b.x(), b.y(), b.z());
    common_pv = common_ev ? common_ev->exact() : Point_3(common_v.x(), common_v.y(), common_v.z());

    Triangle_3 tri0(pp, pq, common_pv);
    Segment_3 seg0(pa, pb);
    Triangle_3 tri1(pa, pb, common_pv);
    Segment_3 seg1(pp, pq);
    return CGAL::do_intersect(tri0, seg0) || CGAL::do_intersect(tri1, seg1);
  }
  else
  {
    // all points are float point, we gonna use fast filter.
    return do_intersect(p, q, common_v, a, b) || do_intersect(a, b, common_v, p, q);
  }
}

}// namespace Geometry
}// namespace Cage