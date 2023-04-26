#include "ExactPoint.h"

namespace Cage
{
namespace Geometry
{

ExactPoint::ExactPoint(double x, double y, double z)
{
  rp = Point_3(x, y, z);
  same = true;
}

ExactPoint::ExactPoint(const Vec3d& src)
{
  rp = Point_3(src.x(), src.y(), src.z());
  same = true;
}

ExactPoint::ExactPoint(const Point_3& src)
{
  rp = src;

  Vec3d fp = approx();
  same = ET(fp.x()) == rp.x() && ET(fp.y()) == rp.y() && ET(fp.z()) == rp.z();
}

ExactPoint::ExactPoint(const Vector_3& src)
{
  rp = Point_3(src.x(), src.y(), src.z());

  Vec3d fp = approx();
  same = ET(fp.x()) == rp.x() && ET(fp.y()) == rp.y() && ET(fp.z()) == rp.z();
}

Vector_3 ExactPoint::exactVec3()
{
  return Vector_3(rp.x(), rp.y(), rp.z());
}

Vec3d ExactPoint::approx() const
{
  return Vec3d(
    CGAL::to_double(rp.x()),
    CGAL::to_double(rp.y()),
    CGAL::to_double(rp.z()));
}

Point_3 ExactPoint::round()const
{
  return Point_3(
    CGAL::to_double(rp.x()),
    CGAL::to_double(rp.y()),
    CGAL::to_double(rp.z()));
}

void ExactPoint::rounded(Vec3d& fp)
{
  fp = approx();
  if (!same)
  {
    rp = Point_3(fp.x(), fp.y(), fp.z());
    same = true;
  }
}

void ExactPoint::fromImplicitPoint(genericPoint* gp)
{
  switch (gp->getType())
  {
  case Point_Type::EXPLICIT3D:
  {
    const explicitPoint3D* expp = static_cast<const explicitPoint3D*>(gp);
    rp = Point_3(expp->X(), expp->Y(), expp->Z());
    same = true;
  }
  break;
  case Point_Type::LPI:
  {
    const implicitPoint3D_LPI* ip = static_cast<const implicitPoint3D_LPI*>(gp);
    // segment
    Vector_3 p(ip->P().X(), ip->P().Y(), ip->P().Z());
    Vector_3 q(ip->Q().X(), ip->Q().Y(), ip->Q().Z());
    // plane
    Vector_3 r(ip->R().X(), ip->R().Y(), ip->R().Z());
    Vector_3 s(ip->S().X(), ip->S().Y(), ip->S().Z());
    Vector_3 t(ip->T().X(), ip->T().Y(), ip->T().Z());

    // refer to lambda3d_LPI_exact or "Indirect Predicates"
    // for compute procedure.
    Vector_3 pq = p - q;
    Vector_3 a23 = CGAL::cross_product(s - r, t - r);
    ET lambda_d = pq * a23;
    ET n = (p - r) * a23;
    Vector_3 lambda = p - pq * (n / lambda_d);
    rp = Point_3(lambda.x(), lambda.y(), lambda.z());
    Vec3d fp = approx();
    same = ET(fp.x()) == rp.x() && ET(fp.y()) == rp.y() && ET(fp.z()) == rp.z();
  }
  break;
  case Point_Type::TPI:
  {
    const implicitPoint3D_TPI* ip = static_cast<const implicitPoint3D_TPI*>(gp);
    // plane 1
    Vector_3 v1(ip->V1().X(), ip->V1().Y(), ip->V1().Z());
    Vector_3 v2(ip->V2().X(), ip->V2().Y(), ip->V2().Z());
    Vector_3 v3(ip->V3().X(), ip->V3().Y(), ip->V3().Z());
    // plane 2
    Vector_3 w1(ip->W1().X(), ip->W1().Y(), ip->W1().Z());
    Vector_3 w2(ip->W2().X(), ip->W2().Y(), ip->W2().Z());
    Vector_3 w3(ip->W3().X(), ip->W3().Y(), ip->W3().Z());
    // plane 3
    Vector_3 u1(ip->U1().X(), ip->U1().Y(), ip->U1().Z());
    Vector_3 u2(ip->U2().X(), ip->U2().Y(), ip->U2().Z());
    Vector_3 u3(ip->U3().X(), ip->U3().Y(), ip->U3().Z());

    // refer to "Exact and Efficient Polyhedral Envelope Containment Check"
    // for compute procedure.
    Vector_3 v32 = v3 - v2;
    Vector_3 v21 = v2 - v1;
    Vector_3 w32 = w3 - w2;
    Vector_3 w21 = w2 - w1;
    Vector_3 u32 = u3 - u2;
    Vector_3 u21 = u2 - u1;

    Vector_3 nv = CGAL::cross_product(v21, v32);
    Vector_3 nw = CGAL::cross_product(w21, w32);
    Vector_3 nu = CGAL::cross_product(u21, u32);

    ET pv = nv * v1, pw = nw * w1, pu = nu * u1;
    ET beta = nv * CGAL::cross_product(nw, nu);

    ET nx = Vector_3(pv, nv.y(), nv.z()) * CGAL::cross_product(
      Vector_3(pw, nw.y(), nw.z()), Vector_3(pu, nu.y(), nu.z()));
    ET ny = Vector_3(nv.x(), pv, nv.z()) * CGAL::cross_product(
      Vector_3(nw.x(), pw, nw.z()), Vector_3(nu.x(), pu, nu.z()));
    ET nz = Vector_3(nv.x(), nv.y(), pv) * CGAL::cross_product(
      Vector_3(nw.x(), nw.y(), pw), Vector_3(nu.x(), nu.y(), pu));
    rp = Point_3(nx / beta, ny / beta, nz / beta);
    Vec3d fp = approx();
    same = ET(fp.x()) == rp.x() && ET(fp.y()) == rp.y() && ET(fp.z()) == rp.z();
  }
  break;
  default:
    ASSERT(false, "wrong point type.");
  }
}

void minimize(Point_3& lhs, const Point_3& rhs)
{
  decltype(lhs.x()) x, y, z;
  x = lhs.x() < rhs.x() ? lhs.x() : rhs.x();
  y = lhs.y() < rhs.y() ? lhs.y() : rhs.y();
  z = lhs.z() < rhs.z() ? lhs.z() : rhs.z();
  lhs = Point_3(x, y, z);
}

void maximize(Point_3& lhs, const Point_3& rhs)
{
  decltype(lhs.x()) x, y, z;
  x = lhs.x() > rhs.x() ? lhs.x() : rhs.x();
  y = lhs.y() > rhs.y() ? lhs.y() : rhs.y();
  z = lhs.z() > rhs.z() ? lhs.z() : rhs.z();
  lhs = Point_3(x, y, z);
}

}// namespace Geometry
}// namespace Cage