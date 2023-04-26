#pragma once
#include "Basic/Types.h"
#include <cmath>

namespace Cage
{
namespace Geometry
{
typedef const Vec3d& CR_Vec3d;

/**********************/
/* Geometry threshold */
/**********************/

struct GeomThr
{
  /// used in:
  /// * check_wrinkle
  /// * flip_ok
  static double area_de_thr;
  /// used in:
  /// * eliminate_almost_degeneration
  static double cos_de_thr;
  /// used in:
  /// * eliminate_almost_degeneration
  static double edge_length_ratio_de_thr;
  /// used in:
  /// * tangential relaxation
  static double relocate_length_thr;
};

/*****************/
/*     Area      */
/*****************/

/// calculate tri area
double triangle_area(CR_Vec3d a, CR_Vec3d b, CR_Vec3d c);

/*****************/
/*     Volume    */
/*****************/

/// calculate tet volume
double tet_volume(CR_Vec3d a, CR_Vec3d b, CR_Vec3d c, CR_Vec3d d);

/******************************/
/*     Coordinate system      */
/******************************/

void make_coordinate_system(const Vec3d& axis_z, Vec3d& axis_x, Vec3d& axis_y);

}// namespace Geometry
}// namespace Cage