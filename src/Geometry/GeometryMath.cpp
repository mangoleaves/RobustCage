#include "GeometryMath.h"
#include "indirect_predicates.h"

namespace Cage
{
namespace Geometry
{

/**********************/
/* Geometry threshold */
/**********************/

double GeomThr::area_de_thr = 1e-8;
double GeomThr::cos_de_thr = 0.965;
double GeomThr::edge_length_ratio_de_thr = 0.2;
double GeomThr::relocate_length_thr = 1e-4;

/*****************/
/*     Area      */
/*****************/

double triangle_area(CR_Vec3d a, CR_Vec3d b, CR_Vec3d c)
{
	return 0.5 * ((a - b) % (c - b)).length();
}

/*****************/
/*     Volume    */
/*****************/

/// calculate tet volume
double tet_volume(CR_Vec3d a, CR_Vec3d b, CR_Vec3d c, CR_Vec3d d)
{
	// detect d is on the tri, or below the tri, or in the tri
	double temp[3][3];
	for (int i = 0; i < 3; ++i)
	{
		temp[0][i] = a[i] - d[i];
		temp[1][i] = b[i] - d[i];
		temp[2][i] = c[i] - d[i];
	}
	return std::abs(temp[0][0] * temp[1][1] * temp[2][2]
		+ temp[0][1] * temp[1][2] * temp[2][0]
		+ temp[0][2] * temp[1][0] * temp[2][1]
		- temp[0][2] * temp[1][1] * temp[2][0]
		- temp[0][1] * temp[1][0] * temp[2][2]
		- temp[0][0] * temp[1][2] * temp[2][1]);
}

/******************************/
/*     Coordinate system      */
/******************************/

void make_coordinate_system(const Vec3d& axis_z, Vec3d& axis_x, Vec3d& axis_y)
{
	Vec3d temp_vec;
	if (abs(axis_z.x()) < 0.9)
		temp_vec = Vec3d(1., 0., 0.);
	else if (abs(axis_z.y()) < 0.9)
		temp_vec = Vec3d(0., 1., 0.);
	else if (abs(axis_z.z()) < 0.9)
		temp_vec = Vec3d(0., 0., 1.);

	axis_x = axis_z.cross(temp_vec).normalized();
	axis_y = axis_z.cross(axis_x).normalized();
}
}// namespace Geometry
}// namespace Cage