#include "SurfaceMeshDefinition.h"
#include "Geometry/GeometryMath.h"

namespace Cage
{
namespace SurfaceMesh
{
/// @brief pre-calculate edge length for all edges
void pre_calculate_edge_length(SMeshT* mesh)
{
	for (EdgeHandle eh : mesh->edges())
	{
		mesh->data(eh).edge_length = mesh->calc_edge_length(eh);
	}
}

/// @brief pre-calculate face area for all faces. 
/// Assume edge length is calculated.
void pre_calculate_face_area(SMeshT* mesh)
{
	for (FaceHandle fh : mesh->faces())
	{
		mesh->data(fh).face_area = calc_face_area(mesh, fh);
	}
}

/// @brief Assume edge length exists, calculate three sector angle of triangle face.
/// @return halfedge: point to sector. double: cos(angle)
std::array<std::pair<HalfedgeHandle, double>, 3>
face_sector_angle(SMeshT* mesh, FaceHandle fh)
{
	auto fhe_it = mesh->fh_begin(fh);
	HalfedgeHandle ha = *fhe_it;fhe_it++;
	HalfedgeHandle hb = *fhe_it;fhe_it++;
	HalfedgeHandle hc = *fhe_it;

	double la = mesh->data(mesh->edge_handle(ha)).edge_length;
	double lb = mesh->data(mesh->edge_handle(hb)).edge_length;
	double lc = mesh->data(mesh->edge_handle(hc)).edge_length;

	double la2 = la * la, lb2 = lb * lb, lc2 = lc * lc;

	double sa = (la2 + lb2 - lc2) / (2. * la * lb);
	double sb = (lb2 + lc2 - la2) / (2. * lb * lc);
	double sc = (la2 + lc2 - lb2) / (2. * la * lc);

	std::array<std::pair<HalfedgeHandle, double>, 3> result;
	result[0] = { ha, sa }; result[1] = { hb, sb }; result[2] = { hc, sc };
	return result;
}

/// @brief Assume edge length exists, use Heron formula to calculate triangle area.
double calc_face_area(SMeshT* mesh, FaceHandle fh)
{
	auto fe_it = mesh->fe_begin(fh);

	double a = mesh->data(*fe_it).edge_length; fe_it++;
	double b = mesh->data(*fe_it).edge_length; fe_it++;
	double c = mesh->data(*fe_it).edge_length;

	double p = 0.5 * (a + b + c);
	double sqr_area = p * (p - a) * (p - b) * (p - c);
	return sqr_area < 0.0 ? 0. : std::sqrt(sqr_area);
}

/// @brief Assume face area exists, sum all of face area.
double calc_mesh_area(SMeshT* mesh)
{
	double total_area = 0.0;
	for (FaceHandle fh : mesh->faces())
	{
		total_area += mesh->data(fh).face_area;
	}
	return total_area;
}

/*****************/
/*     Graph     */
/*****************/

void graph_color_vertex(
	const SMeshT& mesh,
	std::vector<std::vector<int>>& phase)
{
	size_t v_size = mesh.n_vertices();
	std::vector<int> color(v_size, -1);
	std::vector<int> possible_color;   //represent every vertex's posibility of rendering, if the c's element is zero, that mean the posibility of the vertex is c is zero
	int ncolors = 0;
	for (size_t i = 0; i < v_size; i++)
	{
		std::fill(possible_color.begin(), possible_color.end(), 1);
		auto vj_h = mesh.vertex_handle((unsigned)i);
		for (auto itvv = mesh.cvv_begin(vj_h); itvv != mesh.cvv_end(vj_h); itvv++)
		{
			int c = color[itvv->idx()];
			if (c >= 0)
				possible_color[c] = 0;
		}
		int color_id = -1;
		for (auto j = 0; j < possible_color.size(); j++)
		{
			if (possible_color[j] != 0)
			{
				color_id = j;
				break;
			}
		}
		if (color_id < 0)
		{
			color_id = ncolors++;
			possible_color.resize(ncolors);
		}
		color[i] = color_id;
	}
	phase.resize(ncolors);
	for (int i = 0; i < v_size; i++)
		phase[color[i]].push_back(i);
}
}// namespace SurfaceMesh
}// namespace Cage