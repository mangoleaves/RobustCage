#include "GeometryCheck.h"

namespace Cage
{
namespace CageSimp
{

static size_t samples_per_area = 10;
static size_t max_samples_per_face = 20;
static size_t min_samples_per_face = 10;
static size_t max_samples_per_edge = 10;
static size_t min_samples_per_edge = 5;
double cage_infinite_fp = DBL_MAX;

/// @brief calculate size of samples on faces.
/// size of samples is determined by area of the face and its one ring faces.
void get_nb_samples_on_faces(SMeshT* mesh, const std::vector<FaceHandle>& faces)
{
	if (calc_mesh_area(mesh) < GeomThr::area_de_thr) for (FaceHandle f : faces)	// a weird code format lol.
	{
		mesh->data(f).samples_num = min_samples_per_face;
	}
	else for (FaceHandle f : faces)
	{
		double sum_area = 0.0;
		size_t incident_size = 0;
		for (FaceHandle ff : mesh->ff_range(f))
		{
			incident_size++;
			sum_area += mesh->data(ff).face_area;
		}
		size_t nb_face_samples;
		if (sum_area < GeomThr::area_de_thr)
		{
			nb_face_samples = min_samples_per_face;
		}
		else
		{
			nb_face_samples = static_cast<size_t>(
				incident_size * samples_per_area *
				mesh->data(f).face_area / sum_area);

			nb_face_samples = std::min(nb_face_samples, max_samples_per_face);
			nb_face_samples = std::max(nb_face_samples, min_samples_per_face);
		}
		mesh->data(f).samples_num = nb_face_samples;
	}
}

/// @brief calculate size of samples on edges.
/// size of samples is determined by length of edge and area of two adjacent faces.
void get_nb_samples_on_edges(SMeshT* mesh, const std::vector<EdgeHandle>& edges)
{
	for (EdgeHandle e : edges)
	{
		HalfedgeHandle he = mesh->halfedge_handle(e, 0);
		HalfedgeHandle he_oppo = mesh->halfedge_handle(e, 1);

		double face_area = 0.0;
		size_t nb_edge_out_links = 0;
		if (!mesh->is_boundary(he))
		{
			FaceHandle f_he = mesh->face_handle(he);
			if (mesh->data(f_he).samples_num > 0)
			{
				face_area += mesh->data(f_he).face_area;
				nb_edge_out_links += mesh->data(f_he).samples_num;
			}
		}
		if (!mesh->is_boundary(he_oppo))
		{
			FaceHandle f_he_oppo = mesh->face_handle(he_oppo);
			if (mesh->data(f_he_oppo).samples_num > 0)
			{
				face_area += mesh->data(f_he_oppo).face_area;
				nb_edge_out_links += mesh->data(f_he_oppo).samples_num;
			}
		}

		if (nb_edge_out_links > 0)
		{
			double area_per_sample = face_area / nb_edge_out_links;
			double diameter = 2 * std::sqrt(area_per_sample / M_PI);

			size_t nb_samples = static_cast<size_t>(mesh->data(e).edge_length / diameter) / 3;

			nb_samples = std::min(nb_samples, max_samples_per_edge);
			nb_samples = std::max(nb_samples, min_samples_per_edge);

			mesh->data(e).samples_num = nb_samples;
		}
		else
		{
			mesh->data(e).samples_num = min_samples_per_edge;
		}
	}
}

void clear_nb_samples_on_faces(SMeshT* mesh, const std::vector<FaceHandle>& faces)
{
	for (FaceHandle f : faces)
		mesh->data(f).samples_num = 0;
}

void clear_nb_samples_on_edges(SMeshT* mesh, const std::vector<EdgeHandle>& edges)
{
	for (EdgeHandle e : edges)
		mesh->data(e).samples_num = 0;
}

/// @brief generate samples on face depending on size of samples.
std::vector<Vec3d> generate_samples_on_face(SMeshT* mesh, FaceHandle f)
{
	size_t nb_samples = mesh->data(f).samples_num;
	std::vector<Vec3d> inner_samples; inner_samples.reserve(nb_samples);

	HalfedgeHandle ha = mesh->halfedge_handle(f);
	HalfedgeHandle hb = mesh->next_halfedge_handle(ha);
	HalfedgeHandle hc = mesh->next_halfedge_handle(hb);
	Vec3d& b = mesh->point(mesh->to_vertex_handle(ha));
	Vec3d& c = mesh->point(mesh->to_vertex_handle(hb));
	Vec3d& a = mesh->point(mesh->to_vertex_handle(hc));

	if (nb_samples == 1)
	{
		inner_samples.push_back(mesh->calc_face_centroid(f));
	}
	else if (nb_samples == 2)
	{
		double ab = mesh->data(mesh->edge_handle(ha)).edge_length;
		double bc = mesh->data(mesh->edge_handle(hb)).edge_length;
		double ca = mesh->data(mesh->edge_handle(hc)).edge_length;

		size_t longest_edge_index;
		if (bc > ca)
			longest_edge_index = bc > ab ? 0 : 2;
		else
			longest_edge_index = ab < ca ? 1 : 2;

		Vec3d d;
		constexpr double _1_div_3 = 1. / 3.;
		switch (longest_edge_index)
		{
		case 0:
			d = (b + c) * 0.5;
			inner_samples.push_back((a + b + d) * _1_div_3);
			inner_samples.push_back((a + d + c) * _1_div_3);
			break;
		case 1:
			d = (a + c) * 0.5;
			inner_samples.push_back((b + c + d) * _1_div_3);
			inner_samples.push_back((b + d + a) * _1_div_3);
			break;
		case 2:
			d = (a + b) * 0.5;
			inner_samples.push_back((c + a + d) * _1_div_3);
			inner_samples.push_back((c + d + b) * _1_div_3);
			break;
		default:
			break;
		}
	}
	else
	{
		Vec3d ab = b - a, ac = c - a;
		while (inner_samples.size() < nb_samples)
		{
			double u = (double)rand() / (double)RAND_MAX;
			double v = (double)rand() / (double)RAND_MAX;
			if (u + v > 1.0)
			{
				u = 1.0 - u, v = 1.0 - v;
			}
			inner_samples.push_back(a + u * ab + v * ac);
		}
	}
	return inner_samples;
}

/// @brief generate samples on edge depending on size of samples.
std::vector<Vec3d> generate_samples_on_edge(SMeshT* mesh, EdgeHandle e)
{
	HalfedgeHandle he = mesh->halfedge_handle(e, 0);
	VertexHandle from_v = mesh->from_vertex_handle(he), to_v = mesh->to_vertex_handle(he);
	const Vec3d& from_p = mesh->point(from_v), & to_p = mesh->point(to_v);

	Vec3d cur_p = from_p;
	Vec3d step = (to_p - from_p) / (double)(mesh->data(e).samples_num + 1u);
	std::vector<Vec3d> samples;		samples.reserve(mesh->data(e).samples_num);
	for (size_t i = 0; i < mesh->data(e).samples_num; ++i)
	{
		cur_p += step;
		samples.push_back(cur_p);
	}
	return samples;
}

/// @brief out mesh ---link---> in mesh.
/// @param [in] out_mesh links are out-links of out_mesh.
/// @param [in] in_mesh links are in-links of in_mesh.
/// @param [in] in_tree dynamic AABB tree on in_mesh.
/// @param [in] faces faces on out_mesh.
void generate_face_out_links(SMeshT* out_mesh, SMeshT* in_mesh, FaceTree* in_tree, const std::vector<FaceHandle>& faces)
{
	for (FaceHandle f : faces)
	{
		if (out_mesh->data(f).samples_num == 0)
			continue;

		auto samples = generate_samples_on_face(out_mesh, f);
		in_tree->set_hint(*samples.begin());
		for (auto& p : samples)
		{
			auto cp = in_tree->closest_point_and_face_handle(p);
			in_mesh->data(cp.second).face_in_links.push_back({ p, cp.first });
		}
		in_tree->unset_hint();
	}
}

/// @brief out mesh ---link---> in mesh.
void generate_edge_out_links(SMeshT* out_mesh, SMeshT* in_mesh, FaceTree* in_tree, const std::vector<EdgeHandle>& edges)
{
	for (EdgeHandle e : edges)
	{
		if (out_mesh->data(e).samples_num == 0)
			continue;

		auto samples = generate_samples_on_edge(out_mesh, e);
		in_tree->set_hint(*samples.begin());
		for (auto& p : samples)
		{
			auto cp = in_tree->closest_point_and_face_handle(p);
			in_mesh->data(cp.second).face_in_links.push_back({ p, cp.first });
		}
		in_tree->unset_hint();
	}
}

/// @brief out mesh ---link---> in mesh.
void generate_vertex_out_links(SMeshT* out_mesh, SMeshT* in_mesh, FaceTree* in_tree, VertexHandle v)
{
	Vec3d& p = out_mesh->point(v);
	auto cp = in_tree->closest_point_and_face_handle(p);
	in_mesh->data(cp.second).face_in_links.push_back({ p, cp.first });
}

double calc_face_in_error(SMeshT* in_mesh, const std::vector<FaceHandle>& faces)
{
	double distance = 0.0;
	for (FaceHandle fh : faces)
	{
		double face_distance = 0.0;
		for (const Link& link : in_mesh->data(fh).face_in_links)
		{
			double d = (link.first - link.second).length();
			if (d > face_distance)
				face_distance = d;
		}
		in_mesh->data(fh).in_error = face_distance;
		if (face_distance > distance)
			distance = face_distance;
	}
	return distance;
}

#ifdef USE_TREE_SEARCH
/// @brief out_mesh is always a local mesh.
/// So, hint is set before calling this function.
double calc_face_out_error(
	SMeshT* out_mesh, SMeshT* in_mesh,
	DFaceTree* in_tree, const std::vector<FaceHandle>& faces,
	double& threshold)
{
	double distance = 0.0;
	for (FaceHandle f : faces)
	{
		if (out_mesh->data(f).samples_num == 0)
			continue;

		double face_distance = 0.0;
		auto samples = generate_samples_on_face(out_mesh, f);
		for (auto& p : samples)
		{
			auto cd = in_tree->closest_distance_and_face_handle(p);
			if (cd.first >= threshold)
				return DBL_MAX;
			if (cd.first > face_distance)
				face_distance = cd.first;
		}
		out_mesh->data(f).out_error = face_distance;
		if (face_distance > distance)
			distance = face_distance;
	}
	return distance;
}

/// @brief out_mesh is always a local mesh.
/// So, hint is set before calling this function.
double calc_edge_out_error(
	SMeshT* out_mesh, SMeshT* in_mesh,
	DFaceTree* in_tree, const std::vector<EdgeHandle>& edges,
	double& threshold)
{
	double distance = 0.0;
	for (EdgeHandle e : edges)
	{
		if (out_mesh->data(e).samples_num == 0)
			continue;

		double edge_distance = 0.0;
		auto samples = generate_samples_on_edge(out_mesh, e);
		for (auto& p : samples)
		{
			auto cd = in_tree->closest_distance_and_face_handle(p);
			if (cd.first >= threshold)
				return DBL_MAX;
			if (cd.first > edge_distance)
				edge_distance = cd.first;
		}
		out_mesh->data(e).out_error = edge_distance;
		if (edge_distance > distance)
			distance = edge_distance;
	}
	return distance;
}

/// @brief out_mesh is always a local mesh.
/// So, hint is set before calling this function.
double calc_vertex_out_error(
	SMeshT* out_mesh, SMeshT* in_mesh, DFaceTree* in_tree, VertexHandle v)
{
	double distance = 0.0;
	Vec3d& p = out_mesh->point(v);
	auto cd = in_tree->closest_distance_and_face_handle(p);
	out_mesh->data(v).out_error = cd.first;
	if (cd.first > distance)
		distance = cd.first;
	return distance;
}
#else

double calc_face_out_error(
	SMeshT* out_mesh, SMeshT* in_mesh,
	FaceGrid* in_grid, const std::vector<FaceHandle>& faces,
	double& threshold)
{
	double distance = 0.0;
	for (FaceHandle f : faces)
	{
		if (out_mesh->data(f).samples_num == 0)
			continue;

		double face_distance = 0.0;
		auto samples = generate_samples_on_face(out_mesh, f);
		for (auto& p : samples)
		{
			auto cd = in_grid->closest_point(p);
			if (cd.first.second >= threshold)
				return DBL_MAX;
			if (cd.first.second > face_distance)
				face_distance = cd.first.second;
		}
		out_mesh->data(f).out_error = face_distance;
		if (face_distance > distance)
			distance = face_distance;
	}
	return distance;
}

double calc_edge_out_error(
	SMeshT* out_mesh, SMeshT* in_mesh,
	FaceGrid* in_grid, const std::vector<EdgeHandle>& edges,
	double& threshold)
{
	double distance = 0.0;
	for (EdgeHandle e : edges)
	{
		if (out_mesh->data(e).samples_num == 0)
			continue;

		double edge_distance = 0.0;
		auto samples = generate_samples_on_edge(out_mesh, e);
		for (auto& p : samples)
		{
			auto cd = in_grid->closest_point(p);
			if (cd.first.second >= threshold)
				return DBL_MAX;
			if (cd.first.second > edge_distance)
				edge_distance = cd.first.second;
		}
		out_mesh->data(e).out_error = edge_distance;
		if (edge_distance > distance)
			distance = edge_distance;
	}
	return distance;
}

double calc_vertex_out_error(
	SMeshT* out_mesh, SMeshT* in_mesh, FaceGrid* in_grid, VertexHandle v)
{
	double distance = 0.0;
	Vec3d& p = out_mesh->point(v);
	auto cd = in_grid->closest_point(p);
	out_mesh->data(v).out_error = cd.first.second;
	if (cd.first.second > distance)
		distance = cd.first.second;
	return distance;
}
#endif

void calc_out_surround_error(SMeshT* out_mesh, VertexHandle v)
{
	double distance = out_mesh->data(v).out_error;
	for (FaceHandle vf : out_mesh->vf_range(v))
	{
		if (distance < out_mesh->data(vf).in_error)
			distance = out_mesh->data(vf).in_error;
		if (distance < out_mesh->data(vf).out_error)
			distance = out_mesh->data(vf).out_error;
	}
	for (EdgeHandle ve : out_mesh->ve_range(v))
	{
		if (distance < out_mesh->data(ve).out_error)
			distance = out_mesh->data(ve).out_error;
	}
	for (VertexHandle vv : out_mesh->vv_range(v))
	{
		if (distance < out_mesh->data(vv).out_error)
			distance = out_mesh->data(vv).out_error;
	}
	out_mesh->data(v).out_surround_error = distance;
}

void generate_out_links(SMeshT* out_mesh, SMeshT* in_mesh, FaceTree* in_tree)
{
	std::vector<OpenMesh::FaceHandle> faces; faces.reserve(out_mesh->n_faces());
	std::vector<OpenMesh::EdgeHandle> edges; edges.reserve(out_mesh->n_edges());
	for (FaceHandle fh : out_mesh->faces()) faces.push_back(fh);
	for (EdgeHandle eh : out_mesh->edges()) edges.push_back(eh);

	get_nb_samples_on_faces(out_mesh, faces);
	get_nb_samples_on_edges(out_mesh, edges);
	generate_face_out_links(out_mesh, in_mesh, in_tree, faces);
	generate_edge_out_links(out_mesh, in_mesh, in_tree, edges);
	for (VertexHandle v : out_mesh->vertices())
		generate_vertex_out_links(out_mesh, in_mesh, in_tree, v);
}

#ifdef USE_TREE_SEARCH
double calc_out_error(SMeshT* out_mesh, SMeshT* in_mesh, DFaceTree* in_tree, double& threshold)
{
	std::vector<OpenMesh::FaceHandle> faces; faces.reserve(out_mesh->n_faces());
	std::vector<OpenMesh::EdgeHandle> edges; edges.reserve(out_mesh->n_edges());
	for (FaceHandle fh : out_mesh->faces()) /*if (!out_mesh->status(fh).deleted())*/ faces.push_back(fh);
	for (EdgeHandle eh : out_mesh->edges()) /*if (!out_mesh->status(eh).deleted())*/ edges.push_back(eh);

	get_nb_samples_on_faces(out_mesh, faces);
	get_nb_samples_on_edges(out_mesh, edges);
	double distance = 0.0;
	distance = std::max(distance, calc_face_out_error(out_mesh, in_mesh, in_tree, faces, threshold));
	if (distance >= threshold)	return DBL_MAX;
	distance = std::max(distance, calc_edge_out_error(out_mesh, in_mesh, in_tree, edges, threshold));
	if (distance >= threshold)	return DBL_MAX;
	for (VertexHandle v : out_mesh->vertices())
	{
		/*if (!out_mesh->status(v).deleted())*/
		{
			distance = std::max(distance, calc_vertex_out_error(out_mesh, in_mesh, in_tree, v));
			if (distance >= threshold)	return DBL_MAX;
		}
	}
	for (VertexHandle v : out_mesh->vertices())
		calc_out_surround_error(out_mesh, v);
	return distance;
}
#else
double calc_out_error(SMeshT* out_mesh, SMeshT* in_mesh, FaceGrid* in_grid, double& threshold)
{
	std::vector<OpenMesh::FaceHandle> faces; faces.reserve(out_mesh->n_faces());
	std::vector<OpenMesh::EdgeHandle> edges; edges.reserve(out_mesh->n_edges());
	for (FaceHandle fh : out_mesh->faces()) /*if (!out_mesh->status(fh).deleted())*/ faces.push_back(fh);
	for (EdgeHandle eh : out_mesh->edges()) /*if (!out_mesh->status(eh).deleted())*/ edges.push_back(eh);

	get_nb_samples_on_faces(out_mesh, faces);
	get_nb_samples_on_edges(out_mesh, edges);
	double distance = 0.0;
	distance = std::max(distance, calc_face_out_error(out_mesh, in_mesh, in_grid, faces, threshold));
	if (distance >= threshold)	return DBL_MAX;
	distance = std::max(distance, calc_edge_out_error(out_mesh, in_mesh, in_grid, edges, threshold));
	if (distance >= threshold)	return DBL_MAX;
	for (VertexHandle v : out_mesh->vertices())
	{
		/*if (!out_mesh->status(v).deleted())*/
		{
			distance = std::max(distance, calc_vertex_out_error(out_mesh, in_mesh, in_grid, v));
			if (distance >= threshold)	return DBL_MAX;
		}
	}
	for (VertexHandle v : out_mesh->vertices())
		calc_out_surround_error(out_mesh, v);
	return distance;
}
#endif

/// @brief backup local in-links on remeshing mesh.
LinkVector backup_local_in_links(SMeshT* mesh, const std::vector<FaceHandle>& faces)
{
	// count links size
	size_t links_size = 0;
	for (FaceHandle f : faces)
		links_size += mesh->data(f).face_in_links.size();

	LinkVector faces_in_links;
	faces_in_links.reserve(links_size);

	for (FaceHandle f : faces) if (!mesh->data(f).face_in_links.empty())
	{
		faces_in_links.insert(faces_in_links.end(),
			mesh->data(f).face_in_links.begin(), mesh->data(f).face_in_links.end());
		mesh->data(f).face_in_links.clear();
	}
	return faces_in_links;
}

size_t count_face_out_links(SMeshT* mesh)
{
	return std::accumulate(
		mesh->faces().begin(), mesh->faces().end(),
		(size_t)0,
		[&](const size_t& count, const FaceHandle& f) { return count + mesh->data(f).samples_num; }
	);
}

size_t count_edge_out_links(SMeshT* mesh)
{
	return std::accumulate(
		mesh->edges().begin(), mesh->edges().end(),
		(size_t)0,
		[&](const size_t& count, const EdgeHandle& e) { return count + mesh->data(e).samples_num; }
	);
}

void clear_links(SMeshT* mesh)
{
	std::for_each(mesh->faces().begin(), mesh->faces().end(),
		[&](const FaceHandle& f) {mesh->data(f).face_in_links.clear();});
}

/********************/
/*     Predicate    */
/********************/

/// @brief given triangle abc and abp, query whether p and c at same side of ref-plane.
/// ref-plane is formed by a_b_(abc.normal + a).
/// @return 1: same side.  0: on the plane.  -1: different side.
int same_side_wrapped(CR_Vec3d a, CR_Vec3d b, CR_Vec3d c, CR_Vec3d p)
{
	return same_side(
		a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], p[0], p[1], p[2]);
}

/*******************/
/*    Constraint   */
/*******************/

/// @brief check all triangles formed by halfedges and new point
/// will cause wrinkle with respect to original point.
/// @param [in] mesh triangle mesh
/// @param [in] halfedges form a circle(interior) or semi-circle(boundary).
/// @param [in] new_point new position of original center point located on center of halfedges.
/// @return true if new position will cause wrinkle.
bool check_wrinkle(
	SMeshT* mesh,
	const std::vector<HalfedgeHandle>& halfedges,
	CR_Vec3d new_point)
{
	for (HalfedgeHandle h : halfedges)
	{
		if (mesh->data(mesh->edge_handle(h)).edge_length == 0.0)
			continue;

		auto& start_point = mesh->point(mesh->from_vertex_handle(h));
		auto& end_point = mesh->point(mesh->to_vertex_handle(h));
		auto& old_point = mesh->point(mesh->opposite_vh(h));
		/*
		auto start_end = end_point - start_point, start_new = new_point - start_point;
		if (0.5 * (start_end % start_new).length() < GeomThr::area_de_thr)
			return true;
		double radian = get_radian(end_point, new_point, start_point);
		if (radian < GeomThr::radian_de_thr
			|| M_PI - radian < GeomThr::radian_de_thr)
			return true;
		*/
		if (same_side_wrapped(start_point, end_point, old_point, new_point) <= 0)
			return true;
	}
	return false;
}

/// @brief check all triangles formed by halfedges and new point
/// will cause long edge with respect to original point.
/// @param [in] mesh triangle mesh
/// @param [in] halfedges form a circle(interior) or semi-circle(boundary).
/// @param [in] new_point new position of original center point located on center of halfedges.
/// @param [in] new_target_length of new_point.
/// @return true if new position will cause long edge.
bool check_long_edge(
	SMeshT* mesh,
	const std::vector<HalfedgeHandle>& halfedges,
	CR_Vec3d new_point,
	double new_target_length
)
{
	for (HalfedgeHandle h : halfedges)
	{
		VertexHandle vv = mesh->to_vertex_handle(h);
		auto& pp = mesh->point(vv);
		double length = (new_point - pp).sqrnorm();
		double high_edge_len = 4.0 / 3.0 * std::min(mesh->data(vv).target_length, new_target_length);
		if (length > high_edge_len)
			return true;
	}
	return false;
}

/// @brief check all triangles formed by halfedges and new point
/// will be degenerate.
/// @param [in] mesh triangle mesh
/// @param [in] halfedges form a circle(interior) or semi-circle(boundary).
/// @param [in] new_point new position of original center point located on center of halfedges.
/// @param [in] new_ep new exact point, pass nullptr if it is float point.
bool check_degenerate(
	SMeshT* mesh,
	const std::vector<HalfedgeHandle>& halfedges,
	CR_Vec3d new_point,
	const ExactPoint* new_ep)
{
	for (HalfedgeHandle h : halfedges)
	{
		// (1) find other two points to construct new triangle
		Vec3d& from_point = mesh->point(mesh->from_vertex_handle(h));
		Vec3d& to_point = mesh->point(mesh->to_vertex_handle(h));
		const ExactPoint* from_ep = mesh->data(mesh->from_vertex_handle(h)).ep.get();
		const ExactPoint* to_ep = mesh->data(mesh->to_vertex_handle(h)).ep.get();
		// (2) check degenerate
		if (are_points_colinear(from_point, to_point, new_point, from_ep, to_ep, new_ep))
		{
			return true;
		}
	}
	return false;
}

/// @brief check all triangles formed by halfedges and new point
/// will intersect with original mesh or self-intersect.
/// @param [in] mesh triangle mesh
/// @param [in] halfedges form a circle(interior) or semi-circle(boundary).
/// @param [in] new_point new position of original center point located on center of halfedges.
/// @param [in] new_ep new exact point, pass nullptr if it is float point.
/// @param [in] intersect_tree tree for intersection query.
/// @param [in] intersect_grid grid for intersection query.
/// @param [in] self_intersect_tree tree for self-intersection query.
/// @param [in] ignored_faces ignored faces while find self-intersection.
/// @return true if new position will cause intersection.
bool check_intersection(
	SMeshT* mesh,
	const std::vector<HalfedgeHandle>& halfedges,
	CR_Vec3d new_point,
	const ExactPoint* new_ep,
	DFaceTree* intersect_tree,
	FaceGrid* intersect_grid,
	LightDFaceTree* self_intersect_tree,
	const std::set<FaceHandle>& ignored_faces,
	bool check_selfinter,
	bool check_inter
)
{
	if (check_inter)
	{
		// 1. check intersect with original mesh
		for (HalfedgeHandle h : halfedges)
		{
			// (1) find other two points to construct new triangle
			Vec3d& from_point = mesh->point(mesh->from_vertex_handle(h));
			Vec3d& to_point = mesh->point(mesh->to_vertex_handle(h));
			const ExactPoint* from_ep = mesh->data(mesh->from_vertex_handle(h)).ep.get();
			const ExactPoint* to_ep = mesh->data(mesh->to_vertex_handle(h)).ep.get();
			// (2) query intersection
			if (intersect_tree->do_intersect(from_point, to_point, new_point, from_ep, to_ep, new_ep))
			{
				return true;
			}
		}
	}

	if (!check_selfinter)
		return false;

	// 2. check overlap triangle with a common edge
	for (HalfedgeHandle h : halfedges)
	{
		// check two faces adjacent to halfedge h
		Vec3d& start_point = mesh->point(mesh->from_vertex_handle(h));
		Vec3d& end_point = mesh->point(mesh->to_vertex_handle(h));
		ExactPoint* start_ep = mesh->data(mesh->from_vertex_handle(h)).ep.get();
		ExactPoint* end_ep = mesh->data(mesh->to_vertex_handle(h)).ep.get();
		HalfedgeHandle opp_h = mesh->opposite_halfedge_handle(h);
		if (mesh->is_boundary(opp_h))
			continue;
		Vec3d& opp_point = mesh->point(mesh->opposite_vh(opp_h));
		ExactPoint* opp_ep = mesh->data(mesh->opposite_vh(opp_h)).ep.get();
		if (Geometry::triangle_do_overlap(
			new_point, start_point, end_point, opp_point, new_ep, start_ep, end_ep, opp_ep))
		{
			return true;
		}
	}

	// 3. check self-intersect triangle with a common vertex
	for (HalfedgeHandle h : halfedges)
	{
		FaceHandle fh = mesh->face_handle(h);
		// (1) extend face to one ring faces
		const std::set<FaceHandle>& one_ring_faces = mesh->data(fh).one_ring_faces;
		// (2) one_ring_faces - ignored_faces
		std::vector<FaceHandle> faces_with_common_vertex(one_ring_faces.size());
		auto it = std::set_difference(
			one_ring_faces.begin(), one_ring_faces.end(),
			ignored_faces.begin(), ignored_faces.end(),
			faces_with_common_vertex.begin());
		faces_with_common_vertex.resize(it - faces_with_common_vertex.begin());
		// (3) faces_with_common_vertex - face_opposite_halfedge
		if (!mesh->is_boundary(mesh->opposite_halfedge_handle(h)))
		{
			it = std::find(faces_with_common_vertex.begin(), faces_with_common_vertex.end(),
				mesh->opposite_face_handle(h));
			if (it != faces_with_common_vertex.end())
				faces_with_common_vertex.erase(it);
		}
		// (4) find common-vertex intersection
		VertexHandle from_v = mesh->from_vertex_handle(h);
		VertexHandle to_v = mesh->to_vertex_handle(h);
		Vec3d& from_p = mesh->point(from_v);
		Vec3d& to_p = mesh->point(to_v);
		ExactPoint* from_ep = mesh->data(from_v).ep.get();
		ExactPoint* to_ep = mesh->data(to_v).ep.get();

		for (HalfedgeHandle voh : mesh->voh_range(from_v))
		{
			if (std::binary_search(faces_with_common_vertex.begin(), faces_with_common_vertex.end(),
				mesh->face_handle(voh)))
			{
				HalfedgeHandle next_voh = mesh->next_halfedge_handle(voh);
				Vec3d& opp_from_p = mesh->point(mesh->from_vertex_handle(next_voh));
				Vec3d& opp_to_p = mesh->point(mesh->to_vertex_handle(next_voh));
				ExactPoint* opp_from_ep = mesh->data(mesh->from_vertex_handle(next_voh)).ep.get();
				ExactPoint* opp_to_ep = mesh->data(mesh->to_vertex_handle(next_voh)).ep.get();

				if (Geometry::triangle_do_intersect(
					to_p, new_point, from_p, opp_from_p, opp_to_p,
					to_ep, new_ep, from_ep, opp_from_ep, opp_to_ep))
				{
					return true;
				}
			}
		}
		for (HalfedgeHandle voh : mesh->voh_range(to_v))
		{
			if (std::binary_search(faces_with_common_vertex.begin(), faces_with_common_vertex.end(),
				mesh->face_handle(voh)))
			{
				HalfedgeHandle next_voh = mesh->next_halfedge_handle(voh);
				Vec3d& opp_from_p = mesh->point(mesh->from_vertex_handle(next_voh));
				Vec3d& opp_to_p = mesh->point(mesh->to_vertex_handle(next_voh));
				ExactPoint* opp_from_ep = mesh->data(mesh->from_vertex_handle(next_voh)).ep.get();
				ExactPoint* opp_to_ep = mesh->data(mesh->to_vertex_handle(next_voh)).ep.get();

				if (Geometry::triangle_do_intersect(
					new_point, from_p, to_p, opp_from_p, opp_to_p,
					new_ep, from_ep, to_ep, opp_from_ep, opp_to_ep))
				{
					return true;
				}
			}
		}
	}

	// 4. check self-intersection
	for (HalfedgeHandle h : halfedges)
	{
		FaceHandle fh = mesh->face_handle(h);
		// (1) extend ignored faces to adjacent faces of vfh
		const std::set<FaceHandle>& one_ring_faces = mesh->data(fh).one_ring_faces;
		// (2) union two sets to get all ignored faces.
		std::vector<FaceHandle> all_ignored_faces(one_ring_faces.size() + ignored_faces.size());
		auto it = std::set_union(one_ring_faces.begin(), one_ring_faces.end(), ignored_faces.begin(), ignored_faces.end(), all_ignored_faces.begin());
		all_ignored_faces.resize(it - all_ignored_faces.begin());
		// (3) find other two points to construct new triangle
		auto& from_point = mesh->point(mesh->from_vertex_handle(h));
		const auto* from_ep = mesh->data(mesh->from_vertex_handle(h)).ep.get();
		auto& to_point = mesh->point(mesh->to_vertex_handle(h));
		const auto* to_ep = mesh->data(mesh->to_vertex_handle(h)).ep.get();
		// (4) query intersection
		std::vector<int> _all_ignored_faces(all_ignored_faces.size());
		std::transform(all_ignored_faces.begin(), all_ignored_faces.end(),
			_all_ignored_faces.begin(), [&](FaceHandle fh) {return fh.idx();});
		if (self_intersect_tree->do_intersect(
			from_point, to_point, new_point,
			from_ep, to_ep, new_ep, std::move(_all_ignored_faces)))
		{
			return true;
		}
	}
	return false;
}

}// namespace CageSimp
}// namespace Cage