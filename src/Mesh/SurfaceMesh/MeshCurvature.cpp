#include "MeshCurvature.h"
// Eigen
#include "Dense"

namespace Cage
{
namespace SurfaceMesh
{

/******************************/
/*   Principal Curvature      */
/******************************/

PrincipalCurvature::PrincipalCurvature(SMeshT* _mesh)
	:mesh(_mesh)
{}

void PrincipalCurvature::compute_principal_curvature()
{
	size_t nv = mesh->n_vertices();

	//compute vertex normal
	vertex_normal.clear(); vertex_normal.resize(nv, Vec3d(0., 0., 0.));
	for (SMeshT::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it)
	{
		SMeshT::FaceVertexIter fv_it = mesh->fv_iter(*f_it);
		const Vec3d& v0 = mesh->point(*fv_it); size_t vid0 = fv_it->idx(); fv_it++;
		const Vec3d& v1 = mesh->point(*fv_it); size_t vid1 = fv_it->idx(); fv_it++;
		const Vec3d& v2 = mesh->point(*fv_it); size_t vid2 = fv_it->idx();

		Vec3d a = v0 - v1;
		Vec3d b = v1 - v2;
		Vec3d c = v2 - v0;
		double l2a = a.sqrnorm(); double l2b = b.sqrnorm(); double l2c = c.sqrnorm();
		if (!l2a || !l2b || !l2c)
			continue;

		Vec3d face_normal = (a % b);	// weighted with face area
		vertex_normal[vid0] += face_normal * (1.0 / (l2a * l2c));	// why?
		vertex_normal[vid1] += face_normal * (1.0 / (l2b * l2a));
		vertex_normal[vid2] += face_normal * (1.0 / (l2c * l2b));
	}

	principal_K1.clear(); principal_K1.resize(nv, 0.0);
	principal_K2.clear(); principal_K2.resize(nv, 0.0);
	principal_dir1.clear(); principal_dir1.resize(nv, Vec3d(0, 0, 0));
	principal_dir2.clear(); principal_dir2.resize(nv, Vec3d(0, 0, 0));
	std::vector<double> k12(nv, 0.0);

	for (FaceHandle fh : mesh->faces())
	{
		for (HalfedgeHandle fhe : mesh->fh_range(fh))
		{
			VertexHandle from_v = mesh->from_vertex_handle(fhe);
			VertexHandle to_v = mesh->to_vertex_handle(fhe);

			principal_dir1[from_v.idx()] = mesh->point(to_v) - mesh->point(from_v);
		}
	}

	for (VertexHandle v : mesh->vertices())
	{
		int vid = v.idx();
		vertex_normal[vid].normalize();
		principal_dir1[vid] = (principal_dir1[vid] % vertex_normal[vid]).normalized();
		principal_dir2[vid] = (vertex_normal[vid] % principal_dir1[vid]).normalized();
	}

	compute_point_area();

	std::array<Vec3d, 3> ev;
	Vec3d t; Vec3d n; Vec3d b;
	for (SMeshT::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it)
	{
		// Edges
		SMeshT::HalfedgeHandle heh = *mesh->fh_iter(*f_it);
		SMeshT::HalfedgeHandle tmp_heh = mesh->next_halfedge_handle(heh);
		ev[0] = mesh->point(mesh->to_vertex_handle(tmp_heh)) - mesh->point(mesh->from_vertex_handle(tmp_heh));
		tmp_heh = mesh->prev_halfedge_handle(heh);
		ev[1] = mesh->point(mesh->to_vertex_handle(tmp_heh)) - mesh->point(mesh->from_vertex_handle(tmp_heh));
		ev[2] = mesh->point(mesh->to_vertex_handle(heh)) - mesh->point(mesh->from_vertex_handle(heh));

		// N-T-B coordinate system per face
		t = ev[0]; t.normalize();
		n = (ev[0] % ev[1]);
		b = (n % t); b.normalize();

		// Estimate curvature based on variation of normals along edges
		tmp_heh = mesh->next_halfedge_handle(heh);
		Eigen::Vector3d m(0.0, 0.0, 0.0);
		Eigen::Vector3d x;
		Eigen::Matrix3d w; w.setZero();
		for (int j = 0; j < 3; ++j)
		{
			double u = (ev[j] | t);
			double v = (ev[j] | b);
			w(0, 0) += u * u;
			w(0, 1) += u * v;
			//w[1][1] += v*v + u*u; 
			//w[1][2] += u*v; 
			w(2, 2) += v * v;
			/*Vec3d dn = mesh->normal(mesh->to_vertex_handle(temp_heh))
			- mesh->normal(mesh->from_vertex_handle(temp_heh));*/
			Vec3d dn = vertex_normal[mesh->to_vertex_handle(tmp_heh).idx()]
				- vertex_normal[mesh->from_vertex_handle(tmp_heh).idx()];
			double dnu = (dn | t);
			double dnv = (dn | b);
			m(0) += dnu * u;
			m(1) += dnu * v + dnv * u;
			m(2) += dnv * v;
			tmp_heh = mesh->next_halfedge_handle(tmp_heh);
		}
		w(1, 1) = w(0, 0) + w(2, 2);
		w(1, 2) = w(0, 1);
		w(2, 1) = w(1, 2);
		w(1, 0) = w(0, 1);
		//std::cout << w;
		x = w.fullPivHouseholderQr().solve(m);

		tmp_heh = heh;
		for (int j = 0; j < 3; ++j)
		{
			int vertex_id = mesh->from_vertex_handle(tmp_heh).idx();
			double c1, c12, c2;
			proj_curv(t, b, x(0), x(1), x(2),
				principal_dir1[vertex_id], principal_dir2[vertex_id], c1, c12, c2);
			double wt = cornerArea[(*f_it).idx()][vertex_id] / pointArea[vertex_id];
			principal_K1[vertex_id] += wt * c1;
			k12[vertex_id] += wt * c12;
			principal_K2[vertex_id] += wt * c2;
			tmp_heh = mesh->next_halfedge_handle(tmp_heh);
		}
	}

	for (auto v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); ++v_it)
	{
		int vertex_id = (*v_it).idx();
		/*diagonalize_curv(principal_dir1[vertex_id], principal_dir2[vertex_id],
		k1[vertex_id], k12[vertex_id], k2[vertex_id],
		mesh->normal(v_it), principal_dir1[vertex_id], principal_dir2[vertex_id],
		k1[vertex_id], k2[vertex_id]);*/
		diagonalize_curv(principal_dir1[vertex_id], principal_dir2[vertex_id],
			principal_K1[vertex_id], k12[vertex_id], principal_K2[vertex_id],
			vertex_normal[vertex_id], principal_dir1[vertex_id], principal_dir2[vertex_id],
			principal_K1[vertex_id], principal_K2[vertex_id]);
	}
}

void PrincipalCurvature::compute_point_area()
{
	pointArea.resize(mesh->n_vertices(), 0.0);
	cornerArea.resize(mesh->n_faces());

	std::array<SMeshT::Point, 3> e;
	std::array<int, 3> v;

	//#pragma omp parallel for
	for (SMeshT::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); ++f_it)
	{
		int face_index = (*f_it).idx();
		SMeshT::FaceHalfedgeIter fhe_it = mesh->fh_iter(*f_it);

		e[0] = mesh->point(mesh->to_vertex_handle(*fhe_it)) - mesh->point(mesh->from_vertex_handle(*fhe_it));
		v[2] = mesh->to_vertex_handle(*fhe_it).idx();
		++fhe_it;

		e[1] = mesh->point(mesh->to_vertex_handle(*fhe_it)) - mesh->point(mesh->from_vertex_handle(*fhe_it));
		v[0] = mesh->to_vertex_handle(*fhe_it).idx();
		++fhe_it;

		e[2] = mesh->point(mesh->to_vertex_handle(*fhe_it)) - mesh->point(mesh->from_vertex_handle(*fhe_it));
		v[1] = mesh->to_vertex_handle(*fhe_it).idx();

		double area = 0.5f * (e[0] % e[1]).norm();
		double l2[3] = { e[0].sqrnorm(), e[1].sqrnorm(), e[2].sqrnorm() };
		double ew[3] = { l2[0] * (l2[1] + l2[2] - l2[0]),
			l2[1] * (l2[2] + l2[0] - l2[1]),
			l2[2] * (l2[0] + l2[1] - l2[2]) };

		if (ew[0] <= 0.0)
		{
			cornerArea[face_index][v[1]] = -0.25 * l2[2] * area / (e[0] | e[2]);
			cornerArea[face_index][v[2]] = -0.25 * l2[1] * area / (e[0] | e[1]);
			cornerArea[face_index][v[0]] = area - cornerArea[face_index][v[1]] - cornerArea[face_index][v[2]];
		}
		else if (ew[1] <= 0.0)
		{
			cornerArea[face_index][v[2]] = -0.25 * l2[0] * area / (e[1] | e[0]);
			cornerArea[face_index][v[0]] = -0.25 * l2[2] * area / (e[1] | e[2]);
			cornerArea[face_index][v[1]] = area - cornerArea[face_index][v[2]] - cornerArea[face_index][v[0]];
		}
		else if (ew[2] <= 0.0f)
		{
			cornerArea[face_index][v[0]] = -0.25 * l2[1] * area / (e[2] | e[1]);
			cornerArea[face_index][v[1]] = -0.25 * l2[0] * area / (e[2] | e[0]);
			cornerArea[face_index][v[2]] = area - cornerArea[face_index][v[0]] - cornerArea[face_index][v[1]];
		}
		else
		{
			double ewscale = 0.5 * area / (ew[0] + ew[1] + ew[2]);
			for (int j = 0; j < 3; j++)
			{
				cornerArea[face_index][v[j]] = ewscale * (ew[(j + 1) % 3] + ew[(j + 2) % 3]);
			}
		}

		//#pragma omp atomic
		pointArea[v[0]] += cornerArea[face_index][v[0]];
		//#pragma omp atomic
		pointArea[v[1]] += cornerArea[face_index][v[1]];
		//#pragma omp atomic
		pointArea[v[2]] += cornerArea[face_index][v[2]];
	}
}

void PrincipalCurvature::rot_coord_sys(const Vec3d& old_u, const Vec3d& old_v, const Vec3d& new_norm, Vec3d& new_u, Vec3d& new_v)
{
	new_u = old_u;
	new_v = old_v;
	Vec3d old_norm = (old_u % old_v);
	double ndot = (old_norm | new_norm);
	if (ndot <= -1.0)
	{
		new_u = -new_u;
		new_v = -new_v;
		return;
	}
	Vec3d perp_old = new_norm - ndot * old_norm;
	Vec3d dperp = 1.0f / (1 + ndot) * (old_norm + new_norm);
	new_u -= dperp * (new_u | perp_old);
	new_v -= dperp * (new_v | perp_old);
}

void PrincipalCurvature::proj_curv(const Vec3d& old_u, const Vec3d& old_v, double old_ku, double old_kuv, double old_kv, const Vec3d& new_u, const Vec3d& new_v, double& new_ku, double& new_kuv, double& new_kv)
{
	Vec3d r_new_u; Vec3d r_new_v;
	rot_coord_sys(new_u, new_v, (old_u % old_v), r_new_u, r_new_v);

	double u1 = (r_new_u | old_u);
	double v1 = (r_new_u | old_v);
	double u2 = (r_new_v | old_u);
	double v2 = (r_new_v | old_v);
	new_ku = old_ku * u1 * u1 + old_kuv * (2.0f * u1 * v1) + old_kv * v1 * v1;
	new_kuv = old_ku * u1 * u2 + old_kuv * (u1 * v2 + u2 * v1) + old_kv * v1 * v2;
	new_kv = old_ku * u2 * u2 + old_kuv * (2.0f * u2 * v2) + old_kv * v2 * v2;
}

void PrincipalCurvature::diagonalize_curv(const Vec3d& old_u, const Vec3d& old_v, double ku, double kuv, double kv, const Vec3d& new_norm, Vec3d& pdir1, Vec3d& pdir2, double& vk1, double& vk2)
{
	Vec3d r_old_u, r_old_v;
	rot_coord_sys(old_u, old_v, new_norm, r_old_u, r_old_v);

	double c = 1.0, s = 0.0, tt = 0.0;
	if (kuv != 0.0)
	{
		// Jacobi rotation to diagonalize
		double h = 0.5 * (kv - ku) / kuv;
		tt = (h < 0.0) ?
			1.0 / (h - std::sqrt(1.0 + h * h)) :
			1.0 / (h + std::sqrt(1.0 + h * h));
		c = 1.0 / std::sqrt(1.0 + tt * tt);
		s = tt * c;
	}

	vk1 = ku - tt * kuv;
	vk2 = kv + tt * kuv;

	if (std::abs(vk1) >= std::abs(vk2))
	{
		pdir1 = c * r_old_u - s * r_old_v;
	}
	else
	{
		std::swap(vk1, vk2);
		pdir1 = s * r_old_u + c * r_old_v;
	}
	pdir2 = (new_norm % pdir1);
}

/******************************/
/*     Gauss Curvature        */
/******************************/

GaussCurvature::GaussCurvature(SMeshT* _mesh, bool _update)
	:mesh(_mesh)
{
	if (_update)
		update();
}

void GaussCurvature::update()
{
	edge_normal.resize(mesh->n_edges());
	f_vId.resize(mesh->n_faces());
	f_alpha.resize(mesh->n_faces());

	calcFaceAlpha();
	calcEdgeNormal();
}

std::vector<double> GaussCurvature::vertK()
{
	std::vector<double> K; K.resize(mesh->n_vertices());
	for (auto vh : mesh->vertices())
	{
		if (mesh->is_boundary(vh))
			K[vh.idx()] = M_PI;
		else
			K[vh.idx()] = 2 * M_PI;
	}

	for (auto fh : mesh->faces())
	{
		K[f_vId[fh.idx()][0]] -= f_alpha[fh.idx()][0];
		K[f_vId[fh.idx()][1]] -= f_alpha[fh.idx()][1];
		K[f_vId[fh.idx()][2]] -= f_alpha[fh.idx()][2];
	}
	return K;
}

std::vector<double> GaussCurvature::edgeK()
{
	std::vector<double> K; K.resize(mesh->n_edges());
	for (auto eh : mesh->edges())
	{
		if (mesh->is_boundary(eh)) continue;

		auto hh0 = mesh->halfedge_handle(eh, 0);
		auto hh1 = mesh->halfedge_handle(eh, 1);

		const Vec3d& en = edge_normal[eh.idx()];

		auto hh = mesh->next_halfedge_handle(hh0);
		const Vec3d& en1 = edge_normal[mesh->edge_handle(hh).idx()];

		hh = mesh->next_halfedge_handle(hh);
		const Vec3d& en2 = edge_normal[mesh->edge_handle(hh).idx()];

		hh = mesh->next_halfedge_handle(hh1);
		const Vec3d& en3 = edge_normal[mesh->edge_handle(hh).idx()];

		hh = mesh->next_halfedge_handle(hh);
		const Vec3d& en4 = edge_normal[mesh->edge_handle(hh).idx()];

		double i12 = fabs(en1 % en2 | en);
		double i34 = fabs(en3 % en4 | en);
		double i23 = fabs(en2 % en3 | en);
		double i14 = fabs(en4 % en1 | en);

		K[eh.idx()] = i12 + i34 + i23 + i14;
	}
	return K;
}

std::vector<double> GaussCurvature::faceK()
{
	std::vector<double> K;
	for (auto fh : mesh->faces())
	{
		auto fhe_it(mesh->cfh_begin(fh));

		const Vec3d& en0 = edge_normal[mesh->edge_handle(*fhe_it).idx()];		++fhe_it;
		const Vec3d& en1 = edge_normal[mesh->edge_handle(*fhe_it).idx()];		++fhe_it;
		const Vec3d& en2 = edge_normal[mesh->edge_handle(*fhe_it).idx()];

		K[fh.idx()] = fabs(en0 % en1 | en2);
	}
	return K;
}

void GaussCurvature::calcMaxPrincpalCurvature(
	Vec3d dp1, Vec3d dp2,
	Vec3d dn1, Vec3d dn2, Vec3d& KMax)
{
	double ldp1 = dp1.norm();
	dp1 /= ldp1;
	double ldp2 = dp2.norm();
	dp2 /= ldp2;

	double ldn1 = dn1.norm();
	dn1 /= ldn1;
	double ldn2 = dn2.norm();
	dn2 /= ldn2;

	Eigen::Matrix2d DP, DN, Lambda;

	double cosTheta = dp1 | dp2;
	DP(0, 0) = ldp1;	DP(1, 0) = 0;
	DP(0, 1) = ldp2 * cosTheta;
	DP(1, 1) = ldp2 * sqrt(1 - cosTheta * cosTheta);

	cosTheta = dn1 | dn2;
	DN(0, 0) = ldn1;	DN(1, 0) = 0;
	DN(0, 1) = ldn2 * cosTheta;
	DN(1, 1) = ldn2 * sqrt(1 - cosTheta * cosTheta);

	Lambda = DN * DP.inverse();
	Eigen::EigenSolver<Eigen::Matrix2d> eigen_solver(Lambda);
	auto eigenValues = eigen_solver.eigenvalues();
	auto eigenVectors = eigen_solver.eigenvectors();

	assert(eigenValues(0).imag() == 0 && eigenValues(1).imag() == 0);

	int i = fabs(eigenValues(0).real()) > fabs(eigenValues(1).real()) ? 0 : 1;

	if (fabs(eigenValues(i).real()) <= 1e-10)
	{
		KMax = Vec3d(0);
	}
	else
	{
		KMax = eigenVectors.col(i).real()(0) * dp1 + eigenVectors.col(i).real()(1) * (dp2 - (dp2 | dp1) * dp1).normalize();
	}

	/*std::cout << eigenValues(0).real() << " " << eigenValues(1).real() << std::endl;
	std::cout << eigenVectors.col(0).real() << std::endl;
	std::cout << eigenVectors.col(1).real() << std::endl;*/
}

void GaussCurvature::faceKMax(Vec3d* KMax)
{
	for (auto fh : mesh->faces())
	{
		auto fhe_it(mesh->cfh_begin(fh));

		Vec3d p0 = (mesh->point(mesh->from_vertex_handle(*fhe_it)) + mesh->point(mesh->to_vertex_handle(*fhe_it))) * 0.5;
		const Vec3d& en0 = edge_normal[mesh->edge_handle(*fhe_it).idx()];		++fhe_it;

		Vec3d p1 = (mesh->point(mesh->from_vertex_handle(*fhe_it)) + mesh->point(mesh->to_vertex_handle(*fhe_it))) * 0.5;
		const Vec3d& en1 = edge_normal[mesh->edge_handle(*fhe_it).idx()];		++fhe_it;

		Vec3d p2 = (mesh->point(mesh->from_vertex_handle(*fhe_it)) + mesh->point(mesh->to_vertex_handle(*fhe_it))) * 0.5;
		const Vec3d& en2 = edge_normal[mesh->edge_handle(*fhe_it).idx()];

		calcMaxPrincpalCurvature(p1 - p0, p2 - p0, en0 - en1, en0 - en2, KMax[fh.idx()]);
	}
}

void GaussCurvature::calcFaceAlpha()
{
	for (auto fh : mesh->faces())
	{
		auto sector_angle = face_sector_angle(mesh, fh);

		for (size_t i = 0;i < 3;i++)
		{
			f_vId[fh.idx()][i] = mesh->to_vertex_handle(sector_angle[i].first).idx();
			f_alpha[fh.idx()][i] = std::acos(sector_angle[i].second);
		}
	}
}

void GaussCurvature::calcEdgeNormal()
{
	for (auto eh : mesh->edges())
	{
		auto hh0 = mesh->halfedge_handle(eh, 0);
		auto hh1 = mesh->halfedge_handle(eh, 1);
		auto fromvh = mesh->from_vertex_handle(hh0);
		auto tovh = mesh->to_vertex_handle(hh0);

		Vec3d evec = (mesh->point(tovh) - mesh->point(fromvh)).normalized();

		if (mesh->is_boundary(eh))
		{
			if (mesh->is_boundary(hh0))
			{
				edge_normal[eh.idx()] = mesh->normal(mesh->face_handle(hh1));
			}
			else
			{
				edge_normal[eh.idx()] = mesh->normal(mesh->face_handle(hh0));
			}
		}
		else
		{
			const Vec3d& fn0 = mesh->normal(mesh->face_handle(hh0));
			const Vec3d& fn1 = mesh->normal(mesh->face_handle(hh1));
			Vec3d fv0 = fn0 % evec;
			Vec3d fv1 = evec % fn1;
			Vec3d tmp1 = fn0 + fn1;
			Vec3d tmp2 = fv0 + fv1;

			if (tmp1.norm() >= tmp2.norm())
			{
				tmp1.normalize();
				edge_normal[eh.idx()] = tmp1;
			}
			else
			{
				if ((tmp1 | tmp2) < 0)
				{
					tmp2 = -tmp2;
				}

				tmp2.normalize();
				edge_normal[eh.idx()] = tmp2;
			}
		}
	}
}
}// namespace SurfaceMesh
}// namespace Cage