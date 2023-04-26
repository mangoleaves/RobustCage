#include "LatticePointsGenerator.hh"

namespace Cage
{
namespace CageInit
{

LatticePointsGenerator::LatticePointsGenerator(SM::SMeshT* _SMesh, ParamLatticePointsGenerator* _param)
{
  SMesh = _SMesh;
  param = _param;
}

/// @brief generate lattice points with given parameters.
const std::vector<Vec3d>& LatticePointsGenerator::generate()
{
  getBoundingBox();
  buildAABBTree();
  uniformlyDistributed();
  offsetVertices();
  releaseAABBTree();
  return points;
}


///@brief get bounding box for tri_mesh
void LatticePointsGenerator::getBoundingBox()
{
  bbox = BoundingBox(SMesh->point(SM::VertexHandle(0)), SMesh->point(SM::VertexHandle(0)));
  for (SM::VertexHandle vh : SMesh->vertices())
    bbox += SMesh->point(vh);
  double bbox_diag_length = (bbox.max() - bbox.min()).norm();
  length_threshold = bbox_diag_length * param->minDistanceFactor;
  area_threshold = param->areaThresholdRate * sqrt(3) * length_threshold * length_threshold;
}

/// @brief Add lattice points into triangle mesh.
void LatticePointsGenerator::uniformlyDistributed()
{
  REP_FUNC;
  // check input
  ASSERT(SMesh->n_vertices() > 0, "triMesh has no points.");
  ASSERT(param->bboxScale > 1.0, "bboxScale should be a number >= 1.");
  ASSERT(param->pointNumAlongAxis >= 2, "pointNumAlongAxis should be a number >= 2.");

  // scale bbox
  BoundingBox larger_bbox;
  Vec3d center = (bbox.min() + bbox.max()) * 0.5;
  larger_bbox.min() = center + (bbox.min() - center) * param->bboxScale;
  larger_bbox.max() = center + (bbox.max() - center) * param->bboxScale;

  // compute distance on one axis
  uint32_t num = param->pointNumAlongAxis;
  double dis_x = (larger_bbox.max()[0] - larger_bbox.min()[0]) / (num - 1);
  double dis_y = (larger_bbox.max()[1] - larger_bbox.min()[1]) / (num - 1);
  double dis_z = (larger_bbox.max()[2] - larger_bbox.min()[2]) / (num - 1);

  // reserve memory
  points.reserve(points.size() + num * num * num);

  // add lattice points
  Vec3d current_point(larger_bbox.min());
  std::pair<Vec3d, AABBClosestSearch::HeavyTriPtr> hint = aabbTree->best_hint(current_point);

  for (uint32_t cnt_x = 0;cnt_x < num;cnt_x++)
  {
    for (uint32_t cnt_y = 0;cnt_y < num;cnt_y++)
    {
      for (uint32_t cnt_z = 0;cnt_z < num;cnt_z++)
      {
        // calculate distance from current point to mesh surface.
        // we store result in "hint", which means
        // we use result of this time as hint of next time.
        hint = aabbTree->closest_point(current_point, hint);
        double distance = (hint.first - current_point).squaredNorm();
        // reserve lattice points that are far from mesh surface.
        if (distance > length_threshold * length_threshold)
        {
          points.push_back(current_point);
        }
        current_point[2] += dis_z;
      }
      current_point[2] = larger_bbox.min()[2];
      current_point[1] += dis_y;
    }
    current_point[1] = larger_bbox.min()[1];
    current_point[0] += dis_x;
  }
  {Logger::user_logger->info("adding {} lattice points.", points.size());}
}

void LatticePointsGenerator::offsetVertices()
{
  REP_FUNC;
  // check input
  ASSERT(SMesh->n_vertices() > 0, "tri_mesh has no points.");

#ifdef CLUSTER
  if (param->doCluster)
  {
    // cluster
    SM::LloydCluster cluster(SMesh);
    cluster.run(param->clusterSize, param->clusterIter);
    std::vector<SM::FaceHandle> faces = cluster.centerFaces();
    const std::vector<Vec3d>& centroids = cluster.getCentroid();

    for (SM::FaceHandle fh : faces)
    {
      // 1. calculate normal for each face.
      const Normal& normal = SMesh->normal(fh);
      const Point& center = centroids[fh.idx()];

      double avg_length = 0.0;
      for (auto eh : SMesh->fe_range(fh))
        avg_length += SMesh->data(eh).edge_length;
      avg_length /= 3;
      // 2. offset each face's barycenter up and down.
      // NOTE: we only accep new points whose closest point is on current face.
      // we will try offseting barycenter at distance 1.0/0.75/0.5 * avg_edge_length.
      // if all fails, discard this point.
      double offset_length_scale = 1.0;
      auto hint = aabbTree->best_hint(center);
      while (offset_length_scale >= 0.5)
      {
        Vec3d offset = normal * offset_length_scale * avg_length;
        Vec3d up = center + offset;
        auto closest = aabbTree->closest_point(up, hint);
        if (closest.second->index == fh.idx())
        {
          points.push_back(up);
          break;
        }
        else
        {
          offset_length_scale -= 0.25;
        }
      }
    }
  }
  else
  #endif
  {
    SMesh->request_face_normals();
    SMesh->update_face_normals();

    points.reserve(points.size() + SMesh->n_faces());
    for (SM::FaceHandle fh : SMesh->faces())
    {
      const Normal& normal = SMesh->normal(fh);
      auto fv_it = SMesh->fv_iter(fh);
      const Point& a = SMesh->point(*fv_it);fv_it++;
      const Point& b = SMesh->point(*fv_it);fv_it++;
      const Point& c = SMesh->point(*fv_it);

      offsetPoint(fh, a, b, c, normal);
    }
  }

  {Logger::user_logger->info("adding {} offset points.", points.size());}
}

void LatticePointsGenerator::offsetPoint(
  SM::FaceHandle fh, const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& n)
{
  double face_area = 0.5 * ((a - b).cross(c - b)).length();
  if (face_area > area_threshold)
  {
    // subdivide triangle and offset
    offsetPoint(fh, a, 0.5 * (a + b), 0.5 * (a + c), n);
    offsetPoint(fh, b, 0.5 * (b + c), 0.5 * (a + b), n);
    offsetPoint(fh, c, 0.5 * (a + c), 0.5 * (b + c), n);
    offsetPoint(fh, 0.5 * (a + b), 0.5 * (b + c), 0.5 * (a + c), n);
  }
  else
  {
    // offset each face's barycenter up and down.
    // NOTE: we only accep new points whose closest point is on current face.
    Vec3d center = (a + b + c) / 3.0;
    double avg_length = ((a - b).length() + (a - c).length() + (b - c).length()) / 3.0;
    double offset_length_scale = 1.0;
    auto hint = aabbTree->best_hint(center);

    Vec3d offset = n * avg_length;
    Vec3d up = center + offset;
    Vec3d down = center - offset;
    if (aabbTree->closest_point(up, hint).second->index == fh.idx())
      points.push_back(up);
    if (aabbTree->closest_point(down, hint).second->index == fh.idx())
      points.push_back(down);
  }
}

void LatticePointsGenerator::buildAABBTree()
{
  // build AABB tree on mesh
  std::vector<HeavyTriangle> aabb_triangles;
  aabb_triangles.reserve(SMesh->n_faces());
  for (SM::FaceHandle fh : SMesh->faces())
  {
    auto fv_it = SMesh->fv_begin(fh);
    auto& p0 = SMesh->point(*fv_it); fv_it++;
    auto& p1 = SMesh->point(*fv_it); fv_it++;
    auto& p2 = SMesh->point(*fv_it); fv_it++;
    aabb_triangles.emplace_back(p0, p1, p2, fh.idx());
  }
  aabbTree = std::make_unique<AABBClosestSearch>();
  aabbTree->insert(std::move(aabb_triangles));
  aabbTree->build();
  aabb_triangles.clear();
}

void LatticePointsGenerator::releaseAABBTree()
{
  if (aabbTree)
  {
    aabbTree->clear();
    aabbTree.reset();
  }
}
}// namespace CageInit
}// namespace Cage