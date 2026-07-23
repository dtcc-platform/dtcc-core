// Copyright (C) 2022 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_BUILDING_PROCESSOR_H
#define DTCC_BUILDING_PROCESSOR_H

#include "BoundingBox.h"
#include "KDTreeVectorOfVectorsAdaptor.h"
#include "Logging.h"
#include "PointCloudProcessor.h"
#include "Polyfix.h"
#include "Timer.h"
#include "model/GridField.h"
#include "model/Polygon.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

class BuildingProcessor
{
public:

  /// Extract roof points from point cloud.
  ///
  ///
  /// The roof points of a building are defined as all points
  /// of class 6 (Building) that fall within the building
  /// footprint. However, since that classification seems to
  /// be missing in the data from LM, we are currently using
  /// all points (except class 2 and 9).
  ///
  /// @param city The city
  /// @param point_cloud Point cloud (unfiltered)
  static std::vector<std::vector<Vector3D>>
  extract_building_points(const std::vector<Polygon> &footprints,
                          const std::vector<Vector3D> &points)
  {
    info("Computing building points...");
    Timer timer("compute_building_points");

    // Check that point cloud is not empty
    if (points.empty())
      error("empty point cloud");

    // auto tile_city_timer = Timer("Tile City");
    // auto tiles_city = CityProcessor::tile_citymodel(
    //     city, point_cloud, point_cloud.bounding_box, 4, 4);
    // tile_city_timer.stop();
    auto kdt_timer = Timer("ExtractBuildingPoints: BuildKDTree");
    // build a kd-tree for radius search
    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Vector3D>, double,
                                         2 /* dims */>
        my_kd_tree_t;
    my_kd_tree_t pc_index(2, points, 20 /* max leaf */);
    kdt_timer.stop();

    // Iterate over buildings
    std::vector<std::vector<Vector3D>> building_points;
    for (auto &footprint : footprints)
    {
      std::vector<Vector3D> roof_points;
      auto centerPoint = Geometry::polygon_center_2d(footprint);
      double radius = Geometry::polygon_radius_2d(footprint, centerPoint);
      radius *= radius;

      std::vector<double> query_pt{centerPoint.x, centerPoint.y};
      auto radius_t = Timer("RadiusQuery");
      auto indices_dists = pc_index.radius_query(&query_pt[0], radius);
      radius_t.stop();
      for (auto const &ind_pt : indices_dists)
      {
        size_t idx = ind_pt.first;
        const Vector3D &p_3d = points[idx];
        const Vector2D p_2d{p_3d.x, p_3d.y};
        if (Geometry::polygon_contains_2d(footprint, p_2d))
        {
          roof_points.push_back(p_3d);
        }
      }
      building_points.push_back(roof_points);
    }
    return building_points;
  }
};

} // namespace DTCC_BUILDER

#endif
