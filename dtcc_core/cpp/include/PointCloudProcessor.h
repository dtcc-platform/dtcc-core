// Copyright (C) 2020-2022 Dag Wästberg
// Licensed under the MIT License

#ifndef DTCC_POINT_CLOUD_PROCESSOR_H
#define DTCC_POINT_CLOUD_PROCESSOR_H

#include <fstream>
#include <iso646.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "KDTreeVectorOfVectorsAdaptor.h"
#include "nanoflann.hpp"

#include "Geometry.h"
#include "Timer.h"
#include "model/Vector.h"

namespace DTCC_BUILDER
{

class PointCloudProcessor
{
public:

  static std::vector<std::vector<double>>
  knn_nearest_neighbours(std::vector<Vector3D> &points, size_t neighbours)
  {
    size_t pc_size = points.size();
    std::vector<std::vector<double>> neighbour_dist(pc_size);

    if (neighbours <= 0 or neighbours > pc_size)
    {
      neighbours = pc_size;
    }
    neighbours++; // N neighbours other than ourselves

    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Vector3D>, double,
                                         3 /* dims */>
        my_kd_tree_t;
    my_kd_tree_t pc_index(3 /*dim*/, points, 10 /* max leaf */);
    pc_index.index->buildIndex();
    std::vector<size_t> ret_indexes(neighbours);
    std::vector<double> out_dists_sqr(neighbours);
    nanoflann::KNNResultSet<double> resultSet(neighbours);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);

    size_t idx = 0;
    for (auto const &pt : points)
    {
      std::vector<double> query_pt{pt.x, pt.y, pt.z};
      pc_index.query(&query_pt[0], neighbours, &ret_indexes[0],
                     &out_dists_sqr[0]);
      for (size_t i = 1; i < neighbours;
           i++) // start from 1 since 0 is the query point
      {
        neighbour_dist[idx].push_back(std::sqrt(out_dists_sqr[i]));
      }
      idx++;
    }

    return neighbour_dist;
  }

  static std::vector<std::vector<size_t>> points_in_polygons(const std::vector<Vector3D> &points, const std::vector<Polygon> &polygons)
  {
    if (points.empty())
    {
      warning("empty point cloud");
      return std::vector<std::vector<size_t>>();
    }

    typedef KDTreeVectorOfVectorsAdaptor<std::vector<Vector3D>, double,
                                         2 /* dims */> my_kd_tree_t;
    my_kd_tree_t pc_index(2, points, 20 /* max leaf */);
    std::vector<std::vector<size_t>> pip_indices;

    for (auto &polygon : polygons)
    {

      auto centerPoint = Geometry::polygon_center_2d(polygon);
      double radius = Geometry::polygon_radius_2d(polygon, centerPoint);
      radius *= radius;
      std::vector<double> query_pt{centerPoint.x, centerPoint.y};
      auto indices_dists = pc_index.radius_query(&query_pt[0], radius);
      std::vector<size_t> poly_points;

      for (auto const &ind_pt : indices_dists)
      {
        size_t idx = ind_pt.first;
        const Vector3D &p_3d = points[idx];
        const Vector2D p_2d{p_3d.x, p_3d.y};
        if (Geometry::polygon_contains_2d(polygon, p_2d))
        {
         poly_points.push_back(idx);
        }
      }
      pip_indices.push_back(poly_points);
    }
    return pip_indices;
  }

  /// Finds outliers from vector of points by removing all points more than a
  /// given number of standard deviations from the mean distance to their N
  /// nearest neighbours
  ///
  /// @param points The vector of points
  /// @param neighbours Number of neighbours to consider. If less than 1 or
  /// greater than the number of points in the point cloud use all points
  /// @param outlier_margin Number of standard deviations
  /// @return Vector of indices of outlier points
  static std::vector<size_t>
  statistical_outlier_finder(std::vector<Vector3D> &points,
                             size_t neighbours,
                             double outlier_margin,
                             bool verbose = false)
  {
    Timer("StatisticalOurtierFinder");
    // Check that we have enough points
    if (points.size() <= neighbours)
      return std::vector<size_t>();

    std::vector<size_t> outliers;

    auto neighbour_dist = knn_nearest_neighbours(points, neighbours);
    std::vector<double> u_dist_i;

    for (size_t i = 0; i < points.size(); i++)
    {
      double dsum = 0;
      for (auto &d : neighbour_dist[i])
      {
        dsum += d;
      }
      u_dist_i.push_back(dsum / neighbours);
    }

    // Compute mean
    double mean{0};
    for (auto p : u_dist_i)
      mean += p;
    mean /= u_dist_i.size();

    // Compute standard deviation
    double std{0};
    for (auto p : u_dist_i)
      std += (p - mean) * (p - mean);
    std /= u_dist_i.size() - 1;
    std = std::sqrt(std);

    double T = mean + outlier_margin * std;

    // info("T: " + str(T));
    for (size_t i = 0; i < u_dist_i.size(); i++)
    {
      if (u_dist_i[i] > T)
        outliers.push_back(i);
    }

    return outliers;
  }

  /// Remove outliers from Vector<Point3d> using Statistical Outlier algorithm
  ///
  /// @param points vector of points to filter
  /// @param neighbours Number of neighbours to consider. If less than 1 or
  /// greater than the number of points in the point cloud use all points
  /// @param outlier_margin Number of standard deviations
  /// @param verbose give verbose detail
  static void statistical_outlier_remover(std::vector<Vector3D> &points,
                                          size_t neighbours,
                                          double outlier_margin,
                                          bool verbose = false)
  {
    Timer("StatisticalOurtierRemover");
    std::vector<size_t> outliers =
        statistical_outlier_finder(points, neighbours, outlier_margin, verbose);
    if (outliers.size() == 0)
      return;
    std::vector<Vector3D> new_points;
    size_t k = 0;
    for (size_t i = 0; i < points.size(); i++)
    {
      if (k >= outliers.size() || i != outliers[k])
      {
        new_points.push_back(points[i]);
      }
      else
      {
        k++;
      }
    }
    points = new_points;
  }
};

} // namespace DTCC_BUILDER

#endif
