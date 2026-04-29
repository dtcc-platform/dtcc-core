// Copyright (C) 2023 Dag Wästberg
// Licensed under the MIT License
//
// Modified by Anders Logg 2023

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "BuildingProcessor.h"
#include "Intersection.h"
#include "MeshBuilder.h"
#include "MeshProcessor.h"
#include "Smoother.h"
#include "VertexSmoother.h"
#include "VolumeMeshBuilder.h"
#include "model/GridField.h"
#include "model/Mesh.h"
#include "model/Polygon.h"
#include "model/Simplices.h"
#include "model/Vector.h"
#include "model/VolumeMesh.h"
#include "terrain-mesher/zemlya.hpp"

namespace py = pybind11;

namespace DTCC_BUILDER
{

namespace
{

struct BoundaryPoint
{
  double x = 0.0;
  double y = 0.0;
};

struct BoundaryBBox
{
  double min_x = std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();
};

struct BoundaryPolygonData
{
  std::vector<BoundaryPoint> shell;
  std::vector<std::vector<BoundaryPoint>> holes;
  BoundaryBBox bbox;
};

struct BoundaryStatsData
{
  size_t edge_count = 0;
  size_t short_edge_count = 0;
  size_t vertex_count = 0;
  double min_edge_length = std::numeric_limits<double>::infinity();
};

enum class PairRelation
{
  none,
  point_touch,
  close_pair,
  shared_overlap,
};

double squared_distance(const BoundaryPoint &a, const BoundaryPoint &b)
{
  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  return dx * dx + dy * dy;
}

double segment_length(const BoundaryPoint &a, const BoundaryPoint &b)
{
  return std::sqrt(squared_distance(a, b));
}

bool points_equal(const BoundaryPoint &a, const BoundaryPoint &b, double tolerance)
{
  return squared_distance(a, b) <= tolerance * tolerance;
}

double cross(const BoundaryPoint &a, const BoundaryPoint &b, const BoundaryPoint &c)
{
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

double point_to_segment_distance(const BoundaryPoint &point, const BoundaryPoint &start,
                                 const BoundaryPoint &end)
{
  const double dx = end.x - start.x;
  const double dy = end.y - start.y;
  const double base = std::hypot(dx, dy);
  if (base <= 1e-12)
    return std::hypot(point.x - start.x, point.y - start.y);
  return std::abs((point.x - start.x) * dy - (point.y - start.y) * dx) / base;
}

bool bbox_is_valid(const BoundaryBBox &bbox)
{
  return std::isfinite(bbox.min_x) && std::isfinite(bbox.min_y) && std::isfinite(bbox.max_x) &&
         std::isfinite(bbox.max_y);
}

void update_bbox(BoundaryBBox &bbox, const BoundaryPoint &point)
{
  bbox.min_x = std::min(bbox.min_x, point.x);
  bbox.min_y = std::min(bbox.min_y, point.y);
  bbox.max_x = std::max(bbox.max_x, point.x);
  bbox.max_y = std::max(bbox.max_y, point.y);
}

double bbox_distance(const BoundaryBBox &a, const BoundaryBBox &b)
{
  if (!bbox_is_valid(a) || !bbox_is_valid(b))
    return 0.0;

  const double dx = std::max({a.min_x - b.max_x, b.min_x - a.max_x, 0.0});
  const double dy = std::max({a.min_y - b.max_y, b.min_y - a.max_y, 0.0});
  return std::hypot(dx, dy);
}

std::vector<BoundaryPoint> parse_ring(py::handle ring_object)
{
  std::vector<BoundaryPoint> ring;
  const py::sequence ring_sequence = py::reinterpret_borrow<py::sequence>(ring_object);
  ring.reserve(ring_sequence.size());
  for (py::handle point_object : ring_sequence)
  {
    const py::sequence point_sequence = py::reinterpret_borrow<py::sequence>(point_object);
    if (point_sequence.size() < 2)
      continue;
    ring.push_back({point_sequence[0].cast<double>(), point_sequence[1].cast<double>()});
  }

  if (ring.size() >= 2 && points_equal(ring.front(), ring.back(), 0.0))
    ring.pop_back();

  return ring;
}

BoundaryPolygonData parse_boundary_polygon(py::handle polygon_object)
{
  BoundaryPolygonData polygon;
  const py::tuple polygon_tuple = polygon_object.cast<py::tuple>();
  polygon.shell = parse_ring(polygon_tuple[0]);
  for (const auto &point : polygon.shell)
    update_bbox(polygon.bbox, point);

  const py::sequence holes_sequence = py::reinterpret_borrow<py::sequence>(polygon_tuple[1]);
  polygon.holes.reserve(holes_sequence.size());
  for (py::handle hole_object : holes_sequence)
  {
    auto hole = parse_ring(hole_object);
    for (const auto &point : hole)
      update_bbox(polygon.bbox, point);
    polygon.holes.push_back(std::move(hole));
  }
  return polygon;
}

std::vector<BoundaryPolygonData> parse_boundary_polygons(const py::list &polygons_xy)
{
  std::vector<BoundaryPolygonData> polygons;
  polygons.reserve(polygons_xy.size());
  for (py::handle polygon_object : polygons_xy)
    polygons.push_back(parse_boundary_polygon(polygon_object));
  return polygons;
}

void accumulate_ring_stats(const std::vector<BoundaryPoint> &ring, double short_edge_threshold,
                           BoundaryStatsData &stats)
{
  if (ring.size() < 3)
    return;

  stats.vertex_count += ring.size();
  for (size_t i = 0; i < ring.size(); ++i)
  {
    const BoundaryPoint &a = ring[i];
    const BoundaryPoint &b = ring[(i + 1) % ring.size()];
    const double length = segment_length(a, b);
    stats.edge_count += 1;
    stats.min_edge_length = std::min(stats.min_edge_length, length);
    if (short_edge_threshold > 0.0 && length + 1e-12 < short_edge_threshold)
      stats.short_edge_count += 1;
  }
}

BoundaryStatsData compute_boundary_stats(const BoundaryPolygonData &polygon, double threshold)
{
  BoundaryStatsData stats;
  accumulate_ring_stats(polygon.shell, threshold, stats);
  for (const auto &hole : polygon.holes)
    accumulate_ring_stats(hole, threshold, stats);
  return stats;
}

double dot(const BoundaryPoint &a, const BoundaryPoint &b, const BoundaryPoint &c)
{
  return (b.x - a.x) * (c.x - a.x) + (b.y - a.y) * (c.y - a.y);
}

bool on_segment(const BoundaryPoint &a, const BoundaryPoint &b, const BoundaryPoint &c,
                double tolerance)
{
  if (std::abs(cross(a, b, c)) > tolerance)
    return false;
  return c.x >= std::min(a.x, b.x) - tolerance && c.x <= std::max(a.x, b.x) + tolerance &&
         c.y >= std::min(a.y, b.y) - tolerance && c.y <= std::max(a.y, b.y) + tolerance;
}

bool segments_intersect(const BoundaryPoint &a, const BoundaryPoint &b, const BoundaryPoint &c,
                        const BoundaryPoint &d, double tolerance)
{
  const double ab_c = cross(a, b, c);
  const double ab_d = cross(a, b, d);
  const double cd_a = cross(c, d, a);
  const double cd_b = cross(c, d, b);

  if (((ab_c > tolerance && ab_d < -tolerance) || (ab_c < -tolerance && ab_d > tolerance)) &&
      ((cd_a > tolerance && cd_b < -tolerance) || (cd_a < -tolerance && cd_b > tolerance)))
    return true;

  return on_segment(a, b, c, tolerance) || on_segment(a, b, d, tolerance) ||
         on_segment(c, d, a, tolerance) || on_segment(c, d, b, tolerance);
}

double collinear_overlap_length(const BoundaryPoint &a, const BoundaryPoint &b,
                                const BoundaryPoint &c, const BoundaryPoint &d,
                                double tolerance)
{
  if (segment_length(a, b) <= tolerance || segment_length(c, d) <= tolerance)
    return 0.0;
  if (std::abs(cross(a, b, c)) > tolerance || std::abs(cross(a, b, d)) > tolerance)
    return 0.0;

  const bool use_x = std::abs(b.x - a.x) >= std::abs(b.y - a.y);
  const double a0 = use_x ? a.x : a.y;
  const double a1 = use_x ? b.x : b.y;
  const double c0 = use_x ? c.x : c.y;
  const double c1 = use_x ? d.x : d.y;
  const double min_a = std::min(a0, a1);
  const double max_a = std::max(a0, a1);
  const double min_c = std::min(c0, c1);
  const double max_c = std::max(c0, c1);
  const double overlap = std::min(max_a, max_c) - std::max(min_a, min_c);
  return overlap > tolerance ? overlap : 0.0;
}

double segment_segment_distance(const BoundaryPoint &a, const BoundaryPoint &b,
                                const BoundaryPoint &c, const BoundaryPoint &d,
                                double tolerance)
{
  if (segments_intersect(a, b, c, d, tolerance))
    return 0.0;

  return std::min(
      {point_to_segment_distance(a, c, d), point_to_segment_distance(b, c, d),
       point_to_segment_distance(c, a, b), point_to_segment_distance(d, a, b)});
}

void accumulate_ring_pair_relation(const std::vector<BoundaryPoint> &left,
                                   const std::vector<BoundaryPoint> &right, double tolerance,
                                   bool &has_shared_overlap, double &min_distance)
{
  if (left.size() < 2 || right.size() < 2)
    return;

  for (size_t i = 0; i < left.size(); ++i)
  {
    const BoundaryPoint &a = left[i];
    const BoundaryPoint &b = left[(i + 1) % left.size()];
    for (size_t j = 0; j < right.size(); ++j)
    {
      const BoundaryPoint &c = right[j];
      const BoundaryPoint &d = right[(j + 1) % right.size()];
      if (collinear_overlap_length(a, b, c, d, tolerance) > tolerance)
      {
        has_shared_overlap = true;
        return;
      }
      min_distance = std::min(min_distance, segment_segment_distance(a, b, c, d, tolerance));
      if (min_distance <= tolerance)
        return;
    }
  }
}

PairRelation classify_pair_relation(const BoundaryPolygonData &left,
                                    const BoundaryPolygonData &right, double target_scale,
                                    double tolerance)
{
  if (bbox_distance(left.bbox, right.bbox) > target_scale + tolerance)
    return PairRelation::none;

  bool has_shared_overlap = false;
  double min_distance = std::numeric_limits<double>::infinity();
  accumulate_ring_pair_relation(left.shell, right.shell, tolerance, has_shared_overlap, min_distance);
  if (!has_shared_overlap)
  {
    for (const auto &hole : left.holes)
      accumulate_ring_pair_relation(hole, right.shell, tolerance, has_shared_overlap, min_distance);
    for (const auto &hole : right.holes)
      accumulate_ring_pair_relation(left.shell, hole, tolerance, has_shared_overlap, min_distance);
    for (const auto &left_hole : left.holes)
    {
      for (const auto &right_hole : right.holes)
        accumulate_ring_pair_relation(left_hole, right_hole, tolerance, has_shared_overlap,
                                      min_distance);
    }
  }

  if (has_shared_overlap)
    return PairRelation::shared_overlap;
  if (min_distance > target_scale + tolerance)
    return PairRelation::none;
  if (min_distance <= tolerance)
    return PairRelation::point_touch;
  return PairRelation::close_pair;
}

int find_short_collinear_vertex(const std::vector<BoundaryPoint> &points, double target_scale,
                                double grid)
{
  const size_t count = points.size();
  if (count < 3)
    return -1;

  const double line_tolerance = std::max({2.0 * grid, 0.05 * target_scale, 1e-9});
  for (size_t index = 0; index < count; ++index)
  {
    const BoundaryPoint &prev_point = points[(index + count - 1) % count];
    const BoundaryPoint &point = points[index];
    const BoundaryPoint &next_point = points[(index + 1) % count];
    const double prev_length = segment_length(prev_point, point);
    const double next_length = segment_length(point, next_point);
    if (std::min(prev_length, next_length) + 1e-12 >= target_scale)
      continue;
    const double offset = point_to_segment_distance(point, prev_point, next_point);
    if (offset <= line_tolerance)
      return static_cast<int>(index);
  }
  return -1;
}

std::pair<int, int> find_short_step_pair(const std::vector<BoundaryPoint> &points,
                                         double target_scale, double grid)
{
  const size_t count = points.size();
  if (count < 4)
    return {-1, -1};

  const double step_width_tolerance = std::max({4.0 * grid, 0.5 * target_scale, 1e-9});
  const double parallel_tolerance = 0.1;
  for (size_t index = 0; index < count; ++index)
  {
    const BoundaryPoint &a = points[(index + count - 1) % count];
    const BoundaryPoint &b = points[index];
    const BoundaryPoint &c = points[(index + 1) % count];
    const BoundaryPoint &d = points[(index + 2) % count];

    const BoundaryPoint ab{b.x - a.x, b.y - a.y};
    const BoundaryPoint bc{c.x - b.x, c.y - b.y};
    const BoundaryPoint cd{d.x - c.x, d.y - c.y};
    const BoundaryPoint ad{d.x - a.x, d.y - a.y};

    const double len_ab = std::hypot(ab.x, ab.y);
    const double len_bc = std::hypot(bc.x, bc.y);
    const double len_cd = std::hypot(cd.x, cd.y);
    const double len_ad = std::hypot(ad.x, ad.y);
    if (len_ab + 1e-12 >= target_scale || len_cd + 1e-12 >= target_scale)
      continue;
    if (len_bc <= std::max(grid, 1e-9) || len_ad <= std::max(grid, 1e-9))
      continue;

    const double ab_cd_cross = std::abs(ab.x * cd.y - ab.y * cd.x) / std::max(len_ab * len_cd, 1e-12);
    const double ab_cd_dot = (ab.x * cd.x + ab.y * cd.y) / std::max(len_ab * len_cd, 1e-12);
    const double ad_bc_cross = std::abs(ad.x * bc.y - ad.y * bc.x) / std::max(len_ad * len_bc, 1e-12);
    if (ab_cd_cross > parallel_tolerance || ab_cd_dot > -0.5)
      continue;
    if (ad_bc_cross > parallel_tolerance)
      continue;

    const double step_width = std::max(point_to_segment_distance(b, a, d),
                                       point_to_segment_distance(c, a, d));
    if (step_width > step_width_tolerance)
      continue;
    return {static_cast<int>(index), static_cast<int>((index + 1) % count)};
  }

  return {-1, -1};
}

std::pair<std::vector<BoundaryPoint>, int> clean_ring_short_edge_chains(
    const std::vector<BoundaryPoint> &input, double target_scale, double grid)
{
  std::vector<BoundaryPoint> unique;
  unique.reserve(input.size());
  for (const auto &point : input)
  {
    if (!unique.empty() && points_equal(point, unique.back(), 0.0))
      continue;
    unique.push_back(point);
  }
  if (unique.size() >= 2 && points_equal(unique.front(), unique.back(), 0.0))
    unique.pop_back();
  if (unique.size() < 3)
    return {{}, 0};

  int removed_count = 0;
  while (unique.size() >= 3)
  {
    const auto pair = find_short_step_pair(unique, target_scale, grid);
    if (pair.first >= 0)
    {
      std::vector<BoundaryPoint> filtered;
      filtered.reserve(unique.size());
      for (size_t index = 0; index < unique.size(); ++index)
      {
        if (static_cast<int>(index) == pair.first || static_cast<int>(index) == pair.second)
          continue;
        filtered.push_back(unique[index]);
      }
      unique.swap(filtered);
      removed_count += 2;
      continue;
    }

    const int vertex = find_short_collinear_vertex(unique, target_scale, grid);
    if (vertex < 0)
      break;
    unique.erase(unique.begin() + vertex);
    removed_count += 1;
  }

  if (unique.size() < 3)
    return {{}, removed_count};
  return {unique, removed_count};
}

py::list boundary_stats(py::list polygons_xy, double short_edge_threshold)
{
  const auto polygons = parse_boundary_polygons(polygons_xy);
  py::list stats_list;
  for (const auto &polygon : polygons)
  {
    const auto stats = compute_boundary_stats(polygon, short_edge_threshold);
    py::dict item;
    item["edge_count"] = py::int_(stats.edge_count);
    item["short_edge_count"] = py::int_(stats.short_edge_count);
    item["vertex_count"] = py::int_(stats.vertex_count);
    if (std::isfinite(stats.min_edge_length))
      item["min_edge_length"] = py::float_(stats.min_edge_length);
    else
      item["min_edge_length"] = py::none();
    if (bbox_is_valid(polygon.bbox))
      item["bbox"] =
          py::make_tuple(polygon.bbox.min_x, polygon.bbox.min_y, polygon.bbox.max_x, polygon.bbox.max_y);
    else
      item["bbox"] = py::none();
    stats_list.append(item);
  }
  return stats_list;
}

py::list boundary_defect_clusters(py::list polygons_xy, double target_scale, double pair_tolerance)
{
  const auto polygons = parse_boundary_polygons(polygons_xy);
  const double tolerance = std::max(target_scale * 1e-6, 1e-9);
  const size_t polygon_count = polygons.size();
  std::vector<BoundaryStatsData> stats;
  stats.reserve(polygon_count);
  for (const auto &polygon : polygons)
    stats.push_back(compute_boundary_stats(polygon, target_scale));

  std::vector<std::vector<int>> adjacency(polygon_count);
  std::vector<std::vector<PairRelation>> pair_relations(
      polygon_count, std::vector<PairRelation>(polygon_count, PairRelation::none));
  std::set<int> defect_indices;

  for (size_t index = 0; index < polygon_count; ++index)
  {
    if (stats[index].short_edge_count > 0)
      defect_indices.insert(static_cast<int>(index));
  }

  for (size_t left = 0; left < polygon_count; ++left)
  {
    for (size_t right = left + 1; right < polygon_count; ++right)
    {
      if (bbox_distance(polygons[left].bbox, polygons[right].bbox) > std::max(pair_tolerance, target_scale) + tolerance)
        continue;
      const auto relation = classify_pair_relation(polygons[left], polygons[right], target_scale, tolerance);
      pair_relations[left][right] = relation;
      pair_relations[right][left] = relation;
      if (relation == PairRelation::point_touch || relation == PairRelation::close_pair)
      {
        defect_indices.insert(static_cast<int>(left));
        defect_indices.insert(static_cast<int>(right));
        adjacency[left].push_back(static_cast<int>(right));
        adjacency[right].push_back(static_cast<int>(left));
      }
    }
  }

  if (defect_indices.empty())
    return py::list();

  std::set<int> remaining = defect_indices;
  std::vector<std::set<int>> components;
  while (!remaining.empty())
  {
    const int seed = *remaining.begin();
    std::vector<int> stack{seed};
    std::set<int> component;
    while (!stack.empty())
    {
      const int index = stack.back();
      stack.pop_back();
      if (component.count(index))
        continue;
      component.insert(index);
      remaining.erase(index);
      for (int neighbor : adjacency[index])
      {
        if (!component.count(neighbor))
          stack.push_back(neighbor);
      }
    }
    components.push_back(std::move(component));
  }

  py::list clusters;
  for (auto &component : components)
  {
    std::set<int> expanded = component;
    for (int index : component)
    {
      for (size_t other = 0; other < polygon_count; ++other)
      {
        if (bbox_distance(polygons[index].bbox, polygons[other].bbox) <= pair_tolerance + tolerance)
          expanded.insert(static_cast<int>(other));
      }
    }

    bool has_short_edges = false;
    bool has_point_touch = false;
    bool has_close_pair = false;
    int short_edge_count = 0;
    for (int index : expanded)
    {
      short_edge_count += static_cast<int>(stats[index].short_edge_count);
      has_short_edges = has_short_edges || stats[index].short_edge_count > 0;
    }
    for (int left : expanded)
    {
      for (int right : expanded)
      {
        if (right <= left)
          continue;
        has_point_touch = has_point_touch || pair_relations[left][right] == PairRelation::point_touch;
        has_close_pair = has_close_pair || pair_relations[left][right] == PairRelation::close_pair;
      }
    }

    std::string kind = "short_edge_only";
    if ((has_point_touch || has_close_pair) && has_short_edges)
      kind = "mixed_pair_short_edge";
    else if (has_point_touch && !has_close_pair)
      kind = "point_touch_pair";
    else if (has_point_touch || has_close_pair)
      kind = "close_pair";

    py::dict item;
    py::list indices;
    for (int index : expanded)
      indices.append(index);
    item["indices"] = indices;
    item["kind"] = kind;
    item["short_edge_count"] = short_edge_count;
    item["pair_issue_count"] = static_cast<int>(has_point_touch) + static_cast<int>(has_close_pair);
    clusters.append(item);
  }
  return clusters;
}

py::dict rewrite_defect_cluster(py::list polygons_xy, py::dict cluster, double target_scale,
                                double grid)
{
  const auto polygons = parse_boundary_polygons(polygons_xy);
  py::dict result;
  result["supported"] = py::bool_(false);

  if (!cluster.contains("indices") || !cluster.contains("kind"))
    return result;

  const std::string kind = cluster["kind"].cast<std::string>();
  if (kind != "short_edge_only" && kind != "mixed_pair_short_edge")
    return result;

  const py::sequence indices_sequence = py::reinterpret_borrow<py::sequence>(cluster["indices"]);
  py::list changed_indices;
  py::list changed_polygons;
  int changed_count = 0;
  size_t resulting_short_edges = 0;
  double min_edge_length = std::numeric_limits<double>::infinity();

  for (py::handle index_object : indices_sequence)
  {
    const int index = index_object.cast<int>();
    if (index < 0 || static_cast<size_t>(index) >= polygons.size())
      return result;

    const auto &polygon = polygons[static_cast<size_t>(index)];
    const auto cleaned_shell = clean_ring_short_edge_chains(polygon.shell, target_scale, grid);
    if (cleaned_shell.first.empty())
      return result;

    py::list holes;
    bool polygon_changed = cleaned_shell.second > 0;
    std::vector<std::vector<BoundaryPoint>> cleaned_holes;
    cleaned_holes.reserve(polygon.holes.size());
    for (const auto &hole : polygon.holes)
    {
      const auto cleaned_hole = clean_ring_short_edge_chains(hole, target_scale, grid);
      polygon_changed = polygon_changed || cleaned_hole.second > 0;
      if (cleaned_hole.first.empty())
        continue;
      cleaned_holes.push_back(cleaned_hole.first);
    }

    if (!polygon_changed)
      continue;

    BoundaryPolygonData cleaned_polygon;
    cleaned_polygon.shell = cleaned_shell.first;
    cleaned_polygon.holes = cleaned_holes;
    const auto stats = compute_boundary_stats(cleaned_polygon, target_scale);
    resulting_short_edges += stats.short_edge_count;
    min_edge_length = std::min(min_edge_length, stats.min_edge_length);

    py::list shell;
    for (const auto &point : cleaned_polygon.shell)
      shell.append(py::make_tuple(point.x, point.y));
    shell.append(py::make_tuple(cleaned_polygon.shell.front().x, cleaned_polygon.shell.front().y));
    for (const auto &hole : cleaned_polygon.holes)
    {
      py::list hole_ring;
      for (const auto &point : hole)
        hole_ring.append(py::make_tuple(point.x, point.y));
      hole_ring.append(py::make_tuple(hole.front().x, hole.front().y));
      holes.append(hole_ring);
    }

    changed_indices.append(index);
    changed_polygons.append(py::make_tuple(shell, holes));
    changed_count += 1;
  }

  if (changed_count == 0)
    return result;

  result["supported"] = py::bool_(true);
  result["changed_indices"] = changed_indices;
  result["polygons"] = changed_polygons;
  result["short_edge_count"] = resulting_short_edges;
  if (std::isfinite(min_edge_length))
    result["min_edge_length"] = py::float_(min_edge_length);
  else
    result["min_edge_length"] = py::none();
  return result;
}

} // namespace

Polygon create_polygon(py::list vertices, py::list holes)
{
  Polygon poly;
  for (size_t i = 0; i < vertices.size(); i++)
  {
    auto pt = vertices[i].cast<py::tuple>();
    poly.vertices.push_back(Vector2D(pt[0].cast<double>(), pt[1].cast<double>()));
  }
  if (holes.size() > 0)
  {
    for (size_t i = 0; i < holes.size(); i++)
    {
      auto hl = holes[i].cast<py::list>();
      std::vector<Vector2D> hole;
      for (size_t j = 0; j < hl.size(); j++)
      {
        auto pt = hl[j].cast<py::tuple>();
        hole.push_back(Vector2D(pt[0].cast<double>(), pt[1].cast<double>()));
      }
      poly.holes.push_back(hole);
    }
  }
  return poly;
}

Mesh create_mesh(py::array_t<double> vertices, py::array_t<size_t> faces, py::array_t<int> markers)
{
  Mesh mesh;
  auto verts_r = vertices.unchecked<2>();
  auto faces_r = faces.unchecked<2>();
  auto markers_r = markers.unchecked<1>();
  size_t num_vertices = verts_r.shape(0);
  size_t num_faces = faces_r.shape(0);
  size_t num_markers = markers_r.size();

  for (size_t i = 0; i < num_vertices; i++)
  {
    mesh.vertices.push_back(Vector3D(verts_r(i, 0), verts_r(i, 1), verts_r(i, 2)));
  }

  for (size_t i = 0; i < num_faces; i++)
  {
    mesh.faces.push_back(Simplex2D(faces_r(i, 0), faces_r(i, 1), faces_r(i, 2)));
  }

  for (size_t i = 0; i < num_markers; i++)
  {
    mesh.markers.push_back(markers_r(i));
  }

  return mesh;
}

VolumeMesh create_volume_mesh(py::array_t<double> vertices, py::array_t<size_t> cells,
                              py::array_t<int> markers)
{
  VolumeMesh mesh;
  auto verts_r = vertices.unchecked<2>();
  auto cells_r = cells.unchecked<2>();
  auto markers_r = markers.unchecked<1>();
  size_t num_vertices = verts_r.shape(0);
  size_t num_cells = cells_r.shape(0);
  size_t num_markers = markers_r.size();

  for (size_t i = 0; i < num_vertices; i++)
  {
    mesh.vertices.push_back(Vector3D(verts_r(i, 0), verts_r(i, 1), verts_r(i, 2)));
  }

  for (size_t i = 0; i < num_cells; i++)
  {
    mesh.cells.push_back(Simplex3D(cells_r(i, 0), cells_r(i, 1), cells_r(i, 2), cells_r(i, 3)));
  }

  for (size_t i = 0; i < num_markers; i++)
  {
    mesh.markers.push_back(markers_r(i));
  }

  return mesh;
}

py::tuple mesh_as_arrays(const Mesh &mesh)
{
  py::array_t<double> py_vertices(mesh.vertices.size() * 3);
  py::array_t<size_t> py_faces(mesh.faces.size() * 3);
  py::array_t<int> py_markers(mesh.markers.size());
  for (size_t i = 0; i < mesh.vertices.size(); i++)
  {
    py_vertices.mutable_at(i * 3) = mesh.vertices[i].x;
    py_vertices.mutable_at(i * 3 + 1) = mesh.vertices[i].y;
    py_vertices.mutable_at(i * 3 + 2) = mesh.vertices[i].z;
  }
  for (size_t i = 0; i < mesh.faces.size(); i++)
  {
    py_faces.mutable_at(i * 3) = mesh.faces[i].v0;
    py_faces.mutable_at(i * 3 + 1) = mesh.faces[i].v1;
    py_faces.mutable_at(i * 3 + 2) = mesh.faces[i].v2;
  }
  for (size_t i = 0; i < mesh.markers.size(); i++)
  {
    py_markers.mutable_at(i) = mesh.markers[i];
  }
  return py::make_tuple(py_vertices, py_faces, py_markers);
}

Surface create_surface(py::array_t<double> vertices, py::list holes)
{
  Surface surface;
  auto verts_r = vertices.unchecked<2>();
  size_t num_vertices = verts_r.shape(0);
  for (size_t i = 0; i < num_vertices; i++)
  {
    surface.vertices.push_back(Vector3D(verts_r(i, 0), verts_r(i, 1), verts_r(i, 2)));
  }
  for (size_t i = 0; i < holes.size(); i++)
  {
    auto hole = holes[i].cast<py::array_t<double>>();
    auto hole_r = hole.unchecked<2>();
    std::vector<Vector3D> hole_vertices;
    size_t num_hole_vertices = hole_r.shape(0);
    for (size_t j = 0; j < num_hole_vertices; j++)
    {
      hole_vertices.push_back(Vector3D(hole_r(j, 0), hole_r(j, 1), hole_r(j, 2)));
    }
    surface.holes.push_back(hole_vertices);
  }
  return surface;
}

py::list points_in_polygons(const py::array_t<double> &pts, const std::vector<Polygon> &polygons)
{
  py::list in_polygons;
  auto pts_r = pts.unchecked<2>();
  size_t pt_count = pts_r.shape(0);
  std::vector<Vector3D> pc;
  for (size_t i = 0; i < pt_count; i++)
  {
    pc.push_back(Vector3D(pts_r(i, 0), pts_r(i, 1), pts_r(i, 2)));
  }
  auto pips = PointCloudProcessor::points_in_polygons(pc, polygons);
  py::list in_polygons_list;
  for (auto const &pip : pips)
  {
    py::array_t<size_t> indices(pip.size());
    for (size_t i = 0; i < pip.size(); i++)
    {
      indices.mutable_at(i) = pip[i];
    }
    in_polygons_list.append(indices);
  }

  return in_polygons_list;
}

py::list extract_building_points(std::vector<Polygon> &buildings, const py::array_t<double> &pts,
                                 bool statistical_outlier_remover, size_t neighbors,
                                 double outlier_margin)
{
  py::list roof_points;
  auto pts_r = pts.unchecked<2>();
  size_t pt_count = pts_r.shape(0);
  std::vector<Vector3D> pc;
  for (size_t i = 0; i < pt_count; i++)
  {
    pc.push_back(Vector3D(pts_r(i, 0), pts_r(i, 1), pts_r(i, 2)));
  }

  auto _roof_points = BuildingProcessor::extract_building_points(buildings, pc);
  if (statistical_outlier_remover)
  {
    for (auto &rp : _roof_points)
    {
      PointCloudProcessor::statistical_outlier_remover(rp, neighbors, outlier_margin);
    }
  }
  for (auto const &rp : _roof_points)
  {
    py::array_t<double> pts(rp.size() * 3);
    for (size_t i = 0; i < rp.size(); i++)
    {
      pts.mutable_at(i * 3) = rp[i].x;
      pts.mutable_at(i * 3 + 1) = rp[i].y;
      pts.mutable_at(i * 3 + 2) = rp[i].z;
    }
    pts = pts.reshape(std::vector<long>{static_cast<long>(rp.size()), 3});
    roof_points.append(pts);
  }
  return roof_points;
}

MultiSurface create_multisurface(py::list surfaces)
{
  MultiSurface multi_surface;
  for (size_t i = 0; i < surfaces.size(); i++)
  {
    auto surface = surfaces[i].cast<Surface>();
    multi_surface.surfaces.push_back(surface);
  }
  return multi_surface;
}

GridField create_gridfield(py::array_t<double> data, py::tuple bounds, size_t xsize, size_t ysize)
{
  GridField grid_field;
  double px = bounds[0].cast<double>();
  double py = bounds[1].cast<double>();
  double qx = bounds[2].cast<double>();
  double qy = bounds[3].cast<double>();
  auto bbox = BoundingBox2D(Vector2D(px, py), Vector2D(qx, qy));

  grid_field.grid.bounding_box = bbox;
  grid_field.grid.xstep = (qx - px) / xsize;
  grid_field.grid.ystep = (qy - py) / ysize;

  grid_field.grid.xsize = xsize;
  grid_field.grid.ysize = ysize;

  auto data_r = data.unchecked<1>();
  size_t data_count = data_r.size();

  for (size_t i = 0; i < data_count; i++)
  {
    grid_field.values.push_back(data_r(i));
  }

  return grid_field;
}

terrain_mesher::core::RasterDouble gridfield_to_raster(const GridField &grid_field)
{
  if (grid_field.grid.xsize == 0 || grid_field.grid.ysize == 0)
    error("build_terrain_mesh_zemlya: GridField has empty grid dimensions");

  const size_t expected_size = grid_field.grid.xsize * grid_field.grid.ysize;
  if (grid_field.values.size() != expected_size)
  {
    error("build_terrain_mesh_zemlya: GridField values size (" + str(grid_field.values.size()) +
          ") does not match grid dimensions (" + str(expected_size) + ")");
  }

  const double xstep = grid_field.grid.xstep;
  const double ystep = grid_field.grid.ystep;
  const double tolerance = std::max({1.0, std::fabs(xstep), std::fabs(ystep)}) * 1e-12;
  if (std::fabs(xstep - ystep) > tolerance)
  {
    error("build_terrain_mesh_zemlya: Zemlya requires equal x/y spacing, got xstep=" +
          str(xstep) + " and ystep=" + str(ystep));
  }

  terrain_mesher::core::RasterDouble raster(grid_field.grid.xsize, grid_field.grid.ysize);
  raster.set_cell_size(xstep);
  raster.set_pos_x(grid_field.grid.bounding_box.P.x - 0.5 * xstep);
  raster.set_pos_y(grid_field.grid.bounding_box.P.y - 0.5 * ystep);

  for (size_t row = 0; row < grid_field.grid.ysize; row++)
  {
    const size_t source_row = grid_field.grid.ysize - 1 - row;
    for (size_t col = 0; col < grid_field.grid.xsize; col++)
    {
      raster.value(row, col) = grid_field.values[source_row * grid_field.grid.xsize + col];
    }
  }

  return raster;
}

Mesh build_terrain_mesh_zemlya(const GridField &grid_field, double max_error,size_t smoothing_iterations)
{
  auto raster = gridfield_to_raster(grid_field);
  return terrain_mesher::core::generate_zemlya_mesh(std::move(raster), max_error, smoothing_iterations);
}

py::array_t<double> ray_surface_intersection(const Surface &surface,
                                             const py::array_t<double> &py_ray_origin,
                                             const py::array_t<double> &py_ray_vector)
{

  Vector3D ray_origin(py_ray_origin.at(0), py_ray_origin.at(1), py_ray_origin.at(2));
  Vector3D ray_vector(py_ray_vector.at(0), py_ray_vector.at(1), py_ray_vector.at(2));

  auto intersection = Intersection::ray_surface_intersection(surface, ray_origin, ray_vector);

  return py::array_t<double>(3, &intersection.x);
}

py::array_t<double> ray_multisurface_intersection(const MultiSurface &surface,
                                                  const py::array_t<double> &py_ray_origin,
                                                  const py::array_t<double> &py_ray_vector)
{

  Vector3D ray_origin(py_ray_origin.at(0), py_ray_origin.at(1), py_ray_origin.at(2));
  Vector3D ray_vector(py_ray_vector.at(0), py_ray_vector.at(1), py_ray_vector.at(2));

  auto intersection = Intersection::ray_multisurface_intersection(surface, ray_origin, ray_vector);

  return py::array_t<double>(3, &intersection.x);
}

py::array_t<size_t> statistical_outlier_finder(py::array_t<double> &points, size_t neighbors,
                                               double outlier_margin)
{
  auto points_r = points.unchecked<2>();
  size_t num_points = points_r.shape(0);
  std::vector<Vector3D> pc;
  for (size_t i = 0; i < num_points; i++)
  {
    pc.push_back(Vector3D(points_r(i, 0), points_r(i, 1), points_r(i, 2)));
  }

  auto outliers = PointCloudProcessor::statistical_outlier_finder(pc, neighbors, outlier_margin);
  return py::array_t<size_t>(outliers.size(), outliers.data());
}

py::dict compute_boundary_face_markers(const VolumeMesh &mesh)
{
  auto data = MeshProcessor::compute_boundary_facet_markers(mesh);
  py::dict out;
  // out.reserve(data.size());
  for (auto const &kv : data)
  {
    const Simplex2D &f = kv.first;
    int marker = kv.second.first;
    // auto  &n     = kv.second.second;

    py::tuple key = py::make_tuple(f.v0, f.v1, f.v2);
    // py::tuple normal = py::make_tuple(n.x, n.y, n.z);
    // py::tuple val    = py::make_tuple(marker, normal);
    out[key] = marker;
  }
  return out;
}

} // namespace DTCC_BUILDER

PYBIND11_MODULE(_dtcc_builder, m)
{

#ifdef DTCC_HAVE_TRIANGLE
  constexpr bool have_triangle = true;
#else
  constexpr bool have_triangle = false;
#endif

#ifdef DTCC_HAVE_SPADE
  constexpr bool have_spade = true;
#else
  constexpr bool have_spade = false;
#endif

  m.attr("HAVE_TRIANGLE") = py::bool_(have_triangle);
  m.attr("HAVE_SPADE") = py::bool_(have_spade);

  m.def(
      "triangulation_backends",
      []()
      {
        std::vector<std::string> backends;
        if (have_triangle)
          backends.emplace_back("triangle");
        if (have_spade)
          backends.emplace_back("spade");
        if (backends.empty())
          backends.emplace_back("earcut");
        return backends;
      },
      R"pbdoc(
            Report available triangulation backends compiled into the extension.
        )pbdoc");

  py::class_<DTCC_BUILDER::Vector2D>(m, "Vector2D")
      .def(py::init<>())
      .def(
          "__repr__", [](const DTCC_BUILDER::Vector3D &p)
          { return "<Vector3D (" + DTCC_BUILDER::str(p.x) + ", " + DTCC_BUILDER::str(p.y) + ")>"; })
      .def_readonly("x", &DTCC_BUILDER::Vector2D::x)
      .def_readonly("y", &DTCC_BUILDER::Vector2D::y);

  py::class_<DTCC_BUILDER::Vector3D>(m, "Vector3D")
      .def(py::init<>())
      .def("__repr__",
           [](const DTCC_BUILDER::Vector3D &p)
           {
             return "<Vector3D (" + DTCC_BUILDER::str(p.x) + ", " + DTCC_BUILDER::str(p.y) + ", " +
                    DTCC_BUILDER::str(p.z) + ")>";
           })
      .def_readonly("x", &DTCC_BUILDER::Vector3D::x)
      .def_readonly("y", &DTCC_BUILDER::Vector3D::y)
      .def_readonly("z", &DTCC_BUILDER::Vector3D::z);

  // py::class_<DTCC_BUILDER::Vector3D>(m, "Vector3D")
  //     .def(py::init<>())
  //     .def_readonly("x", &DTCC_BUILDER::Vector3D::x)
  //     .def_readonly("y", &DTCC_BUILDER::Vector3D::y)
  //     .def_readonly("z", &DTCC_BUILDER::Vector3D::z);

  py::class_<DTCC_BUILDER::BoundingBox2D>(m, "bounding_box")
      .def(py::init<>())
      .def_readonly("P", &DTCC_BUILDER::BoundingBox2D::P)
      .def_readonly("Q", &DTCC_BUILDER::BoundingBox2D::Q);

  py::class_<DTCC_BUILDER::Polygon>(m, "Polygon")
      .def(py::init<>())
      .def_readonly("vertices", &DTCC_BUILDER::Polygon::vertices)
      .def_readonly("holes", &DTCC_BUILDER::Polygon::holes);

  py::class_<DTCC_BUILDER::GridField>(m, "GridField")
      .def(py::init<>())
      .def_readonly("grid", &DTCC_BUILDER::GridField::grid)
      .def_readonly("values", &DTCC_BUILDER::GridField::values);

  py::class_<DTCC_BUILDER::Grid>(m, "Grid")
      .def(py::init<>())
      .def_readonly("xsize", &DTCC_BUILDER::Grid::xsize)
      .def_readonly("ysize", &DTCC_BUILDER::Grid::ysize)
      .def_readonly("xstep", &DTCC_BUILDER::Grid::xstep)
      .def_readonly("ystep", &DTCC_BUILDER::Grid::ystep);

  py::class_<DTCC_BUILDER::Simplex2D>(m, "Simplex2D")
      .def(py::init<>())
      .def_readonly("v0", &DTCC_BUILDER::Simplex2D::v0)
      .def_readonly("v1", &DTCC_BUILDER::Simplex2D::v1)
      .def_readonly("v2", &DTCC_BUILDER::Simplex2D::v2);

  py::class_<DTCC_BUILDER::Simplex3D>(m, "Simplex3D")
      .def(py::init<>())
      .def_readonly("v0", &DTCC_BUILDER::Simplex3D::v0)
      .def_readonly("v1", &DTCC_BUILDER::Simplex3D::v1)
      .def_readonly("v2", &DTCC_BUILDER::Simplex3D::v2)
      .def_readonly("v3", &DTCC_BUILDER::Simplex3D::v3);

  py::class_<DTCC_BUILDER::Mesh>(m, "Mesh")
      .def(py::init<>())
      .def_readonly("vertices", &DTCC_BUILDER::Mesh::vertices)
      .def_readonly("faces", &DTCC_BUILDER::Mesh::faces)
      .def_readonly("normals", &DTCC_BUILDER::Mesh::normals)
      .def_readonly("markers", &DTCC_BUILDER::Mesh::markers)
      .def(
          "from_cpp",
          [](const DTCC_BUILDER::Mesh &m)
          {
            py::object conv = py::module::import("dtcc_core.builder.model_conversion")
                                  .attr("builder_mesh_to_mesh");
            return conv(m);
          },
          R"pbdoc(
                 Convert this C++ Mesh back into a Python model.Mesh.
             )pbdoc");
  ;

  py::class_<DTCC_BUILDER::VolumeMesh>(m, "VolumeMesh")
      .def(py::init<>())
      .def_readonly("num_layers", &DTCC_BUILDER::VolumeMesh::num_layers)
      .def_readonly("vertices", &DTCC_BUILDER::VolumeMesh::vertices)
      .def_readonly("cells", &DTCC_BUILDER::VolumeMesh::cells)
      .def_readonly("markers", &DTCC_BUILDER::VolumeMesh::markers)
      .def(
          "from_cpp",
          [](const DTCC_BUILDER::VolumeMesh &m)
          {
            py::object conv = py::module::import("dtcc_core.builder.model_conversion")
                                  .attr("builder_volume_mesh_to_volume_mesh");
            return conv(m);
          },
          R"pbdoc(
                 Convert this C++ Volume Mesh back into a Python model.VolumeMesh.
             )pbdoc");

  py::class_<DTCC_BUILDER::Surface>(m, "Surface")
      .def(py::init<>())
      .def_readonly("vertices", &DTCC_BUILDER::Surface::vertices)
      .def_readonly("holes", &DTCC_BUILDER::Surface::holes);

  py::class_<DTCC_BUILDER::MultiSurface>(m, "MultiSurface")
      .def(py::init<>())
      .def_readonly("surfaces", &DTCC_BUILDER::MultiSurface::surfaces);

  m.def("create_polygon", &DTCC_BUILDER::create_polygon, "Create C++ polygon");

  m.def("create_mesh", &DTCC_BUILDER::create_mesh, "Create C++ mesh");

  m.def("create_volume_mesh", &DTCC_BUILDER::create_volume_mesh, "Create C++ volume mesh");

  m.def("mesh_as_arrays", &DTCC_BUILDER::mesh_as_arrays, "Create C++ mesh");

  m.def("create_gridfield", &DTCC_BUILDER::create_gridfield, "Create C++ grid field");

  m.def("build_terrain_mesh_zemlya", &DTCC_BUILDER::build_terrain_mesh_zemlya,
        "Build a terrain mesh from a GridField using the Zemlya terrain mesher");

  m.def("extract_building_points", &DTCC_BUILDER::extract_building_points,
        "Compute building points from point cloud");

  m.def("points_in_polygons", &DTCC_BUILDER::points_in_polygons, "Find points inside polygons");

  m.def("boundary_stats", &DTCC_BUILDER::boundary_stats,
        "Compute lightweight polygon boundary statistics for cleaner runtime shortcuts");

  m.def("boundary_defect_clusters", &DTCC_BUILDER::boundary_defect_clusters,
        "Detect short-edge and pair-defect clusters for the cleaner");

  m.def("rewrite_defect_cluster", &DTCC_BUILDER::rewrite_defect_cluster,
        "Apply deterministic short-edge boundary rewrites for simple cleaner clusters");

  m.def("smooth_field", &DTCC_BUILDER::VertexSmoother::smooth_field, "Smooth grid field");

  // m.def("build_mesh", &DTCC_BUILDER::MeshBuilder::build_mesh,
  //       "build mesh for city, returning a list of meshes");

  m.def("build_city_flat_mesh", &DTCC_BUILDER::MeshBuilder::build_city_flat_mesh,
        "build city flat mesh");

  m.def("build_terrain_surface_mesh", &DTCC_BUILDER::MeshBuilder::build_terrain_surface_mesh,
        "build terrain surface mesh");
  m.def("build_terrain_surface_mesh_from_ground_mesh",
        &DTCC_BUILDER::MeshBuilder::build_terrain_surface_mesh_from_ground_mesh,
        "build terrain surface mesh from a prebuilt ground mesh");

  m.def("build_city_surface_mesh", &DTCC_BUILDER::MeshBuilder::build_city_surface_mesh,
        "build city surface mesh");
  m.def("build_city_surface_mesh_from_terrain_mesh",
        &DTCC_BUILDER::MeshBuilder::build_city_surface_mesh_from_terrain_mesh,
        "build city surface mesh from a prebuilt terrain mesh");

  m.def("layer_ground_mesh", &DTCC_BUILDER::MeshBuilder::layer_ground_mesh, "Layer ground mesh");

  m.def("smooth_volume_mesh", &DTCC_BUILDER::Smoother::smooth_volume_mesh, "Smooth volume mesh");

  m.def("trim_volume_mesh", &DTCC_BUILDER::MeshBuilder::trim_volume_mesh,
        "Trim volume mesh by removing cells inside buildings");

  // m.def("extrude_footprint", &DTCC_BUILDER::MeshBuilder::extrude_footprint,
  //       "Extrude footprint to a mesh");

  m.def("compute_boundary_mesh", &DTCC_BUILDER::MeshProcessor::compute_boundary_mesh,
        "Compute boundary mesh from volume mesh");

  m.def("compute_boundary_face_markers", &DTCC_BUILDER::compute_boundary_face_markers,
        "Compute markers and outward normals for volume mesh boundary faces");

  m.def("compute_boundary_mesh", &DTCC_BUILDER::MeshProcessor::compute_boundary_mesh,
        "Compute boundary mesh from volume mesh");

  m.def("compute_open_mesh", &DTCC_BUILDER::MeshProcessor::compute_open_mesh,
        "Compute open mesh from boundary, excluding top and sides");

  m.def("merge_meshes", &DTCC_BUILDER::MeshProcessor::merge_meshes,
        "Merge meshes into a single mesh");

  m.def("snap_mesh_vertices", &DTCC_BUILDER::MeshProcessor::snap_vertices, "Snap mesh vertices");

  m.def("create_surface", &DTCC_BUILDER::create_surface, "Create C++ surface");

  m.def("create_multisurface", &DTCC_BUILDER::create_multisurface, "Create C++ multisurface");

  m.def("mesh_surface", &DTCC_BUILDER::MeshBuilder::mesh_surface,
        "Create triangulated mesh from surface");

  m.def("mesh_multisurface", &DTCC_BUILDER::MeshBuilder::mesh_multisurface,
        "Create triangulated mesh from multisurface");
  m.def("mesh_multisurfaces", &DTCC_BUILDER::MeshBuilder::mesh_multisurfaces,
        "Create a lits of triangulated meshes from a list of multisurfaces");

  m.def("ray_surface_intersection", &DTCC_BUILDER::ray_surface_intersection,
        "Compute ray-surface intersection");

  m.def("ray_multisurface_intersection", &DTCC_BUILDER::ray_multisurface_intersection,
        "Compute ray-multisurface intersection");

  m.def("statistical_outlier_finder", &DTCC_BUILDER::statistical_outlier_finder,
        "Find statistical outliers in point cloud");

  py::class_<DTCC_BUILDER::VolumeMeshBuilder>(m, "VolumeMeshBuilder")
      .def(py::init<const std::vector<DTCC_BUILDER::Surface> &, const DTCC_BUILDER::GridField &,
                    DTCC_BUILDER::Mesh &, double>(),
           py::arg("buildings"), py::arg("dem"), py::arg("ground_mesh"), py::arg("domain_height"),
           "Constructor for VolumeMeshBuilder taking city, dem, ground_mesh, "
           "and domain_height as arguments")
      .def("build", &DTCC_BUILDER::VolumeMeshBuilder::build,
           "Layers the ground mesh and returns a VolumeMesh")
      // Expose public variables directly
      .def_readwrite("domain_height", &DTCC_BUILDER::VolumeMeshBuilder::domain_height)
      .def_readwrite("top_height", &DTCC_BUILDER::VolumeMeshBuilder::top_height)
      // If you need to expose std::vectors or similar, pybind11/stl.h header
      // takes care of this. For custom types like City, GridField, Mesh, ensure
      // you've also provided bindings for them.
      ;
}
