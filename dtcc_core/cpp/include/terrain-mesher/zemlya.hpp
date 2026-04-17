#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "model/Mesh.h"
#include "terrain-mesher/geometry.hpp"
#include "terrain-mesher/raster.hpp"
#include "terrain-mesher/topology.hpp"

namespace terrain_mesher::core {

namespace detail {

struct Candidate {
    int x = 0;
    int y = 0;
    double z = 0.0;
    double importance = -std::numeric_limits<double>::max();
    int token = 0;
    dt_ptr triangle;

    void consider(int sample_x, int sample_y, double sample_z, double sample_importance) noexcept {
        if (sample_importance > importance) {
            x = sample_x;
            y = sample_y;
            z = sample_z;
            importance = sample_importance;
        }
    }

    bool operator<(const Candidate& other) const noexcept {
        return importance < other.importance;
    }
}; 

class CandidateList {
  public:
    void push_back(const Candidate& candidate) { candidates_.push(candidate); }
    bool empty() const noexcept { return candidates_.empty(); }

    Candidate grab_greatest() {
        if (candidates_.empty()) {
            return Candidate{};
        }
        Candidate result = candidates_.top();
        candidates_.pop();
        return result;
    }

  private:
    std::priority_queue<Candidate> candidates_;
};

void order_triangle_points(std::array<Point2, 3>& points) noexcept;
void compute_plane(Plane& plane, dt_ptr triangle, const RasterDouble& raster);

class TerraBaseMesh : protected DelaunayMesh {
  public:
    void load_raster(RasterDouble raster);

  protected:
    RasterDouble raster_;
    void repair_point(int px, int py);
};

inline int pow2_int(int exponent) {
    return 1 << exponent;
}

inline double average_of(double d1, double d2, double d3, double d4, double no_data_value) {
    int count = 0;
    double sum = 0.0;
    for (double value : {d1, d2, d3, d4}) {
        if (is_no_data(value, no_data_value)) {
            continue;
        }
        ++count;
        sum += value;
    }
    return count > 0 ? sum / static_cast<double>(count) : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace detail

class ZemlyaMesh : public detail::TerraBaseMesh {
  public:
    void greedy_insert(double max_error);
    void scan_triangle(dt_ptr triangle) override;
    DTCC_BUILDER::Mesh convert_to_mesh() const;

  private:
    void scan_triangle_line(const Plane& plane,
                            int y,
                            double x1,
                            double x2,
                            detail::Candidate& candidate,
                            double no_data_value);

    RasterDouble sample_;
    RasterDouble insert_;
    RasterDouble result_;
    Raster<char> used_;
    Raster<int> token_;

    detail::CandidateList candidates_;
    double max_error_ = 0.0;
    int counter_ = 0;
    int current_level_ = 0;
    int max_level_ = 0;
};

DTCC_BUILDER::Mesh generate_zemlya_mesh(RasterDouble raster, double max_error);

inline void detail::order_triangle_points(std::array<Point2, 3>& points) noexcept {
    if (points[0].y > points[1].y) {
        std::swap(points[0], points[1]);
    }
    if (points[1].y > points[2].y) {
        std::swap(points[1], points[2]);
    }
    if (points[0].y > points[1].y) {
        std::swap(points[0], points[1]);
    }
}

inline void detail::compute_plane(Plane& plane, dt_ptr triangle, const RasterDouble& raster) {
    const Point2 p1 = triangle->point1();
    const Point2 p2 = triangle->point2();
    const Point2 p3 = triangle->point3();

    Point3 v1{
        p1.x,
        p1.y,
        raster.value(static_cast<std::size_t>(p1.y), static_cast<std::size_t>(p1.x)),
    };
    Point3 v2{
        p2.x,
        p2.y,
        raster.value(static_cast<std::size_t>(p2.y), static_cast<std::size_t>(p2.x)),
    };
    Point3 v3{
        p3.x,
        p3.y,
        raster.value(static_cast<std::size_t>(p3.y), static_cast<std::size_t>(p3.x)),
    };
    plane.init(v1, v2, v3);
}

inline void detail::TerraBaseMesh::load_raster(RasterDouble raster) {
    raster_ = std::move(raster);
}

inline void detail::TerraBaseMesh::repair_point(int px, int py) {
    double& point = raster_.value(static_cast<std::size_t>(py), static_cast<std::size_t>(px));
    const double value = sample_nearest_valid_avg(
        raster_,
        static_cast<std::size_t>(py),
        static_cast<std::size_t>(px));
    if (is_no_data(value, raster_.get_no_data_value())) {
        point = 0.0;
    } else {
        point = value;
    }
}

inline void ZemlyaMesh::greedy_insert(double max_error) {
    max_error_ = max_error;
    counter_ = 0;
    const int width = static_cast<int>(raster_.get_width());
    const int height = static_cast<int>(raster_.get_height());
    max_level_ = static_cast<int>(std::ceil(std::log2(static_cast<double>(std::max(width, height)))));

    sample_.allocate(raster_.get_width(), raster_.get_height());
    sample_.set_all(std::numeric_limits<double>::quiet_NaN());

    const double no_data_value = raster_.get_no_data_value();

    for (int level = max_level_ - 1; level >= 1; --level) {
        const int step = max_level_ - level;
        const int stride = detail::pow2_int(step);

        for (int y = 0; y < height; y += stride) {
            for (int x = 0; x < width; x += stride) {
                if (step == 1) {
                    const double v1 = (y < height && x < width)
                        ? raster_.value(y, x)
                        : std::numeric_limits<double>::quiet_NaN();
                    const double v2 = (y < height && x + 1 < width)
                        ? raster_.value(y, x + 1)
                        : std::numeric_limits<double>::quiet_NaN();
                    const double v3 = (y + 1 < height && x < width)
                        ? raster_.value(y + 1, x)
                        : std::numeric_limits<double>::quiet_NaN();
                    const double v4 = (y + 1 < height && x + 1 < width)
                        ? raster_.value(y + 1, x + 1)
                        : std::numeric_limits<double>::quiet_NaN();

                    if (y + 1 < height && x + 1 < width) {
                        sample_.value(y + 1, x + 1) =
                            detail::average_of(v1, v2, v3, v4, no_data_value);
                    }
                } else {
                    const int offset = detail::pow2_int(step - 1);
                    const int delta = detail::pow2_int(step - 2);

                    auto sample_value = [&](int row, int column) {
                        if (row < 0 || column < 0 || row >= height || column >= width) {
                            return std::numeric_limits<double>::quiet_NaN();
                        }
                        return sample_.value(
                            static_cast<std::size_t>(row),
                            static_cast<std::size_t>(column));
                    };

                    const double v1 = sample_value(y + offset - delta, x + offset - delta);
                    const double v2 = sample_value(y + offset - delta, x + offset + delta);
                    const double v3 = sample_value(y + offset + delta, x + offset - delta);
                    const double v4 = sample_value(y + offset + delta, x + offset + delta);

                    if (y + offset < height && x + offset < width) {
                        sample_.value(y + offset, x + offset) =
                            detail::average_of(v1, v2, v3, v4, no_data_value);
                    }
                }
            }
        }
    }

    repair_point(0, 0);
    repair_point(0, height - 1);
    repair_point(width - 1, height - 1);
    repair_point(width - 1, 0);

    result_.allocate(raster_.get_width(), raster_.get_height());
    result_.set_all(std::numeric_limits<double>::quiet_NaN());
    result_.value(0, 0) = raster_.value(0, 0);
    result_.value(height - 1, 0) = raster_.value(height - 1, 0);
    result_.value(height - 1, width - 1) = raster_.value(height - 1, width - 1);
    result_.value(0, width - 1) = raster_.value(0, width - 1);

    insert_.allocate(raster_.get_width(), raster_.get_height());
    insert_.set_all(std::numeric_limits<double>::quiet_NaN());

    used_.allocate(raster_.get_width(), raster_.get_height());

    token_.allocate(raster_.get_width(), raster_.get_height());
    token_.set_all(0);

    init_mesh(
        Point2{0.0, 0.0},
        Point2{0.0, static_cast<double>(height - 1)},
        Point2{static_cast<double>(width - 1), static_cast<double>(height - 1)},
        Point2{static_cast<double>(width - 1), 0.0});

    for (int level = 1; level <= max_level_; ++level) {
        current_level_ = level;
        used_.set_all(0);

        if (level >= 5 && level <= max_level_ - 1) {
            const int step = max_level_ - level;

            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    double& z = insert_.value(y, x);
                    if (is_no_data(z, no_data_value)) {
                        continue;
                    }
                    z = raster_.value(y, x);
                }
            }

            const int stride = detail::pow2_int(step);
            const int offset = detail::pow2_int(step - 1);
            for (int y = 0; y < height; y += stride) {
                for (int x = 0; x < width; x += stride) {
                    if (y + offset < height && x + offset < width) {
                        insert_.value(y + offset, x + offset) =
                            raster_.value(y + offset, x + offset);
                    }
                }
            }
        } else if (level < max_level_) {
            const int step = max_level_ - level;

            if (step >= 3) {
                const int delta = detail::pow2_int(step - 3);

                for (int y = 0; y < height; ++y) {
                    for (int x = 0; x < width; ++x) {
                        double& z = insert_.value(y, x);
                        if (is_no_data(z, no_data_value)) {
                            continue;
                        }

                        auto sample_value = [&](int row, int column) {
                            if (row < 0 || column < 0 || row >= height || column >= width) {
                                return std::numeric_limits<double>::quiet_NaN();
                            }
                            return sample_.value(
                                static_cast<std::size_t>(row),
                                static_cast<std::size_t>(column));
                        };

                        const double v1 = sample_value(y - delta, x - delta);
                        const double v2 = sample_value(y - delta, x + delta);
                        const double v3 = sample_value(y + delta, x - delta);
                        const double v4 = sample_value(y + delta, x + delta);
                        const double average = detail::average_of(v1, v2, v3, v4, no_data_value);
                        if (!is_no_data(average, no_data_value)) {
                            z = average;
                        }
                    }
                }
            }

            const int stride = detail::pow2_int(step);
            const int offset = detail::pow2_int(step - 1);
            for (int y = 0; y < height; y += stride) {
                for (int x = 0; x < width; x += stride) {
                    if (y + offset < height && x + offset < width) {
                        insert_.value(y + offset, x + offset) =
                            sample_.value(y + offset, x + offset);
                    }
                }
            }
        }

        dt_ptr triangle = first_face_;
        while (triangle) {
            scan_triangle(triangle);
            triangle = triangle->get_link();
        }

        while (!candidates_.empty()) {
            detail::Candidate candidate = candidates_.grab_greatest();
            if (candidate.importance < max_error_) {
                continue;
            }
            if (token_.value(candidate.y, candidate.x) != candidate.token) {
                continue;
            }

            result_.value(candidate.y, candidate.x) = candidate.z;
            used_.value(candidate.y, candidate.x) = 1;
            insert(
                Point2{static_cast<double>(candidate.x), static_cast<double>(candidate.y)},
                candidate.triangle);
        }
    }
}

inline void ZemlyaMesh::scan_triangle_line(const Plane& plane,
                                           int y,
                                           double x1,
                                           double x2,
                                           detail::Candidate& candidate,
                                           double no_data_value) {
    const int start_x = static_cast<int>(std::ceil(std::min(x1, x2)));
    const int end_x = static_cast<int>(std::floor(std::max(x1, x2)));
    if (start_x > end_x) {
        return;
    }

    double z0 = plane.eval(static_cast<double>(start_x), static_cast<double>(y));
    const double dz = plane.a;

    for (int x = start_x; x <= end_x; ++x) {
        if (!used_.value(y, x)) {
            const double z =
                current_level_ == max_level_ ? raster_.value(y, x) : insert_.value(y, x);
            if (!is_no_data(z, no_data_value)) {
                candidate.consider(x, y, z, std::fabs(z - z0));
            }
        }
        z0 += dz;
    }
}

inline void ZemlyaMesh::scan_triangle(dt_ptr triangle) {
    Plane plane;
    detail::compute_plane(plane, triangle, result_);

    std::array<Point2, 3> by_y = {triangle->point1(), triangle->point2(), triangle->point3()};
    detail::order_triangle_points(by_y);

    const double v0_x = by_y[0].x;
    const double v0_y = by_y[0].y;
    const double v1_x = by_y[1].x;
    const double v1_y = by_y[1].y;
    const double v2_x = by_y[2].x;
    const double v2_y = by_y[2].y;

    detail::Candidate candidate;
    candidate.importance = -std::numeric_limits<double>::max();
    candidate.token = counter_++;
    candidate.triangle = triangle;

    const double dx2 = (v2_x - v0_x) / (v2_y - v0_y);
    const double no_data_value = raster_.get_no_data_value();

    if (v1_y != v0_y) {
        const double dx1 = (v1_x - v0_x) / (v1_y - v0_y);
        double x1 = v0_x;
        double x2 = v0_x;
        for (int y = static_cast<int>(v0_y); y < static_cast<int>(v1_y); ++y) {
            scan_triangle_line(plane, y, x1, x2, candidate, no_data_value);
            x1 += dx1;
            x2 += dx2;
        }
    }

    if (v2_y != v1_y) {
        const double dx1 = (v2_x - v1_x) / (v2_y - v1_y);
        double x1 = v1_x;
        double x2 = v0_x;
        for (int y = static_cast<int>(v1_y); y <= static_cast<int>(v2_y); ++y) {
            scan_triangle_line(plane, y, x1, x2, candidate, no_data_value);
            x1 += dx1;
            x2 += dx2;
        }
    }

    token_.value(candidate.y, candidate.x) = candidate.token;
    candidates_.push_back(candidate);
}

inline DTCC_BUILDER::Mesh ZemlyaMesh::convert_to_mesh() const {
    DTCC_BUILDER::Mesh mesh;
    Raster<std::int32_t> vertex_id(raster_.get_width(), raster_.get_height());
    vertex_id.set_all(-1);

    mesh.vertices.reserve(raster_.get_width() * raster_.get_height());

    for (std::size_t y = 0; y < raster_.get_height(); ++y) {
        for (std::size_t x = 0; x < raster_.get_width(); ++x) {
            const double z = result_.value(y, x);
            if (is_no_data(z, raster_.get_no_data_value())) {
                continue;
            }

            mesh.vertices.emplace_back(raster_.col2x(x), raster_.row2y(y), z);
            vertex_id.value(y, x) = static_cast<std::int32_t>(mesh.vertices.size() - 1);
        }
    }

    dt_ptr triangle = first_face_;
    while (triangle) {
        const Point2 p1 = triangle->point1();
        const Point2 p2 = triangle->point2();
        const Point2 p3 = triangle->point3();

        auto push_index = [&](const Point2& point) {
            return static_cast<std::size_t>(vertex_id.value(
                static_cast<std::size_t>(point.y),
                static_cast<std::size_t>(point.x)));
        };

        if (!ccw(p1, p2, p3)) {
            mesh.faces.emplace_back(push_index(p1), push_index(p2), push_index(p3));
        } else {
            mesh.faces.emplace_back(push_index(p3), push_index(p2), push_index(p1));
        }

        triangle = triangle->get_link();
    }

    return mesh;
}

inline DTCC_BUILDER::Mesh generate_zemlya_mesh(RasterDouble raster, double max_error) {
    ZemlyaMesh mesh;
    mesh.load_raster(std::move(raster));
    mesh.greedy_insert(max_error);
    return mesh.convert_to_mesh();
}

}  // namespace terrain_mesher::core
