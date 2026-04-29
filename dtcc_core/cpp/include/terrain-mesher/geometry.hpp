#pragma once

#include <cmath>
#include <cstddef>

namespace terrain_mesher::core {

constexpr double kEpsilon = 1e-6;

struct Point2 {
    double x = 0.0;
    double y = 0.0;

    double operator[](std::size_t index) const noexcept {
        return index == 0 ? x : y;
    }

    Point2 operator-(const Point2& other) const noexcept {
        return {x - other.x, y - other.y};
    }

    bool operator==(const Point2& other) const noexcept {
        return x == other.x && y == other.y;
    }

    bool operator!=(const Point2& other) const noexcept { return !(*this == other); }

    double length() const noexcept { return std::sqrt(x * x + y * y); }
};

struct Point3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Plane {
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;

    void init(const Point3& p, const Point3& q, const Point3& r) noexcept;
    double eval(double x_value, double y_value) const noexcept {
        return a * x_value + b * y_value + c;
    }
};

double tri_area(const Point2& a, const Point2& b, const Point2& c) noexcept;
bool ccw(const Point2& a, const Point2& b, const Point2& c) noexcept;
bool right_of(const Point2& point, const Point2& origin, const Point2& destination) noexcept;
bool left_of(const Point2& point, const Point2& origin, const Point2& destination) noexcept;
bool in_circle(const Point2& a, const Point2& b, const Point2& c, const Point2& d) noexcept;

class Line {
  public:
    Line(const Point2& p, const Point2& q);
    double eval(const Point2& p) const noexcept { return a_ * p.x + b_ * p.y + c_; }

  private:
    double a_ = 0.0;
    double b_ = 0.0;
    double c_ = 0.0;
};

inline void Plane::init(const Point3& p, const Point3& q, const Point3& r) noexcept {
    const double ux = q.x - p.x;
    const double uy = q.y - p.y;
    const double uz = q.z - p.z;

    const double vx = r.x - p.x;
    const double vy = r.y - p.y;
    const double vz = r.z - p.z;

    const double denominator = ux * vy - uy * vx;
    a = (uz * vy - uy * vz) / denominator;
    b = (ux * vz - uz * vx) / denominator;
    c = p.z - a * p.x - b * p.y;
}

inline double tri_area(const Point2& a, const Point2& b, const Point2& c) noexcept {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

inline bool ccw(const Point2& a, const Point2& b, const Point2& c) noexcept {
    return tri_area(a, b, c) > 0.0;
}

inline bool right_of(const Point2& point,
                     const Point2& origin,
                     const Point2& destination) noexcept {
    return ccw(point, destination, origin);
}

inline bool left_of(const Point2& point,
                    const Point2& origin,
                    const Point2& destination) noexcept {
    return ccw(point, origin, destination);
}

inline bool in_circle(const Point2& a,
                      const Point2& b,
                      const Point2& c,
                      const Point2& d) noexcept {
    return (a[0] * a[0] + a[1] * a[1]) * tri_area(b, c, d)
            - (b[0] * b[0] + b[1] * b[1]) * tri_area(a, c, d)
            + (c[0] * c[0] + c[1] * c[1]) * tri_area(a, b, d)
            - (d[0] * d[0] + d[1] * d[1]) * tri_area(a, b, c)
        > kEpsilon;
}

inline Line::Line(const Point2& p, const Point2& q) {
    const Point2 t = q - p;
    const double length = t.length();
    a_ = t.y / length;
    b_ = -t.x / length;
    c_ = -(a_ * p.x + b_ * p.y);
}

}  // namespace terrain_mesher::core
