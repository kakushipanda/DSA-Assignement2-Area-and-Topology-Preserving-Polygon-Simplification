// geometry.hpp
// Shared low-level geometry primitives used by the simplifier.
//
// This header intentionally stays lightweight: it defines the Point/Segment/Line
// types plus the small collection of vector, area, intersection, and coordinate
// frame helpers that the higher-level simplification code depends on.
#pragma once

#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace simplify {

constexpr double kEpsilon = 1e-9;

struct Point {
    double x = 0.0;
    double y = 0.0;
};

inline Point operator+(const Point& a, const Point& b) { return {a.x + b.x, a.y + b.y}; }
inline Point operator-(const Point& a, const Point& b) { return {a.x - b.x, a.y - b.y}; }
inline Point operator*(const Point& p, double s) { return {p.x * s, p.y * s}; }
inline Point operator/(const Point& p, double s) { return {p.x / s, p.y / s}; }

inline double dot(const Point& a, const Point& b) { return a.x * b.x + a.y * b.y; }
inline double cross(const Point& a, const Point& b) { return a.x * b.y - a.y * b.x; }
inline double cross(const Point& a, const Point& b, const Point& c) { return cross(b - a, c - a); }
inline double squared_length(const Point& p) { return dot(p, p); }
inline double length(const Point& p) { return std::sqrt(squared_length(p)); }
inline bool nearly_equal(double a, double b, double eps = kEpsilon) { return std::abs(a - b) <= eps; }
inline bool same_point(const Point& a, const Point& b, double eps = kEpsilon) {
    return std::abs(a.x - b.x) <= eps && std::abs(a.y - b.y) <= eps;
}

struct Segment {
    Point a;
    Point b;
};

struct Line {
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
};

// Returns the signed shoelace area of a ring. Counter-clockwise rings are
// positive and clockwise rings are negative.
inline double signed_ring_area(const std::vector<Point>& ring) {
    if (ring.size() < 3) {
        return 0.0;
    }
    double area = 0.0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const Point& p = ring[i];
        const Point& q = ring[(i + 1) % ring.size()];
        area += cross(p, q);
    }
    return 0.5 * area;
}

// Standard orientation predicate used for segment-intersection tests.
inline double orientation(const Point& a, const Point& b, const Point& c) {
    return cross(a, b, c);
}

// Returns true when p lies on the closed segment ab, within tolerance.
inline bool on_segment(const Point& a, const Point& b, const Point& p, double eps = kEpsilon) {
    if (std::abs(orientation(a, b, p)) > eps) {
        return false;
    }
    return p.x <= std::max(a.x, b.x) + eps && p.x + eps >= std::min(a.x, b.x) &&
           p.y <= std::max(a.y, b.y) + eps && p.y + eps >= std::min(a.y, b.y);
}

// General-purpose segment intersection test that accepts both proper crossings
// and endpoint/collinear touches.
inline bool segments_intersect(const Segment& s1, const Segment& s2, double eps = kEpsilon) {
    const double o1 = orientation(s1.a, s1.b, s2.a);
    const double o2 = orientation(s1.a, s1.b, s2.b);
    const double o3 = orientation(s2.a, s2.b, s1.a);
    const double o4 = orientation(s2.a, s2.b, s1.b);

    const bool proper = ((o1 > eps && o2 < -eps) || (o1 < -eps && o2 > eps)) &&
                        ((o3 > eps && o4 < -eps) || (o3 < -eps && o4 > eps));
    if (proper) {
        return true;
    }

    if (std::abs(o1) <= eps && on_segment(s1.a, s1.b, s2.a, eps)) return true;
    if (std::abs(o2) <= eps && on_segment(s1.a, s1.b, s2.b, eps)) return true;
    if (std::abs(o3) <= eps && on_segment(s2.a, s2.b, s1.a, eps)) return true;
    if (std::abs(o4) <= eps && on_segment(s2.a, s2.b, s1.b, eps)) return true;
    return false;
}

// Computes the intersection of two infinite lines given by point pairs.
inline std::optional<Point> line_intersection(const Point& p1, const Point& p2, const Point& q1, const Point& q2) {
    const Point r = p2 - p1;
    const Point s = q2 - q1;
    const double denom = cross(r, s);
    if (std::abs(denom) <= kEpsilon) {
        return std::nullopt;
    }
    const double t = cross(q1 - p1, s) / denom;
    return p1 + r * t;
}

// Computes the intersection of two infinite lines in ax + by + c = 0 form.
inline std::optional<Point> line_intersection(const Line& l1, const Line& l2) {
    const double det = l1.a * l2.b - l2.a * l1.b;
    if (std::abs(det) <= kEpsilon) {
        return std::nullopt;
    }
    const double x = (l2.c * l1.b - l1.c * l2.b) / det;
    const double y = (l1.c * l2.a - l2.c * l1.a) / det;
    return Point{x, y};
}

// Builds the implicit line representation passing through p and q.
inline Line through_points(const Point& p, const Point& q) {
    return {p.y - q.y, q.x - p.x, p.x * q.y - p.y * q.x};
}

// Returns the signed perpendicular distance from p to the directed line a->b.
inline double point_line_signed_distance(const Point& p, const Point& a, const Point& b) {
    const Point ab = b - a;
    const double len = length(ab);
    if (len <= kEpsilon) {
        return 0.0;
    }
    return cross(ab, p - a) / len;
}

// Rotates/translates p into a local frame whose origin is `origin` and whose
// x-axis is aligned with a caller-supplied angle.
inline Point rotate_into_frame(const Point& p, const Point& origin, double cos_theta, double sin_theta) {
    const double dx = p.x - origin.x;
    const double dy = p.y - origin.y;
    return {dx * cos_theta + dy * sin_theta, -dx * sin_theta + dy * cos_theta};
}

}  // namespace simplify
