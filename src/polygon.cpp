// polygon.cpp
// Core implementation of the polygon simplifier.
//
// The code below handles four main jobs:
//   1. reading CSV polygon input into linked-list rings,
//   2. building and validating local collapse candidates,
//   3. greedily applying the best valid collapses with topology checks, and
//   4. exporting/reporting the final simplified polygon.
#include "polygon.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace simplify {
namespace {

struct LocalSegment {
    Point p;
    Point q;
};

struct IntervalSegment {
    double x0 = 0.0;
    double x1 = 0.0;
    double y0 = 0.0;
    double y1 = 0.0;
};

// Computes the line of all possible replacement points E that preserve the
// signed area when the chain A-B-C-D is replaced by A-E-D.
Line area_preserving_line(const Point& a, const Point& b, const Point& c, const Point& d) {
    return {a.y - d.y, d.x - a.x, cross(a, b) + cross(b, c) + cross(c, d)};
}

double local_areal_displacement(const Point& a, const Point& b, const Point& c, const Point& d, const Point& e);

// Chooses the best placement for E by intersecting the area-preserving line
// with the neighboring support lines AB and CD, then picking the lower local
// displacement option.
std::optional<Point> choose_replacement(const Point& a, const Point& b, const Point& c, const Point& d) {
    const Line line = area_preserving_line(a, b, c, d);
    const Line ab = through_points(a, b);
    const Line cd = through_points(c, d);

    const auto on_ab = line_intersection(line, ab);
    const auto on_cd = line_intersection(line, cd);

    if (!on_ab && !on_cd) {
        return std::nullopt;
    }
    if (!on_cd) {
        return on_ab;
    }
    if (!on_ab) {
        return on_cd;
    }

    const double displacement_ab = local_areal_displacement(a, b, c, d, *on_ab);
    const double displacement_cd = local_areal_displacement(a, b, c, d, *on_cd);
    if (displacement_cd + kEpsilon < displacement_ab) {
        return on_cd;
    }
    return on_ab;
}

// Converts a polyline into x-ordered segments for the local displacement
// integration performed in a frame aligned with A-D.
std::vector<IntervalSegment> chain_to_intervals(const std::vector<Point>& chain) {
    std::vector<IntervalSegment> result;
    for (std::size_t i = 0; i + 1 < chain.size(); ++i) {
        Point p = chain[i];
        Point q = chain[i + 1];
        if (std::abs(p.x - q.x) <= kEpsilon) {
            continue;
        }
        if (p.x > q.x) {
            std::swap(p, q);
        }
        result.push_back({p.x, q.x, p.y, q.y});
    }
    return result;
}

// Interpolates the y value of a segment at a given x location.
std::optional<double> y_at_x(const IntervalSegment& segment, double x) {
    if (x + kEpsilon < segment.x0 || x - kEpsilon > segment.x1) {
        return std::nullopt;
    }
    const double t = (x - segment.x0) / (segment.x1 - segment.x0);
    return segment.y0 + (segment.y1 - segment.y0) * t;
}

// Approximates the area between the original chain A-B-C-D and the replacement
// chain A-E-D by integrating their separation in a local coordinate frame.
double local_areal_displacement(const Point& a, const Point& b, const Point& c, const Point& d, const Point& e) {
    Point axis = d - a;
    const double len = length(axis);
    if (len <= kEpsilon) {
        return std::numeric_limits<double>::infinity();
    }
    const double cos_theta = axis.x / len;
    const double sin_theta = axis.y / len;

    std::vector<Point> old_chain = {
        rotate_into_frame(a, a, cos_theta, sin_theta),
        rotate_into_frame(b, a, cos_theta, sin_theta),
        rotate_into_frame(c, a, cos_theta, sin_theta),
        rotate_into_frame(d, a, cos_theta, sin_theta)};
    std::vector<Point> new_chain = {
        rotate_into_frame(a, a, cos_theta, sin_theta),
        rotate_into_frame(e, a, cos_theta, sin_theta),
        rotate_into_frame(d, a, cos_theta, sin_theta)};

    std::vector<double> xs;
    xs.reserve(16);
    for (const Point& p : old_chain) xs.push_back(p.x);
    for (const Point& p : new_chain) xs.push_back(p.x);

    for (std::size_t i = 0; i + 1 < old_chain.size(); ++i) {
        for (std::size_t j = 0; j + 1 < new_chain.size(); ++j) {
            if (const auto inter = line_intersection(old_chain[i], old_chain[i + 1], new_chain[j], new_chain[j + 1])) {
                xs.push_back(inter->x);
            }
        }
    }

    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end(), [](double lhs, double rhs) {
                 return std::abs(lhs - rhs) <= kEpsilon;
             }),
             xs.end());

    const auto old_intervals = chain_to_intervals(old_chain);
    const auto new_intervals = chain_to_intervals(new_chain);

    double area = 0.0;
    for (std::size_t i = 0; i + 1 < xs.size(); ++i) {
        const double x0 = xs[i];
        const double x1 = xs[i + 1];
        if (x1 - x0 <= kEpsilon) {
            continue;
        }
        const double mid = 0.5 * (x0 + x1);

        std::optional<double> old_y0;
        std::optional<double> old_y1;
        std::optional<double> new_y0;
        std::optional<double> new_y1;

        for (const auto& seg : old_intervals) {
            if (auto y = y_at_x(seg, mid)) {
                old_y0 = y_at_x(seg, x0);
                old_y1 = y_at_x(seg, x1);
                break;
            }
        }
        for (const auto& seg : new_intervals) {
            if (auto y = y_at_x(seg, mid)) {
                new_y0 = y_at_x(seg, x0);
                new_y1 = y_at_x(seg, x1);
                break;
            }
        }

        if (!old_y0 || !old_y1 || !new_y0 || !new_y1) {
            continue;
        }

        area += 0.5 * (std::abs(*old_y0 - *new_y0) + std::abs(*old_y1 - *new_y1)) * (x1 - x0);
    }

    return area;
}


    } // namespace
} // namespace simplify