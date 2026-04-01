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

// Verifies that a priority-queue candidate still matches the live linked-list
// topology after earlier collapses may have modified neighboring vertices.
bool candidate_is_current(const Candidate& candidate) {
    if (!candidate.a || !candidate.b || !candidate.c || !candidate.d) {
        return false;
    }
    if (!candidate.a->alive || !candidate.b->alive || !candidate.c->alive || !candidate.d->alive) {
        return false;
    }
    if (candidate.a->next != candidate.b || candidate.b->next != candidate.c || candidate.c->next != candidate.d) {
        return false;
    }
    if (candidate.b->prev != candidate.a || candidate.c->prev != candidate.b || candidate.d->prev != candidate.c) {
        return false;
    }
    return candidate.a->version == candidate.a_version && candidate.b->version == candidate.b_version &&
           candidate.c->version == candidate.c_version && candidate.d->version == candidate.d_version;
}

// Builds a fresh candidate collapse around the current vertex B.
std::optional<Candidate> make_candidate(Vertex* b) {
    if (!b || !b->alive) {
        return std::nullopt;
    }
    Vertex* a = b->prev;
    Vertex* c = b->next;
    if (!a || !c) {
        return std::nullopt;
    }
    Vertex* d = c->next;
    if (!d || !a->alive || !c->alive || !d->alive) {
        return std::nullopt;
    }
    if (a == b || b == c || c == d || d == a) {
        return std::nullopt;
    }

    const auto replacement = choose_replacement(a->point, b->point, c->point, d->point);
    if (!replacement) {
        return std::nullopt;
    }

    const double displacement = local_areal_displacement(a->point, b->point, c->point, d->point, *replacement);
    if (!std::isfinite(displacement)) {
        return std::nullopt;
    }

    return Candidate{b, a, c, d, *replacement, displacement, a->version, b->version, c->version, d->version};
}

// Shared-endpoint contacts with adjacent segments are allowed during topology
// checks, so we explicitly identify those benign cases.
bool segments_share_allowed_endpoint(const Segment& first, const Segment& second) {
    const bool share_a = same_point(first.a, second.a) || same_point(first.a, second.b);
    const bool share_b = same_point(first.b, second.a) || same_point(first.b, second.b);
    return share_a || share_b;
}

// Collects every currently live ring edge. The implementation is simple and
// favors readability over aggressive spatial indexing.
std::vector<std::pair<Vertex*, Vertex*>> collect_segments(const PolygonData& polygon) {
    std::vector<std::pair<Vertex*, Vertex*>> segments;
    for (const Ring& ring : polygon.rings) {
        if (!ring.head || ring.size < 2) {
            continue;
        }
        Vertex* current = ring.head;
        for (std::size_t i = 0; i < ring.size; ++i) {
            segments.emplace_back(current, current->next);
            current = current->next;
        }
    }
    return segments;
}

// Rejects a collapse if either new edge A-E or E-D would intersect an existing
// edge outside the locally replaced neighborhood.
bool topology_ok(const PolygonData& polygon, const Candidate& candidate) {
    Segment new_first{candidate.a->point, candidate.replacement};
    Segment new_second{candidate.replacement, candidate.d->point};

    const std::set<int> ignored_ids = {
        candidate.a->id,
        candidate.b->id,
        candidate.c->id,
        candidate.d->id,
        candidate.a->prev->id,
        candidate.d->next->id,
    };

    for (const auto& [u, v] : collect_segments(polygon)) {
        if (!u->alive || !v->alive) {
            continue;
        }
        if (ignored_ids.count(u->id) && ignored_ids.count(v->id)) {
            continue;
        }
        Segment existing{u->point, v->point};
        if (segments_intersect(existing, new_first) && !segments_share_allowed_endpoint(existing, new_first)) {
            return false;
        }
        if (segments_intersect(existing, new_second) && !segments_share_allowed_endpoint(existing, new_second)) {
            return false;
        }
    }
    return true;
}


    } // namespace
} // namespace simplify