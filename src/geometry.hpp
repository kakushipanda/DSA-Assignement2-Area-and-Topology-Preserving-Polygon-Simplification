// geometry.hpp
// Low-level geometry primitives for the APSC polygon simplifier.
// Matches the architecture of the reference single-file implementation:
// plain structs, free inline functions, no std::optional, no heap allocation
// in hot paths.
#pragma once

#include <cmath>
#include <vector>

namespace simplify {

constexpr double kEps = 1e-9;

// ---------------------------------------------------------------------------
// Core point type
// ---------------------------------------------------------------------------
struct Point {
    double x = 0.0;
    double y = 0.0;
    Point() = default;
    Point(double x_, double y_) : x(x_), y(y_) {}
};

inline Point operator+(const Point& a, const Point& b) { return {a.x + b.x, a.y + b.y}; }
inline Point operator-(const Point& a, const Point& b) { return {a.x - b.x, a.y - b.y}; }
inline Point operator*(const Point& p, double s)       { return {p.x * s,   p.y * s};   }

inline double dot(const Point& a, const Point& b)   { return a.x * b.x + a.y * b.y; }
inline double cross2(const Point& a, const Point& b) { return a.x * b.y - a.y * b.x; }

// Signed cross product of vectors (b-a) and (c-a) — the core orientation predicate.
inline double cross(const Point& a, const Point& b, const Point& c) {
    return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
}

inline double length(const Point& p) { return std::sqrt(dot(p, p)); }

// ---------------------------------------------------------------------------
// Segment / line helpers
// ---------------------------------------------------------------------------

// Returns true when p lies on closed segment [a,b], within a loose tolerance.
inline bool on_segment(const Point& a, const Point& b, const Point& p) {
    return std::min(a.x, b.x) - 1e-12 <= p.x && p.x <= std::max(a.x, b.x) + 1e-12 &&
           std::min(a.y, b.y) - 1e-12 <= p.y && p.y <= std::max(a.y, b.y) + 1e-12;
}

// Proper segment intersection: endpoint touches do count (collinear overlap
// is treated as an intersection, matching the reference implementation).
inline bool segments_properly_intersect(const Point& p1, const Point& p2,
                                        const Point& p3, const Point& p4) {
    const double d1 = cross(p3, p4, p1);
    const double d2 = cross(p3, p4, p2);
    const double d3 = cross(p1, p2, p3);
    const double d4 = cross(p1, p2, p4);

    if (((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) &&
        ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0)))
        return true;

    if (std::fabs(d1) < 1e-12 && on_segment(p3, p4, p1)) return true;
    if (std::fabs(d2) < 1e-12 && on_segment(p3, p4, p2)) return true;
    if (std::fabs(d3) < 1e-12 && on_segment(p1, p2, p3)) return true;
    if (std::fabs(d4) < 1e-12 && on_segment(p1, p2, p4)) return true;
    return false;
}

// Signed shoelace area of a closed ring (CCW → positive).
inline double signed_ring_area(const std::vector<Point>& ring) {
    const int n = static_cast<int>(ring.size());
    if (n < 3) return 0.0;
    double area = 0.0;
    for (int i = 0; i < n; ++i) {
        const Point& p = ring[i];
        const Point& q = ring[(i + 1) % n];
        area += p.x * q.y - q.x * p.y;
    }
    return 0.5 * area;
}

// ---------------------------------------------------------------------------
// Line-E math (APSC placement function)
// ax + by + c = 0 form.
// ---------------------------------------------------------------------------
struct LineEq {
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
};

// Computes the area-preserving line E that passes through A and D such that
// triangles ABE and CDE have equal and opposite signed areas to AEB and ECD.
inline LineEq compute_line_e(const Point& A, const Point& B,
                              const Point& C, const Point& D) {
    LineEq eq;
    eq.a = D.y - A.y;
    eq.b = A.x - D.x;
    eq.c = -B.y * A.x + (A.y - C.y) * B.x + (B.y - D.y) * C.x + C.y * D.x;
    return eq;
}

// True when line E is parallel to (coincident with) segment AD — fallback case.
inline bool line_e_coincident_with_ad(const LineEq& eq, const Point& A) {
    const double c_ad = -(eq.a * A.x + eq.b * A.y);
    const double norm = std::sqrt(eq.a * eq.a + eq.b * eq.b);
    if (norm < 1e-30) return true;
    return std::fabs(eq.c - c_ad) / norm < 1e-12;
}

// Intersects line E (implicit) with the line through P and Q (parametric).
// Returns false if they are parallel.
inline bool intersect_line_e_with_segment(const LineEq& eq,
                                          const Point& P, const Point& Q,
                                          Point& result) {
    const double a2  = Q.y - P.y;
    const double b2  = P.x - Q.x;
    const double c2  = a2 * P.x + b2 * P.y;
    const double det = eq.a * b2 - a2 * eq.b;
    if (std::fabs(det) < 1e-20) return false;
    result.x = ((-eq.c) * b2 - c2 * eq.b) / det;
    result.y = (eq.a  * c2 - a2 * (-eq.c)) / det;
    return true;
}

// Full APSC placement: given A-B-C-D, compute replacement vertex E on line E
// following the tie-breaking rules from Kronenfeld's Fig 6
inline Point compute_replacement_point(const Point& A, const Point& B,
                                       const Point& C, const Point& D) {
    const LineEq eq = compute_line_e(A, B, C, D);

    // Fig 6a: line E is the same as AD → E = D.
    if (line_e_coincident_with_ad(eq, A)) return D;

    const double c_ad   = -(eq.a * A.x + eq.b * A.y);
    const double fB     = eq.a * B.x + eq.b * B.y + c_ad;
    const double fC     = eq.a * C.x + eq.b * C.y + c_ad;
    const double norm   = std::sqrt(eq.a * eq.a + eq.b * eq.b);
    const double dB     = std::fabs(fB) / norm;
    const double dC     = std::fabs(fC) / norm;

    // B and C on same side of AD?
    const bool same_side = (dB < 1e-12 || dC < 1e-12) ? true : (fB * fC > 0.0);

    // Fig 6c tie-breaker: prefer the nearer of B and C as anchor segment.
    bool use_ab;
    if (same_side) {
        if (std::fabs(dB - dC) < 1e-12 * std::max(1.0, dB + dC))
            use_ab = (dB <= dC + 1e-15);
        else
            use_ab = (dB > dC);
    } else {
        const double fE_on_ad = -(eq.c - c_ad);
        use_ab = (fB * fE_on_ad > 0.0);
    }

    Point result;
    if (use_ab) {
        if (intersect_line_e_with_segment(eq, A, B, result)) return result;
        if (intersect_line_e_with_segment(eq, C, D, result)) return result;
        return B;
    } else {
        if (intersect_line_e_with_segment(eq, C, D, result)) return result;
        if (intersect_line_e_with_segment(eq, A, B, result)) return result;
        return C;
    }
}

// ---------------------------------------------------------------------------
// Areal displacement: area between polyline A-B-C-D and polyline A-E-D.
// Uses fixed-size stack arrays — zero heap allocation.
// ---------------------------------------------------------------------------
inline double shoelace_area(const Point* pts, int n) {
    double area = 0.0;
    for (int i = 0; i < n; ++i) {
        const int j = (i + 1) % n;
        area += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
    }
    return 0.5 * area;
}

inline double compute_displacement(const Point& A, const Point& B,
                                   const Point& C, const Point& D,
                                   const Point& E) {
    // Fast path (Green's theorem): when the two chains A→B→C→D and A→E→D
    // don't properly cross, the displacement equals the absolute shoelace
    // area of the closed polygon A→B→C→D→E→A — five cross products, O(1).
    {
        Point segs[5][2] = {{A,B},{B,C},{C,D},{A,E},{E,D}};
        const int pairs[4][2] = {{0,4},{1,3},{1,4},{2,3}};
        bool has_crossing = false;
        for (int p = 0; p < 4; ++p) {
            const int i = pairs[p][0], j = pairs[p][1];
            const double dx1 = segs[i][1].x - segs[i][0].x;
            const double dy1 = segs[i][1].y - segs[i][0].y;
            const double dx2 = segs[j][1].x - segs[j][0].x;
            const double dy2 = segs[j][1].y - segs[j][0].y;
            const double denom = dx1 * dy2 - dy1 * dx2;
            if (std::fabs(denom) < 1e-20) continue;
            const double dx3 = segs[j][0].x - segs[i][0].x;
            const double dy3 = segs[j][0].y - segs[i][0].y;
            const double t = (dx3 * dy2 - dy3 * dx2) / denom;
            const double u = (dx3 * dy1 - dy3 * dx1) / denom;
            if (t > 1e-12 && t < 1.0 - 1e-12 && u > 1e-12 && u < 1.0 - 1e-12) {
                has_crossing = true;
                break;
            }
        }
        if (!has_crossing) {
            const double area = cross2(A, B) + cross2(B, C) + cross2(C, D)
                              + cross2(D, E) + cross2(E, A);
            return std::fabs(0.5 * area);
        }
    }

    // Slow path: chains cross — full lobe decomposition needed.
    struct XPt { double gp1, gp2; Point pt; };
    XPt xpts[16];
    int nxp = 0;

    Point segs[5][2] = {{A,B},{B,C},{C,D},{A,E},{E,D}};
    const int pairs[4][2] = {{0,4},{1,3},{1,4},{2,3}};
    for (int p = 0; p < 4; ++p) {
        const int i = pairs[p][0], j = pairs[p][1];
        const double dx1 = segs[i][1].x - segs[i][0].x;
        const double dy1 = segs[i][1].y - segs[i][0].y;
        const double dx2 = segs[j][1].x - segs[j][0].x;
        const double dy2 = segs[j][1].y - segs[j][0].y;
        const double denom = dx1 * dy2 - dy1 * dx2;
        if (std::fabs(denom) < 1e-20) continue;
        const double dx3 = segs[j][0].x - segs[i][0].x;
        const double dy3 = segs[j][0].y - segs[i][0].y;
        const double t = (dx3 * dy2 - dy3 * dx2) / denom;
        const double u = (dx3 * dy1 - dy3 * dx1) / denom;
        if (t > 1e-12 && t < 1.0 - 1e-12 && u > 1e-12 && u < 1.0 - 1e-12)
            xpts[nxp++] = {static_cast<double>(i) + t,
                           static_cast<double>(j - 3) + u,
                           {segs[i][0].x + t * dx1, segs[i][0].y + t * dy1}};
    }

    // Endpoint-on-segment degenerate cases.
    auto pt_on_seg = [](const Point& P, const Point& S0, const Point& S1, double& t) -> bool {
        const double dx = S1.x - S0.x, dy = S1.y - S0.y;
        const double len2 = dx * dx + dy * dy;
        if (len2 < 1e-30) return false;
        t = ((P.x - S0.x) * dx + (P.y - S0.y) * dy) / len2;
        if (t < 1e-9 || t > 1.0 - 1e-9) return false;
        const double cx = (P.x - S0.x) * dy - (P.y - S0.y) * dx;
        return std::fabs(cx) < std::sqrt(len2) * 1e-9;
    };
    double tv;
    if (pt_on_seg(B, A, E, tv)) xpts[nxp++] = {1.0,       tv,       B};
    if (pt_on_seg(B, E, D, tv)) xpts[nxp++] = {1.0,       1.0 + tv, B};
    if (pt_on_seg(C, A, E, tv)) xpts[nxp++] = {2.0,       tv,       C};
    if (pt_on_seg(C, E, D, tv)) xpts[nxp++] = {2.0,       1.0 + tv, C};
    if (pt_on_seg(E, A, B, tv)) xpts[nxp++] = {tv,        1.0,      E};
    if (pt_on_seg(E, B, C, tv)) xpts[nxp++] = {1.0 + tv,  1.0,      E};
    if (pt_on_seg(E, C, D, tv)) xpts[nxp++] = {2.0 + tv,  1.0,      E};

    // Build ordered chains.
    struct OPt { double gp; Point pt; };
    OPt chain1[20], chain2[20];
    int n1 = 0, n2 = 0;
    chain1[n1++] = {0.0, A}; chain1[n1++] = {1.0, B};
    chain1[n1++] = {2.0, C}; chain1[n1++] = {3.0, D};
    chain2[n2++] = {0.0, A}; chain2[n2++] = {1.0, E}; chain2[n2++] = {2.0, D};
    for (int i = 0; i < nxp; ++i) {
        chain1[n1++] = {xpts[i].gp1, xpts[i].pt};
        chain2[n2++] = {xpts[i].gp2, xpts[i].pt};
    }
    auto cmp = [](const OPt& a, const OPt& b) { return a.gp < b.gp; };
    std::sort(chain1, chain1 + n1, cmp);
    std::sort(chain2, chain2 + n2, cmp);

    auto dedup = [](OPt* arr, int& n) {
        int j = 0;
        for (int i = 0; i < n; ++i)
            if (j == 0 || std::fabs(arr[i].gp - arr[j-1].gp) > 1e-15)
                arr[j++] = arr[i];
        n = j;
    };
    dedup(chain1, n1);
    dedup(chain2, n2);

    // Shared boundary points.
    struct SharedPoints { int i1, i2; };
    SharedPoints shared[20];
    int nsh = 0;
    shared[nsh++] = {0, 0};
    for (int x = 0; x < nxp; ++x) {
        int i1 = -1, i2 = -1;
        for (int k = 0; k < n1; ++k)
            if (std::fabs(chain1[k].gp - xpts[x].gp1) < 1e-12) { i1 = k; break; }
        for (int k = 0; k < n2; ++k)
            if (std::fabs(chain2[k].gp - xpts[x].gp2) < 1e-12) { i2 = k; break; }
        if (i1 >= 0 && i2 >= 0) shared[nsh++] = {i1, i2};
    }
    shared[nsh++] = {n1 - 1, n2 - 1};
    std::sort(shared, shared + nsh, [](const SharedPoints& a, const SharedPoints& b) { return a.i1 < b.i1; });
    { int j = 0; for (int i = 0; i < nsh; ++i) { if (j==0||shared[i].i1!=shared[j-1].i1) shared[j++]=shared[i]; } nsh=j; }

    double total = 0.0;
    for (int s = 0; s + 1 < nsh; ++s) {
        Point sub[20];
        int ns = 0;
        for (int k = shared[s].i1; k <= shared[s+1].i1; ++k) sub[ns++] = chain1[k].pt;
        for (int k = shared[s+1].i2 - 1; k > shared[s].i2; --k) sub[ns++] = chain2[k].pt;
        if (ns >= 3) total += std::fabs(shoelace_area(sub, ns));
    }
    return total;
}

}  // namespace simplify
