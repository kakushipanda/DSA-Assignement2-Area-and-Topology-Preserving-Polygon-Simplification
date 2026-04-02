// polygon.cpp
// APSC polygon simplifier — implementation.
//
// Follows the reference single-file design closely:
// - one flat global vertex pool (glb_pool), one rings array (glb_rings),
//   one global spatial grid (glb_grid);
// - in-place vertex collapse: B ← E, C marked dead;
// - lazy priority-queue invalidation via per-vertex version counters;
// - fixed-size stack arrays in every hot-path function (zero heap allocation
//   per collapse iteration beyond the priority-queue node itself);
// - generation-counter deduplication in the spatial grid (no set).
#include "polygon.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>

namespace simplify {

// ---------------------------------------------------------------------------
// Global state
// ---------------------------------------------------------------------------
std::vector<Vertex> glb_pool;
std::vector<Ring>   glb_rings;
SpatialGrid         glb_grid;

// ---------------------------------------------------------------------------
// SpatialGrid method bodies
// ---------------------------------------------------------------------------
void SpatialGrid::init(double cell_size, int pool_size) {
    inv_cell_size = 1.0 / cell_size;
    cells.clear();
    cells.reserve(static_cast<std::size_t>(pool_size) * 2);
    query_gen.assign(static_cast<std::size_t>(pool_size), 0);
    cur_gen = 0;
}

void SpatialGrid::add_edge(const std::vector<Vertex>& pool, int vi) {
    const int nxt = pool[static_cast<std::size_t>(vi)].next;
    if (nxt < 0 || static_cast<std::size_t>(nxt) >= pool.size()) return;
    if (!pool[static_cast<std::size_t>(nxt)].alive) return;

    const Point& a = pool[static_cast<std::size_t>(vi)].pt;
    const Point& b = pool[static_cast<std::size_t>(nxt)].pt;
    const int x0 = to_grid(std::min(a.x, b.x));
    const int x1 = to_grid(std::max(a.x, b.x));
    const int y0 = to_grid(std::min(a.y, b.y));
    const int y1 = to_grid(std::max(a.y, b.y));
    for (int ix = x0; ix <= x1; ++ix)
        for (int iy = y0; iy <= y1; ++iy)
            cells[cell_key(ix, iy)].push_back(vi);
}

void SpatialGrid::remove_edge(int vi, const Point& a, const Point& b) {
    const int x0 = to_grid(std::min(a.x, b.x));
    const int x1 = to_grid(std::max(a.x, b.x));
    const int y0 = to_grid(std::min(a.y, b.y));
    const int y1 = to_grid(std::max(a.y, b.y));
    for (int ix = x0; ix <= x1; ++ix) {
        for (int iy = y0; iy <= y1; ++iy) {
            auto it = cells.find(cell_key(ix, iy));
            if (it == cells.end()) continue;
            auto& v = it->second;
            for (int k = static_cast<int>(v.size()) - 1; k >= 0; --k) {
                if (v[k] == vi) { v[k] = v.back(); v.pop_back(); break; }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------
int read_polygon_csv(const std::string& path) {
    std::ifstream fin(path);
    if (!fin.is_open())
        throw std::runtime_error("Cannot open input file: " + path);

    std::string line;
    if (!std::getline(fin, line))
        throw std::runtime_error("Input file is empty: " + path);

    std::map<int, std::vector<std::pair<int, Point>>> grouped;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string tok;
        std::getline(ss, tok, ','); const int ring_id   = std::stoi(tok);
        std::getline(ss, tok, ','); const int vertex_id = std::stoi(tok);
        std::getline(ss, tok, ','); const double x      = std::stod(tok);
        std::getline(ss, tok, ','); const double y      = std::stod(tok);
        grouped[ring_id].push_back({vertex_id, {x, y}});
    }

    glb_pool.clear();
    glb_rings.clear();

    int total = 0;
    for (auto& [ring_id, entries] : grouped) {
        std::sort(entries.begin(), entries.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });

        const int n    = static_cast<int>(entries.size());
        const int base = static_cast<int>(glb_pool.size());
        const int ridx = static_cast<int>(glb_rings.size());

        glb_pool.resize(glb_pool.size() + static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            Vertex& v  = glb_pool[static_cast<std::size_t>(base + i)];
            v.pt       = entries[static_cast<std::size_t>(i)].second;
            v.ring_id  = ring_id;
            v.ring_idx = ridx;
            v.alive    = true;
            v.version  = 0;
            v.prev     = base + ((i - 1 + n) % n);
            v.next     = base + ((i + 1)     % n);
        }

        Ring r;
        r.ring_id = ring_id;
        r.head    = base;
        r.size    = n;
        glb_rings.push_back(r);
        total += n;
    }
    return total;
}

// ---------------------------------------------------------------------------
// Area helpers
// ---------------------------------------------------------------------------
double compute_total_signed_area() {
    double area = 0.0;
    for (const Ring& r : glb_rings) {
        if (r.size < 3 || r.head < 0) continue;
        int v = r.head;
        do {
            const int nxt = glb_pool[static_cast<std::size_t>(v)].next;
            area += glb_pool[static_cast<std::size_t>(v)].pt.x *
                        glb_pool[static_cast<std::size_t>(nxt)].pt.y
                  - glb_pool[static_cast<std::size_t>(nxt)].pt.x *
                        glb_pool[static_cast<std::size_t>(v)].pt.y;
            v = nxt;
        } while (v != r.head);
    }
    return 0.5 * area;
}

// ---------------------------------------------------------------------------
// Intersection check via grid
// ---------------------------------------------------------------------------
static bool seg_intersects_any(const Point& P, const Point& Q,
                                int exA, int exB, int exC, int exD) {
    bool found = false;
    glb_grid.query_for_each(P, Q, [&](int vi) -> bool {
        if (!glb_pool[static_cast<std::size_t>(vi)].alive) return false;
        const int nxt = glb_pool[static_cast<std::size_t>(vi)].next;
        if (vi == exA || vi == exB || vi == exC || vi == exD ||
            nxt == exA || nxt == exB || nxt == exC || nxt == exD)
            return false;
        if (!glb_pool[static_cast<std::size_t>(nxt)].alive) return false;
        if (segments_properly_intersect(P, Q,
                glb_pool[static_cast<std::size_t>(vi)].pt,
                glb_pool[static_cast<std::size_t>(nxt)].pt)) {
            found = true;
            return true;
        }
        return false;
    });
    return found;
}

// ---------------------------------------------------------------------------
// Main simplification loop
// ---------------------------------------------------------------------------
double simplify_polygon(int target_vertices) {
    const int num_rings = static_cast<int>(glb_rings.size());

    // ------------------------------------------------------------------
    // Spatial grid initialisation
    // ------------------------------------------------------------------
    {
        double min_x =  1e300, max_x = -1e300;
        double min_y =  1e300, max_y = -1e300;
        for (const Vertex& v : glb_pool) {
            if (!v.alive) continue;
            min_x = std::min(min_x, v.pt.x); max_x = std::max(max_x, v.pt.x);
            min_y = std::min(min_y, v.pt.y); max_y = std::max(max_y, v.pt.y);
        }
        double span = std::max(max_x - min_x, max_y - min_y);
        double cs   = span / (4.0 * std::sqrt(static_cast<double>(glb_pool.size())));
        if (cs < 1e-12) cs = 1.0;

        glb_grid.init(cs, static_cast<int>(glb_pool.size()));

        for (int r = 0; r < num_rings; ++r) {
            if (glb_rings[r].size < 2 || glb_rings[r].head < 0) continue;
            int v = glb_rings[r].head;
            do {
                glb_grid.add_edge(glb_pool, v);
                v = glb_pool[static_cast<std::size_t>(v)].next;
            } while (v != glb_rings[r].head);
        }
    }

    // ------------------------------------------------------------------
    // Candidate factory
    // ------------------------------------------------------------------
    int seq_counter = 0;
    auto make_collapse = [&](int ai, int bi, int ci, int di) -> Collapse {
        Collapse col;
        col.A = ai; col.B = bi; col.C = ci; col.D = di;
        col.vA = glb_pool[static_cast<std::size_t>(ai)].version;
        col.vB = glb_pool[static_cast<std::size_t>(bi)].version;
        col.vC = glb_pool[static_cast<std::size_t>(ci)].version;
        col.vD = glb_pool[static_cast<std::size_t>(di)].version;
        col.E  = compute_replacement_point(
                     glb_pool[static_cast<std::size_t>(ai)].pt,
                     glb_pool[static_cast<std::size_t>(bi)].pt,
                     glb_pool[static_cast<std::size_t>(ci)].pt,
                     glb_pool[static_cast<std::size_t>(di)].pt);
        col.displacement = compute_displacement(
                     glb_pool[static_cast<std::size_t>(ai)].pt,
                     glb_pool[static_cast<std::size_t>(bi)].pt,
                     glb_pool[static_cast<std::size_t>(ci)].pt,
                     glb_pool[static_cast<std::size_t>(di)].pt,
                     col.E);
        col.seq_id = seq_counter++;
        return col;
    };

    // ------------------------------------------------------------------
    // Helper: validate adjacency before pushing
    // ------------------------------------------------------------------
    auto try_push = [&](std::priority_queue<Collapse,
                                            std::vector<Collapse>,
                                            std::greater<Collapse>>& pq,
                        int ai, int bi, int ci, int di) {
        const std::size_t sz = glb_pool.size();
        if (ai < 0 || bi < 0 || ci < 0 || di < 0) return;
        if (static_cast<std::size_t>(ai) >= sz || static_cast<std::size_t>(bi) >= sz ||
            static_cast<std::size_t>(ci) >= sz || static_cast<std::size_t>(di) >= sz) return;
        if (!glb_pool[static_cast<std::size_t>(ai)].alive ||
            !glb_pool[static_cast<std::size_t>(bi)].alive ||
            !glb_pool[static_cast<std::size_t>(ci)].alive ||
            !glb_pool[static_cast<std::size_t>(di)].alive) return;
        if (glb_pool[static_cast<std::size_t>(ai)].next != bi ||
            glb_pool[static_cast<std::size_t>(bi)].next != ci ||
            glb_pool[static_cast<std::size_t>(ci)].next != di) return;
        pq.push(make_collapse(ai, bi, ci, di));
    };

    // ------------------------------------------------------------------
    // Build initial priority queue (pre-allocated to avoid rehashing)
    // ------------------------------------------------------------------
    std::priority_queue<Collapse, std::vector<Collapse>, std::greater<Collapse>> pq;
    {
        std::vector<Collapse> initial;
        initial.reserve(static_cast<std::size_t>(glb_pool.size()));
        for (int r = 0; r < num_rings; ++r) {
            if (glb_rings[r].size < 4 || glb_rings[r].head < 0) continue;
            int ai = glb_rings[r].head;
            do {
                int bi = glb_pool[static_cast<std::size_t>(ai)].next;
                int ci = glb_pool[static_cast<std::size_t>(bi)].next;
                int di = glb_pool[static_cast<std::size_t>(ci)].next;
                initial.push_back(make_collapse(ai, bi, ci, di));
                ai = glb_pool[static_cast<std::size_t>(ai)].next;
            } while (ai != glb_rings[r].head);
        }
        pq = std::priority_queue<Collapse, std::vector<Collapse>,
                                 std::greater<Collapse>>(
                 std::greater<Collapse>(), std::move(initial));
    }

    // ------------------------------------------------------------------
    // APSC main loop
    // ------------------------------------------------------------------
    int total_vertices = 0;
    for (const Ring& r : glb_rings) total_vertices += r.size;

    double total_displacement = 0.0;

    while (total_vertices > target_vertices && !pq.empty()) {
        Collapse col = pq.top();
        pq.pop();

        // Bounds check
        {
            const std::size_t sz = glb_pool.size();
            if (col.A < 0 || col.B < 0 || col.C < 0 || col.D < 0) continue;
            if (static_cast<std::size_t>(col.A) >= sz ||
                static_cast<std::size_t>(col.B) >= sz ||
                static_cast<std::size_t>(col.C) >= sz ||
                static_cast<std::size_t>(col.D) >= sz) continue;
        }

        const Vertex& vA = glb_pool[static_cast<std::size_t>(col.A)];
        const Vertex& vB = glb_pool[static_cast<std::size_t>(col.B)];
        const Vertex& vC = glb_pool[static_cast<std::size_t>(col.C)];
        const Vertex& vD = glb_pool[static_cast<std::size_t>(col.D)];

        // Lazy invalidation via version + alive checks
        if (!vA.alive || !vB.alive || !vC.alive || !vD.alive) continue;
        if (vA.version != col.vA || vB.version != col.vB ||
            vC.version != col.vC || vD.version != col.vD) continue;
        if (vB.prev != col.A || vB.next != col.C || vC.next != col.D) continue;

        Ring& ring = glb_rings[static_cast<std::size_t>(vB.ring_idx)];
        if (ring.size <= 3) continue;

        const int prev_a = vA.prev;
        const int next_d = vD.next;

        // Topology: new edges A→E and E→D must not cross any existing edge.
        if (seg_intersects_any(vA.pt, col.E, prev_a, col.A, col.B, col.C))
            continue;
        if (seg_intersects_any(col.E, vD.pt, col.B, col.C, col.D, next_d))
            continue;

        // ------------------------------------------------------------------
        // Commit
        // ------------------------------------------------------------------
        const Point old_b = vB.pt;
        const Point old_c = vC.pt;

        glb_grid.remove_edge(col.A, vA.pt, old_b);
        glb_grid.remove_edge(col.B, old_b, old_c);
        glb_grid.remove_edge(col.C, old_c, vD.pt);

        glb_pool[static_cast<std::size_t>(col.B)].pt      = col.E;
        glb_pool[static_cast<std::size_t>(col.B)].version++;
        glb_pool[static_cast<std::size_t>(col.C)].alive   = false;
        glb_pool[static_cast<std::size_t>(col.B)].next    = col.D;
        glb_pool[static_cast<std::size_t>(col.D)].prev    = col.B;

        if (ring.head == col.C) ring.head = col.B;

        glb_grid.add_edge(glb_pool, col.A);
        glb_grid.add_edge(glb_pool, col.B);

        ring.size--;
        total_vertices--;
        total_displacement += col.displacement;

        // ------------------------------------------------------------------
        // Enqueue neighbours (E = col.B after mutation)
        // ------------------------------------------------------------------
        if (ring.size >= 4) {
            const int Ei = col.B;
            const int Ai = col.A;
            const int Di = col.D;

            try_push(pq, glb_pool[static_cast<std::size_t>(Ai)].prev, Ai, Ei, Di);
            try_push(pq, Ai, Ei, Di, glb_pool[static_cast<std::size_t>(Di)].next);

            const int prev_ai   = glb_pool[static_cast<std::size_t>(Ai)].prev;
            const int prev_prev = glb_pool[static_cast<std::size_t>(prev_ai)].prev;
            if (prev_prev != Ei)
                try_push(pq, prev_prev, prev_ai, Ai, Ei);

            const int next_di   = glb_pool[static_cast<std::size_t>(Di)].next;
            const int next_next = glb_pool[static_cast<std::size_t>(next_di)].next;
            if (next_next != Ei)
                try_push(pq, Ei, Di, next_di, next_next);
        }
    }

    return total_displacement;
}

// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------
void write_output() {
    std::cout << "ring_id,vertex_id,x,y\n";
    for (const Ring& r : glb_rings) {
        if (r.head < 0 || r.size <= 0) continue;
        int vid = 0;
        int v   = r.head;
        do {
            const Vertex& vtx = glb_pool[static_cast<std::size_t>(v)];
            std::ostringstream xs, ys;
            xs << std::setprecision(10) << vtx.pt.x;
            ys << std::setprecision(10) << vtx.pt.y;
            std::cout << r.ring_id << ',' << vid++ << ','
                      << xs.str() << ',' << ys.str() << '\n';
            v = vtx.next;
        } while (v != r.head);
    }
}

}  // namespace simplify
