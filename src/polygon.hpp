// polygon.hpp
// Data structures and function declarations for the APSC polygon simplifier.
//
// Architecture mirrors the reference implementation:
// - One flat global vertex pool (std::vector<Vertex>) — cache-friendly,
//   no pointer chasing, trivially resizable for 400k+ vertices.
// - Integer prev/next links and a version counter for O(1) lazy invalidation.
// - Collapse candidates carry snapshotted versions; stale entries are skipped
//   without a secondary data structure.
// - SpatialGrid is a plain struct with a generation counter for O(1)
//   per-query deduplication — no std::set or std::unordered_set in the hot path.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "geometry.hpp"

namespace simplify {

// ---------------------------------------------------------------------------
// Vertex pool node
// ---------------------------------------------------------------------------
struct Vertex {
    Point  pt;
    int    prev     = -1;
    int    next     = -1;
    int    ring_id  = -1;   // original ring identifier (from CSV)
    int    ring_idx = -1;   // index into the rings[] array
    bool   alive    = true;
    int    version  = 0;    // bumped on every in-place coordinate update
};

// ---------------------------------------------------------------------------
// Ring book-keeping (one entry per input ring)
// ---------------------------------------------------------------------------
struct Ring {
    int ring_id = -1;   // original CSV ring_id
    int head    = -1;   // pool index of any one live vertex in the ring
    int size    = 0;    // current live vertex count
};

// ---------------------------------------------------------------------------
// Collapse candidate — carries a full snapshot so the priority queue never
// needs to be mutated; stale entries are detected in O(1) by version check.
// ---------------------------------------------------------------------------
struct Collapse {
    int    A = -1, B = -1, C = -1, D = -1;
    int    vA = 0, vB = 0, vC = 0, vD = 0;  // version snapshots
    Point  E;                                  // pre-computed replacement point
    double displacement = 0.0;
    double total_cost   = 0.0;
    int    seq_id       = 0;                   // FIFO tiebreaker

    bool operator>(const Collapse& o) const {
        if (total_cost != o.total_cost) return total_cost > o.total_cost;
        return seq_id > o.seq_id;
    }
};

// ---------------------------------------------------------------------------
// Spatial grid — flat hash map from cell key to list of edge-start indices.
// queryGen + curGen give O(1) per-vertex deduplication without any set.
// Sized once at startup; queryGen is indexed directly by pool index so it
// must be as large as pool.size().
// ---------------------------------------------------------------------------
struct SpatialGrid {
    double inv_cell_size = 1.0;
    std::unordered_map<int64_t, std::vector<int>> cells;
    std::vector<int> query_gen;
    int cur_gen = 0;

    static int64_t cell_key(int ix, int iy) {
        return (static_cast<int64_t>(ix) << 32) |
               static_cast<int64_t>(static_cast<unsigned int>(iy));
    }

    // Must be called once before any add/remove/query.  pool_size must equal
    // the total number of slots in the vertex pool (pool.size()).
    void init(double cell_size, int pool_size);

    int to_grid(double v) const {
        return static_cast<int>(std::floor(v * inv_cell_size));
    }

    // Register the edge from pool[vi] to pool[pool[vi].next].
    // Caller must ensure pool[vi].next is valid and alive.
    void add_edge(const std::vector<Vertex>& pool, int vi);

    // Remove the edge that used to run from vi through old endpoints a→b.
    void remove_edge(int vi, const Point& a, const Point& b);

    template <typename F>
    bool traverse_cells(const Point& A, const Point& B, F&& func) const {
        if (inv_cell_size <= 0.0) return false;
        const double cs = 1.0 / inv_cell_size;

        int ix = to_grid(A.x);
        int iy = to_grid(A.y);
        const int ex = to_grid(B.x);
        const int ey = to_grid(B.y);

        if (func(cell_key(ix, iy))) return true;
        if (ix == ex && iy == ey) return false;

        const double dx = B.x - A.x;
        const double dy = B.y - A.y;

        const int step_x = (dx > 0.0) ? 1 : (dx < 0.0) ? -1 : 0;
        const int step_y = (dy > 0.0) ? 1 : (dy < 0.0) ? -1 : 0;

        const double inf = std::numeric_limits<double>::infinity();

        double t_max_x = inf;
        double t_max_y = inf;
        double t_delta_x = inf;
        double t_delta_y = inf;

        if (step_x != 0) {
            const double next_x = (step_x > 0) ? (static_cast<double>(ix + 1) * cs)
                                               : (static_cast<double>(ix) * cs);
            t_max_x = (next_x - A.x) / dx;
            t_delta_x = cs / std::abs(dx);
        }

        if (step_y != 0) {
            const double next_y = (step_y > 0) ? (static_cast<double>(iy + 1) * cs)
                                               : (static_cast<double>(iy) * cs);
            t_max_y = (next_y - A.y) / dy;
            t_delta_y = cs / std::abs(dy);
        }

        const int max_steps = std::abs(ex - ix) + std::abs(ey - iy) + 4;
        for (int steps = 0; steps < max_steps && !(ix == ex && iy == ey); ++steps) {
            if (t_max_x < t_max_y) {
                ix += step_x;
                t_max_x += t_delta_x;
            } else if (t_max_y < t_max_x) {
                iy += step_y;
                t_max_y += t_delta_y;
            } else {
                ix += step_x;
                iy += step_y;
                t_max_x += t_delta_x;
                t_max_y += t_delta_y;
            }
            if (func(cell_key(ix, iy))) return true;
        }
        return false;
    }

    // Iterate over every unique edge whose bounding box overlaps segment P→Q,
    // calling func(edge_start_index).  func returns true to stop early.
    template <typename F>
    void query_for_each(const Point& P, const Point& Q, F&& func) {
        if (++cur_gen == std::numeric_limits<int>::max()) {
            std::fill(query_gen.begin(), query_gen.end(), 0);
            cur_gen = 1;
        }

        traverse_cells(P, Q, [&](int64_t key) -> bool {
            auto it = cells.find(key);
            if (it == cells.end()) return false;
            for (int vi : it->second) {
                if (vi < 0 || static_cast<std::size_t>(vi) >= query_gen.size())
                    continue;
                if (query_gen[vi] == cur_gen) continue;
                query_gen[vi] = cur_gen;
                if (func(vi)) return true;
            }
            return false;
        });
    }
};

// ---------------------------------------------------------------------------
// Top-level I/O and simplification interface
// ---------------------------------------------------------------------------

// Shared mutable state — global pools mirroring the reference design so that
// all functions share one allocation without passing large structs by value.
extern std::vector<Vertex> glb_pool;
extern std::vector<Ring>   glb_rings;
extern SpatialGrid         glb_grid;

// Read ring_id,vertex_id,x,y CSV into glb_pool / glb_rings.
// Returns total vertex count.
int read_polygon_csv(const std::string& path);

// Run the APSC simplification loop until total live vertices <= target.
// Returns total areal displacement accumulated.
double simplify_polygon(int target_vertices);

// Compute signed shoelace area over all rings in glb_rings / glb_pool.
double compute_total_signed_area();

// Emit ring_id,vertex_id,x,y rows to stdout for all live vertices.
void write_output();

}  // namespace simplify
