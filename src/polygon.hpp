// polygon.hpp
// Higher-level polygon and simplification data structures.
//
// The simplifier operates on rings represented as circular doubly-linked lists
// of vertices so that local collapses can update neighbors cheaply. This header
// defines those structures plus the public entry points used by main.cpp.
#pragma once

#include <memory>
#include <optional>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

#include "geometry.hpp"

namespace simplify {

struct Vertex {
    int id = -1;
    int ring_id = -1;
    Point point;
    Vertex* prev = nullptr;
    Vertex* next = nullptr;
    bool alive = true;
    std::size_t version = 0;
};

// Each input/output ring keeps a stable identifier and a pointer to the current
// traversal head in the circular linked list.
struct Ring {
    int ring_id = -1;
    Vertex* head = nullptr;
    std::size_t size = 0;
};

// Candidate collapse record stored in the priority queue. It captures the
// current four-vertex window A-B-C-D, the computed replacement point E, and
// version numbers used to discard stale queue entries after local edits.
struct Candidate {
    Vertex* b = nullptr;
    Vertex* a = nullptr;
    Vertex* c = nullptr;
    Vertex* d = nullptr;
    Point replacement;
    double displacement = 0.0;
    std::size_t a_version = 0;
    std::size_t b_version = 0;
    std::size_t c_version = 0;
    std::size_t d_version = 0;
};

// Min-heap ordering for collapse candidates: lower areal displacement wins,
// with a stable tie-breaker on the candidate's B vertex id.
struct CandidateCompare {
    bool operator()(const Candidate& lhs, const Candidate& rhs) const {
        if (!nearly_equal(lhs.displacement, rhs.displacement)) {
            return lhs.displacement > rhs.displacement;
        }
        return lhs.b->id > rhs.b->id;
    }
};

// In-memory polygon representation containing every ring plus owning storage
// for all vertices, including newly created Steiner points.
struct PolygonData {
    std::vector<Ring> rings;
    std::vector<std::unique_ptr<Vertex>> storage;
    std::size_t total_vertices = 0;
};

// Plain exported ring representation used for final output/reporting.
struct OutputRing {
    int ring_id = -1;
    std::vector<Point> vertices;
};

// Final program output bundle returned by the simplifier core.
struct ProgramOutput {
    std::vector<OutputRing> rings;
    double input_signed_area = 0.0;
    double output_signed_area = 0.0;
    double total_areal_displacement = 0.0;
};

// Parses the input CSV file into linked-list rings.
PolygonData read_polygon_csv(const std::string& path);

// Runs the simplification algorithm until the requested target is reached or no
// more valid collapses remain.
ProgramOutput simplify_polygon(PolygonData& polygon, std::size_t target_vertices);

// Converts the linked-list representation back into plain vectors for printing.
std::vector<OutputRing> export_rings(const PolygonData& polygon);

// Computes the algebraic sum of ring areas across the polygon.
double compute_total_signed_area(const std::vector<OutputRing>& rings);

// Estimates the total areal displacement between the original and simplified
// ring sets by comparing vertical cross-sections.
double compute_total_areal_displacement(const std::vector<OutputRing>& original,
                                        const std::vector<OutputRing>& simplified);

}  // namespace