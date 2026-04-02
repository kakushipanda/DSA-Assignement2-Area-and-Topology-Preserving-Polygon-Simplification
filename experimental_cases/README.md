# Additional Experimental Edge Cases (20 total)

These datasets are analytically constructed to cover boundary conditions and geometric/topological stress patterns that are either absent or underrepresented in the original provided suite.

**Important:** expected outputs in `output_*.txt` were derived by manual geometric calculation (shoelace signed area and unchanged ring geometry), **not** by running `./simplify` to auto-generate golden files.

## Case inventory

| # | Case | Input | Expected output | Target | Vertices | Holes | Primary stress area |
|---:|---|---|---|---:|---:|---:|---|
| 1 | `axis_aligned_square` | `input_axis_aligned_square.csv` | `output_axis_aligned_square.txt` | 4 | 4 | 0 | Tiny convex baseline |
| 2 | `four_holes_compact` | `input_four_holes_compact.csv` | `output_four_holes_compact.txt` | 16 | 16 | 4 | Dense multi-hole topology |
| 3 | `large_with_hole` | `input_large_with_hole.csv` | `output_large_with_hole.txt` | 8 | 8 | 1 | Very large coordinates + hole subtraction |
| 4 | `minimal_triangle` | `input_minimal_triangle.csv` | `output_minimal_triangle.txt` | 3 | 3 | 0 | Minimum legal ring size |
| 5 | `narrow_corridor` | `input_narrow_corridor.csv` | `output_narrow_corridor.txt` | 8 | 8 | 0 | Narrow-gap geometry |
| 6 | `near_collinear_large` | `input_near_collinear_large.csv` | `output_near_collinear_large.txt` | 5 | 5 | 0 | Near-degenerate edges at large scale |
| 7 | `regular128` | `input_regular128.csv` | `output_regular128.txt` | 128 | 128 | 0 | High-n smooth ring baseline |
| 8 | `regular32` | `input_regular32.csv` | `output_regular32.txt` | 32 | 32 | 0 | Medium-n smooth ring baseline |
| 9 | `star20` | `input_star20.csv` | `output_star20.txt` | 20 | 20 | 0 | Alternating curvature / spikes |
| 10 | `triangle_hole_minimal` | `input_triangle_hole_minimal.csv` | `output_triangle_hole_minimal.txt` | 7 | 7 | 1 | Minimal interior ring (3 vertices) |
| 11 | `skinny_rectangle_6` | `input_skinny_rectangle_6.csv` | `output_skinny_rectangle_6.txt` | 6 | 6 | 0 | Very high aspect-ratio region |
| 12 | `tiny_hole_clearance` | `input_tiny_hole_clearance.csv` | `output_tiny_hole_clearance.txt` | 7 | 7 | 1 | Tiny hole area and tight floating-point margin |
| 13 | `eight_holes_field` | `input_eight_holes_field.csv` | `output_eight_holes_field.txt` | 28 | 28 | 8 | Larger hole count than provided suite |
| 14 | `regular256` | `input_regular256.csv` | `output_regular256.txt` | 256 | 256 | 0 | Very high-n smooth ring scaling point |
| 15 | `negative_coords_mixed` | `input_negative_coords_mixed.csv` | `output_negative_coords_mixed.txt` | 5 | 5 | 0 | Negative-coordinate handling |
| 16 | `spike24` | `input_spike24.csv` | `output_spike24.txt` | 24 | 24 | 0 | Frequent sharp turns |
| 17 | `dual_hole_narrow_gap` | `input_dual_hole_narrow_gap.csv` | `output_dual_hole_narrow_gap.txt` | 12 | 12 | 2 | Two holes separated by a small gap |
| 18 | `offset_large_small` | `input_offset_large_small.csv` | `output_offset_large_small.txt` | 8 | 8 | 1 | Large base coordinates + small local feature |
| 19 | `zigzag40_band` | `input_zigzag40_band.csv` | `output_zigzag40_band.txt` | 42 | 42 | 0 | Long oscillating boundary |
| 20 | `concave_notch_10` | `input_concave_notch_10.csv` | `output_concave_notch_10.txt` | 10 | 10 | 0 | Concave notch / re-entrant angles |

## Why these complement the original provided suite

The instructor suite already covers moderate-hole examples and large natural lake boundaries. These 20 cases extend coverage for:
- strict lower-bound ring cardinality,
- tiny/skinny geometry,
- large-coordinate numeric sensitivity,
- higher hole counts,
- regular-shape high-n scaling anchors,
- negative-coordinate domains,
- repeated sharp-curvature patterns,
- narrow clearances between holes and around boundaries.

All expected outputs follow the assignment output format and include signed-area values in scientific notation.
