# Makefile
# Builds the simplifier and provides test/compare targets for both
# test_cases/ and experimental_cases/.
#
# Output naming convention:
#   input_foo.csv  ->  my_output_foo.txt   (same directory as input)
#
# Targets:
#   make                  build simplify + alias binary
#   make test             run test_cases/, write my_output_*.txt, diff vs output_*.txt
#   make test-extra       run experimental_cases/, write my_output_*.txt, diff vs output_*.txt
#   make test-all         both of the above
#   make clean            remove object files, binaries, and all my_output_*.txt files

# ---------------------------------------------------------------------------
# Compiler
# ---------------------------------------------------------------------------
CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pedantic -Isrc
LDFLAGS  :=

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------
SRC := $(filter-out src/main.cpp src/benchmark.cpp, $(wildcard src/*.cpp))

OBJ := $(SRC:.cpp=.o)
OBJ_MAIN := $(OBJ) src/main.o
OBJ_BENCHMARK := $(OBJ) src/benchmark.o

# ---------------------------------------------------------------------------
# Test cases  format: input_path:target_vertices
# The my_output and expected-output paths are derived automatically:
#   input dir  + "my_output_" + stem + ".txt"
#   input dir  + "output_"    + stem + ".txt"
# ---------------------------------------------------------------------------
TEST_CASES := \
	test_cases/input_rectangle_with_two_holes.csv:7 \
	test_cases/input_cushion_with_hexagonal_hole.csv:13 \
	test_cases/input_blob_with_two_holes.csv:17 \
	test_cases/input_wavy_with_three_holes.csv:21 \
	test_cases/input_lake_with_two_islands.csv:17 \
	test_cases/input_original_01.csv:99 \
	test_cases/input_original_02.csv:99 \
	test_cases/input_original_03.csv:99 \
	test_cases/input_original_04.csv:99 \
	test_cases/input_original_05.csv:99 \
	test_cases/input_original_06.csv:99 \
	test_cases/input_original_07.csv:99 \
	test_cases/input_original_08.csv:99 \
	test_cases/input_original_09.csv:99 \
	test_cases/input_original_10.csv:99

EXTRA_TEST_CASES := \
	experimental_cases/input_minimal_triangle.csv:3 \
	experimental_cases/input_axis_aligned_square.csv:4 \
	experimental_cases/input_triangle_hole_minimal.csv:7 \
	experimental_cases/input_four_holes_compact.csv:16 \
	experimental_cases/input_narrow_corridor.csv:8 \
	experimental_cases/input_near_collinear_large.csv:5 \
	experimental_cases/input_large_with_hole.csv:8 \
	experimental_cases/input_regular32.csv:32 \
	experimental_cases/input_regular128.csv:128 \
	experimental_cases/input_star20.csv:20 \
	experimental_cases/input_skinny_rectangle_6.csv:6 \
	experimental_cases/input_tiny_hole_clearance.csv:7 \
	experimental_cases/input_eight_holes_field.csv:28 \
	experimental_cases/input_regular256.csv:256 \
	experimental_cases/input_negative_coords_mixed.csv:5 \
	experimental_cases/input_spike24.csv:24 \
	experimental_cases/input_dual_hole_narrow_gap.csv:12 \
	experimental_cases/input_offset_large_small.csv:8 \
	experimental_cases/input_zigzag40_band.csv:42 \
	experimental_cases/input_concave_notch_10.csv:10

# ---------------------------------------------------------------------------
# Helper macro: run one test suite and diff against expected outputs.
#
# Usage:  $(call RUN_SUITE, <case-list>, <suite-label>)
#
# For each "dir/input_STEM.csv:TARGET" entry it:
#   1. Derives output path:   dir/my_output_STEM.txt
#   2. Derives expected path: dir/output_STEM.txt
#   3. Runs ./simplify and captures stdout to the output path.
#   4. Diffs output vs expected; prints PASS/FAIL per case.
#   5. Prints a summary line at the end.
# ---------------------------------------------------------------------------
define RUN_SUITE
	@echo ""; \
	echo "=== $(2) ==="; \
	passed=0; failed=0; \
	for entry in $(1); do \
		input_file=$${entry%%:*}; \
		target=$${entry##*:}; \
		dir=$$(dirname "$$input_file"); \
		basename=$$(basename "$$input_file"); \
		stem=$${basename#input_}; \
		stem=$${stem%.csv}; \
		my_output="$$dir/my_output_$$stem.txt"; \
		expected="$$dir/output_$$stem.txt"; \
		./simplify "$$input_file" "$$target" > "$$my_output"; \
		if [ ! -f "$$expected" ]; then \
			echo "  [SKIP]  $$input_file  (no expected output at $$expected)"; \
		elif diff -q --strip-trailing-cr "$$expected" "$$my_output" > /dev/null 2>&1; then \
			echo "  [PASS]  $$input_file -> $$my_output"; \
			passed=$$((passed + 1)); \
		else \
			echo "  [FAIL]  $$input_file -> $$my_output"; \
			diff --strip-trailing-cr "$$expected" "$$my_output" | head -30; \
			failed=$$((failed + 1)); \
		fi; \
	done; \
	echo ""; \
	echo "  $(2) results: $$passed passed, $$failed failed."; \
	if [ $$failed -ne 0 ]; then exit 1; fi
endef

# ---------------------------------------------------------------------------
# Build targets
# ---------------------------------------------------------------------------
.PHONY: all clean test test-extra test-all evaluate

all: simplify simplify_benchmark area_and_topology_preserving_polygon_simplification

simplify: $(OBJ_MAIN)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ_MAIN) $(LDFLAGS)

simplify_benchmark: $(OBJ_BENCHMARK)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ_BENCHMARK) $(LDFLAGS)

area_and_topology_preserving_polygon_simplification: simplify
	cp simplify area_and_topology_preserving_polygon_simplification
	cp simplify_benchmark benchmark

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ---------------------------------------------------------------------------
# Test targets
# ---------------------------------------------------------------------------

# Standard test suite only
test: simplify
	$(call RUN_SUITE,$(TEST_CASES),Standard test cases)

# Experimental test suite only
test-extra: simplify
	$(call RUN_SUITE,$(EXTRA_TEST_CASES),Experimental test cases)

# Both suites; exits non-zero if either has any failure
test-all: simplify
	$(call RUN_SUITE,$(TEST_CASES),Standard test cases)
	$(call RUN_SUITE,$(EXTRA_TEST_CASES),Experimental test cases)

# ---------------------------------------------------------------------------
# Evaluation (unchanged from original)
# ---------------------------------------------------------------------------
evaluate: simplify
	python3 scripts/run_experimental_evaluation.py

# ---------------------------------------------------------------------------
# Clean
# ---------------------------------------------------------------------------
clean:
	rm -f src/*.o \
	      simplify \
	      area_and_topology_preserving_polygon_simplification \
	      test_cases/my_output_*.txt \
	      experimental_cases/my_output_*.txt\
				simplify_benchmark \
				benchmark
