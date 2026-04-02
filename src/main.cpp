// main.cpp
// Command-line front-end for the APSC polygon simplifier.
// Argument parsing and final output formatting only — all logic lives in
// polygon.cpp / geometry.hpp.
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>

#include "polygon.hpp"

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            std::cerr << "Usage: ./simplify <input_file.csv> <target_vertices>\n";
            return 1;
        }

        const std::string path      = argv[1];
        const int         target    = std::atoi(argv[2]);

        // Read input into global pool.
        const int total_in = simplify::read_polygon_csv(path);
        const double input_area = simplify::compute_total_signed_area();

        // Early exit: already at or below target.
        if (total_in <= target) {
            simplify::write_output();
            std::printf("Total signed area in input: %e\n",  input_area);
            std::printf("Total signed area in output: %e\n", input_area);
            std::printf("Total areal displacement: %e\n",    0.0);
            return 0;
        }

        // Run APSC simplification.
        const double displacement = simplify::simplify_polygon(target);
        const double output_area  = simplify::compute_total_signed_area();

        simplify::write_output();
        std::printf("Total signed area in input: %e\n",  input_area);
        std::printf("Total signed area in output: %e\n", output_area);
        std::printf("Total areal displacement: %e\n",    displacement);

        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
