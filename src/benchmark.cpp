#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#include <sys/time.h>
#elif defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <psapi.h>
#endif

#include "polygon.hpp"

namespace fs = std::filesystem;

namespace {

struct Options {
    fs::path input_dir;
    fs::path output_csv = "results.csv";
    int iterations = 5;
    int warmup_runs = 1;
    int jobs = 1;
    int verbose = 1;
};

static std::string os_string() {
#if defined(_WIN32)
    return "Windows";
#elif defined(__APPLE__)
    return "macOS";
#elif defined(__linux__)
    return "Linux";
#elif defined(__unix__)
    return "Unix";
#else
    return "UnknownOS";
#endif
}

static std::string now_utc_iso8601() {
    using clock = std::chrono::system_clock;
    const auto now = clock::now();
    const std::time_t t = clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
    return oss.str();
}

static void log_msg(int verbose, int level, const std::string& s) {
    if (verbose >= level) std::cerr << s << '\n';
}

static bool is_csv_file(const fs::path& p) {
    if (!p.has_extension()) return false;
    auto ext = p.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return ext == ".csv";
}

static std::string csv_escape(const std::string& s) {
    bool need_quotes = false;
    for (char c : s) {
        if (c == '"' || c == ',' || c == '\n' || c == '\r') { need_quotes = true; break; }
    }
    if (!need_quotes) return s;
    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('"');
    for (char c : s) {
        if (c == '"') out.push_back('"');
        out.push_back(c);
    }
    out.push_back('"');
    return out;
}

static bool is_finite(double v) {
    return std::isfinite(v) != 0;
}

static double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    if (n % 2 == 1) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

static double mean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    const double sum = std::accumulate(v.begin(), v.end(), 0.0);
    return sum / static_cast<double>(v.size());
}

static double stdev_sample(const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    const double m = mean(v);
    double acc = 0.0;
    for (double x : v) {
        const double d = x - m;
        acc += d * d;
    }
    return std::sqrt(acc / static_cast<double>(v.size() - 1));
}

static long long peak_rss_kb() {
#if defined(__unix__) || defined(__APPLE__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) return -1;
#if defined(__APPLE__)
    return static_cast<long long>(usage.ru_maxrss / 1024);
#else
    return static_cast<long long>(usage.ru_maxrss);
#endif
#elif defined(_WIN32)
    PROCESS_MEMORY_COUNTERS_EX pmc{};
    if (!GetProcessMemoryInfo(GetCurrentProcess(),
                              reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&pmc),
                              sizeof(pmc))) {
        return -1;
    }
    return static_cast<long long>(pmc.PeakWorkingSetSize / 1024);
#else
    return -1;
#endif
}

static int clamp_target(int original, int target) {
    if (original <= 0) return 0;
    if (target < 1) target = 1;
    if (target > original) target = original;
    return target;
}

static int pct_target_round(int original, double pct) {
    const double raw = pct * static_cast<double>(original);
    if (!is_finite(raw)) return clamp_target(original, 1);
    const long long rounded = static_cast<long long>(std::llround(raw));
    return clamp_target(original, static_cast<int>(rounded));
}

struct RunStats {
    std::vector<double> times_ms;
    std::vector<long long> peak_kb_delta;
    std::vector<double> displacement;
};

struct AggregatedRow {
    std::string filename;
    int original_vertices = 0;
    int target_vertices = 0;
    double target_percentage = 0.0;
    double mean_time_ms = 0.0;
    double median_time_ms = 0.0;
    double stdev_time_ms = 0.0;
    long long peak_mem_kb = 0;
    double mean_displacement = 0.0;
    bool ok = true;
    std::string error;
};

static std::mutex g_simplify_mutex;

static std::optional<RunStats> benchmark_one_target(
    const fs::path& csv_path,
    int target_vertices,
    int warmup_runs,
    int iterations,
    int verbose,
    std::string& out_error
) {
    RunStats stats;
    stats.times_ms.reserve(static_cast<std::size_t>(iterations));
    stats.peak_kb_delta.reserve(static_cast<std::size_t>(iterations));
    stats.displacement.reserve(static_cast<std::size_t>(iterations));

    auto do_one = [&](bool record) -> bool {
        try {
            std::lock_guard<std::mutex> lock(g_simplify_mutex);

            const int total_in = simplify::read_polygon_csv(csv_path.string());
            if (total_in <= 0) {
                out_error = "read_polygon_csv returned non-positive vertex count";
                return false;
            }

            const long long base_peak = peak_rss_kb();
            const auto t0 = std::chrono::high_resolution_clock::now();
            const double disp = simplify::simplify_polygon(target_vertices);
            const auto t1 = std::chrono::high_resolution_clock::now();
            const long long after_peak = peak_rss_kb();

            const std::chrono::duration<double, std::milli> dt = t1 - t0;
            const double ms = dt.count();

            if (!is_finite(ms) || ms < 0.0) {
                out_error = "Non-finite or negative timing observed";
                return false;
            }
            if (!is_finite(disp) || disp < 0.0) {
                out_error = "Non-finite or negative areal displacement observed";
                return false;
            }

            long long delta_kb = -1;
            if (base_peak >= 0 && after_peak >= 0) {
                delta_kb = after_peak - base_peak;
                if (delta_kb < 0) delta_kb = 0;
            }

            if (record) {
                stats.times_ms.push_back(ms);
                stats.displacement.push_back(disp);
                stats.peak_kb_delta.push_back(delta_kb);
            }

            return true;
        } catch (const std::exception& e) {
            out_error = e.what();
            return false;
        } catch (...) {
            out_error = "Unknown exception";
            return false;
        }
    };

    for (int i = 0; i < warmup_runs; ++i) {
        if (!do_one(false)) return std::nullopt;
    }

    for (int i = 0; i < iterations; ++i) {
        if (!do_one(true)) return std::nullopt;
        if (verbose >= 2) {
            std::ostringstream oss;
            oss << "  iter " << (i + 1) << "/" << iterations
                << " time_ms=" << std::fixed << std::setprecision(3) << stats.times_ms.back()
                << " disp=" << std::scientific << std::setprecision(6) << stats.displacement.back()
                << " peak_delta_kb=" << stats.peak_kb_delta.back();
            log_msg(verbose, 2, oss.str());
        }
    }

    return stats;
}

static AggregatedRow aggregate(
    const fs::path& csv_path,
    int original_vertices,
    int target_vertices,
    double pct,
    const RunStats& stats
) {
    AggregatedRow row;
    row.filename = csv_path.string();
    row.original_vertices = original_vertices;
    row.target_vertices = target_vertices;
    row.target_percentage = pct;

    row.mean_time_ms = mean(stats.times_ms);
    row.median_time_ms = median(stats.times_ms);
    row.stdev_time_ms = stdev_sample(stats.times_ms);
    row.mean_displacement = mean(stats.displacement);

    long long peak = 0;
    for (long long v : stats.peak_kb_delta) {
        if (v > peak) peak = v;
    }
    row.peak_mem_kb = peak;

    if (!is_finite(row.mean_time_ms) || !is_finite(row.stdev_time_ms) ||
        !is_finite(row.mean_displacement) || !is_finite(row.target_percentage)) {
        row.ok = false;
        row.error = "Non-finite aggregated metric";
    }
    return row;
}

static std::vector<fs::path> discover_csv_files(const fs::path& root, int verbose) {
    std::vector<fs::path> files;
    std::error_code ec;

    if (!fs::exists(root, ec)) {
        throw std::runtime_error("Input directory does not exist: " + root.string());
    }
    if (!fs::is_directory(root, ec)) {
        throw std::runtime_error("Input path is not a directory: " + root.string());
    }

    fs::directory_options opts = fs::directory_options::skip_permission_denied;
    try {
        for (fs::recursive_directory_iterator it(root, opts, ec), end; it != end; it.increment(ec)) {
            if (ec) {
                log_msg(verbose, 1, "Filesystem warning: " + ec.message());
                ec.clear();
                continue;
            }
            const auto& entry = *it;
            if (!entry.is_regular_file(ec)) { ec.clear(); continue; }
            const fs::path p = entry.path();
            if (!is_csv_file(p)) continue;
            files.push_back(p);
        }
    } catch (const fs::filesystem_error& e) {
        throw std::runtime_error(std::string("Filesystem error: ") + e.what());
    }

    std::sort(files.begin(), files.end());
    return files;
}

static Options parse_args(int argc, char** argv) {
    Options opt;
    if (argc < 2) {
        throw std::runtime_error(
            "Usage: benchmark <input_dir> [--out results.csv] [--iters N] [--warmup N] [--jobs N] [--verbose 0|1|2]");
    }
    opt.input_dir = fs::path(argv[1]);

    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        auto need_value = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + flag);
            return argv[++i];
        };

        if (a == "--out") {
            opt.output_csv = fs::path(need_value(a));
        } else if (a == "--iters") {
            opt.iterations = std::stoi(need_value(a));
        } else if (a == "--warmup") {
            opt.warmup_runs = std::stoi(need_value(a));
        } else if (a == "--jobs") {
            opt.jobs = std::stoi(need_value(a));
        } else if (a == "--verbose") {
            opt.verbose = std::stoi(need_value(a));
        } else {
            throw std::runtime_error("Unknown argument: " + a);
        }
    }

    if (opt.iterations < 5) opt.iterations = 5;
    if (opt.warmup_runs < 0) opt.warmup_runs = 0;
    if (opt.jobs < 1) opt.jobs = 1;
    if (opt.verbose < 0) opt.verbose = 0;
    if (opt.verbose > 2) opt.verbose = 2;

    return opt;
}

static void write_results_csv(
    const fs::path& out_path,
    const std::vector<AggregatedRow>& rows,
    const Options& opt,
    const std::vector<double>& target_pcts
) {
    std::ofstream fout(out_path, std::ios::binary);
    if (!fout.is_open()) {
        throw std::runtime_error("Cannot open output file for writing: " + out_path.string());
    }

    fout << "# benchmark.cpp results\n";
    fout << "# timestamp_utc=" << now_utc_iso8601() << "\n";
    fout << "# os=" << os_string() << "\n";
    fout << "# hw_threads=" << std::thread::hardware_concurrency() << "\n";
    fout << "# input_dir=" << csv_escape(opt.input_dir.string()) << "\n";
    fout << "# iters=" << opt.iterations << ", warmup=" << opt.warmup_runs << ", jobs=" << opt.jobs << "\n";
    fout << "# targets=";
    for (std::size_t i = 0; i < target_pcts.size(); ++i) {
        if (i) fout << ",";
        fout << std::fixed << std::setprecision(2) << (target_pcts[i] * 100.0);
    }
    fout << "\n";

    fout << "Filename,Original_Vertices,Target_Vertices,Target_Percentage,Execution_Time_MS,Peak_Memory_KB,Areal_Displacement,Std_Deviation_MS\n";

    std::size_t ok_count = 0;
    double sum_time = 0.0;
    double sum_peak = 0.0;
    double sum_disp = 0.0;

    for (const auto& r : rows) {
        fout
            << csv_escape(r.filename) << ','
            << r.original_vertices << ','
            << r.target_vertices << ','
            << std::fixed << std::setprecision(2) << (r.target_percentage * 100.0) << ','
            << std::fixed << std::setprecision(3) << r.mean_time_ms << ','
            << r.peak_mem_kb << ','
            << std::scientific << std::setprecision(10) << r.mean_displacement << ','
            << std::fixed << std::setprecision(3) << r.stdev_time_ms
            << "\n";

        if (r.ok) {
            ++ok_count;
            sum_time += r.mean_time_ms;
            if (r.peak_mem_kb >= 0) sum_peak += static_cast<double>(r.peak_mem_kb);
            sum_disp += r.mean_displacement;
        }
    }

    fout << "# summary_total_rows=" << rows.size() << "\n";
    fout << "# summary_ok_rows=" << ok_count << "\n";
    if (ok_count > 0) {
        fout << "# summary_avg_mean_time_ms=" << std::fixed << std::setprecision(3) << (sum_time / ok_count) << "\n";
        fout << "# summary_avg_peak_mem_kb=" << std::fixed << std::setprecision(3) << (sum_peak / ok_count) << "\n";
        fout << "# summary_avg_areal_displacement=" << std::scientific << std::setprecision(10) << (sum_disp / ok_count) << "\n";
        fout << "# summary_mem_efficiency_ratio_kb_per_ms="
             << std::fixed << std::setprecision(6) << ((sum_time > 0.0) ? (sum_peak / sum_time) : 0.0) << "\n";
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options opt = parse_args(argc, argv);
        const std::vector<double> target_pcts = {0.49, 0.22, 0.15, 0.07};

        log_msg(opt.verbose, 1, "Discovering CSV files under: " + opt.input_dir.string());
        const std::vector<fs::path> files = discover_csv_files(opt.input_dir, opt.verbose);
        if (files.empty()) {
            throw std::runtime_error("No .csv files found under: " + opt.input_dir.string());
        }
        log_msg(opt.verbose, 1, "Found " + std::to_string(files.size()) + " CSV file(s).");

        std::vector<AggregatedRow> all_rows;
        all_rows.reserve(files.size() * target_pcts.size());
        std::mutex rows_mutex;

        std::atomic<std::size_t> next_idx{0};

        auto worker = [&](int tid) {
            (void)tid;
            for (;;) {
                const std::size_t idx = next_idx.fetch_add(1);
                if (idx >= files.size()) break;

                const fs::path& csv_path = files[idx];
                if (opt.verbose >= 1) {
                    log_msg(opt.verbose, 1, "Processing: " + csv_path.string());
                }

                int total_in = 0;
                {
                    try {
                        std::lock_guard<std::mutex> lock(g_simplify_mutex);
                        total_in = simplify::read_polygon_csv(csv_path.string());
                    } catch (const std::exception& e) {
                        AggregatedRow err;
                        err.filename = csv_path.string();
                        err.ok = false;
                        err.error = e.what();
                        std::lock_guard<std::mutex> lk(rows_mutex);
                        all_rows.push_back(err);
                        continue;
                    } catch (...) {
                        AggregatedRow err;
                        err.filename = csv_path.string();
                        err.ok = false;
                        err.error = "Unknown exception during read";
                        std::lock_guard<std::mutex> lk(rows_mutex);
                        all_rows.push_back(err);
                        continue;
                    }
                }

                if (total_in <= 0) {
                    AggregatedRow err;
                    err.filename = csv_path.string();
                    err.ok = false;
                    err.error = "Invalid vertex count from read_polygon_csv";
                    std::lock_guard<std::mutex> lk(rows_mutex);
                    all_rows.push_back(err);
                    continue;
                }

                for (double pct : target_pcts) {
                    const int target = pct_target_round(total_in, pct);

                    std::string err;
                    auto stats_opt = benchmark_one_target(
                        csv_path,
                        target,
                        opt.warmup_runs,
                        opt.iterations,
                        opt.verbose,
                        err
                    );

                    AggregatedRow row;
                    if (!stats_opt.has_value()) {
                        row.filename = csv_path.string();
                        row.original_vertices = total_in;
                        row.target_vertices = target;
                        row.target_percentage = pct;
                        row.ok = false;
                        row.error = err.empty() ? "Benchmark run failed" : err;
                    } else {
                        row = aggregate(csv_path, total_in, target, pct, *stats_opt);
                    }

                    std::lock_guard<std::mutex> lk(rows_mutex);
                    all_rows.push_back(std::move(row));
                }
            }
        };

        const int jobs = opt.jobs;
        std::vector<std::thread> threads;
        threads.reserve(static_cast<std::size_t>(jobs));
        for (int t = 0; t < jobs; ++t) threads.emplace_back(worker, t);
        for (auto& th : threads) th.join();

        for (const auto& r : all_rows) {
            if (!r.ok && opt.verbose >= 1) {
                log_msg(opt.verbose, 1, "Error: " + r.filename + " : " + r.error);
            }
        }

        write_results_csv(opt.output_csv, all_rows, opt, target_pcts);
        log_msg(opt.verbose, 1, "Wrote: " + opt.output_csv.string());
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}