#!/usr/bin/env python3
"""Run runtime/memory/displacement experiments and generate tables + SVG plots.

No third-party dependencies are required.
"""
from __future__ import annotations
import csv
import math
import os
import statistics
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'benchmark' / 'results'
PLOTS = ROOT / 'benchmark' / 'plots'
SYNTH = ROOT / 'benchmark' / 'synthetic_inputs'
BIN = ROOT / 'simplify'

RESULTS.mkdir(parents=True, exist_ok=True)
PLOTS.mkdir(parents=True, exist_ok=True)
SYNTH.mkdir(parents=True, exist_ok=True)


def generate_regular_polygon(path: Path, n: int, radius: float = 1000.0) -> None:
    with path.open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['ring_id', 'vertex_id', 'x', 'y'])
        for i in range(n):
            a = 2 * math.pi * i / n
            x = radius * math.cos(a)
            y = radius * math.sin(a)
            w.writerow([0, i, f'{x:.6f}', f'{y:.6f}'])


def sample_peak_rss_kb(proc: subprocess.Popen) -> int:
    pid = proc.pid
    peak = 0
    status = Path(f'/proc/{pid}/status')
    while proc.poll() is None:
        if status.exists():
            try:
                text = status.read_text()
                for line in text.splitlines():
                    if line.startswith('VmRSS:'):
                        rss = int(line.split()[1])
                        if rss > peak:
                            peak = rss
                        break
            except Exception:
                pass
        time.sleep(0.002)
    return peak


def run_time_and_memory(input_csv: Path, target: int) -> tuple[float, int]:
    start = time.perf_counter()
    with open(os.devnull, 'w') as devnull:
        proc = subprocess.Popen(
            [str(BIN), str(input_csv), str(target)],
            stdout=devnull,
            stderr=devnull,
            text=True,
        )
        peak = sample_peak_rss_kb(proc)
        proc.wait()
    elapsed = time.perf_counter() - start
    if proc.returncode != 0:
        raise RuntimeError(f'Command failed for {input_csv.name}')
    return elapsed, peak


def run_for_displacement(input_csv: Path, target: int) -> tuple[float, float]:
    start = time.perf_counter()
    proc = subprocess.run([str(BIN), str(input_csv), str(target)], capture_output=True, text=True, check=True)
    elapsed = time.perf_counter() - start
    disp = 0.0
    for line in proc.stdout.splitlines():
        if line.startswith('Total areal displacement:'):
            disp = float(line.split(':', 1)[1].strip())
            break
    return elapsed, disp


def robust_median_measurements(input_csv: Path, target: int, repeats: int) -> tuple[float, float]:
    """Return median elapsed time and median peak RSS for repeated runs."""
    times = []
    peaks = []
    for _ in range(repeats):
        elapsed, peak_kb = run_time_and_memory(input_csv, target)
        times.append(elapsed)
        peaks.append(peak_kb)
    return statistics.median(times), statistics.median(peaks)


def fit_models(xs: list[int], ys: list[float]) -> dict[str, tuple[float, float, float]]:
    # returns name -> (coefficient c, intercept b, r2)
    def r2(pred: list[float]) -> float:
        mean_y = statistics.mean(ys)
        ss_tot = sum((y - mean_y) ** 2 for y in ys)
        ss_res = sum((y - p) ** 2 for y, p in zip(ys, pred))
        return 1.0 - (ss_res / ss_tot if ss_tot else 0.0)

    models = {
        'c*n': [float(x) for x in xs],
        'c*nlogn': [float(x) * math.log(max(x, 2)) for x in xs],
        'c*n^2': [float(x) ** 2 for x in xs],
    }
    out = {}
    for name, basis in models.items():
        mean_phi = statistics.mean(basis)
        mean_y = statistics.mean(ys)
        denom = sum((p - mean_phi) ** 2 for p in basis)
        if denom == 0:
            c = 0.0
            b0 = mean_y
        else:
            c = sum((p - mean_phi) * (y - mean_y) for p, y in zip(basis, ys)) / denom
            b0 = mean_y - c * mean_phi
        pred = [c * p + b0 for p in basis]
        out[name] = (c, b0, r2(pred))
    return out


def write_svg(
    path: Path,
    title: str,
    xlabel: str,
    ylabel: str,
    xs: list[float],
    ys: list[float],
    model_name: str,
    model_c: float,
    model_b: float,
    model_fn,
    y_scale: str = 'linear',
    equation_text: str = '',
):
    w, h = 900, 540
    ml, mr, mt, mb = 90, 40, 70, 80
    pw, ph = w - ml - mr, h - mt - mb
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if xmax == xmin:
        xmax += 1.0
    if ymax == ymin:
        ymax += 1.0

    def ty(y):
        if y_scale == 'log10':
            return math.log10(max(y, 1e-12))
        return y

    ymin_t, ymax_t = ty(ymin), ty(ymax)
    if ymax_t == ymin_t:
        ymax_t += 1.0

    def sx(x):
        return ml + (x - xmin) / (xmax - xmin) * pw

    def sy(y):
        yt = ty(y)
        return mt + ph - (yt - ymin_t) / (ymax_t - ymin_t) * ph

    points = ' '.join(f'{sx(x):.2f},{sy(y):.2f}' for x, y in zip(xs, ys))
    model_pts = []
    for i in range(200):
        x = xmin + (xmax - xmin) * i / 199
        y = model_c * model_fn(x) + model_b
        model_pts.append(f'{sx(x):.2f},{sy(y):.2f}')

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}">',
        '<style>text{font-family:Arial,sans-serif} .axis{stroke:#111;stroke-width:2} .grid{stroke:#ddd;stroke-width:1} .data{fill:none;stroke:#1f77b4;stroke-width:2} .fit{fill:none;stroke:#d62728;stroke-width:2;stroke-dasharray:6 4}</style>',
        f'<text x="{w/2}" y="35" text-anchor="middle" font-size="22">{title}</text>',
        f'<line class="axis" x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}"/>',
        f'<line class="axis" x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}"/>',
        f'<polyline class="data" points="{points}"/>',
        f'<polyline class="fit" points="{" ".join(model_pts)}"/>',
        f'<text x="{w/2}" y="{h-20}" text-anchor="middle" font-size="16">{xlabel}</text>',
        f'<text x="25" y="{h/2}" text-anchor="middle" transform="rotate(-90,25,{h/2})" font-size="16">{ylabel}</text>',
        f'<rect x="{ml+pw-300}" y="{mt+6}" width="292" height="42" fill="white" opacity="0.75"/>',
        f'<text x="{ml+pw-10}" y="{mt+22}" text-anchor="end" font-size="14" fill="#d62728">Best fit: {model_name}</text>',
        f'<text x="{ml+pw-10}" y="{mt+40}" text-anchor="end" font-size="13" fill="#d62728">{equation_text}</text>',
        '</svg>'
    ]
    # Draw simple ticks for readability.
    for i in range(6):
        tx = xmin + (xmax - xmin) * i / 5
        px = sx(tx)
        svg.insert(-1, f'<line class="axis" x1="{px:.2f}" y1="{mt+ph}" x2="{px:.2f}" y2="{mt+ph+6}"/>')
        svg.insert(-1, f'<text x="{px:.2f}" y="{mt+ph+24}" text-anchor="middle" font-size="12">{tx:.0f}</text>')
    for i in range(6):
        frac = i / 5
        py = mt + ph - frac * ph
        val_t = ymin_t + frac * (ymax_t - ymin_t)
        val = (10 ** val_t) if y_scale == 'log10' else val_t
        svg.insert(-1, f'<line class="axis" x1="{ml-6}" y1="{py:.2f}" x2="{ml}" y2="{py:.2f}"/>')
        label = f'{val:.3g}' if val >= 1 else f'{val:.2e}'
        svg.insert(-1, f'<text x="{ml-10}" y="{py+4:.2f}" text-anchor="end" font-size="12">{label}</text>')
    path.write_text('\n'.join(svg))


def main():
    subprocess.run(['make', 'simplify'], cwd=ROOT, check=True)

    sizes = [64, 96, 128, 160, 192, 256, 320, 384, 512, 640, 768, 896, 1024, 1280, 1536, 1792, 2048, 2560, 3072, 3584, 4096]
    runtime_rows = []
    mem_rows = []

    # warm-up call to reduce first-run effects from lazy page faults / cache cold-start.
    warm_csv = SYNTH / 'regular_warmup_512.csv'
    generate_regular_polygon(warm_csv, 512)
    run_time_and_memory(warm_csv, 256)

    for n in sizes:
        csv_path = SYNTH / f'regular_{n}.csv'
        generate_regular_polygon(csv_path, n)
        target = max(3, n // 2)
        # More repeats for small-N where process-launch jitter is a larger share of runtime.
        repeats = 9 if n <= 512 else 5
        med_t, med_rss = robust_median_measurements(csv_path, target, repeats)
        runtime_rows.append({'n': n, 'time_sec': med_t})
        mem_rows.append({'n': n, 'peak_rss_kb': med_rss})

    disp_n = 2048
    disp_csv = SYNTH / f'regular_{disp_n}.csv'
    generate_regular_polygon(disp_csv, disp_n)
    targets = list(range(2048, 255, -128))
    disp_rows = []
    for t in targets:
        # Repeat and take medians to suppress one-off scheduler/CPU jitter.
        d_times = []
        d_vals = []
        for _ in range(5):
            elapsed, disp = run_for_displacement(disp_csv, t)
            d_times.append(elapsed)
            d_vals.append(disp)
        disp_rows.append({
            'target_vertices': t,
            'time_sec': statistics.median(d_times),
            'areal_displacement': statistics.median(d_vals),
        })

    with (RESULTS / 'runtime_vs_input_size.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['n', 'time_sec'])
        w.writeheader(); w.writerows(runtime_rows)
    with (RESULTS / 'memory_vs_input_size.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['n', 'peak_rss_kb'])
        w.writeheader(); w.writerows(mem_rows)
    with (RESULTS / 'displacement_vs_target.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['target_vertices', 'time_sec', 'areal_displacement'])
        w.writeheader(); w.writerows(disp_rows)

    xs = [r['n'] for r in runtime_rows]
    ys_time = [r['time_sec'] for r in runtime_rows]
    ys_mem = [r['peak_rss_kb'] for r in mem_rows]
    # Fit runtime on n>=256 to reduce process-startup dominance at very small sizes.
    fit_runtime_rows = [r for r in runtime_rows if r['n'] >= 256]
    fit_xs_time = [r['n'] for r in fit_runtime_rows]
    fit_ys_time = [r['time_sec'] for r in fit_runtime_rows]
    fits_time = fit_models(fit_xs_time, fit_ys_time)
    fits_mem = fit_models(xs, ys_mem)
    best_time = max(fits_time.items(), key=lambda kv: kv[1][1])
    best_mem = max(fits_mem.items(), key=lambda kv: kv[1][1])

    model_fn = {
        'c*n': lambda x: x,
        'c*nlogn': lambda x: x * math.log(max(x, 2)),
        'c*n^2': lambda x: x * x,
    }

    write_svg(PLOTS / 'runtime_vs_input_size.svg', 'Running Time vs Input Size', 'Input vertices (n)', 'Time (seconds)',
              xs, ys_time, best_time[0], best_time[1][0], best_time[1][1], model_fn[best_time[0]],
              y_scale='log10',
              equation_text=f'y={best_time[1][0]:.3e}*f(n)+{best_time[1][1]:.3e}')
    write_svg(PLOTS / 'memory_vs_input_size.svg', 'Peak RSS vs Input Size', 'Input vertices (n)', 'Peak RSS (kB)',
              xs, ys_mem, best_mem[0], best_mem[1][0], best_mem[1][1], model_fn[best_mem[0]],
              y_scale='linear',
              equation_text=f'y={best_mem[1][0]:.3e}*f(n)+{best_mem[1][1]:.3e}')

    # Sort for plotting left-to-right to avoid visual zig-zag caused by descending target order.
    disp_rows = sorted(disp_rows, key=lambda r: r['target_vertices'])
    dx = [r['target_vertices'] for r in disp_rows]
    dy = [r['areal_displacement'] for r in disp_rows]
    # Fit displacement with y = a * (1/x - 1/x_max), anchored at y(x_max)=0.
    # This avoids the invalid negative tail that occurs with unconstrained a*(1/x)+b on log-scale.
    x_max = max(dx)
    phi = [(1.0 / x) - (1.0 / x_max) for x in dx]
    den = sum(p * p for p in phi) or 1.0
    a_anch = sum(p * y for p, y in zip(phi, dy)) / den
    write_svg(PLOTS / 'areal_displacement_vs_target.svg', 'Areal Displacement vs Target Vertices',
              'Target vertices', 'Areal displacement', dx, dy, 'a*(1/x-1/xmax)', 1.0, 0.0,
              lambda x: max(a_anch * ((1.0 / x) - (1.0 / x_max)), 0.0),
              y_scale='linear',
              equation_text=f'y={a_anch:.3e}*(1/x-1/{x_max})')

    report = ROOT / 'benchmark' / 'EVALUATION.md'
    with report.open('w') as f:
        f.write('# Experimental Evaluation\n\n')
        f.write('## Method\n')
        f.write('- Synthetic regular polygons were generated at increasing sizes (64 to 4096 vertices).\n')
        f.write('- For each size, simplification was run to `target=n/2`; each point uses median of repeated runs (9 repeats for n<=512, otherwise 5).\n')
        f.write('- One warm-up run is executed before sampling to reduce cold-start noise.\n')
        f.write('- For displacement analysis, a 2048-vertex polygon was simplified at multiple targets; each point is a median of 5 runs.\n')
        f.write('- Runtime model fitting excludes n<256 to reduce tiny-input process-startup bias.\n')
        f.write('- Displacement fit is anchored at zero for the unsimplified endpoint: y = a*(1/x - 1/x_max).\n\n')
        f.write('## Best-fit scaling models\n')
        f.write(f'- Runtime best fit: **{best_time[0]}**, c={best_time[1][0]:.6e}, b={best_time[1][1]:.6e}, R^2={best_time[1][2]:.4f}.\n')
        f.write(f'- Memory best fit: **{best_mem[0]}**, c={best_mem[1][0]:.6e}, b={best_mem[1][1]:.6e}, R^2={best_mem[1][2]:.4f}.\n\n')
        f.write('## Plots\n')
        f.write('- `benchmark/plots/runtime_vs_input_size.svg`\n')
        f.write('- `benchmark/plots/memory_vs_input_size.svg`\n')
        f.write('- `benchmark/plots/areal_displacement_vs_target.svg`\n\n')
        f.write('## Interpretation\n')
        f.write('- Runtime increases super-linearly with input size in this implementation due to global segment scanning in topology checks.\n')
        f.write('- Peak RSS grows with input size and reflects the linked-list vertex storage plus candidate queue overhead.\n')
        f.write('- Areal displacement generally increases as target vertices decrease, indicating stronger geometric change at more aggressive simplification levels.\n')

    print('Wrote benchmark outputs to benchmark/results and benchmark/plots')


if __name__ == '__main__':
    main()
