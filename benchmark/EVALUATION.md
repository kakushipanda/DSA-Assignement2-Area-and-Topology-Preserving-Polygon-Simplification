# Experimental Evaluation

## Method
- Synthetic regular polygons were generated at increasing sizes (64 to 4096 vertices).
- For each size, simplification was run to `target=n/2`; each point uses median of repeated runs (9 repeats for n<=512, otherwise 5).
- One warm-up run is executed before sampling to reduce cold-start noise.
- For displacement analysis, a 2048-vertex polygon was simplified at multiple targets; each point is a median of 5 runs.
- Runtime model fitting excludes n<256 to reduce tiny-input process-startup bias.
- Displacement fit is anchored at zero for the unsimplified endpoint: y = a*(1/x² - 1/x_max²).

## Best-fit scaling models
- Runtime best fit: **c*n**, c=3.708718e-06, b=1.694256e-03, R^2=0.9764.
- Memory best fit: **c*sqrt(n)**, c=9.771478e+01, b=-5.310725e+02, R^2=0.6790.

## Plots
- `benchmark/plots/runtime_vs_input_size.svg`
- `benchmark/plots/memory_vs_input_size.svg`
- `benchmark/plots/areal_displacement_vs_target.svg`

## Interpretation
- Runtime increases super-linearly with input size in this implementation due to global segment scanning in topology checks.
- Peak RSS grows with input size and reflects the linked-list vertex storage plus candidate queue overhead.
- Areal displacement generally increases as target vertices decrease, indicating stronger geometric change at more aggressive simplification levels.
