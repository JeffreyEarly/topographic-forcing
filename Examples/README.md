# Examples

## Exact sinusoidal-ridge wave generation

Run the validated prescribed-generation example with:

```matlab
benchmark = SinusoidalRidgeWaveGenerationBenchmark;
```

The example drives a $50$ m sinusoidal ridge with a spatially uniform, $5$ cm/s M2 current for one tidal period. It compares adaptive WaveVortexModel integrations with the exact forced interaction coefficients and returns all numerical and scientific diagnostics in `benchmark`.

The first figure shows:

- convergence toward the exact modal coefficients;
- wave energy and cumulative bottom pressure work;
- energy in the two wave branches; and
- the final vertical-mode energy distribution.

The second figure shows the terrain, local bottom work density, vertical velocity, and displacement at the time of maximum exact wave energy.

The interactive default uses resolution `[32 4 17]`. Use `resolution=[8 4 5]` for the inexpensive automated-reference case, `shouldMakeFigures=false` to skip graphics, or set `relativeTolerances` to choose the adaptive integrations. The example creates no output files.
