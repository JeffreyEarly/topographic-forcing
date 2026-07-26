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

## Goff abyssal-hill wave generation

Run the broadband variable-stratification example with:

```matlab
result = GoffAbyssalHillWaveGenerationExample;
```

The example uses a deterministic two-dimensional Goff terrain realization with $100$ m RMS height, an exponential $N^2(z)$ profile, and a prescribed $5$ cm/s M2 current. Its four panels show the terrain, realized and target terrain spectra, wave energy and cumulative bottom work, and the final energy distribution across horizontal wavenumber and vertical mode.

The default resolution is `[64 64 17]`, and the shortest terrain wavelength is $40$ km so that it is resolved with the default antialias layout. Set `resolution`, `minimumWavelength`, `numberOfOutputTimes`, `relativeTolerance`, or `shouldMakeFigures` to configure the run. No output files are created.

Generate the same terrain independently with:

```matlab
[virtualDepth,h,diagnostics] = ...
    WVBottomWaveGenerationForcing.goffAbyssalHillTopography( ...
        wvt,rmsHeight=100,minimumWavelength=40e3,randomSeed=2023);
```

## Performance benchmark

Profile construction and ordinary forcing calls with:

```matlab
benchmark = BottomWaveGenerationPerformanceBenchmark;
```

The returned table reports construction time, application time, active wave coefficients, and estimated stored response bytes at three resolutions. Timings are intended for comparisons on the same machine and are not fixed acceptance thresholds.

## Uniform-depth wave scattering

Run the autonomous-scattering frequency oracle with:

```matlab
result = UniformDepthWaveScatteringBenchmark;
```

A uniform upward-positive bottom offset changes the physical depth from $D$ to $D-h_0$. The example constructs the first-order physical coefficient generator, compares its eigenfrequency with the exact constant-depth dispersion relation at three values of $h_0$, and confirms the expected $O(h_0^2)$ difference. It also integrates a coupled eigenvector with adaptive `ode78` and compares the resulting wave coefficients with a matrix exponential.

The default figure shows the frequency correction and its second-order error. Use `shouldMakeFigures=false`, `resolution`, `topographicHeights`, or `relativeTolerance` to configure the benchmark.

## Sinusoidal-ridge wave scattering

Run a complete autonomous calculation with:

```matlab
result = SinusoidalRidgeWaveScatteringExample;
```

The example initializes one $+$-branch internal wave, scatters it from a sinusoidal ridge, and shows the incident mode, the $k\pm k_h$ sidebands, and the flat and first-order physical energy histories. A second figure shows the final bottom velocity and bottom displacement.

The model uses the adaptive integrator and standard WaveVortexModel NetCDF output. After the run, `WVBottomWaveScatteringForcing.bottomDisplacementFromFile` reconstructs the bottom displacement from the saved `wave-vortex` coefficients; no extra prognostic displacement is added to the model. Temporary output is deleted automatically. Pass `outputPath` to retain the file, or use `resolution`, `numberOfWavePeriods`, `numberOfOutputTimes`, `relativeTolerance`, `shouldAntialias`, and `shouldMakeFigures` to configure the run.
