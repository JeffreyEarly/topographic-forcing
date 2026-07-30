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

## Eddy-tide validation experiment

Run the 30-day constant-stratification eddy experiment and its matched no-eddy control with:

```matlab
[~,~,eddyFile] = EddyTideTopographicForcingSimulation;
[~,~,noEddyFile] = EddyTideTopographicForcingSimulation(includeEddy=false);
```

The setup retains the domain, shallow eddy, nonlinear advection, and damping of the JPO2026 minimal simulation. It removes the initialized fixed-amplitude wave beam. Both runs instead start with zero wave coefficients and use a one-M2-period ramp of a 5 cm/s zonal barotropic current over the same seed-2023 Goff terrain. `WVBottomWaveGenerationForcing` generates the tide and automatically excludes the exact nonzero support of `WVAdaptiveDamping`. Autonomous `WVBottomWaveScatteringForcing` is opt-in with `shouldUseScattering=true` and is not part of the default validation pair.

The defaults use `Nxy=128`, hourly output, $100$ m RMS terrain, corner wavenumber $10^{-4}\ \mathrm{m^{-1}}$, and a 20 km minimum terrain wavelength. The cutoff lies above the 17.84 km antialiased resolution limit for this grid. Hourly records resolve the twice-tidal terms used by the exact energy and APV-enstrophy diagnostics. Use `horizontalDomainSize` to replace the default four-mode-one-wavelength square domain; the vertical resolution is then selected by `WVStratification.verticalResolutionForHorizontalResolution`. Use `maxT`, `Nxy`, `outputInterval`, `minimumWavelength`, and the other documented options for shorter tests or longer continuations. An explicit `outputFilename` overrides the descriptive default name.

After both runs finish, build standard WaveVortexModelDiagnostics files and create the normalized energy figure with:

```matlab
[figureHandle,energy,diagnosticsFiles] = ...
    AnalyzeEddyTideEnergy(eddyFile,noEddyFile);
```

The upper panel compares wave energy and the lower panel compares geostrophic energy. Eddy and no-eddy curves share the eddy simulation's initial geostrophic energy as their normalization. The returned `energy` structure contains the raw and normalized time series, and the default PNG is saved beside the model output. Set `shouldExport=false` for an in-memory figure or specify `exportPath`, `exportResolution`, and `shouldOverwriteExisting`.

Create the complete basic figure set after a paired run with:

```matlab
[figures,summary] = AnalyzeEddyTideBasicFigures(eddyFile,noEddyFile);
```

This produces paired energy and maximum-wave-speed histories, final horizontal/vertical-mode wave spectra, and final surface wave/geostrophic vertical-vorticity maps. The MAT summary records the plotted values, extrema, integrated generation, spectral peaks, and output paths.

For forcing-resolved energy and potential-enstrophy pathways, run:

```matlab
[figures,budget,diagnosticsFiles] = AnalyzeEddyTideBudgets(eddyFile,noEddyFile);
```

The three figures compare reservoirs, cumulative forcing budgets, and nonlinear-triad pathways. The returned structure contains raw rates, cumulative integrals, integrated contributions, closure metrics, and eddy-minus-control differences. The implementation uses the lower-level WVDiagnostics flux APIs directly; it does not depend on the mirror-triad variables required by `summarizeSourcesSinksReservoirs`.

Create the Shakespeare-style radius-depth comparison of eddy velocity and adaptive-damping energy removal with:

```matlab
[figureHandle,analysis,diagnosticsFile] = AnalyzeEddyTideAdaptiveDamping(eddyFile);
```

The default averages days 25--30, uses the domain-centered eddy, and overlays the eddy-aligned \(u_\theta=0.02\ \mathrm{m\,s^{-1}}\) contour. The damping panel is the exact signed work of the spectral closure, not a pointwise non-negative viscous dissipation rate; see [the accompanying note](AnalyzeEddyTideAdaptiveDamping.md).

After the generation-only control is available, assess whether autonomous scattering is perturbative with:

```matlab
[benchmark,figureHandle] = GoffScatteringPerturbationBenchmark(noEddyFile);
```

The benchmark scales the same Goff realization from 100 down to 3.125 m RMS height by factors of two, verifies the expected first- and second-order exponents, and measures first-order physical-energy drift in short scattering-only integrations. Scattering remains opt-in even when a reduced terrain amplitude passes the reported gates.

The model output can be extended without rebuilding the experiment:

```matlab
model = WVModel.modelFromFile(char(eddyFile));
model.integrateToTime(60*86400);
model.closeNetCDFFile();
```

Rerun `AnalyzeEddyTideEnergy` to append diagnostics for the new records. WaveVortexModel 4.1.1 or newer is required for the simulation, and WaveVortexModelDiagnostics 1.0.6 or newer is required for the analysis.

This 30-day run is intended to validate the new infrastructure before longer experiments. It follows the JPO2026 WaveVortexModel configuration rather than literally reproducing the $500$ m MITgCM grid in [Shakespeare (2023)](https://doi.org/10.1175/JPO-D-23-0127.1).

The repository-level [README](../README.md#500-km-shakespeare-comparison-pair) records the exact matched $500$ km, $256^2\times45$ Shakespeare-comparison pair that was run first to day 30 and then restarted to day 50. It includes the two process-ready simulation calls, exact output filenames, reduced $0.05/\sqrt{10}\ \mathrm{m\,s^{-1}}$ tidal amplitude, 6 km terrain cutoff, diagnostics command, and day-50 cross-machine reference values.

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

## Gaussian-ridge wave scattering and movie

Run a localized wave-packet calculation with:

```matlab
result = GaussianRidgeWaveScatteringExample;
```

The example retains the constant-$N$, mode-one M2 configuration and six-wavelength domain of the earlier Gaussian-ridge benchmark, but uses a gentle first-order ridge with default height and criticality

```math
h_0/D=0.05,
\qquad
\max |h_x|/\mu=0.05.
```

It initializes a rightward mode-one packet, removes nonlinear advection, registers only `WVBottomWaveScatteringForcing`, and advances the model with adaptive `ode78`. The returned result contains every sampled `Ap`, `Am`, and `A0`, modal and energy diagnostics, and selected $x$--$z$ snapshots. The two figures show the ridge and packet evolution, reflected and higher-mode energy, and the flat versus first-order physical energy histories.

No file is written by default. To retain standard restartable WaveVortexModel output at the same cadence, supply a path:

```matlab
result = GaussianRidgeWaveScatteringExample( ...
    outputPath="gaussian-ridge.nc", ...
    shouldOverwriteExisting=false);
```

Render the saved records independently with:

```matlab
movie = GaussianRidgeWaveScatteringMovie( ...
    "gaussian-ridge.nc", ...
    field="u", ...
    frameRate=15);
```

The renderer reconstructs `"u"`, `"w"`, or `"eta"` from the saved wave--vortex coefficients. It masks values beneath the physical bottom and determines its fixed symmetric color scale from wet points only. The large $x$--$z$ section uses a labelled 20-times vertical exaggeration by default; set `verticalExaggeration` to another positive value when needed. The companion panels show the first-order depth-integrated wave-energy profile, including its bottom correction, and the modal energy histories. The MPEG-4 and PNG poster names are derived from the NetCDF path unless `videoPath` or `posterPath` is supplied. Use `iTime` for a strictly increasing subset of records and `forcingName` when the file contains multiple scattering forcings.

The poster is selected at the saved time with maximum reflected-plus-higher-mode energy. The returned movie metadata includes all output paths, selected time indices, frame count and dimensions, duration, field limits, dominant period, first-mode wavelength, packet group velocity, and poster index.
