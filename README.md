# Mean-depth bottom wave generation

`topographic-forcing` provides a fast, first-order bottom wave generator for WaveVortexModel. The formulation retains the ordinary rigid-lid wave--vortex basis and represents weak topography through the mean-depth bottom condition.

For a prescribed, horizontally uniform barotropic current, the bottom velocity is

```math
g_b(\boldsymbol{x},t)
=
\boldsymbol{U}_{\mathrm{bt}}(t)\boldsymbol{\cdot}\nabla_H h(\boldsymbol{x}).
```

The wave forcing is obtained directly from the bottom pressure of each complete phase-inclusive mode:

```math
\dot A_\alpha
=
\frac{1}{A E_\alpha}
\int_A p_{\alpha,d}^*(\boldsymbol{x},t)g_b(\boldsymbol{x},t)\,dA,
\qquad
\alpha\in\{+,-\}.
```

This construction adds only to the wave coefficients $A_+$ and $A_-$. It therefore produces no direct linear interior QGPV tendency and leaves the balanced forcing coefficient $F_0$ unchanged. It requires no finite-terrain pressure solve, modal terrain matrix, or artificial bottom-localized vertical envelope.

The prescribed barotropic current is an external energy reservoir. Wave energy is not conserved by itself; its rate of increase must equal the bottom pressure work supplied by the prescribed current.

WaveVortexModel 4.1.1 or newer is required. The forcing classes themselves have no diagnostics dependency. The eddy-tide analysis example additionally requires WaveVortexModelDiagnostics 1.0.6 or newer.

## Quick start

Construct a supported wave-bearing transform, prescribe the terrain and complex barotropic velocity amplitude, and register the spectral forcing:

```matlab
forcing = WVBottomWaveGenerationForcing(wvt, ...
    topographicHeight=h, ...
    barotropicVelocityAmplitude=[0.05; 0]);
wvt.removeAllForcing();
wvt.addForcing(forcing);
```

The velocity amplitude is the complex two-component vector $\widehat{\boldsymbol U}_{\mathrm{bt}}$ in meters per second. The default frequency is M2; `frequency`, `rampDuration`, `startTime`, and `name` are optional constructor arguments.

The constructor precomputes the bottom-pressure projection on the transform's native spectral layout. It supports `WVTransformBoussinesq`, both hydrostatic and nonhydrostatic `WVTransformConstantStratification` configurations, and `WVTransformHydrostatic`. These transforms provide their optimized modal endpoint factors through WaveVortexModel's `waveModeVerticalStructureAtIndex` API. Each subsequent forcing call evaluates the prescribed current, combines two response arrays per wave branch, and applies WaveVortexModel's interaction phases. There is no runtime pressure solve, FFT, spatial projection, or modal coupling matrix. Transforms with either value of `shouldAntialias` are supported.

`forcingWithResolutionOfTransform` spectrally transfers the terrain and rebuilds all modal responses for the new transform. The physical forcing configuration is also included when its parent transform or model is written to NetCDF; transform-derived response arrays are rebuilt after restoration.

Milestones 1--6, 8, and 9 of the development [roadmap](milestones.md) are implemented on the `mean-depth-wave-generator` branch. The optional comparison in Milestone 7 was intentionally skipped.

A separate [second-order roadmap](second-order-milestones.md) specifies the strict second Born hierarchy, its endpoint-convergence gate, and the tests required before it can be presented as a physical $O(h^2)$ scattering approximation. That hierarchy is planned work and is not part of the current forcing classes.

The scientific and computational status of the earlier repositories and branches is summarized in [PRIOR_APPROACHES.md](PRIOR_APPROACHES.md).

## Exact sinusoidal-ridge example

Run the known-solution benchmark and create its two diagnostic figures with:

```matlab
benchmark = SinusoidalRidgeWaveGenerationBenchmark;
```

The example uses a uniform M2 current over one sinusoidal terrain component. It compares three adaptive `ode78` integrations with the analytically integrated interaction coefficients, verifies that wave-energy growth equals the work done by bottom pressure, and shows the generated vertical velocity and displacement fields. Pass `shouldMakeFigures=false` for diagnostics without graphics or set `resolution` and `relativeTolerances` explicitly.

At the automated reference resolution, the coefficient errors for tolerances $10^{-6}$, $10^{-8}$, and $10^{-10}$ are approximately $6.4\times10^{-8}$, $5.6\times10^{-10}$, and $5.8\times10^{-12}$. The tight-run energy/work error is approximately $2.5\times10^{-12}$. The balanced tendency is exactly zero, while an independent physical-space calculation confirms the linear-QGPV source vanishes to the accuracy of the discrete derivative transforms.

## Goff abyssal-hill example

Generate deterministic periodic Goff topography with:

```matlab
[virtualDepth,h,diagnostics] = ...
    WVBottomWaveGenerationForcing.goffAbyssalHillTopography( ...
        wvt,minimumWavelength=40e3);
```

The generator uses a local random stream, returns upward-positive zero-mean terrain with the requested post-filter RMS, and reports realized slope and radial-spectrum diagnostics. It does not change MATLAB's global random state.

Run the variable-stratification broadband example with:

```matlab
result = GoffAbyssalHillWaveGenerationExample;
```

The example drives a seed-2023, $100$ m RMS Goff field with a uniform M2 current and shows the terrain spectrum, wave-energy distribution, and bottom-work budget. Use `BottomWaveGenerationPerformanceBenchmark` to report construction, storage, and forcing-application costs at several resolutions.

## Autonomous first-order wave scattering

`WVBottomWaveScatteringForcing` is a separate forcing for scattering an existing wave field. It evaluates the wave-only mean-depth bottom velocity

```math
g_b
=
\boldsymbol u_{H,d}\boldsymbol{\cdot}\nabla_Hh
-
h\,\partial_z w_d
```

and applies the same bottom-pressure projection to the two wave branches:

```matlab
scattering = WVBottomWaveScatteringForcing( ...
    wvt,topographicHeight=h);
wvt.removeAllForcing();
wvt.addForcing(scattering);
```

The ordinary forcing call reconstructs only the three required bottom fields, performs horizontal pseudospectral products, and projects the resulting bottom velocity. It supports the same Boussinesq, constant-stratification, and hydrostatic transforms as the generation forcing through WaveVortexModel's endpoint-factor API. It does not evaluate full three-dimensional fields, solve for pressure, assemble a modal terrain matrix, or modify the balanced tendency.

Bottom displacement is a postprocessed diagnostic rather than a model state. Save ordinary WaveVortexModel output and integrate the saved bottom velocity with:

```matlab
diagnostics = ...
    WVBottomWaveScatteringForcing.bottomDisplacementFromFile( ...
        "scattering-output.nc");
```

The method reads `t`, `Ap`, and `Am` from the standard `wave-vortex` group, reconstructs $g_b$ at the selected accepted output times, and applies cumulative trapezoidal quadrature to $\partial_t\eta_d=g_b$. Use `iTime` to select output records, `forcingName` when the file contains more than one scattering forcing, and `initialBottomDisplacement` for a nonzero initial condition.

Two examples exercise the autonomous formulation:

```matlab
uniform = UniformDepthWaveScatteringBenchmark;
ridge = SinusoidalRidgeWaveScatteringExample;
gaussian = GaussianRidgeWaveScatteringExample;
```

The uniform-depth benchmark verifies the first-order frequency correction against the exact depth-$D-h_0$ dispersion relation. The sinusoidal-ridge example writes standard output, displays the scattered sidebands and energy exchange, and reconstructs the bottom displacement from that file.

The Gaussian-ridge example follows a localized rightward mode-one M2 wave packet as it crosses a gentle subcritical ridge. Its default terrain has $h_0/D=0.05$ and $\max|h_x|/\mu=0.05$, keeping the calculation in the controlled first-order regime. It uses adaptive `ode78`, returns the sampled `Ap`, `Am`, and `A0` trajectory and diagnostics in memory, and creates no output file unless one is requested:

```matlab
result = GaussianRidgeWaveScatteringExample( ...
    outputPath="gaussian-ridge.nc");
```

The output is standard restartable WaveVortexModel NetCDF. A separate renderer reads only that saved file and creates an $x$--$z$ movie and PNG poster:

```matlab
movie = GaussianRidgeWaveScatteringMovie( ...
    "gaussian-ridge.nc",field="u");
```

The movie masks the reconstructed field beneath the physical bottom, uses wet points only for its fixed color scale, and shows the first-order depth-integrated wave-energy profile alongside the evolving rightward, reflected, and higher-mode energy. Its default section uses a labelled 20-times vertical exaggeration, configurable with `verticalExaggeration`. The renderer also supports `field="w"` and `field="eta"`, a subset of saved records through `iTime`, and explicit `videoPath` and `posterPath` values. Existing files are never replaced unless `shouldOverwriteExisting=true`.

## Thirty-day eddy-tide validation

`EddyTideTopographicForcingSimulation` replaces the initialized fixed-amplitude wave beam in the JPO2026 minimal eddy-tide calculation with waves generated by a prescribed barotropic tide over Goff topography. The validation calculation is generation-only by default. Autonomous first-order scattering remains available with `shouldUseScattering=true`, but is kept outside the production pair until its $O(h^2)$ residual is demonstrably small. The default transform retains the constant $N^2=2\times10^{-5}\ \mathrm{s^{-2}}$ JPO configuration, the four-mode-one-wavelength horizontal domain, the shallow $10$ cm/s eddy, nonlinear advection, adaptive damping, and antialiasing.

The initial validation runs use `Nxy=128`, a 30-day final time, and hourly output. The terrain has $100$ m RMS height, corner wavenumber $10^{-4}\ \mathrm{m^{-1}}$, a 20 km short-wavelength cutoff, and seed 2023. The 20 km cutoff is the nearest round value above the 17.84 km antialiased resolution limit. A zonal 5 cm/s M2 current ramps up over one tidal period. Both simulations start without waves; the matched control differs only by omitting the initial eddy. Hourly records resolve the twice-tidal terms needed for forcing-resolved energy and APV-enstrophy budgets.

Set `horizontalDomainSize` to use an explicit square domain instead of the default four-mode-one-wavelength domain. The simulation passes that size to `WVStratification.verticalResolutionForHorizontalResolution`, so the vertical grid remains consistent with the requested horizontal resolution and stratification.

Run the pair and build their energy comparison with:

```matlab
[~,~,eddyFile] = EddyTideTopographicForcingSimulation;
[~,~,noEddyFile] = EddyTideTopographicForcingSimulation(includeEddy=false);
[figureHandle,energy,diagnosticsFiles] = AnalyzeEddyTideEnergy(eddyFile,noEddyFile);
[basicFigures,summary] = AnalyzeEddyTideBasicFigures(eddyFile,noEddyFile);
[budgetFigures,budget] = AnalyzeEddyTideBudgets(eddyFile,noEddyFile);
```

The lean energy analysis creates or updates standard WaveVortexModelDiagnostics files and plots wave and geostrophic energy. `AnalyzeEddyTideBasicFigures` adds maximum wave-speed histories, final wave spectra, and final surface vertical-vorticity maps while saving the plotted diagnostics in a MAT summary. `AnalyzeEddyTideBudgets` additionally diagnoses geostrophic, mean-density-anomaly, wave, total-energy, quadratic potential-enstrophy, and exact APV-enstrophy reservoirs. It integrates forcing and nonlinear-triad rates, reports closure residuals, forms eddy-minus-control differences, and writes reservoir, forcing-budget, and triad-budget figures.

Test autonomous scattering against the generated control without altering its output:

```matlab
[benchmark,figureHandle] = GoffScatteringPerturbationBenchmark(noEddyFile);
```

The benchmark holds the terrain realization fixed while scaling its RMS height from 100 down to 3.125 m by factors of two. It verifies $O(h)$ first-order work and $O(h^2)$ autonomous residuals, runs five-M2-period scattering-only integrations, and reports the largest height satisfying the residual and energy-drift gates.

The model files are ordinary restartable WaveVortexModel output. Extend a completed run to 60 days with:

```matlab
model = WVModel.modelFromFile(char(eddyFile));
model.integrateToTime(60*86400);
model.closeNetCDFFile();
```

Rerunning `AnalyzeEddyTideEnergy` after an extension appends diagnostics for the new saved times. Existing simulation and figure files are never replaced unless their corresponding `shouldOverwriteExisting` option is true.

This is initially an infrastructure-validation version of the experiment, not a literal reproduction of the $500$ m-grid MITgCM calculation in [Shakespeare (2023)](https://doi.org/10.1175/JPO-D-23-0127.1). The 20 km terrain cutoff is resolved by the default antialiased WVM grid. “No eddy” means no initial eddy; nonlinear wave rectification may subsequently generate balanced flow. Longer integrations and manuscript analysis can reuse the same restart and diagnostics workflow after the paired 30-day runs are verified.

### 500 km Shakespeare-comparison pair

The first Shakespeare-comparison attempt consists of two matched generation-only simulations, one initialized with the shallow eddy and one initialized without it. The intent is to move toward the key eddy-tide calculation in [Shakespeare (2023)](https://doi.org/10.1175/JPO-D-23-0127.1), not to claim a literal reproduction. Relative to the default infrastructure test, this pair uses the paper's $500$ km square domain and increases the horizontal spectral resolution to $256^2$. It remains coarser than the paper's $500$ m MITgCM grid, but a spectral discretization makes a direct grid-spacing comparison conservative. The automatically selected vertical grid has 45 points and retains 29 wave modes.

The pair uses constant $N^2=2\times10^{-5}\ \mathrm{s^{-2}}$, the nonhydrostatic antialiased transform, nonlinear advection, adaptive damping, the $10$ cm/s shallow eddy, and the same generation-only topographic forcing in both cases. The terrain is a seed-2023 Goff realization with $100$ m RMS height, corner wavenumber $10^{-4}\ \mathrm{m^{-1}}$, and a 6 km minimum wavelength. The 6 km cutoff is the nearest round cutoff above the 5.882 km antialiased limit. The zonal M2-current amplitude is $0.05/\sqrt{10}\ \mathrm{m\,s^{-1}}$; this is the reduced-input experiment, with one tenth of the nominal wave-energy input because generation scales quadratically with current amplitude. The tide ramps over one M2 period. Full state is saved every 6 hours.

To repeat the actual calculation, run the following two calls in separate MATLAB processes so that they integrate concurrently:

```matlab
[~,~,eddyFile] = EddyTideTopographicForcingSimulation( ...
    includeEddy=true, ...
    horizontalDomainSize=500e3, ...
    Nxy=256, ...
    maxT=30*86400, ...
    outputInterval=6*3600, ...
    barotropicVelocityAmplitude=[0.05/sqrt(10); 0], ...
    rmsHeight=100, ...
    cornerWavenumber=1e-4, ...
    minimumWavelength=6e3, ...
    randomSeed=2023, ...
    shouldUseScattering=false, ...
    shouldAntialias=true, ...
    outputFilename="eddy-tide-topographic-eddy-generation-only-Lxy500km-Nxy256-Ein10-hrms100m-lmin6km-seed2023.nc");
```

```matlab
[~,~,noEddyFile] = EddyTideTopographicForcingSimulation( ...
    includeEddy=false, ...
    horizontalDomainSize=500e3, ...
    Nxy=256, ...
    maxT=30*86400, ...
    outputInterval=6*3600, ...
    barotropicVelocityAmplitude=[0.05/sqrt(10); 0], ...
    rmsHeight=100, ...
    cornerWavenumber=1e-4, ...
    minimumWavelength=6e3, ...
    randomSeed=2023, ...
    shouldUseScattering=false, ...
    shouldAntialias=true, ...
    outputFilename="eddy-tide-topographic-no-eddy-generation-only-Lxy500km-Nxy256-Ein10-hrms100m-lmin6km-seed2023.nc");
```

The reference files were then restarted independently and extended by 20 days:

```matlab
model = WVModel.modelFromFile(char(eddyFile));
model.integrateToTime(50*86400);
model.closeNetCDFFile();
```

Apply the same restart commands to `noEddyFile`. Each completed file contains 201 records from day 0 through day 50. Reproduce the four day-50 figures and MAT summary with:

```matlab
[figures,summary] = AnalyzeEddyTideBasicFigures( ...
    eddyFile,noEddyFile, ...
    exportPrefix="eddy-tide-Lxy500km-Nxy256-Ein10-lmin6km-day50");
```

As cross-machine reference values, the day-50 eddy run has wave and geostrophic energies $0.7760$ and $1.4376\ \mathrm{m^3\,s^{-2}}$; the no-initial-eddy run has $0.8171$ and $0.1446\ \mathrm{m^3\,s^{-2}}$. Their largest 6-hour-sampled three-dimensional wave speeds are 13.31 and 12.76 cm/s, respectively. Both final wave spectra peak in the 62.5 km, mode-3 bin and remain surprisingly similar. The experiment therefore has not yet reproduced the strong eddy-induced scattering seen by Shakespeare. This is consistent with the current setup being generation-only: autonomous `WVBottomWaveScatteringForcing` is deliberately disabled, so spectral differences can arise only through the resolved nonlinear dynamics and advection.

## Scientific scope

Both forcings are accurate through first order in terrain height. The prescribed generator supports broadband terrain with arbitrary stationary stratification, transform-resolution rebuilding, and restart persistence. The autonomous forcing adds wave--wave scattering with the same production behavior. Its strict perturbation hierarchy conserves the first-order physical energy; when the first-order operator is iterated autonomously, the unresolved energy tendency is $O(h^2)$.

The planned strict second Born model will retain the zeroth-, first-, and second-order coefficients separately and stop before uncontrolled higher-order feedback. See the [second-order milestones](second-order-milestones.md) for its scientific gates and implementation sequence.

The following remain outside the initial proof of concept:

- exact finite-amplitude terrain dynamics;
- nonlinear terrain-aware advection;
- independent bottom-buoyancy dynamics;
- backreaction on a dynamically evolving barotropic tide;
- exact finite-amplitude energy conservation under autonomous first-order scattering.

WaveVortexModel and the analysis-only WaveVortexModelDiagnostics package remain external dependencies and are not modified or vendored by this repository. Run the automated suite with:

```matlab
results = runTests;
```

The runner resolves WaveVortexModel from an explicit `waveVortexModelRoot` option, `WAVE_VORTEX_MODEL_ROOT`, or the sibling `wave-vortex-model` repository. It resolves WaveVortexModelDiagnostics in the same way from `waveVortexModelDiagnosticsRoot`, `WAVE_VORTEX_MODEL_DIAGNOSTICS_ROOT`, or the sibling `wave-vortex-model-diagnostics` repository.
