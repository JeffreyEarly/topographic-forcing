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

## Quick start

Construct a `WVTransformBoussinesq`, prescribe the terrain and complex barotropic velocity amplitude, and register the spectral forcing:

```matlab
forcing = WVBottomWaveGenerationForcing(wvt, ...
    topographicHeight=h, ...
    barotropicVelocityAmplitude=[0.05; 0]);
wvt.removeAllForcing();
wvt.addForcing(forcing);
```

The velocity amplitude is the complex two-component vector $\widehat{\boldsymbol U}_{\mathrm{bt}}$ in meters per second. The default frequency is M2; `frequency`, `rampDuration`, `startTime`, and `name` are optional constructor arguments.

The constructor precomputes the bottom-pressure projection on the transform's native spectral layout. It supports any stationary stratification represented by `WVTransformBoussinesq`. Each subsequent forcing call evaluates the prescribed current, combines two response arrays per wave branch, and applies WaveVortexModel's interaction phases. There is no runtime pressure solve, FFT, spatial projection, or modal coupling matrix. Transforms with either value of `shouldAntialias` are supported.

`forcingWithResolutionOfTransform` spectrally transfers the terrain and rebuilds all modal responses for the new transform. The physical forcing configuration is also included when its parent transform or model is written to NetCDF; transform-derived response arrays are rebuilt after restoration.

Milestones 1--6, 8, and 9 of the development [roadmap](milestones.md) are implemented on the `mean-depth-wave-generator` branch. The optional comparison in Milestone 7 was intentionally skipped.

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

The ordinary forcing call reconstructs only the three required bottom fields, performs horizontal pseudospectral products, and projects the resulting bottom velocity. It does not evaluate full three-dimensional fields, solve for pressure, assemble a modal terrain matrix, or modify the balanced tendency.

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
```

The uniform-depth benchmark verifies the first-order frequency correction against the exact depth-$D-h_0$ dispersion relation. The sinusoidal-ridge example writes standard output, displays the scattered sidebands and energy exchange, and reconstructs the bottom displacement from that file.

## Scientific scope

Both forcings are accurate through first order in terrain height. The prescribed generator supports broadband terrain with arbitrary stationary stratification, transform-resolution rebuilding, and restart persistence. The autonomous forcing adds wave--wave scattering with the same production behavior. Its strict perturbation hierarchy conserves the first-order physical energy; when the first-order operator is iterated autonomously, the unresolved energy tendency is $O(h^2)$.

The following remain outside the initial proof of concept:

- exact finite-amplitude terrain dynamics;
- nonlinear terrain-aware advection;
- independent bottom-buoyancy dynamics;
- backreaction on a dynamically evolving barotropic tide;
- exact finite-amplitude energy conservation under autonomous first-order scattering.

WaveVortexModel remains an external dependency and is not modified or vendored by this repository. Run the automated suite with:

```matlab
results = runTests;
```

The runner resolves WaveVortexModel from an explicit `waveVortexModelRoot` option, `WAVE_VORTEX_MODEL_ROOT`, or the sibling `wave-vortex-model` repository.
