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

Construct a constant-stratification `WVTransformBoussinesq`, prescribe the terrain and complex barotropic velocity amplitude, and register the spectral forcing:

```matlab
forcing = WVBottomWaveGenerationForcing(wvt, ...
    topographicHeight=h, ...
    barotropicVelocityAmplitude=[0.05; 0]);
wvt.removeAllForcing();
wvt.addForcing(forcing);
```

The velocity amplitude is the complex two-component vector $\widehat{\boldsymbol U}_{\mathrm{bt}}$ in meters per second. The default frequency is M2; `frequency`, `rampDuration`, `startTime`, and `name` are optional constructor arguments.

The constructor precomputes the bottom-pressure projection on the transform's native spectral layout. Each subsequent forcing call evaluates the prescribed current, combines two response arrays per wave branch, and applies WaveVortexModel's interaction phases. There is no runtime pressure solve, FFT, spatial projection, or modal coupling matrix. Transforms with either value of `shouldAntialias` are supported.

Milestones 1--6 of the development [roadmap](milestones.md) are implemented on the `mean-depth-wave-generator` branch. Explicit conversion to a different transform resolution remains deferred.

The scientific and computational status of the earlier repositories and branches is summarized in [PRIOR_APPROACHES.md](PRIOR_APPROACHES.md).

## Exact sinusoidal-ridge example

Run the known-solution benchmark and create its two diagnostic figures with:

```matlab
benchmark = SinusoidalRidgeWaveGenerationBenchmark;
```

The example uses a uniform M2 current over one sinusoidal terrain component. It compares three adaptive `ode78` integrations with the analytically integrated interaction coefficients, verifies that wave-energy growth equals the work done by bottom pressure, and shows the generated vertical velocity and displacement fields. Pass `shouldMakeFigures=false` for diagnostics without graphics or set `resolution` and `relativeTolerances` explicitly.

At the automated reference resolution, the coefficient errors for tolerances $10^{-6}$, $10^{-8}$, and $10^{-10}$ are approximately $6.4\times10^{-8}$, $5.6\times10^{-10}$, and $5.8\times10^{-12}$. The tight-run energy/work error is approximately $2.5\times10^{-12}$. The balanced tendency is exactly zero, while an independent physical-space calculation confirms the linear-QGPV source vanishes to the accuracy of the discrete derivative transforms.

## Scientific scope

The initial generator is accurate through first order in terrain height. It is intended to establish bottom generation without direct linear PV forcing before adding autonomous wave scattering, arbitrary stratification, persistence, or dynamic barotropic backreaction.

The following remain outside the initial proof of concept:

- exact finite-amplitude terrain dynamics;
- nonlinear terrain-aware advection;
- independent bottom-buoyancy dynamics;
- backreaction on a dynamically evolving barotropic tide;
- exact conservation of wave energy when the barotropic energy reservoir is omitted.

WaveVortexModel remains an external dependency and is not modified or vendored by this repository. Run the automated suite with:

```matlab
results = runTests;
```

The runner resolves WaveVortexModel from an explicit `waveVortexModelRoot` option, `WAVE_VORTEX_MODEL_ROOT`, or the sibling `wave-vortex-model` repository.
