# Mean-depth bottom wave generation

`topographic-forcing` is being redirected toward a fast, first-order bottom wave generator for WaveVortexModel. The new formulation retains the ordinary rigid-lid wave--vortex basis and represents weak topography through the mean-depth bottom condition rather than through mapped-coordinate volume terms.

For a prescribed, horizontally uniform barotropic current, the bottom velocity is

\[
g_b(\boldsymbol{x},t)
=
\boldsymbol{U}_{\mathrm{bt}}(t)\boldsymbol{\cdot}\nabla_H h(\boldsymbol{x}).
\]

The wave forcing is obtained directly from the bottom pressure of each complete phase-inclusive mode:

\[
\dot A_\alpha
=
\frac{1}{A E_\alpha}
\int_A p_{\alpha,d}^*(\boldsymbol{x},t)g_b(\boldsymbol{x},t)\,dA,
\qquad
\alpha\in\{+,-\}.
\]

This construction adds only to the wave coefficients \(A_+\) and \(A_-\). It therefore produces no direct linear interior QGPV tendency and leaves the balanced forcing coefficient \(F_0\) unchanged. It requires no finite-terrain pressure solve, modal terrain matrix, or artificial bottom-localized vertical envelope.

The prescribed barotropic current is an external energy reservoir. Wave energy is not conserved by itself; its rate of increase must equal the bottom pressure work supplied by the prescribed current.

## Development status

Development of the new formulation will occur on the `mean-depth-wave-generator` branch according to [milestones.md](milestones.md). The existing `WVExactTopographicForcing` source on this branch is the preserved checkpoint of the earlier mapped strong-form attempt. It is not the recommended formulation and will be replaced in Milestone 1.

The initial public API is planned to be

```matlab
forcing = WVBottomWaveGenerationForcing(wvt, ...
    topographicHeight=h, ...
    barotropicVelocityAmplitude=Uhat, ...
    frequency=omega, ...
    rampDuration=rampDuration, ...
    startTime=wvt.t, ...
    name="bottom wave generation");
```

Here `Uhat` is a finite complex two-component vector defining a horizontally uniform harmonic barotropic velocity. The first implementation will target `WVTransformBoussinesq`, constant stratification, stationary periodic terrain, and linear wave generation from a prescribed current.

The scientific and computational status of the earlier repositories and branches is summarized in [PRIOR_APPROACHES.md](PRIOR_APPROACHES.md).

## Scientific scope

The initial generator is accurate through first order in terrain height. It is intended to establish bottom generation without direct linear PV forcing before adding autonomous wave scattering, arbitrary stratification, persistence, or dynamic barotropic backreaction.

The following remain outside the initial proof of concept:

- exact finite-amplitude terrain dynamics;
- nonlinear terrain-aware advection;
- independent bottom-buoyancy dynamics;
- backreaction on a dynamically evolving barotropic tide;
- exact conservation of wave energy when the barotropic energy reservoir is omitted.

WaveVortexModel remains an external dependency and is not modified or vendored by this repository.
