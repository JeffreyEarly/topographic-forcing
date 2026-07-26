# Mean-depth bottom wave-generation milestones

## Objective

Develop `WVBottomWaveGenerationForcing`, a fast spectral forcing for `WVTransformBoussinesq` that generates internal waves from the first-order mean-depth bottom condition without directly forcing linear interior QGPV. The proof of concept will prescribe a horizontally uniform barotropic tide and project its bottom velocity onto the ordinary rigid-lid wave modes using their bottom pressure values.

For stationary topography $h(\boldsymbol{x})$ and a prescribed barotropic current,

```math
\boldsymbol{U}_{\mathrm{bt}}(t)
=
R(t)\mathrm{Re}
\left\{
\widehat{\boldsymbol{U}}_{\mathrm{bt}}
e^{-i\omega(t-t_0)}
\right\},
\qquad
g_b
=
\boldsymbol{U}_{\mathrm{bt}}\boldsymbol{\cdot}\nabla_Hh.
```

The phase-inclusive modal tendency is

```math
\dot A_\alpha
=
\frac{1}{A E_\alpha}
\int_Ap_{\alpha,d}^*(\boldsymbol{x},t)g_b(\boldsymbol{x},t)\,dA,
\qquad
\alpha\in\{+,-\}.
```

Only $F_+$ and $F_-$ are modified. The forcing leaves the incoming $F_0$ unchanged, requires no finite-terrain pressure solve or modal terrain matrix, and is first order in terrain height.

## Fixed conventions and boundaries

- `topographicHeight` is a real, finite, stationary, upward-positive field sampled on the transform's periodic $N_x\times N_y$ horizontal grid.
- `barotropicVelocityAmplitude` is a finite complex two-component vector $\widehat{\boldsymbol U}_{\mathrm{bt}}$ with units of velocity.
- `frequency` is a finite positive angular frequency, `rampDuration` is finite and nonnegative, and `startTime` is finite.
- The transform target is `WVTransformBoussinesq` with arbitrary stationary $N^2(z)$ supported by the transform.
- The forcing uses the transform's native retained wavenumbers and is agnostic to `shouldAntialias`.
- The initial barotropic current is prescribed and horizontally uniform. It is an external energy reservoir rather than a prognostic part of the model state.
- Complete wave modes include their spatial and temporal phases. WaveVortexModel's stored interaction coefficients therefore receive the phase conversion required by the conjugated output-mode pressure.
- The initial generator applies the additive term $g_b=\boldsymbol U_{\mathrm{bt}}\boldsymbol{\cdot}\nabla_Hh$. Autonomous wave scattering is deferred to Milestone 9.
- WaveVortexModel remains an external dependency. This research repository will not initially be released as an MPM package.

## Milestone 1: Repository and scientific contract

- [x] Complete

### Purpose

Replace the mapped strong-form experiment with a narrowly defined public API for prescribed bottom wave generation.

### Dependencies

None.

### Deliverables

- Replace `WVExactTopographicForcing` and its tests with `WVBottomWaveGenerationForcing`.
- Implement the public constructor:

  ```matlab
  forcing = WVBottomWaveGenerationForcing(wvt, ...
      topographicHeight=h, ...
      barotropicVelocityAmplitude=Uhat, ...
      frequency=omega, ...
      rampDuration=rampDuration, ...
      startTime=wvt.t, ...
      name="bottom wave generation");
  ```

- Validate the supported transform, terrain dimensions and values, the two-component velocity amplitude, frequency, ramp duration, start time, and name with structured class-specific errors.
- Store authoritative scientific inputs as read-only state.
- Retain a portable `runTests` entry point that resolves WaveVortexModel from an explicit option, then `WAVE_VORTEX_MODEL_ROOT`, then the sibling repository, and restores the original MATLAB path.
- Document that the forcing is first order in terrain height, wave-only, and energetically open because the barotropic current is prescribed.

### Automated acceptance

- Constructor tests accept flat, uniform-offset, sinusoidal, and broadband real terrain arrays.
- Constructor tests reject unsupported transforms, bad terrain shapes, complex or nonfinite terrain, invalid velocity amplitudes, nonpositive frequency, negative ramp duration, and nonfinite start time.
- The test runner leaves the user's permanent MATLAB path unchanged.
- The old public class name and mapped-coordinate scientific claims no longer appear in active source or tests.

## Milestone 2: Boundary-projection oracle

- [x] Complete

### Purpose

Establish the exact WaveVortexModel normalization and phase convention for the bottom Green-identity projection before optimizing it.

### Dependencies

Milestone 1.

### Deliverables

- Build a low-resolution reference routine that reconstructs the complete phase-inclusive bottom pressure $p_{\alpha,d}(\boldsymbol{x},t)$ of every active $+$ and $-$ mode.
- Evaluate $g_b$ in physical space and compute each modal tendency by direct horizontal quadrature.
- Evaluate the same overlap from the Fourier coefficient of $g_b$ selected by the output mode.
- Treat the transform's half-complex Fourier layout, inactive entries, horizontal zero mode, branch normalization, modal energy $E_\alpha$, and conjugacy explicitly.
- Keep the oracle independent of the production forcing callback.

### Automated acceptance

- Direct quadrature and Fourier selection agree within $10^{-12}$ relative error for deterministic single-mode and random band-limited terrain.
- Results agree at $t=t_0$ and at nonzero model times containing nondegenerate wave phases.
- Reconstructed real boundary fields and coefficient conjugacy close to $10^{-12}$.
- No balanced coefficient is evaluated or retained by the oracle.

## Milestone 3: Fast prescribed-generation kernel

- [x] Complete

### Purpose

Reduce prescribed barotropic generation to precomputed spectral response arrays and inexpensive time-dependent scalar factors.

### Dependencies

Milestone 2.

### Deliverables

- Precompute the spectral terrain gradients and the projections of unit currents in the $x$ and $y$ directions onto every active wave mode.
- Verify the equivalent complex-amplitude identity

```math
\widehat g_b(\boldsymbol K)
=
i\boldsymbol K\boldsymbol{\cdot}
\widehat{\boldsymbol U}_{\mathrm{bt}}\,
\widehat h(\boldsymbol K).
```

- At runtime evaluate the real ramped current, combine the two precomputed responses, apply the output-mode interaction phases componentwise, and return $F_+$ and $F_-$.
- Preserve exact zeros at inactive and excluded coefficient locations.
- Leave the incoming $F_0$ unchanged rather than assigning a new zero array.

### Automated acceptance

- Production tendencies agree with the direct-quadrature oracle within $10^{-12}$ at multiple model times.
- The flat-terrain and zero-current limits vanish exactly apart from signed floating-point zero.
- Tendencies scale linearly with terrain amplitude and barotropic velocity to $10^{-12}$.
- Repeated evaluation is deterministic and does not modify `wvt.t`, `Ap`, `Am`, or `A0`.

## Milestone 4: WaveVortexModel forcing integration

- [x] Complete

### Purpose

Connect the validated kernel to the standard WaveVortexModel spectral forcing and adaptive integration path.

### Dependencies

Milestones 1 through 3.

### Deliverables

- Implement `WVBottomWaveGenerationForcing` as `WVForcingType("Spectral")`.
- In `addSpectralForcing`, add only this forcing's contributions to the incoming `Fp` and `Fm` and return the incoming `F0` unchanged.
- Preserve contributions from every previously registered forcing.
- Run a model from rest after removing nonlinear advection and registering only the bottom generator.
- Use WaveVortexModel's default adaptive `ode78` integrator.
- During this milestone, return a structured unsupported-operation error from resolution conversion; Milestone 8 replaces this temporary gate with spectral rebuilding.

### Automated acceptance

- Direct callback results agree with the independently evaluated pack, phase, and projection calculation at zero and nonzero model times.
- Existing `Fp`, `Fm`, and `F0` arrays are incremented exactly once or left unchanged as appropriate.
- A short adaptive run produces nonzero wave coefficients and identically zero direct balanced forcing.
- An ordinary forcing call performs no pressure solve, wave--vortex transform, terrain matrix construction, or runtime horizontal transform.

## Milestone 5: No-PV and energy-work scientific gates

- [x] Complete

### Purpose

Demonstrate that the implemented source generates waves without directly generating linear interior QGPV and that its energetics reproduce the bottom pressure work.

### Dependencies

Milestone 4.

### Deliverables

- Verify that the isolated forcing has $F_0=0$ in the transform's independent coefficient layout.
- Reconstruct the physical wave tendencies and evaluate their linear QGPV source with the transform's native derivative operators.
- Compare the modal wave-energy tendency with the bottom work in WaveVortexModel's energy-per-density normalization,

```math
P_b
=
\frac1A
\int_A\frac{p_{w,d}}{\rho_0}g_b\,dA,
```

  where $p_{w,d}$ is the bottom pressure reconstructed from the evolving wave state.
- Test deterministic random wave states, states reached during forced evolution, and both wave branches.
- Verify amplitude scaling and the flat-terrain and zero-current limits independently.

### Automated acceptance

- The isolated balanced and modal QGPV tendencies are exactly zero. The independently differentiated physical linear-QGPV tendency is below $10^{-10}$ relative to its vorticity-plus-stretching cancellation scale.
- Modal source power and bottom pressure work agree within $10^{-12}$ relative error for one RHS evaluation.
- Terrain and current rescaling produce the corresponding linear source rescaling within $10^{-12}$.
- Any failure of the PV or work identities blocks Milestone 6 and is resolved without empirical correction factors.

## Milestone 6: Known-solution linear benchmark

- [x] Complete

### Purpose

Validate the complete forcing and adaptive model evolution against an analytically integrated modal response.

### Dependencies

Milestone 5.

### Deliverables

- Use

```math
[L_x,L_y,D]=[20,20,2]\ {\rm km},
\qquad
N^2=2\times10^{-5}\ {\rm s}^{-2},
\qquad
h(x)=50\ {\rm m}\cos(2\pi x/L_x),
```

  at latitude $45^\circ$, with a $5\ {\rm cm\,s^{-1}}$ $x$-directed M2 barotropic velocity and no startup ramp.
- Use resolution `[8 4 5]` for the automated oracle and `[32 4 17]` for the default visual example.
- Initialize from rest, remove nonlinear advection, and register only `WVBottomWaveGenerationForcing`.
- Integrate the phase-inclusive forced coefficient equations analytically for every excited mode and use them as the reference solution.
- Run `WVModel` with adaptive relative and absolute tolerances $10^{-6}$, $10^{-8}$, and $10^{-10}$.
- Report coefficient error, wave-energy growth, integrated bottom work, QGPV, branch energy, and vertical-mode distribution without committing generated output.
- Provide a user-facing example that returns the benchmark diagnostics and, by default, shows:
  - exact-solution convergence, energy-work closure, branch energy, and vertical-mode distribution;
  - the ridge, local bottom work density, vertical velocity, and displacement at the time of maximum exact wave energy.

### Automated acceptance

- Only the terrain's horizontal Fourier pair is directly forced.
- Coefficient errors decrease with adaptive tolerance and are below $10^{-8}$ at the tightest tolerance.
- Integrated bottom work and wave-energy change agree within $10^{-8}$ at the tightest tolerance and converge together as tolerances tighten.
- The balanced coefficient and linear QGPV norms remain below $10^{-10}$ of the wave response.
- Both example figures render with finite plotted data and no generated files are required.

## Milestone 7: Comparison with Pseudo-topography

- [x] Skipped by decision

### Status

This optional comparison was intentionally skipped after the prescribed generator passed its direct scientific gates. It is not a dependency of the production implementation.

No comparison code or dependency on `Pseudo-topography` is added.

## Milestone 8: Generalization and production behavior

- [x] Complete

### Purpose

Extend the validated generator to practical terrain, stratification, resolution, persistence, and performance requirements without changing its scientific definition.

### Dependencies

Milestone 6.

### Deliverables

- Support broadband periodic terrain and arbitrary stationary $N^2(z)$ available to `WVTransformBoussinesq`.
- Add a deterministic periodic Goff abyssal-hill generator by adapting the validated implementation from `Pseudo-topography` without introducing a repository dependency.
- Implement `forcingWithResolutionOfTransform` by resampling the authoritative terrain and rebuilding all modal projection factors for the new transform.
- Verify explicit rebuilding for WaveVortexModel transforms created at alternate antialias resolutions.
- Add restart persistence for the terrain, barotropic amplitude, frequency, ramp duration, start time, and forcing name; rebuild transform-derived response arrays after restoration.
- Add a variable-stratification Goff-terrain model example showing the terrain spectrum, generated wave-energy distribution, and bottom-work budget.
- Profile construction, storage, and runtime application at three resolutions.
- Retain the pressure-free and matrix-free runtime path.

### Automated acceptance

- The optimized response agrees with the direct-quadrature oracle within $10^{-10}$ at reference resolution for constant and variable stratification.
- Resolution conversion preserves forcing phase, conjugacy, and the zero-PV identity.
- Restart round trips reproduce forcing tendencies within $10^{-12}$.
- Adaptive continuation after restart agrees with an uninterrupted control within $10^{-10}$ in the stored wave coefficients.
- The Goff generator reproduces its requested statistics and cutoff, is deterministic without changing the global random state, and its example renders finite diagnostics without generated files.
- Runtime application contains only scalar current evaluation, componentwise phase factors, and additions of precomputed wave arrays.

## Milestone 9: Optional autonomous wave scattering

- [ ] Complete

### Purpose

Extend the validated additive generator to first-order scattering of an evolving zero-interior-APV wave field.

### Dependencies

Milestone 8.

### Deliverables

- Reconstruct the instantaneous bottom fields of the wave state and evaluate

```math
g_b
=
\boldsymbol u_{H,d}\boldsymbol{\cdot}\nabla_Hh
-
h\,\partial_zw_d.
```

- Apply the same pressure-weighted projection only to $A_+$ and $A_-$.
- Advance the separate bottom displacement diagnostic with $\partial_t\eta_d=g_b$.
- Verify the sinusoidal-terrain sideband rule $\boldsymbol K\mapsto\boldsymbol K\pm\boldsymbol q$ and linear leading-order dependence on terrain amplitude.
- Evaluate the first-order physical energy and verify that the autonomously iterated model's residual scales as $O(h^2)$.

### Automated acceptance

- Direct linear interior QGPV remains zero to $10^{-12}$.
- Leading sidebands scale linearly with terrain amplitude, while the first-order energy residual converges quadratically.
- The bottom displacement diagnostic agrees with the time integral of $g_b$.
- Dynamic barotropic backreaction, independent bottom buoyancy, nonlinear terrain dynamics, and finite-amplitude exactness remain explicitly unsupported.

## Planning and implementation cadence

- Plan and implement Milestones 1--4 together as the first executable vertical slice.
- Plan Milestones 5--6 together after the forcing callback passes its algebraic tests. Milestone 5 is a blocking scientific gate.
- Proceed directly from Milestone 6 to Milestone 8; Milestone 7 is intentionally skipped.
- Treat Milestone 9 as a separate extension; it is not required to establish the prescribed bottom wave generator.

## Definition of done

The initial proof of concept is complete when Milestones 1--6 pass: the spectral forcing agrees with the Green-identity oracle, directly forces no balanced coefficient or linear QGPV, reproduces bottom pressure work, and converges to the analytic sinusoidal-terrain response under adaptive integration.

The prescribed generator is production-ready for research use when Milestone 8 passes. Promotion into WaveVortexModel, an MPM release, dynamic barotropic backreaction, independent bottom buoyancy, nonlinear terrain dynamics, and exact finite-amplitude topography remain outside this roadmap.
