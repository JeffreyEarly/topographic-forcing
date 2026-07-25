# Exact topographic forcing milestones

## Objective

Develop `WVExactTopographicForcing`, a research add-on that evaluates the exact flow-linear bottom-topography terms pseudospectrally and projects them onto the existing rigid-lid wave--vortex modes of `WVTransformBoussinesq`. The proof of concept will use the mapped coordinate \(\xi\), retain the finite-terrain quadratic energy and potential-enstrophy laws, and avoid modal terrain matrices, modifications to WaveVortexModel, and preservation of the flat-reference invariants.

The initial roadmap is deliberately linear. Complete terrain-aware advection \(\mathcal N\) and nonlinear buoyancy \(\mathcal B\) will be considered only after the exact linear formulation has passed its analytical, boundary, conservation, resolution, and performance tests.

## Fixed conventions and boundaries

- `topographicHeight` is a real, finite, stationary, upward-positive, horizontally periodic \(N_x\times N_y\) field.
- The mapped-depth factor is \(\gamma=1-h/D\), where \(D=\mathtt{wvt.Lz}\), and every accepted terrain must satisfy \(\gamma>0\) everywhere.
- The first supported transform is `WVTransformBoussinesq`.
- The transform coordinate `wvt.z` is interpreted as \(\xi\). The transform fields `wvt.u`, `wvt.v`, `wvt.w`, `wvt.eta`, and `wvt.p` are interpreted as \(\hat u\), \(\hat v\), \(\hat w\), \(\hat\eta\), and \(\hat p\).
- The existing rigid-lid eigenmodes remain the interior projection basis. The proof of concept must nevertheless test whether a separate bottom-buoyancy degree of freedom is required to represent the sloping-bottom displacement law.
- The forcing is evaluated in physical space and projected once by the existing `WVTransformBoussinesq.nonlinearFlux` path.
- The public forcing is named `WVExactTopographicForcing` to avoid collision with the conservative surrogate named `WVTopographicForcing` in the separate `conservative-terrain-forcing` repository.
- WaveVortexModel remains an external dependency. This repository will not initially be an MPM package and will not modify the WaveVortexModel authoring repository.

## Milestone 1: Repository and scientific contract

- [x] Complete

### Purpose

Establish the research-add-on structure and fix the scientific conventions and public interface before implementing the mapped dynamics.

### Dependencies

None.

### Deliverables

- Add MATLAB source at the repository root, with internal numerical helpers placed in `private/`.
- Add `UnitTests/` for automated tests and `Examples/` for reproducible model runs.
- Define the initial public constructor:

  ```matlab
  forcing = WVExactTopographicForcing(wvt,topographicHeight=h);
  ```

- Validate that `wvt` is a `WVTransformBoussinesq`.
- Validate that `topographicHeight` is real, finite, horizontally periodic, and has size \(N_x\times N_y\).
- Require \(\gamma=1-h/D>0\) everywhere and report structured class-specific errors for invalid inputs.
- Treat the terrain and all derived geometric factors as read-only scientific state after construction.
- Add a portable `runTests` entry point that resolves WaveVortexModel from an explicit option, then `WAVE_VORTEX_MODEL_ROOT`, then the sibling repository, and restores the original MATLAB path afterward.
- Add concise documentation stating that the forcing is exact in prescribed terrain height but linear in flow amplitude.

### Automated acceptance

- Constructor tests accept flat, uniform-height, and sinusoidal terrain.
- Constructor tests reject unsupported transforms, incorrect terrain dimensions, complex or nonfinite terrain, and any terrain for which \(\gamma\leq0\).
- The test suite resolves the external WaveVortexModel dependency without changing the user's permanent MATLAB path.
- No WaveVortexModel source file is modified or copied into this repository.

## Milestone 2: Mapped geometry and physical-field reconstruction

- [x] Complete

### Purpose

Implement and validate the stationary coordinate map independently of the terrain tendency and the WaveVortexModel forcing interface.

### Dependencies

Milestone 1.

### Deliverables

- Precompute \(\gamma\), \(\partial_x\gamma\), \(\partial_y\gamma\), \(\partial_x\ln\gamma\), and \(\partial_y\ln\gamma\) on the transform's horizontal grid.
- Store broadcast-ready representations of the mapped vertical coordinate \(\xi\) and the physical coordinate \(z=\gamma\xi\).
- Reconstruct the physical velocities from the hatted modal fields:

  \[
  u=\frac{\hat u}{\gamma},
  \qquad
  v=\frac{\hat v}{\gamma},
  \qquad
  w=\hat w+\xi\hat{\boldsymbol u}_H\boldsymbol{\cdot}\nabla_H\ln\gamma.
  \]

- Add diagnostic helpers for physical coordinates, physical velocity, mapped volume quadrature, and bottom-normal velocity.
- Use WaveVortexModel's horizontal spectral derivatives so that terrain gradients follow the same truncation convention as the state.

### Automated acceptance

- For \(h=0\), the mapped and physical coordinates and velocities agree to roundoff.
- For uniform \(h\), \(z=\gamma\xi\), \(u=\hat u/\gamma\), \(v=\hat v/\gamma\), and \(w=\hat w\) agree pointwise with their analytic values.
- For sinusoidal terrain and deterministic modal states, the reconstructed bottom velocity satisfies

  \[
  w_b=\boldsymbol u_{H,b}\boldsymbol{\cdot}\nabla_Hh
  \]

  to the accuracy of the horizontal spectral derivative.
- The mapped quadrature reduces to the ordinary transform quadrature when \(\gamma=1\).

## Milestone 3: Exact linear terrain tendency kernel

- [x] Complete

### Purpose

Implement the flow-linear, exact-in-terrain spatial tendency as a separately testable numerical kernel before connecting it to `WVForcing`.

### Dependencies

Milestones 1 and 2.

### Deliverables

- Begin with discretely constant \(N^2\), for which \(N^2(\gamma\xi)-N^2(\xi)=0\).
- Reconstruct \((\hat u,\hat v,\hat w,\hat\eta,\hat p)\) from the transform at the current model time.
- Evaluate the fixed-\(\xi\) pressure derivatives with WaveVortexModel's existing horizontal and \(F\)-to-\(G\) vertical derivative actions.
- Implement

  \[
  \mathcal T_u=
  \frac{1}{\rho_0}
  \left[
  (\gamma-1)\partial_x\hat p
  -\xi\partial_x\gamma\,\partial_\xi\hat p
  \right],
  \qquad
  \mathcal T_v=
  \frac{1}{\rho_0}
  \left[
  (\gamma-1)\partial_y\hat p
  -\xi\partial_y\gamma\,\partial_\xi\hat p
  \right].
  \]

- Compute the finite-terrain horizontal acceleration once:

  \[
  \mathcal H_u^L=f\hat v-\frac{1}{\rho_0}\partial_x\hat p-\mathcal T_u,
  \qquad
  \mathcal H_v^L=-f\hat u-\frac{1}{\rho_0}\partial_y\hat p-\mathcal T_v.
  \]

- Implement

  \[
  \mathcal T_w=
  \xi\boldsymbol{\mathcal H}^L\boldsymbol{\cdot}\nabla_H\ln\gamma
  +\frac{\gamma^{-1}-1}{\rho_0}\partial_\xi\hat p,
  \qquad
  \mathcal T_\eta=
  -\xi\hat{\boldsymbol u}_H\boldsymbol{\cdot}\nabla_H\ln\gamma.
  \]

- Return the right-hand-side tendency

  \[
  (F_u,F_v,F_w,F_\eta)
  =
  -(\mathcal T_u,\mathcal T_v,\mathcal T_w,\mathcal T_\eta)
  \]

  without performing a wave--vortex transform inside the kernel.
- Keep \(\hat p\) reconstruction from `wvt.p` as an explicit, testable scientific assumption rather than hiding it in the implementation.

### Automated acceptance

- Every component of the terrain tendency vanishes for \(h=0\).
- Uniform terrain eliminates every horizontal-slope contribution and leaves only the expected constant-\(\gamma\) pressure rescaling.
- Deterministic real modal states produce real spatial tendencies and conjugate spectral tendencies after projection.
- Repeated evaluation at the same state is deterministic and does not modify `wvt.t`, `Ap`, `Am`, or `A0`.

## Milestone 4: WaveVortexModel forcing integration

- [x] Complete

### Purpose

Connect the validated spatial kernel to the existing WaveVortexModel projection and adaptive integration path without introducing a modal terrain operator.

### Dependencies

Milestones 1 through 3.

### Deliverables

- Implement `WVExactTopographicForcing` as a `WVForcingType("NonhydrostaticSpatial")` subclass.
- In `addNonhydrostaticSpatialForcing`, add only the forcing's spatial contribution to the incoming `Fu`, `Fv`, `Fw`, and `Feta` arrays.
- Let `WVTransformBoussinesq.nonlinearFlux` perform the single call to `transformUVWEtaToWaveVortex` and its componentwise interaction-phase conversion.
- Do not transform fields, remove phases, or manipulate `Fp`, `Fm`, and `F0` inside the forcing.
- Run the linear terrain model through an ordinary `WVModel(wvt)` after calling `wvt.removeAllForcing()` and registering only `WVExactTopographicForcing`.
- Use the default adaptive `ode78` integrator for model tests.
- Leave antialiasing and resolution conversion disabled until Milestone 8.

### Automated acceptance

- A direct call to `wvt.nonlinearFlux` agrees with an independent projection of the terrain kernel through `transformUVWEtaToWaveVortex`.
- The agreement holds at \(t=t_0\) and at nonzero \(t-t_0\), including nondegenerate wave frequencies.
- Pre-existing spatial forcing arrays are preserved and incremented exactly once.
- A flat-terrain forcing produces zero coefficient tendency at every tested model time.
- An ordinary forcing call performs no modal-matrix construction and no explicit pressure solve.

## Milestone 5: Uniform-depth exact-solution oracle

- [ ] Complete

### Purpose

Validate every metric factor, the modal pressure reconstruction, the wave--vortex projection, and the interaction representation against an independently generated exact solution.

### Dependencies

Milestones 1 through 4.

### Deliverables

- Use constant terrain \(h=h_0\), giving physical depth

  \[
  H=D-h_0=\gamma D.
  \]

- Construct a flat `WVTransformBoussinesq` of physical depth \(H\) and a mapped `WVTransformBoussinesq` of reference depth \(D\), using identical horizontal geometry, resolution, latitude, and constant \(N^2\).
- Initialize one internal-wave mode in the depth-\(H\) transform and map its fields into the depth-\(D\) variables using

  \[
  \hat u=\gamma u,
  \qquad
  \hat v=\gamma v,
  \qquad
  \hat w=w,
  \qquad
  \hat\eta=\eta,
  \qquad
  z=\gamma\xi.
  \]

- Project the mapped initial state into the reference-depth transform and integrate it with only `WVExactTopographicForcing`.
- Compare against the independent depth-\(H\) solution and the analytic frequency

  \[
  \omega_H^2=
  \frac{N^2K^2+f^2(j\pi/H)^2}
  {K^2+(j\pi/H)^2}.
  \]

- Test \(\gamma=1\) and at least two nontrivial positive depth ratios.
- Repeat each nontrivial case with successively tighter adaptive relative and absolute tolerances.

### Automated acceptance

- The measured frequency agrees with \(\omega_H\) to the spatial and temporal discretization tolerance.
- At the tightest tolerance, the mapped velocity and displacement fields agree with the independent depth-\(H\) solution to \(10^{-8}\) relative error.
- Errors decrease consistently as adaptive tolerances are tightened.
- The \(\gamma=1\) case reduces to ordinary WaveVortexModel linear evolution.
- Any failure blocks Milestone 6 and is diagnosed as a metric, pressure, projection, initialization, or phase-convention error rather than absorbed into an empirical correction.

## Milestone 6: Sinusoidal-slope and scalar-representation test

- [ ] Complete

### Purpose

Exercise horizontal terrain coupling and determine whether the existing interior \(G\)-space is sufficient for the physical displacement variable at a sloping bottom.

### Dependencies

Milestone 5.

### Deliverables

- Use constant \(N^2\), a single low-wavenumber incident wave, and

  \[
  h(x)=h_a\cos(k_hx).
  \]

- Choose the incident and terrain modes far enough below truncation that the leading products do not alias.
- Verify that \(y\)-independent terrain preserves meridional wavenumber.
- Verify that the \(O(h_a/D)\) response contains the \(k-k_h\) and \(k+k_h\) sidebands.
- Repeat with at least three small values of \(h_a/D\) to distinguish leading linear scattering from higher-order exact-in-height harmonics.
- Verify the physical bottom kinematic relation from the reconstructed velocity.
- Evaluate the raw bottom displacement tendency from the mapped equation and compare it with the bottom tendency reconstructed after wave--vortex projection.
- Treat the displacement comparison as a scientific gate. If the homogeneous \(G\)-space cannot reproduce the required bottom law under resolution refinement, stop before Milestone 7 and revise the state representation to include an explicit bottom-buoyancy degree of freedom.

### Automated acceptance

- Spectral tendency outside the permitted meridional wavenumber is below \(10^{-12}\) relative to the total tendency.
- Leading sideband amplitudes converge linearly with \(h_a/D\).
- Harmonics absent from the first-order terrain expansion decrease quadratically with \(h_a/D\).
- The physical bottom-normal velocity residual converges to the horizontal spectral differentiation error.
- The bottom displacement-tendency comparison either converges with resolution or produces an explicit blocked diagnostic identifying the missing boundary degree of freedom; it must not be silently ignored.

## Milestone 7: Finite-terrain conservation and end-to-end linear run

- [ ] Complete

### Purpose

Verify that the projected spatially truncated system retains the finite-terrain quadratic laws and produces a stable multi-period linear scattering calculation.

### Dependencies

Milestone 6, including successful completion of its scalar-representation gate.

### Deliverables

- Implement the finite-terrain linear energy diagnostic

  \[
  \mathcal E_\gamma^{(2)}
  =
  \frac{\rho_0}{2A}
  \int_A\int_{-D}^{0}
  \left[
  \gamma^{-1}(\hat u^2+\hat v^2)
  +\gamma w^2
  +\gamma N^2(\gamma\xi)\hat\eta^2
  \right]d\xi\,dA.
  \]

- Implement \(q_\gamma\) and

  \[
  \mathcal Z_\gamma^{(2)}
  =
  \frac{1}{2A}
  \int_A\int_{-D}^{0}
  \gamma q_\gamma^2\,d\xi\,dA
  \]

  using the mapped physical-horizontal derivatives and the same spectral differentiation conventions as the forcing.
- Compute instantaneous invariant tendencies from a projected terrain RHS independently of full time integration.
- Run a multi-period sinusoidal-hill calculation with only the exact linear terrain forcing and adaptive `ode78`.
- Repeat the run with tighter integration tolerances and increased spatial resolution.
- Report `wvt.totalEnergy` only as the intentionally nonconserved flat-reference comparison.

### Automated acceptance

- Instantaneous relative tendencies of \(\mathcal E_\gamma^{(2)}\) and \(\mathcal Z_\gamma^{(2)}\) converge to the spatial truncation and quadrature error.
- Integrated invariant drift decreases with adaptive tolerance until it reaches the spatial-error floor.
- Increasing resolution lowers the spatial-error floor.
- The solution preserves coefficient conjugacy and remains finite for the complete benchmark interval.
- The example reports incident, reflected, and higher-mode energy together with finite-terrain invariant drift.

## Milestone 8: Dealiasing, resolution changes, and performance

- [ ] Complete

### Purpose

Make the validated linear forcing compatible with WaveVortexModel's production pseudospectral path and demonstrate that its runtime cost remains comparable to ordinary spatial nonlinear forcing.

### Dependencies

Milestones 1 through 7.

### Deliverables

- Implement `forcingWithResolutionOfTransform` by evaluating or Fourier-resampling the canonical terrain on the target transform and rebuilding every derived mapped factor.
- Enable WaveVortexModel's standard antialiasing path.
- Compare dealiased results with deliberately band-limited undealiased calculations for which aliasing is analytically excluded.
- Reuse reconstructed fields and derivative results within one RHS evaluation to avoid repeated transforms and unnecessary temporary arrays.
- Profile a terrain RHS evaluation against `WVNonlinearAdvection` at three increasing resolutions.
- Confirm that runtime work consists of field reconstruction, pseudospectral derivatives and products, one existing wave--vortex projection, and no modal-matrix construction or diagnostic pressure solve.

### Automated acceptance

- Resolution-converted forcing has terrain arrays and mapped factors compatible with the target transform.
- Band-limited dealiased and undealiased tendencies agree to \(10^{-10}\) relative error at the reference resolution.
- Dealiased sinusoidal runs converge under horizontal and vertical refinement.
- Benchmarks report wall time and peak stored bytes for the terrain forcing and `WVNonlinearAdvection`.
- Instrumentation confirms that an ordinary forcing call constructs no dense modal operator and performs no explicit pressure solve.

## Milestone 9: Arbitrary stationary stratification

- [ ] Complete

### Purpose

Extend the validated constant-stratification implementation to the full stationary \(N^2(z)\) scope without changing its projection architecture.

### Dependencies

Milestone 8.

### Deliverables

- Evaluate \(N^2(\gamma\xi)\) through the transform's stationary stratification profile and validate that every mapped value is real, finite, and positive.
- Add the linear stratification terrain term

  \[
  \left[N^2(\gamma\xi)-N^2(\xi)\right]\hat\eta
  \]

  to \(\mathcal T_w\).
- Route constant and nonconstant profiles through one generalized scientific implementation.
- Use an exponential \(N^2(z)\) profile as the primary nonconstant test case.
- Repeat the uniform-depth oracle, sinusoidal coupling tests, finite-terrain conservation diagnostics, antialiasing checks, and resolution-convergence runs.

### Automated acceptance

- A constant profile passed through the generalized path agrees with the Milestone 8 constant-\(N\) tendency to \(10^{-12}\) relative error.
- The exponential-profile uniform-depth run agrees with an independent flat transform of the corresponding physical depth and profile to \(10^{-8}\) relative error at the tightest tolerance.
- The exponential-profile sinusoidal run satisfies the symmetry, boundary, conjugacy, invariant, and resolution-convergence criteria established in Milestones 6 through 8.
- Invalid mapped stratification values are rejected with structured errors before integration begins.

## Planning and implementation cadence

### Planning batch A: Milestones 1--4

Plan Milestones 1--4 together. Their public interface, field conventions, derivative choices, forcing accumulation, and interaction-phase handling form one vertical slice and must not be designed independently. Implement them sequentially so that geometry and the spatial kernel can be tested before the `WVForcing` wrapper.

### Planning batch B: Milestones 5--7

Plan Milestones 5--7 together after Milestone 4 passes. The uniform-depth oracle, sloping-boundary test, pressure/state validation, and finite-terrain invariant diagnostics must share one state-mapping and comparison harness. Implement Milestone 5 first. Milestone 6 is a required scientific gate, and Milestone 7 must not begin until the displacement representation passes that gate.

### Planning batch C: Milestones 8--9

Plan Milestones 8 and 9 together only after the linear scientific gates pass. Dealiasing, resolution conversion, and arbitrary stratification all require terrain and profile data to be reconstructed on alternate transforms and should share that infrastructure.

### Deferred nonlinear planning

Do not plan the complete \(\mathcal N+\mathcal T+\mathcal B\) implementation until Milestones 1--9 establish the linear formulation. That work will replace, rather than supplement, `WVNonlinearAdvection` and requires a separate roadmap for nonlinear APE, APV, aliasing, and buoyancy validation.

## Definition of done

The exact linear topographic-forcing roadmap is complete when:

- the flat-terrain tendency vanishes;
- the uniform-depth calculation matches its independent exact solution;
- the sinusoidal calculation satisfies its symmetry, sideband, amplitude-scaling, and bottom-kinematic checks;
- the selected displacement representation passes the bottom-boundary gate;
- instantaneous finite-terrain energy and potential-enstrophy tendencies converge to the spatial-discretization error;
- adaptive model errors and invariant drift converge with tolerance and resolution;
- constant and arbitrary stationary stratification pass the same scientific tests;
- dealiased and deliberately alias-free calculations agree; and
- ordinary forcing calls remain matrix-free and use the existing WaveVortexModel projection.

## Out of scope

- Complete nonlinear terrain-aware advection \(\mathcal N\) and nonlinear buoyancy \(\mathcal B\).
- Dense or sparse modal terrain operators.
- The flat-energy- and flat-potential-enstrophy-preserving conservative surrogate.
- Modifications to WaveVortexModel or a new `WVTransform` subclass.
- NetCDF persistence, restart support, MPM packaging, release automation, and generated website documentation.
- A new time integrator or a requirement that an adaptive time step conserve invariants to machine precision after every accepted step.
