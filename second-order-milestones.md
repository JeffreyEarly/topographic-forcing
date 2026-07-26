# Second-order bottom wave-scattering milestones

## Objective

Develop and validate a strict second Born approximation for linear wave scattering by weak stationary topography. The calculation will retain the ordinary rigid-lid wave basis and the pressure-weighted mean-bottom projection already implemented by `WVBottomWaveScatteringForcing`, while keeping the zeroth-, first-, and second-order coefficient fields separate.

With `topographicHeight` equal to the physical terrain height, the stored interaction-representation contributions will satisfy

```math
\dot{\boldsymbol A}^{(0)}=0,
\qquad
\dot{\boldsymbol A}^{(1)}
=
\mathcal T_h(t)\boldsymbol A^{(0)},
\qquad
\dot{\boldsymbol A}^{(2)}
=
\mathcal T_h(t)\boldsymbol A^{(1)}.
```

The reconstructed approximation is

```math
\boldsymbol A^{[2]}
=
\boldsymbol A^{(0)}
+
\boldsymbol A^{(1)}
+
\boldsymbol A^{(2)}.
```

The production calculation will remain pressure-free and matrix-free. One right-hand-side evaluation will apply the existing two-dimensional bottom projection twice and will never feed the second-order result back into itself.

## Fixed scientific and numerical contract

- The target is `WVTransformBoussinesq` with stationary periodic terrain and arbitrary stationary stratification supported by the transform.
- The initial state may contain either wave branch. Balanced coefficients remain constant and do not participate in the scattering calculation.
- Complete wave modes include their spatial and temporal phases. `ApOrder0`, `ApOrder1`, and `ApOrder2`, and the corresponding minus-branch arrays, are interaction-representation coefficients.
- The strict hierarchy is an asymptotic approximation at fixed time. Validation times remain a fixed number of reference-depth wave periods as terrain amplitude tends to zero.
- The exact second-order bottom velocity requires the strong first-order endpoint fields:

```math
g_b^{(2)}
=
\boldsymbol u_{H,d}^{(1)}
\boldsymbol{\cdot}\nabla_H h
-
h\,\partial_z w_d^{(1)}.
```

- Reapplying the existing finite modal endpoint kernel is the proposed second-Born closure. Its agreement with the strong boundary-value problem is a blocking scientific question, not an assumed consequence of interior energy-norm convergence.
- The ordinary autonomously iterated first-order forcing remains available as a separate model. The strict hierarchy will not change its behavior or persistence format.
- The default adaptive `ode78` path will be used. The implementation will be agnostic to `shouldAntialias` and will preserve explicit antialias resolution rebuilding.

## Milestone 1: Shared arbitrary-state boundary projection

- [ ] Complete

### Purpose

Make the existing first-order scattering kernel reusable for specified coefficient arrays without changing the public behavior of `WVBottomWaveScatteringForcing`.

### Dependencies

The completed first-order Milestone 9.

### Deliverables

- Add a documented developer method:

  ```matlab
  [Fp,Fm,gBottom,bottomFields] = ...
      forcing.tendencyFromWaveCoefficients(Ap,Am,t);
  ```

- Validate that `Ap` and `Am` match the originating transform layout and that `t` is finite.
- Return the endpoint fields `u`, `v`, and `dWdz` used to construct

  ```math
  g_b
  =
  u_d\,\partial_x h
  +
  v_d\,\partial_y h
  -
  h\,\partial_z w_d.
  ```

- Make `bottomVelocityFromWaveState` and `addSpectralForcing` delegate to this single calculation.
- Keep terrain transfer, conjugacy reconstruction, phase conversion, pressure projection, and incoming `F0` behavior unchanged.

### Automated acceptance

- The refactored forcing agrees with the pre-refactor test oracle within $10^{-13}$ relative error.
- Explicit-array evaluation agrees with `addSpectralForcing` at nonzero model time within $10^{-13}$.
- Input coefficient arrays and the transform state are not modified.
- Flat terrain gives zero tendency, and spectral conjugacy closes at roundoff.
- Constant and variable stratification and both antialias settings pass the existing suite.

## Milestone 2: Strict second-order adaptive state

- [ ] Complete

### Purpose

Represent the triangular perturbation hierarchy as a pure adaptive right-hand side rather than hidden mutable state inside a forcing callback.

### Dependencies

Milestone 1.

### Deliverables

- Add `WVSecondOrderBottomWaveScattering` as a `WVObservingSystem` subclass with the public construction:

  ```matlab
  model = WVModel(wvt,shouldUseLinearDynamics=true);
  hierarchy = WVSecondOrderBottomWaveScattering( ...
      model,topographicHeight=h);
  model.addFluxedObservingSystem(hierarchy);
  ```

- Capture the transform's initial wave coefficients as read-only `ApOrder0` and `AmOrder0`.
- Integrate four arrays: `ApOrder1`, `AmOrder1`, `ApOrder2`, and `AmOrder2`, initialized to zero unless explicitly restored.
- At every trial time, evaluate

  ```math
  (\dot A_+^{(1)},\dot A_-^{(1)})
  =
  \mathcal T_h(t)
  (A_+^{(0)},A_-^{(0)})
  ```

  and

  ```math
  (\dot A_+^{(2)},\dot A_-^{(2)})
  =
  \mathcal T_h(t)
  (A_+^{(1)},A_-^{(1)}).
  ```

- Update `wvt.Ap` and `wvt.Am` with the sum of all three orders at accepted output times. Leave `wvt.A0` equal to its initial value.
- Reject simultaneous use of the ordinary `WVCoefficients` flux system or an autonomously registered scattering forcing, because either would double-advance the wave state.
- Supply order-appropriate adaptive absolute tolerances using `Apm_TE_factor`.

### Automated acceptance

- A zero first-order correction produces an initially zero second-order tendency.
- Doubling terrain doubles the first-order tendency and quadruples the second-order tendency evaluated on a correspondingly doubled first-order state.
- The transform always reports the exact sum of the three stored orders after accepted updates.
- Trial right-hand-side calls are deterministic and do not depend on rejected adaptive steps.
- The balanced coefficients remain bitwise unchanged.
- Tightening adaptive tolerance gives convergent order-separated and total coefficients.

## Milestone 3: Uniform-depth blocking oracle

- [ ] Complete

### Purpose

Determine whether the projected second Born closure reproduces the physical second-order boundary perturbation before applying it to sloping terrain.

### Dependencies

Milestone 2.

### Deliverables

- Extend the constant-$N$ uniform-depth benchmark to compare the strict hierarchy with the independent exact solution at physical depth $H=D-h_0$.
- Expand the exact frequency and phase through second order about $H=D$.
- Compare the order-separated coefficients with the corresponding derivatives of the exact depth-$H$ solution.
- Evaluate the first-order correction's bottom horizontal velocity and vertical-velocity derivative independently from the exact solution, and compare them with the finite modal endpoint reconstruction used by the second kernel call.
- Repeat the calculation at increasing vertical resolution and several signed uniform offsets. Negative offsets are used only by the centered derivative oracle and represent a deeper comparison domain.

### Automated acceptance

- First-order coefficient error scales as $O(h_0^2)$.
- The total strict second-order coefficient and phase errors scale as $O(h_0^3)$, with an observed exponent between 2.8 and 3.2.
- The finest-amplitude total relative coefficient error is below $10^{-8}$.
- Independently evaluated and modal-reconstructed first-order endpoint fields converge to within $10^{-9}$ relative error.
- Interior linear QGPV remains below $10^{-12}$, and `A0` remains unchanged.
- Failure of endpoint convergence blocks Milestones 4--7 and prevents describing the implementation as physically second order.

## Milestone 4: Sinusoidal-terrain second Born structure

- [ ] Complete

### Purpose

Verify horizontal selection rules, terrain-amplitude ordering, and the relationship between the strict truncation and autonomous feedback.

### Dependencies

Milestone 3.

### Deliverables

- Use one incident wave and

  ```math
  h(\boldsymbol x)
  =
  h_a\cos(\boldsymbol q\boldsymbol{\cdot}\boldsymbol x).
  ```

- Verify that the first-order response is confined to $\boldsymbol K\pm\boldsymbol q$.
- Verify that the second-order response is confined to $\boldsymbol K$ and $\boldsymbol K\pm2\boldsymbol q$.
- Compare the strict total with `WVBottomWaveScatteringForcing` over fixed time while reducing $h_a$.
- Evaluate the translated bottom condition through second order from independently reconstructed order-separated fields.
- Repeat at increasing horizontal and vertical resolution to distinguish amplitude ordering from truncation error.

### Automated acceptance

- First-order coefficient amplitude scales as $h_a$ with exponent $1\pm0.05$.
- Second-order coefficient amplitude scales as $h_a^2$ with exponent $2\pm0.05$.
- Spectral energy outside the permitted first- and second-order supports is below $10^{-12}$ of the corresponding order total.
- The strict second-order and autonomous discrete solutions differ as $O(h_a^3)$ at fixed time.
- The translated-bottom residual scales as $O(h_a^3)$ after spatial convergence.
- Both antialias settings give the same band-limited result at common retained modes.

## Milestone 5: Second-order energy and no-PV gates

- [ ] Complete

### Purpose

Verify the conservation laws of the separated perturbation hierarchy independently of coefficient accuracy.

### Dependencies

Milestone 4.

### Deliverables

- Diagnose

  ```math
  \mathcal E^{(0)}
  =
  \frac12
  \left\langle
  \Psi^{(0)},\Psi^{(0)}
  \right\rangle_{\mathcal E_0},
  ```

  the existing $\mathcal E^{(1)}$, and

  ```math
  \mathcal E^{(2)}
  =
  \operatorname{Re}
  \left\langle
  \Psi^{(0)},\Psi^{(2)}
  \right\rangle_{\mathcal E_0}
  +
  \frac12
  \left\langle
  \Psi^{(1)},\Psi^{(1)}
  \right\rangle_{\mathcal E_0}
  -
  \frac{\rho_0}{A}
  \int_A
  h\,
  \operatorname{Re}
  \left(
  \boldsymbol u_{H,d}^{(0)*}
  \boldsymbol{\cdot}
  \boldsymbol u_{H,d}^{(1)}
  \right)dA.
  ```

- Compute instantaneous energy-coefficient tendencies directly from the two order-separated right-hand sides.
- Reconstruct the linear QGPV tendency at each perturbation order.
- Diagnose first- and second-order bottom displacement afterward from accepted output:

  ```math
  \partial_t\eta_d^{(1)}=g_b^{(1)},
  \qquad
  \partial_t\eta_d^{(2)}=g_b^{(2)}.
  ```

### Automated acceptance

- Instantaneous tendencies of $\mathcal E^{(0)}$, $\mathcal E^{(1)}$, and $\mathcal E^{(2)}$ close within $10^{-11}$ relative to the individual exchange terms.
- Drift of each retained energy coefficient converges with adaptive tolerance.
- The reconstructed strict physical-energy residual scales as $O(h^3)$.
- Both wave-order QGPV tendencies remain below $10^{-12}$.
- Direct balanced tendencies are exactly zero.
- Bottom-displacement quadrature converges at second order with saved output interval.

## Milestone 6: Examples and comparison

- [ ] Complete

### Purpose

Show where the second-order approximation adds useful scattering physics and where the weak-terrain expansion remains inadequate.

### Dependencies

Milestone 5.

### Deliverables

- Add `SecondOrderUniformDepthWaveScatteringBenchmark` showing zeroth-, first-, and second-order phase errors against the exact depth-$H$ solution.
- Add `SecondOrderSinusoidalRidgeScatteringExample` showing the incident mode, first sidebands, second harmonics, energy coefficients, and bottom-condition residual.
- Extend the Gaussian-ridge workflow with an optional strict second-order calculation and side-by-side diagnostics against the autonomous first-order model.
- Keep the controlled default $h_0/D=0.05$. Treat larger ridges as exploratory until amplitude and resolution convergence are demonstrated.
- Reuse the corrected wet-domain movie visualization only after standard total coefficients can be written at accepted output times.

### Automated acceptance

- Both examples render finite figures without requiring generated files.
- The uniform-depth example passes the cubic-error gate at its automated reference resolution.
- The sinusoidal example displays nonzero second harmonics with the expected quadratic amplitude.
- The Gaussian example reports order-separated energy, reflected energy, higher-mode energy, and runtime without claiming finite-amplitude accuracy.
- Default and tighter adaptive tolerances agree to their requested error levels.

## Milestone 7: Persistence, resolution changes, and performance

- [ ] Complete

### Purpose

Make the validated hierarchy usable for restartable research calculations without changing its scientific definition.

### Dependencies

Milestone 6.

### Deliverables

- Persist terrain, fixed zeroth-order coefficients, all four integrated correction arrays, hierarchy version, and construction metadata in a dedicated observing-system group.
- Write the total `Ap`, `Am`, and unchanged `A0` through the standard `wave-vortex` output group.
- Reconstruct the hierarchy during `WVModel.modelFromFile` and continue with the ordinary adaptive integrator.
- Rebuild terrain factors and all order arrays through `observingSystemWithResolutionOfTransform`.
- Profile one right-hand-side evaluation and confirm it contains two boundary-kernel applications and no three-dimensional transform, diagnostic pressure solve, or modal matrix.

### Automated acceptance

- Restart continuation agrees with an uninterrupted control within $10^{-10}$ in every order-separated coefficient array.
- Resolution conversion preserves terrain, conjugacy, phase, and total-state reconstruction.
- Explicit-antialias rebuilding agrees with the equivalent manually constructed transform.
- Runtime is approximately twice the terrain work of one autonomous scattering-forcing call, reported separately from adaptive step-count differences.
- The complete repository suite, `checkcode`, and `git diff --check` pass.

## Planning and implementation cadence

- Plan and implement Milestones 1--2 together. They define one pure reusable kernel and the minimum triangular adaptive state.
- Plan Milestones 3--5 together. Milestone 3 is a blocking physical oracle; Milestones 4--5 may proceed only after it passes.
- Plan Milestones 6--7 together after the scientific gates pass. Examples, persistence, and optimization must not redefine the second-order equations.

## Definition of done

The strict second-order proof of concept is established when Milestones 1--5 pass: the endpoint reconstruction converges to the independent uniform-depth solution, first and second scattering obey their selection and amplitude laws, the total fixed-time error is cubic in terrain amplitude, every retained energy coefficient is conserved to spatial accuracy, and no interior linear QGPV is generated.

The implementation is ready for restartable research examples when Milestones 6--7 also pass. The original prescribed generator and autonomous first-order scattering forcing remain supported and scientifically distinct.

Dynamic barotropic backreaction, independent geostrophic bottom buoyancy, nonlinear terrain dynamics, a uniformly valid long-time resummation, and exact finite-amplitude topography remain outside this roadmap.
