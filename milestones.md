# Terrain-energy Galerkin milestones

> **Checkpoint after Milestone 6.8:** the boundary-complete weak oracle derives APV from terrain-dependent geostrophic test states inside the primitive weak space and closes the trusted tangent identities at roundoff. This resolves the Milestone-6.7 representation question without replacement APV rows or corrections. Milestones 7–13 remain inactive until separately planned.

## Objective

Develop `WVTerrainEnergyGalerkin`, a boundary-dynamical Galerkin system for the linear rotating Boussinesq equations over stationary bottom topography. The formulation is linear in flow amplitude and exact in the resolved terrain. The immediate goal is to construct one adjoint-consistent Fourier–Galerkin projection that conserves physical energy, enforces projected bottom evolution, and either conserves projected APV exactly or converges to stationary APV on a trusted physical band. Pressure remains a Lagrange multiplier during construction and is removed only from the validated reduced evolution.

The augmented state begins with

```math
\boldsymbol x
=
\begin{pmatrix}
\boldsymbol a\\
\boldsymbol b
\end{pmatrix},
\qquad
B\dot{\boldsymbol x}=R_h\boldsymbol x,
```

where $\boldsymbol a$ contains the existing volume coordinates and $\boldsymbol b$ contains independent bottom coordinates. The mathematical specification begins with `terrain-energy-galerkin.tex` at commit `72c967c` in the `ape-apv-bottom-topography` literature repository. The boundary-energy and potential-enstrophy analysis in `boundary-energy-enstrophy.tex` at commit `716c216` supplies the subsequent correction to the numerical state: a nonzero bottom displacement cannot be represented as an isolated scalar lift if its dynamically associated balanced velocity, pressure, and APV structure are omitted.

The implementation will use ordinary hydrostatic modes as economical vertical coordinates. The public coefficient ordering remains fixed, and the complete zero-APV balanced bottom inversion remains available for field reconstruction. The augmented descriptor will nevertheless keep the bottom coefficient independent until the volume equations, pressure constraint, and strong bottom row have been assembled together. The boundary Green identity will then determine whether the reduced problem uses a positive generalized metric, a signed metric, or physical energy alone.

The completed mean-depth generator and scattering implementation is retained as reusable engineering infrastructure. Its scientific roadmap is archived in [mean-depth-wave-generator-milestones.md](mean-depth-wave-generator-milestones.md). Neither existing forcing is used as the terrain-energy evolution operator.

The completed dense finite-terrain forms remain diagnostic oracles. Milestones 4.5–6.3 record the successive bottom-state, pressure-ordering, local-slope, and primitive-coordinate audits. Milestones 6.4–6.5 show that coupled volume and boundary PV provide a closed local and periodic QG representation. Milestone 6.6 shows that simply using projected APV and bottom PV as replacement coordinates yields a complete descriptor but not one equivalent to the primitive weak dynamics. Its row-equivalence and physical-energy failures remain decisive. Milestone 6.7 returns to the unmodified primitive energy and exchange forms and shows that ordinary dealiased projection alone does not recover the APV cancellation. Milestone 6.8 supplies the missing terrain derivative of the geostrophic test space and verifies the resulting weak APV Green identity directly.

## Fixed conventions and boundaries

- The public target is `WVTransformBoussinesq`; WaveVortexModel remains an external dependency and is not modified.
- The initial balanced bottom-inversion construction requires a nonzero, spatially uniform Coriolis parameter $f$. Equatorial and variable-$f$ extensions are outside the present roadmap.
- `topographicHeight` is a real, finite, stationary, upward-positive, horizontally periodic field on the transform grid.
- The mapped geometry is

```math
\gamma=1-\frac{h}{D}>0.
```

- The prognostic state contains volume coordinates for $(\hat u,\hat v,\hat w,\hat\eta)$ and one independent bottom-displacement coordinate for each retained horizontal Fourier coefficient. The complete balanced bottom inversion supplies linked fields for reconstruction without replacing that independent boundary coordinate. Pressure remains a Lagrange multiplier while the augmented constrained system is assembled and is absent from ordinary validated reduced evolution.
- The velocity satisfies mapped continuity and homogeneous normal-flow conditions,

```math
\nabla_\xi\boldsymbol{\cdot}\hat{\boldsymbol u}=0,
\qquad
\hat w(-D)=\hat w(0)=0.
```

- Surface displacement vanishes. The bottom value of $\hat\eta$ is unrestricted and evolves through the displacement equation.
- Independent surface-buoyancy anomalies are excluded. Bottom buoyancy is retained through the zero-frequency balanced bottom-inversion states; it is not represented by a displacement-only coordinate.
- Hydrostatic modes are coordinates only. The projected equations, recovered flat modes, finite-terrain modes, and evolution remain nonhydrostatic.
- Galerkin states use complex horizontal Fourier coefficients with explicit conjugacy maps for real physical fields.
- Terrain products use a common oversampled quadrature and adjoint-consistent reconstruction and projection.
- The formulation is not an additive `WVForcing`; ordinary evolution acts through the approved finite-terrain generator, terrain-mode phases, or an energy-preserving implicit step.

## Milestone 1: Repository and scientific contract

- [x] Complete

### Purpose

Establish the standalone public API and state conventions without changing the existing mean-depth forcing classes or WaveVortexModel.

### Dependencies

None.

### Deliverables

- Add `WVTerrainEnergyGalerkin` as a standalone scientific class with the expensive construction path

  ```matlab
  problem = WVTerrainEnergyGalerkin.fromTopography( ...
      wvt, ...
      topographicHeight=h, ...
      verticalModeIndices=j, ...
      horizontalOversamplingFactor=2);
  ```

- Support only `WVTransformBoussinesq` initially.
- Validate terrain size, reality, finiteness, periodic layout, and $\gamma>0$ with structured class-specific errors.
- Treat `verticalModeIndices` as the retained hydrostatic reference coordinates and require a positive integer `horizontalOversamplingFactor`.
- Store authoritative geometry, terrain, basis settings, originating transform metadata, and construction diagnostics as read-only scientific state.
- Document units, coefficient ordering, bottom-displacement coordinates, real-field conjugacy, energy normalization, and the distinction between the Galerkin system and `WVForcing`.
- Extend the portable `runTests` entry point without changing permanent MATLAB paths.

### Automated acceptance

- Constructor tests accept flat, uniform-offset, sinusoidal, and broadband terrain.
- Constructor tests reject unsupported transforms, bad terrain dimensions, complex or nonfinite terrain, $\gamma\leq0$, invalid mode indices, and invalid oversampling factors.
- The originating transform and authoritative terrain are not modified during construction.
- Existing mean-depth generation and scattering tests remain green.
- No Galerkin class subclasses `WVForcing` or registers an additive terrain forcing.

## Milestone 2: Mixed hydrostatic and bottom-displacement prototype

- [x] Complete — historical prototype superseded by Milestone 4.7

### Purpose

Construct the initial mixed coordinate space used to test the finite-terrain weak equations. This milestone established the reconstruction and adjoint machinery, but its displacement-only bottom coordinate was later shown to be dynamically incomplete.

### Dependencies

Milestone 1.

### Deliverables

- Construct a matched `WVTransformHydrostatic` using the target domain, grid, stratification, Coriolis parameter, density, gravity, and antialias convention.
- Use its ordinary $F$ and $G$ modes for the homogeneous interior fields.
- Add one bottom-displacement coefficient per retained horizontal wavenumber with

```math
\chi_{\kappa b}(\xi)=
\begin{cases}
\sinh(-\kappa\xi)/\sinh(\kappa D),&\kappa>0,\\
-\xi/D,&\kappa=0.
\end{cases}
```

- Represent displacement as

```math
\hat\eta_{\boldsymbol K}(\xi)
=
\eta_{b,\boldsymbol K}\chi_{\kappa b}(\xi)
+
\sum_j\eta_{\boldsymbol K}^jG^j(\xi).
```

- Implement deterministic packing and unpacking for the interior modal coefficients and bottom coefficients.
- Implement reconstruction of $(\hat u,\hat v,\hat w,\hat\eta)$, reconstruction of physical $(u,v,w)$, and the exact mapped volume quadrature.
- Implement the adjoints of every reconstruction map using the same vertical weights and horizontal normalization.
- Implement real-field conjugacy and the discrete finite-terrain APV map.
- Verify directly that writing the linear displacement equation using $\eta_i=(1-\gamma)\xi+\hat\eta$ produces the same interior and bottom evolution as the $\hat\eta$ formulation.

### Automated acceptance

- Pack/unpack and reconstruction/projection round trips close within $10^{-12}$ relative error.
- The bottom function has the requested endpoint values to roundoff and introduces exactly one bottom value per horizontal coefficient.
- Reconstructed velocity satisfies mapped continuity and homogeneous $\hat w$ boundary values to $10^{-12}$.
- Real coefficient sets reconstruct real fields and return with conjugacy defects below $10^{-12}$.
- Numerical reconstruction and projection satisfy their weighted adjoint identity within $10^{-12}$.
- The $\eta_i$ and $\hat\eta$ tendencies agree at every grid point and at the bottom within $10^{-12}$.

## Milestone 3: Flat nonhydrostatic dense-oracle prototype

- [x] Complete — historical prototype to be rerun in the Milestone-4.7 basis

### Purpose

Prove that the initial mixed hydrostatic coordinates recover the flat nonhydrostatic wave–vortex frequencies and identify the state-space defect exposed by finite terrain. Milestone 4.7 must repeat this oracle after replacing the displacement-only bottom coordinate.

### Dependencies

Milestone 2.

### Deliverables

- Form dense $E_0$, $J_0$, and $Q_0$ by probing every mixed-basis coordinate at low resolution.
- Solve the generalized Hermitian problem independently at each horizontal wavenumber:

```math
iJ_0\boldsymbol c=\omega E_0\boldsymbol c.
```

- Recover internal-wave, inertial, geostrophic, mean-density-anomaly, and independent bottom-buoyancy subspaces.
- For constant stratification, compare with

```math
\omega_{\kappa j}^2
=
\frac{N^2\kappa^2+f^2m_j^2}
{\kappa^2+m_j^2}.
```

- For arbitrary stationary stratification, compare frequencies, energy-normalized eigenfunctions, polarization, and APV with `WVTransformBoussinesq`.
- Report convergence with retained hydrostatic vertical modes and the indicator $\kappa/m_j$.
- Record raw structural residuals before any explicit symmetrization.

### Automated acceptance

- Raw matrices satisfy

```math
\frac{\lVert E_0-E_0^*\rVert}{\lVert E_0\rVert}\leq10^{-13},
\qquad
\frac{\lVert J_0+J_0^*\rVert}{\lVert J_0\rVert}\leq10^{-13}.
```

- $E_0$ is positive definite on the retained state space.
- Constant-stratification frequencies agree with the analytic dispersion relation to approximately $10^{-11}$ relative error.
- Variable-stratification frequencies and eigenfunctions converge to the directly computed nonhydrostatic modes as vertical coordinates are added.
- Nonzero-frequency modes have negligible discrete APV and the zero-frequency dimension matches the then-expected balanced, MDA, and displacement-only bottom coordinates.
- At this stage, failure to recover the flat problem or acceptable conditioning blocked Milestone 4. The revised flat oracle in Milestone 4.7 now governs further development.

## Milestone 4: Exact finite-terrain dense-form prototype

- [x] Complete — historical prototype to be rebuilt after Milestone 4.7

### Purpose

Establish the finite-terrain weak forms, exact in the resolved terrain height $h$, in the initial displacement-only basis. These forms and their quadrature implementation remain reusable, but their matrices are not the final oracle because the underlying boundary state was incomplete.

### Dependencies

Milestone 3.

### Deliverables

- Evaluate

```math
\gamma,\qquad
\gamma^{-1},\qquad
\nabla_H\ln\gamma,\qquad
N^2(\gamma\xi)
```

  on a common oversampled grid.
- Form dense $E_\gamma$, $J_\gamma$, and $Q_\gamma$ by direct quadrature of every mixed-basis pair.
- Preserve the raw matrices and their structural defects for diagnosis before any roundoff-level restoration.
- Verify the exact flat limit and random-state identity

```math
\mathrm{Re}
\left(\boldsymbol a^*J_\gamma\boldsymbol a\right)=0.
```

- Use constant $\gamma$ as an exact mapped-coordinate oracle and compare against an independent flat transform of physical depth $H=\gamma D$.
- Separate quadrature, truncation, and eigensolver errors in the diagnostics.

### Automated acceptance

- Raw finite-terrain structural residuals are below $10^{-12}$.
- $E_\gamma$ is positive definite and remains well conditioned at the documented reference terrains.
- Random-state finite-terrain energy tendencies close to roundoff.
- For $h=0$, all finite-terrain matrices reproduce the Milestone-3 flat matrices within $10^{-12}$.
- Uniform-depth frequencies, mapped eigenfunctions, and energy normalization converge to the independent solution of physical depth $H$ within $10^{-9}$.
- Repeated dense construction is deterministic to $10^{-13}$.

## Milestone 4.5: Constrained conservative closure audit

- [x] Complete — historical incompatible exit; superseded by Milestone 4.7

### Purpose

Determine whether the then-current mixed state space admits an energy-preserving evolution that also conserves discrete finite-terrain APV and enforces the strong bottom-displacement evolution exactly. Do not assume that these constraints are jointly feasible.

Milestone 4 guarantees the raw energy identity but does not imply

```math
Q_\gamma E_\gamma^{-1}J_\gamma=0
```

or

```math
B E_\gamma^{-1}J_\gamma=R_h,
```

where $Q_\gamma$ maps state coefficients to resolved APV, $B$ extracts the bottom value of $\hat\eta$, and $R_h$ maps the state to $\boldsymbol u_{H,b}\boldsymbol{\cdot}\nabla_Hh$. The failure of either identity prevents the raw finite-dimensional system from passing Milestone 5 even though it conserves energy.

### Dependencies

Milestone 4.

### Deliverables

- Preserve the raw Milestone-4 matrices and diagnostics unchanged.
- Factor the finite-terrain energy as

  ```math
  E_\gamma=S^*S
  ```

  and introduce energy coordinates $\boldsymbol b=S\boldsymbol a$ with raw generator

  ```math
  K_0=S^{-*}J_\gamma S^{-1},
  \qquad
  K_0=-K_0^*.
  ```

- Form the transformed constraint maps

  ```math
  \overline Q=Q_\gamma S^{-1},
  \qquad
  \overline B=BS^{-1},
  \qquad
  \overline R=R_hS^{-1}.
  ```

- Use rank-revealing dense linear algebra to determine whether there exists a complete-state operator $K_c$ satisfying

  ```math
  K_c=-K_c^*,
  \qquad
  \overline QK_c=0,
  \qquad
  \overline BK_c=\overline R.
  ```

- Diagnose compatibility before constructing a correction. Let $P_Q$ be the orthogonal projector onto $\ker\overline Q$. Because skew-Hermiticity and $\overline QK_c=0$ imply $K_c=P_QK_cP_Q$, report the defects in the necessary APV–bottom compatibility conditions, including

  ```math
  \overline R(I-P_Q)=0
  ```

  and skew-Hermiticity of $\overline R\,\overline B^*$.
- If the constraints are feasible, compute the unique minimum-change closure

  ```math
  \underset{K_c}{\operatorname{minimize}}
  \quad
  \lVert K_c-K_0\rVert_F
  ```

  subject to the three exact constraints. Implement this first as a low-resolution dense oracle with explicit real and imaginary linear constraints.
- Transform a feasible closure back to the original coordinates,

  ```math
  J_\gamma^c=S^*K_cS,
  \qquad
  E_\gamma\dot{\boldsymbol a}=J_\gamma^c\boldsymbol a,
  ```

  while retaining $J_\gamma$ as the unmodified Galerkin reference.
- Report the relative correction $\lVert K_c-K_0\rVert_F/\lVert K_0\rVert_F$, constraint residuals, ranks and nullities, balanced dimension, and dependence on terrain amplitude and horizontal and vertical resolution.
- Test the constraints on the entire then-current state space, including the APV-bearing balanced and displacement-only bottom coordinates. Do not obtain feasibility by silently deleting balanced columns or restricting the audit to a selected wave subspace.

### Automated acceptance and exits

- The audit uses flat, uniform-depth, and sinusoidal-terrain cases and is deterministic under repeated construction.
- For $h=0$, the raw generator is already feasible and the minimum-change closure satisfies $K_c=K_0$ to $10^{-12}$ relative accuracy.
- A feasible closure must satisfy

  ```math
  \frac{\lVert K_c+K_c^*\rVert}
  {\max(1,\lVert K_c\rVert)}
  \leq10^{-13},
  ```

  ```math
  \frac{\lVert\overline QK_c\rVert}
  {\max(1,\lVert K_c\rVert)}
  \leq10^{-12},
  \qquad
  \frac{\lVert\overline BK_c-\overline R\rVert}
  {\max(1,\lVert\overline R\rVert)}
  \leq10^{-12}.
  ```

- Random-state tests verify zero finite-terrain energy tendency, constant discrete APV, and exact discrete bottom evolution.
- Every nonzero-frequency eigenvector of a feasible $K_c$ has negligible APV, its frequencies are real to the eigensolver tolerance, and the dimension of its zero-frequency space is reported against the expected balanced count.
- The correction norm and each raw compatibility defect are reported over terrain-amplitude and resolution sweeps. No acceptance threshold is imposed on the correction norm until its convergence behavior is known.
- The side quest has two scientifically valid exits:
  - **Feasible:** the exact constraints close at the stated tolerances and the correction is reproducible. This was the audit's original criterion for proceeding with $J_\gamma^c$; the closure path is now superseded by Milestone 4.7.
  - **Incompatible:** the rank-revealing solve establishes a nonzero minimum constraint residual, identifies the failed compatibility condition and implicated state subspace, and Milestone 5 remains blocked pending a revised state or constraint formulation.
- Approximate feasibility obtained by relaxed tolerances, empirical correction factors, or removal of physical state coordinates does not pass this gate.

### Outcome

The dense audit preserves the entire then-current mixed state and constructs the resolved bottom row by projecting the oversampled physical product into the retained bottom Fourier space. Flat and uniform-depth cases are admissible: the APV nullity equals the 56-dimensional flat wave space in the automated reference problem, all three constraints close at roundoff, and the minimum-change closure is the raw generator.

The 20 m sinusoidal-terrain reference is incompatible. Under the documented numerical-rank criterion, its APV map has rank 65 and nullity 35, leaving 21 fewer resolved zero-APV directions than the flat wave count. Its minimum bottom-constraint residual is approximately $2.46\times10^{-5}$, compared with the $10^{-12}$ gate; both the APV–bottom right compatibility and skew compatibility conditions fail materially. The audit therefore returns no constrained terrain generator. This establishes the incompatible exit without deleting balanced coordinates or relaxing the APV rank. The original pointwise-APV path stopped here; Milestone 4.6 records the subsequent quadratic-invariant audit.

## Milestone 4.6: Quadratic energy–enstrophy closure audit

- [x] Complete — historical incompatible exit; superseded by Milestone 4.7

### Purpose

Determine whether the then-current mixed state space admits a minimum-change evolution that preserves the finite-terrain energy and total discrete quadratic potential enstrophy while enforcing the resolved strong bottom-displacement evolution. This replaces the pointwise APV constraint that Milestone 4.5 proved incompatible; it does not alter the energy or potential-enstrophy norms.

The quadrature-weighted APV map already stored by the dense oracle defines

```math
Z_\gamma=Q_\gamma^*Q_\gamma,
\qquad
\mathcal Z_\gamma^{(2)}
=\frac12\boldsymbol a^*Z_\gamma\boldsymbol a.
```

With

```math
E_\gamma=S^*S,
\qquad
\boldsymbol b=S\boldsymbol a,
\qquad
K_0=S^{-*}J_\gamma S^{-1},
```

define the potential-enstrophy metric in energy coordinates by

```math
G_\gamma=S^{-*}Z_\gamma S^{-1}.
```

An energy-skew generator $K_c$ conserves both quadratic invariants for every state precisely when

```math
K_c=-K_c^*,
\qquad
[G_\gamma,K_c]=0.
```

Indeed,

```math
\frac{d}{dt}\frac12\boldsymbol b^*\boldsymbol b
=\frac12\boldsymbol b^*(K_c^*+K_c)\boldsymbol b=0,
\qquad
\frac{d}{dt}\frac12\boldsymbol b^*G_\gamma\boldsymbol b
=\frac12\boldsymbol b^*(K_c^*G_\gamma+G_\gamma K_c)\boldsymbol b=0.
```

The resolved strong bottom equation remains

```math
\overline B K_c=\overline R,
\qquad
\overline B=BS^{-1},
\qquad
\overline R=R_hS^{-1}.
```

### Dependencies

Milestone 4 and the completed incompatibility diagnosis from Milestone 4.5.

### Deliverables

- Preserve all raw Milestone-4 forms and Milestone-4.5 audit results unchanged.
- Form $Z_\gamma$ and $G_\gamma$ from the existing quadrature-weighted APV map. Verify directly that $\boldsymbol a^*Z_\gamma\boldsymbol a/2$ equals the discrete volume integral of $|q_\gamma|^2/2$.
- Diagonalize the Hermitian matrix $G_\gamma$ to expose its exact commutant: in this basis, an entry of $K_c$ may be nonzero only between equal eigenvalues. Construct the numerical nullspace of the commutator with a rank criterion derived from the eigendecomposition's backward error, and test the resulting operator against the original unmodified $G_\gamma$.
- Treat eigenvalues as numerically indistinguishable only when their backward-error intervals overlap at machine precision. Record every multiplicity, spectral separation, numerical-nullspace decision, and original commutator residual; do not introduce a tunable physical tolerance that weakens exact potential-enstrophy conservation.
- Impose the bottom equation and real-field conjugacy on the remaining free skew-Hermitian block coefficients with rank-revealing dense linear algebra.
- If feasible, compute the unique minimum-change closure

  ```math
  \underset{K_c}{\operatorname{minimize}}
  \quad
  \lVert K_c-K_0\rVert_F
  ```

  subject to skew-Hermiticity, commutation with $G_\gamma$, the bottom equation, and real-field conjugacy.
- Transform a feasible closure back to the original coordinates,

  ```math
  J_\gamma^c=S^*K_cS,
  \qquad
  E_\gamma\dot{\boldsymbol a}=J_\gamma^c\boldsymbol a,
  ```

  while retaining $J_\gamma$ and $K_0$ as the unmodified Galerkin references.
- Store the raw and constrained generators, their difference, numerical ranks, nullities, eigenspace metadata, invariant defects, bottom residuals, conjugacy defects, and relative correction norm.
- Run terrain-amplitude and horizontal-, vertical-, and oversampling-resolution sweeps. Report separately the operator correction, bottom residual, invariant residuals, resolved APV tendency, and modal APV content.
- Do not require $Q_\gamma K_c=0$. Pointwise APV conservation and negligible APV for every oscillatory eigenvector remain diagnostics rather than defining constraints.

### Automated acceptance and exits

- Test flat, uniform-depth, and sinusoidal-terrain cases, including repeated construction for determinism.
- For flat and uniform-depth cases, require the raw generator to satisfy the quadratic constraints and $K_c=K_0$ within $10^{-12}$ relative accuracy.
- A feasible closure must satisfy

  ```math
  \frac{\lVert K_c+K_c^*\rVert_F}
  {\max(\operatorname{realmin},\lVert K_c\rVert_F)}
  \leq10^{-13},
  ```

  ```math
  \frac{\lVert G_\gamma K_c-K_cG_\gamma\rVert_F}
  {\max(\operatorname{realmin},\lVert G_\gamma\rVert_F\lVert K_c\rVert_F)}
  \leq10^{-12},
  ```

  ```math
  \frac{\lVert\overline B K_c-\overline R\rVert_F}
  {\max(\operatorname{realmin},\lVert\overline B\rVert_F\lVert K_c\rVert_F+\lVert\overline R\rVert_F)}
  \leq10^{-12}.
  ```

- Random complex and conjugate-symmetric real states must have zero instantaneous finite-terrain energy and quadratic potential-enstrophy tendencies at the stated tolerances.
- The constrained generator must preserve real-field conjugacy, have real frequencies to eigensolver tolerance, and satisfy the resolved bottom relation for random states and eigenvectors.
- Report $Q_\gamma K_c$ and the APV content of every eigenvector without requiring either to vanish. Verify instead that APV redistribution leaves the total quadratic potential enstrophy unchanged.
- The audit has three scientifically valid exits:
  - **Feasible and convergent:** all exact constraints close, the correction tends to zero with terrain amplitude and decreases or stabilizes under resolution refinement, and the corrected dynamics remain close to the raw oracle. This was the audit's original criterion for proceeding with $J_\gamma^c$; that closure path is now superseded by Milestone 4.7.
  - **Feasible but dynamically poor:** the algebraic constraints close, but the correction remains order one, fails to converge, or produces modal structure inconsistent with the raw oracle. Retain the result as a diagnostic and keep Milestone 5 blocked.
  - **Incompatible:** the rank-revealing solve establishes a nonzero minimum bottom or conjugacy residual within the exact commutant of $G_\gamma$. Keep Milestone 5 blocked and next investigate revised boundary forms or state coordinates.
- Approximate feasibility obtained by merging spectrally distinct eigenvalues, relaxing invariant tolerances, deleting physical coordinates, or changing the energy or potential-enstrophy norm does not pass this gate.

### Outcome

The audit uses unitary maps from independent real physical coordinates to the full signed Fourier state. In these coordinates the energy and quadratic potential-enstrophy forms are real symmetric, conjugacy is structural, and the closure problem separates into independent real skew-symmetric solves within the numerically resolved equal-eigenvalue subspaces of $G_\gamma$. Flat and uniform-depth reference problems satisfy all constraints with the raw generator and require zero correction.

The sinusoidal-terrain problem is incompatible. At the automated $[4,4,5]$ resolution with horizontal oversampling factor two, a 20 m sinusoid has a minimum bottom residual of approximately $2.46\times10^{-5}$, equal to $88.5\%$ of the resolved bottom target; 23 of its 44 enstrophy eigenspaces fail the bottom constraint. The raw generator's relative enstrophy-commutator defect is approximately $2.70\times10^{-3}$, so the failure cannot be hidden by the absolute scaling of the dimensional enstrophy metric.

The minimum residual scales linearly with terrain amplitude from 1 m through 20 m. The incompatible exit persists for horizontal resolutions $[4,4,5]$ and $[6,4,5]$, vertical resolution $[4,4,7]$, and horizontal oversampling factors one through three; the residual remains approximately $81\%$--$89\%$ of the bottom target. No quadratic-invariant-preserving terrain generator is returned. The closure route stopped here; Milestone 4.7 implements the resulting decision to revise the state coordinates.

## Milestone 4.7: Balanced bottom-inversion basis

- [x] Complete — blocking basis reset passed

### Purpose

Replace the displacement-only bottom function from Milestone 2 with the complete flat balanced state associated with a prescribed bottom displacement. This repairs the state-space omission identified by the boundary energy–enstrophy analysis before any further finite-terrain closure or eigensolve is attempted.

### Dependencies

The completed Milestones 4.5 and 4.6 incompatibility diagnoses, together with the mathematical state-space analysis in `boundary-energy-enstrophy.tex`.

### Deliverables

- Preserve the public meaning of `bottomDisplacement`: its coefficient remains the bottom value $\hat\eta_b$.
- For every retained $\kappa>0$, solve the balanced zero-APV bottom-inversion problem

  ```math
  \partial_{\xi\xi}\eta_B
  -
  \frac{\kappa^2N^2(\xi)}{f^2}\eta_B
  =0,
  \qquad
  \eta_B(-D)=1,
  \qquad
  \eta_B(0)=0.
  ```

- Recover the dynamically linked balanced fields from

  ```math
  \psi_B
  =
  -\frac{f}{\kappa^2}\partial_\xi\eta_B,
  \qquad
  u_B=-i\ell\psi_B,
  \qquad
  v_B=ik\psi_B,
  \qquad
  w_B=0,
  \qquad
  p_B=\rho_0f\psi_B.
  ```

  Verify directly that

  ```math
  q_B
  =
  -\kappa^2\psi_B
  -
  f\partial_\xi\eta_B
  =0.
  ```

- Cache the bottom-inversion solution once per unique $\kappa$ and reconstruct its horizontal phase through the existing full-signed Fourier layout.
- For $\kappa=0$, retain the mean-density-anomaly coordinate with $\eta_B=-\xi/D$ and classify it explicitly as the horizontally uniform balanced state rather than a zero-APV bottom inversion.
- Replace the use of `bottomFunction` inside `constructBasisBlocks` with a private bottom-inversion builder that returns $\eta_B$, $\partial_\xi\eta_B$, $\psi_B$, and the linked horizontal velocity factors. Keep the public constructor unchanged.
- Replace only the scientific content of the existing bottom coordinate. Preserve the public state ordering, packing, unpacking, conjugacy maps, and coefficient normalization.
- Extend reconstruction and adjoint projection so the bottom coefficient contributes all linked fields, not only $\hat\eta$.
- Rebuild the flat dense forms and eigensystem in the revised basis. Identify the ordinary interior geostrophic coordinates, the zero-APV bottom-inversion coordinate, the mean-density anomaly, and the nonzero-frequency wave subspace separately.
- Retain the displacement-only basis and both closure-audit implementations as historical regression diagnostics, but do not use them to construct the revised generator.

### Automated acceptance

- The bottom-inversion function has $\eta_B(-D)=1$ and $\eta_B(0)=0$ to $10^{-12}$ and converges under vertical refinement.
- For constant $N^2$, the numerical bottom-inversion solution and its derivative agree with the analytic hyperbolic-function solution to $10^{-11}$ relative error.
- For arbitrary stationary $N^2$, the differential-equation residual and boundary residuals decrease with vertical resolution.
- The reconstructed bottom-inversion state satisfies geostrophic and hydrostatic balance, has $w_B=0$, and has relative APV norm below $10^{-11}$.
- Its flat exchange column is zero to $10^{-12}$ relative accuracy, while its energy is finite and strictly positive.
- Pack/unpack, reconstruction/projection adjointness, real-field conjugacy, mapped continuity, and bottom-value normalization retain the Milestone-2 tolerances.
- The revised flat eigensystem recovers the WaveVortexModel nonhydrostatic wave frequencies and polarizations at the Milestone-3 tolerances and has the analytically expected balanced dimension.
- Failure of the balanced inversion, flat eigensystem, or conditioning checks blocks Milestone 4.8.

### Result

The displacement-only coordinate was replaced by the complete zero-APV balanced inversion computed with `IMSurfaceGeostrophicModes` and `IMSolverSpectral` from the isolated `InternalModesEVP` checkout. The public coefficient remains the bottom displacement, while reconstruction now supplies its linked streamfunction, horizontal velocity, hydrostatic pressure, and displacement profiles.

Native WaveVortexModel vertical quadrature does not accurately integrate the thin bottom mode at low resolution. The flat oracle therefore assembles its energy row and column from the exact geostrophic boundary-inversion identity. With this correction, the flat energy is positive, the bottom state has exactly zero exchange and APV columns, the ordinary nonhydrostatic waves are unchanged, and the balanced nullspace gains the required bottom-inversion direction for every retained nonzero horizontal coefficient.

All Milestone-4.7 analytic, variable-stratification, balance, conjugacy, adjointness, conditioning, APV, and flat-eigensystem tests pass. The finite-terrain bottom-mode quadrature and compatibility question remain deliberately unresolved until Milestone 4.8.

## Milestone 4.7.1: Common flat energy–enstrophy–dynamical basis

- [x] Complete — common physical basis passed

### Purpose

Resolve the degeneracy of the flat zero-frequency eigenspace without changing the physical state or its invariants. The dynamical eigensolve separates waves from stationary states, while physical potential enstrophy separates the APV-bearing geostrophic modes from the zero-APV bottom inversion.

### Dependencies

Milestone 4.7 and the common-mode analysis in `terrain-energy-galerkin.tex`.

### Deliverables

- Form the physical flat potential-enstrophy matrix from the existing APV map and vertical quadrature:

  ```math
  Z_0=Q_0^*W_\xi Q_0.
  ```

- Retain the flat dynamical eigensolve

  ```math
  iJ_0\boldsymbol c=\omega E_0\boldsymbol c
  ```

  and identify its stationary subspace with the existing scale-aware frequency tolerance.
- Within that subspace, solve

  ```math
  C_0^*Z_0C_0\boldsymbol d
  =
  \lambda C_0^*E_0C_0\boldsymbol d.
  ```

- Rotate only the stationary vectors. Preserve the nonzero-frequency wave frequencies and vectors, apart from immaterial phase or degenerate-subspace choices.
- For every retained $\kappa>0$, identify one stationary zero-APV bottom inversion and the APV-bearing stationary complement.
- Construct the physical-energy projector

  ```math
  P_B=\boldsymbol c_B(\boldsymbol c_B^*E_0).
  ```

- Treat $\kappa=0$ as the separately constrained MDA sector and do not create an independent bottom projector there.
- Preserve the public constructor, state ordering, packing, unpacking, and coefficient normalization. Add only diagnostic flat-block fields for $Z_0$, the enstrophy eigenvalues, mode indices, structural residuals, and the bottom projector.
- Use the physical energy and physical volume potential enstrophy. Do not introduce a generalized bottom energy or a separately conserved boundary invariant.

### Automated acceptance

- $Z_0$ is Hermitian positive semidefinite to roundoff, and the complete flat vectors remain $E_0$-orthonormal.
- The stationary block of the transformed potential-enstrophy matrix is diagonal to $10^{-12}$ relative error.
- Every $\kappa>0$ block contains exactly one stationary zero-APV mode. It agrees with the complete bottom inversion up to phase and normalization.
- All other stationary modes have positive potential enstrophy and are $E_0$-orthogonal to the bottom inversion.
- The projector satisfies

  ```math
  P_B^2=P_B,
  \qquad
  P_B^*E_0=E_0P_B,
  ```

  retains the bottom mode, and annihilates the APV-bearing stationary modes to $10^{-11}$.
- Constant and variable stratification pass with antialiasing enabled and disabled.
- Analytic wave dispersion, wave polarization, zero wave APV, and zero wave bottom coefficient retain the Milestone-3 and Milestone-4.7 tolerances.
- Failure of this common flat basis blocks Milestone 4.8.

### Result

The flat oracle now forms $Z_0$ from the existing APV map and vertical quadrature, then diagonalizes it only inside the energy-orthonormal zero-frequency block. Every retained nonzero horizontal coefficient has exactly one stationary zero-APV vector, which agrees with the complete bottom-inversion coordinate up to phase and normalization; the remaining stationary vectors have positive potential enstrophy.

Each flat block stores the common classification and the physical-energy projector onto the bottom inversion. The projector is idempotent and energy self-adjoint at roundoff, retains the bottom mode, and annihilates the APV-bearing stationary complement. The $\kappa=0$ block retains its separately constrained MDA classification and has no artificial bottom projector. Constant and variable stratification, both antialias settings, the analytic wave oracle, the complete repository suite, and static analysis all pass.

## Milestone 4.8: Boundary-complete finite-terrain compatibility gate

- [x] Complete — incompatible scientific exit; Milestone 5 blocked

### Purpose

Determine whether the unmodified finite-terrain weak forms become mutually compatible once the bottom coefficient carries its complete balanced state. This gate tests the physical discretization itself; it does not construct a minimum-change closure.

### Dependencies

Milestone 4.7.1.

### Deliverables

- Rebuild the raw dense $E_\gamma$, $J_\gamma$, $Q_\gamma$, bottom extractor $B$, and strong bottom row $R_h$ in the revised basis using the existing common oversampled quadrature.
- Reuse `buildFiniteTerrainForms` and `buildBottomConstraintMaps` after updating them to consume the revised basis blocks; do not fork a second finite-terrain quadrature path.
- Define the raw generator

  ```math
  L_\gamma=E_\gamma^{-1}J_\gamma
  ```

  and retain the raw Hermitian and skew-Hermitian defects before any roundoff-level restoration.
- Test the three compatibility identities directly:

  ```math
  Q_\gamma L_\gamma=0,
  \qquad
  BL_\gamma=R_h,
  ```

  ```math
  L_\gamma^*Z_\gamma
  +
  Z_\gamma L_\gamma
  =0,
  \qquad
  Z_\gamma=Q_\gamma^*Q_\gamma.
  ```

- Test the corresponding random-state tendencies of finite-terrain energy and quadratic potential enstrophy, the resolved bottom evolution, and real-field conjugacy.
- Use flat terrain and uniform depth as exact reference cases. Use weak sinusoidal terrain for the first coupled compatibility test, followed by terrain-amplitude, horizontal-resolution, vertical-resolution, and oversampling sweeps.
- Diagnose the identities by state subspace and by equation row so any remaining defect can be assigned to interior waves, ordinary geostrophic modes, bottom inversion, mean density, quadrature, or truncation.
- Preserve Milestones 4.5 and 4.6 as evidence for the rejected displacement-only state. Do not apply either minimum-change closure to the revised basis unless a later, separately approved roadmap establishes a new reason to do so.

### Automated acceptance and exits

- For flat and uniform-depth cases, all three compatibility identities close within $10^{-12}$ relative error and the unmodified raw generator reproduces the corresponding exact linear dynamics.
- For the sinusoidal reference terrain, require

  ```math
  \frac{\lVert Q_\gamma L_\gamma\rVert_F}
  {\max(\operatorname{realmin},\lVert Q_\gamma\rVert_F\lVert L_\gamma\rVert_F)}
  \leq10^{-10},
  ```

  ```math
  \frac{\lVert BL_\gamma-R_h\rVert_F}
  {\max(\operatorname{realmin},\lVert B\rVert_F\lVert L_\gamma\rVert_F+\lVert R_h\rVert_F)}
  \leq10^{-10},
  ```

  ```math
  \frac{\lVert L_\gamma^*Z_\gamma+Z_\gamma L_\gamma\rVert_F}
  {\max(\operatorname{realmin},\lVert Z_\gamma\rVert_F\lVert L_\gamma\rVert_F)}
  \leq10^{-10}.
  ```

- Raw $E_\gamma$ and $J_\gamma$ retain their structural tolerances, $E_\gamma$ remains positive definite, and random-state energy tendencies close to roundoff.
- Bottom-evolution, APV, and enstrophy defects decrease under independent horizontal, vertical, and oversampling refinement rather than approaching a nonzero plateau.
- The gate has two scientifically valid exits:
  - **Compatible and convergent:** the raw forms pass the stated identities and refinement tests. Milestone 5 proceeds with the unmodified $E_\gamma$, $J_\gamma$, and $Q_\gamma$.
  - **Incompatible:** a residual remains above tolerance or converges to a nonzero limit. Stop before Milestone 5 and revise the boundary-complete state, weak forms, or adjoint-consistent discretization. Do not relax the invariants or substitute a constrained surrogate.

### Result

The read-only `auditFiniteTerrainCompatibility` method now evaluates the unmodified generator and retains complete residual matrices, deterministic random-state tests, APV rank, flat-common-subspace diagnostics, equation-row diagnostics, and the pre-restoration forms. It does not call either historical closure audit and does not alter the generator.

Flat terrain and uniform depth pass the energy, APV, bottom, quadratic-enstrophy, and conjugacy identities within $10^{-12}$. Direct common quadrature of the boundary-intensified bottom inversion converges toward the exact flat bottom rows under vertical refinement.

For the 20 m sinusoidal reference terrain at resolution $[4,4,5]$ and horizontal oversampling factor two, the normalized operational defects are:

```math
\frac{\lVert Q_\gamma L_\gamma\rVert_F}
{\lVert Q_\gamma\rVert_F\lVert L_\gamma\rVert_F}
=1.89\times10^{-7},
```

```math
\frac{\lVert BL_\gamma-R_h\rVert_F}
{\lVert B\rVert_F\lVert L_\gamma\rVert_F+\lVert R_h\rVert_F}
=2.99\times10^{-2},
```

```math
\frac{\lVert L_\gamma^*Z_\gamma+Z_\gamma L_\gamma\rVert_F}
{\lVert Z_\gamma\rVert_F\lVert L_\gamma\rVert_F}
=4.78\times10^{-8}.
```

The energy and conjugacy defects are $1.87\times10^{-23}$ and $4.07\times10^{-16}$. The raw form defects remain near roundoff and the scaled energy reciprocal condition number is $6.09\times10^{-2}$, so the failure is not caused by loss of Hermitian structure or an ill-conditioned energy form.

Oversampling factors two and three give the same incompatible residuals. Across the two finest horizontal resolutions, the APV, bottom, and enstrophy defects change by approximately $9.6\%$, $12.5\%$, and $6.2\%$. Vertical refinement reduces the residuals but does not overcome the independent horizontal and oversampling plateaus. At the reference discretization, the energy-scaled APV map has numerical nullity 34 while the raw generator has rank 76, which is already incompatible with mapping the complete generator range into the APV nullspace. Right-action diagnostics show nonzero defects from the wave, APV-bearing geostrophic, and bottom-inversion input subspaces; the mean-density sector has negligible direct APV and bottom defects.

Milestone 4.8 therefore takes its incompatible exit. Milestone 5 must not proceed from these raw forms. No corrected closure, empirical symmetrization, or relaxed invariant has been introduced.

## Milestone 5: Explicit boundary-dynamical descriptor system

- [x] Complete — blocking gate passed

### Purpose

Test whether the Milestone-4.8 incompatibility is resolved when bottom displacement is an independent dynamical coordinate and the bottom equation participates in the constrained evolution problem.

### Dependencies

Completed Milestones 4.7–4.8, including the common flat energy–enstrophy basis and the incompatible volume-only compatibility audit.

### Deliverables

- Construct the finite-dimensional state

```math
\boldsymbol x
=
\begin{pmatrix}
\boldsymbol a\\
\boldsymbol b
\end{pmatrix},
```

where $\boldsymbol a$ contains the existing interior coordinates and $\boldsymbol b$ contains independent bottom-displacement coefficients.
- Preserve the public coefficient ordering. Use the complete balanced bottom inversion only to reconstruct the fields associated with $\boldsymbol b$; do not identify the bottom coordinate with a volume equivalence class.
- Assemble one mixed descriptor problem containing the volume momentum and displacement equations, incompressibility with pressure retained as a Lagrange multiplier, the strong bottom row

```math
B\dot{\boldsymbol x}=R_h\boldsymbol x,
```

and independent volume and bottom test equations.
- Retain pressure and all algebraic constraints until the complete augmented saddle-point system has been assembled.
- Record the descriptor rank, constraint rank, pressure nullspace, stationary nullspace, and number of independent prognostic coordinates.
- Do not eliminate pressure or project onto the divergence-free volume space before the bottom row has been included.

### Automated acceptance

- Flat and uniform-depth systems reproduce the existing nonhydrostatic waves, stationary geostrophic modes, complete bottom inversion, and mean-density-anomaly states.
- Continuity, pressure gauge, surface condition, and bottom evolution close within $10^{-12}$.
- The bottom coordinate remains independent of the volume equivalence class under packing, reconstruction, and constraint elimination.
- The descriptor rank and nullspaces agree with the expected prognostic, diagnostic, and stationary dimensions.
- Failure of any rank, independence, or residual check blocks Milestone 6.

### Outcome

The implementation uses the private primitive ordering

```math
\boldsymbol y
=
(\boldsymbol u_F,\boldsymbol v_F,\boldsymbol w_G,
\boldsymbol\eta_G,\eta_b,\boldsymbol p_F)^T.
```

Pressure has no time derivative. The bottom equation supplies the endpoint displacement row, continuity supplies the algebraic constraint rows, and the null horizontal-mean continuity row is replaced by the pressure gauge. Finite generalized eigenvectors are mapped back to the unchanged public coefficient layout.

At the automated `[4,4,5]` reference resolution, both flat depth and the uniform 150 m depth change have 100 finite prognostic modes and 71 infinite diagnostic modes across the retained horizontal blocks. The largest flat/uniform descriptor residual is `6.32e-14`; continuity, pressure gauge, public mapping, and bottom-coordinate mapping close between `4.77e-35` and `1.74e-15`. The largest frequency error is `5.25e-14`. All descriptor, constraint, pressure, finite-mode, and stationary-mode ranks match their expected values. Milestone 5 therefore passes.

## Milestone 6: Boundary Green identity and invariant discovery

- [x] Complete — continuum branch P; incompatible discrete exit

### Purpose

Derive rather than assume the quadratic structure associated with the active bottom coordinate, and determine which invariant formulation the later implementation must follow.

### Dependencies

Milestone 5.

### Deliverables

- Use constant stratification and a local constant-slope problem, for which horizontal Fourier coefficients remain independent.
- Start from the material bottom relation

```math
-i\omega\eta_b
=
\boldsymbol u_{H,b}\boldsymbol{\cdot}\nabla_H h
```

and derive the complete boundary Green identity before choosing an inner product.
- Seek the most general Hermitian quadratic boundary contribution permitted by the physical boundary invariants. Permit volume–boundary cross terms when a pure $|\eta_b|^2$ contribution is insufficient.
- Compare the derived form with the free-surface generalized-energy construction, Yassin's $L^2\oplus\mathbb C$ endpoint inner product, the hydrostatic or quasigeostrophic constant-slope limit, and the full nonhydrostatic constant-$N$ problem.
- Classify the result into exactly one of the following branches:
  - **H — positive generalized metric:** a positive boundary-complete invariant exists.
  - **K — signed metric:** the derived form is nondegenerate but indefinite and requires Pontryagin or Krein treatment.
  - **P — physical-energy-only:** no separate nontrivial boundary quadratic invariant exists, but the mixed descriptor preserves physical energy.
  - **Incompatible:** the mixed equations fail physical energy, pointwise APV, or bottom compatibility.
- Record the selected branch together with the relevant metric's definiteness, rank, numerical conditioning, and nullspace.

### Automated acceptance

- The boundary form follows from the Green identity or an exact boundary Casimir and never from residual fitting.
- The augmented generator satisfies the applicable skew-adjoint identity within $10^{-11}$.
- Physical energy, pointwise APV, and bottom evolution close independently.
- The flat limit contains stationary boundary modes, while finite slope activates at least one boundary or topographic mode.
- Frequencies and mode counts agree with the applicable Yassin and analytic limits.
- Branch classification is reproducible from definiteness, rank, and conditioning diagnostics.
- An incompatible result stops the roadmap before Milestone 7.

### Outcome

Freezing the mapped coefficients at `h=0` while retaining the local slope `s=grad(h)` gives

```math
w
=
\hat w-\frac{\xi}{D}\hat{\boldsymbol u}_H\boldsymbol{\cdot}\boldsymbol s.
```

The horizontal metric-pressure terms cancel pointwise against the metric contribution inside physical `w`; mapped continuity and homogeneous mapped normal velocity remove the remaining pressure integral. The resulting continuous metric is physical energy. The bottom value instead obeys

```math
\partial_t\hat\eta_b
=
\hat{\boldsymbol u}_{H,b}\boldsymbol{\cdot}\boldsymbol s,
```

so a standalone endpoint form proportional to `abs(eta_b)^2` is not conserved. The permitted volume-boundary cross term is already the first-order cross term in physical vertical kinetic energy. No separate nontrivial boundary Casimir is derived. The continuous local problem is therefore branch **P — physical-energy-only**, rather than H or K.

The discrete test uses constant `N2=2e-5 s^-2`, domain `[24,20,1.2] km`, resolution `[4,4,5]`, latitude 45 degrees, and slope `[0.01,0]`. Its first-order tangent is obtained from the unmodified descriptor about zero slope; halving the differentiation increment changes the tangent by `2.58e-9` relatively. Pressure cancellation and the strong bottom row close to `1.26e-16` and `3.77e-13`, and finite slope activates nonzero-frequency modes with bottom participation. However, the physical-energy, pointwise-APV, and quadratic-potential-enstrophy defects are respectively `1.27e-6`, `1.19e-4`, and `9.75e-5`; the maximum real growth rate is `2.73e-6 s^-1`.

The implemented primitive F–G descriptor therefore takes the **incompatible** discrete exit even though the continuous Green identity selects branch P. No empirical symmetrization, fitted closure, or altered invariant has been introduced. Batch D1 stops here, and Milestone 7 must not begin from this descriptor.

## Milestone 6.1: Branch-P primitive oracle and F–G repair audit

- [x] Complete — independent oracle establishes an APV-compatible-space blocker

### Purpose

Determine whether the Milestone-6 failure came from the physical Branch-P formulation, the hydrostatic F–G coordinates, their embedded balanced bottom inversion, or the eigenvector-based public-coordinate reconstruction.

### Dependencies

Milestones 5–6.

### Deliverables

- Add `auditBranchPDiscreteOracle` without changing the public Galerkin coefficient layout.
- Construct an independent primitive polynomial oracle using Gauss–Legendre quadrature and the vertical spaces

```math
\mathcal F_N=\operatorname{span}\{P_0,\ldots,P_N\},
\qquad
\mathcal G_N^0=(1-r^2)\operatorname{span}\{P_0,\ldots,P_{N-1}\},
```

with one eta-only bottom function

```math
\chi_b=\frac{1-r}{2}.
```

- Retain pressure through the constrained saddle solve and differentiate that saddle system analytically with respect to the local bottom slope.
- Assemble the coupled momentum–displacement rows with the physical terrain-energy weak form. In particular, include the metric-induced physical vertical velocity carried by every horizontal test velocity.
- Audit the first-order coefficient identities for weak evolution, physical energy, pointwise APV, strong bottom evolution, and quadratic potential enstrophy.
- Reduce the existing F–G descriptor directly through its saddle system, compare its native and public-coordinate forms, and do not construct its generator from generalized eigenvectors.
- Permit a repair only if the independent primitive oracle passes every physical identity.

### Automated acceptance

- Endpoint values and the polynomial Green identity close within $10^{-13}$.
- The constrained saddle residual, continuity tangency, pressure work, weak evolution, physical energy, and bottom evolution close within $10^{-11}$.
- The analytic slope derivative agrees with a centered independent calculation within $10^{-9}$.
- The flat mode-one frequency converges to the constant-$N$ nonhydrostatic dispersion relation.
- Every physical identity remains convergent for polynomial degrees 4, 6, 8, and 10 and for zonal, meridional, and oblique slopes.
- The existing public packing, ordering, and Fourier-conjugacy maps remain unchanged.
- If the independent polynomial oracle fails a mandatory physical identity without refinement, record the blocker and apply no F–G repair.

### Outcome

The eta-only polynomial coordinate is independent: its velocity and pressure fields vanish, its bottom displacement is one, and its surface displacement is zero. The polynomial endpoint residual is zero and the quadrature Green-identity defect remains below `3.1e-15`. At polynomial degree eight, the flat mode-one frequency agrees with the analytic constant-$N$ nonhydrostatic frequency to `1.27e-13`.

The first uncoupled primitive projection reproduced the Milestone-6 energy failure. Reassembling the momentum equations as the coupled finite-terrain energy weak form identified and removed that defect: for degree six and slope `[0.01,0]`, the weak-evolution, physical-energy, bottom-evolution, pressure-work, continuity-tangency, and saddle residuals are respectively

```text
6.79e-16, 8.71e-16, 1.53e-15, 3.26e-14, 9.51e-20, 1.17e-16.
```

The analytic slope derivative agrees with centered differencing to `9.82e-12`, and the finite-slope remainder has the expected second-order absolute scaling.

Pointwise APV and quadratic potential enstrophy do not pass. Their normalized first-order defects are `6.97e-1` and `4.90e-1` at degree six. Across degrees 4, 6, 8, and 10, the APV defect ranges from `6.85e-1` to `7.04e-1`, while the enstrophy defect ranges from `4.82e-1` to `4.96e-1`; neither tends toward zero. The same split occurs for zonal, meridional, and oblique slopes.

The direct F–G saddle reduction is worse but consistent with the earlier audit: its native first-order energy and APV defects are approximately `6.6e-3` and `7.8e-1`, while strong bottom evolution remains at roundoff. Mapping that reduction to the unchanged public coordinates does not remove the failure.

Milestone 6.1 therefore takes its blocking exit. The result does **not** invalidate the continuous Branch-P equations: it shows that neither the present F–G descriptor nor the independent finite polynomial trial/test pair is a pointwise-APV-compatible discretization of them. Because the independent oracle fails a mandatory identity, no F–G repair, empirical correction, or public-layout change is applied. Periodic terrain and modal construction remain blocked.

## Milestone 6.2: APV-compatible local primitive investigation

- [x] Complete — the frozen cross-slope problem has a mathematical APV obstruction

### Purpose

Determine whether compatible mixed polynomial spaces, vorticity--divergence variables, or an explicit APV coordinate can simultaneously preserve the frozen Branch-P weak equations, physical energy, stationary APV, and strong bottom evolution without an empirical correction.

### Dependencies

Milestone 6.1.

### Deliverables

- Add `auditAPVCompatiblePrimitive` without changing the public Galerkin coefficient layout.
- Use the compatible polynomial sequence

```math
u,v\in P_N,
\qquad
\hat w\in P_{N+1}^{00},
\qquad
\eta\in P_{N+1}^{0s}.
```

- Verify the rectangular product rule between the surface-zero displacement space and the horizontal-velocity space.
- Test an energy-weak mixed descriptor with one additional pressure polynomial and the complete surface-zero space as the vertical-momentum test space.
- Replace the horizontal momentum rows by vorticity--divergence and APV rows in an independent candidate.
- Transform that candidate to coordinates whose first block is the discrete APV itself.
- Derive the APV tendency of the frozen constant-slope equations analytically before interpreting a numerical residual.
- Retain the strong bottom-displacement row and prohibit empirical symmetrization, APV projection, or fitted closure.
- Stop before periodic terrain or modal construction.

### Automated acceptance

- The rectangular polynomial product rule, continuity rank, admissible dimension, and bottom endpoint values close within $10^{-11}$.
- The energy-weak candidate preserves its weak equation, physical energy, continuity, and bottom evolution within $10^{-11}$.
- The vorticity--divergence candidate preserves APV, potential enstrophy, continuity, and bottom evolution within $10^{-11}$.
- The explicit-APV coordinate map is complete and its APV tendency rows vanish within $10^{-11}$.
- The flat reference recovers the constant-$N$ nonhydrostatic mode-one frequency spectrally.
- Zonal, meridional, oblique, and sign-reversed slopes reproduce the analytic cross-slope source.
- The result is unchanged in classification under vertical-degree and quadrature refinement and remains valid for a smooth nonconstant $N^2(z)$.
- A working oracle is accepted only if one unmodified candidate preserves physical energy, pointwise APV, and strong bottom evolution simultaneously.

### Analytic result

For the frozen equations used by the Branch-P local oracle, the APV is

```math
q_{\boldsymbol s}
=
ikv-i\ell u
+\frac{s_x}{D}\left(v+\xi v_\xi\right)
-\frac{s_y}{D}\left(u+\xi u_\xi\right)
-f\eta_\xi .
```

Substitution of the frozen horizontal momentum, displacement, and mapped-continuity equations gives

```math
\boxed{
\partial_t q_{\boldsymbol s}
=
\frac{i(k s_y-\ell s_x)}{\rho_0D}\,p .
}
```

The Coriolis and stretching contributions cancel through mapped continuity. The two metric pressure-gradient terms cancel only when $k s_y-\ell s_x=0$. Thus a cross-slope Fourier component of the frozen local problem creates APV before any vertical discretization is chosen. This does not contradict stationary APV in the exact terrain equations: freezing the terrain value while retaining its slope removes horizontal coefficient variation needed by the exact APV cancellation.

### Numerical outcome

At polynomial degree six, horizontal mode `[1,0]`, and slope `[0,0.01]`, the energy-weak mixed descriptor gives the normalized defects

```text
weak evolution       4.38e-19
physical energy      6.23e-19
bottom evolution     1.50e-15
pointwise APV        1.11e-05
potential enstrophy  1.86e-06
```

The direct polynomial evaluation of the boxed analytic identity closes to `5.52e-16`. The mixed descriptor's computed APV tendency agrees with that analytic pressure source to `6.14e-3` relatively. The compatible rectangular product-rule defect is `1.01e-13`. Refinement from degrees 4 through 10 and doubled quadrature preserve the same classification.

The vorticity--divergence candidate instead gives

```text
pointwise APV        2.38e-17
potential enstrophy  1.38e-17
bottom evolution     1.45e-15
physical energy      2.31e-07
weak evolution       1.59e-05
```

The explicit-APV coordinate map has full rank and its APV rows vanish to `1.43e-20`, but it is only a coordinate representation of this same APV-enforced, non-energy-conserving candidate. It does not repair the weak equations.

When the wavenumber and slope are parallel, the analytic source vanishes. The energy-weak descriptor then preserves physical energy, APV, potential enstrophy, and bottom evolution at roundoff. This is a valid two-dimensional special oracle but not a general constant-slope solution.

Milestone 6.2 therefore records a **mathematical blocker for the frozen cross-slope Branch-P system**. No finite-space change can preserve both its unmodified weak equations and $QL=0$ when the analytic source is nonzero. The next scientific formulation must retain the globally varying small-terrain coefficients or derive a covariant local/WKB system containing the missing connection terms. Periodic terrain and modal construction remain out of scope.

## Milestone 6.3: Global small-terrain primitive compatibility oracle

- [x] Complete — the global equations remove the frozen-slope source, but the finite primitive representation is not APV-compatible

### Purpose

Determine whether retaining the globally varying first-order terrain coefficients restores the APV cancellation lost by the frozen local model before attempting finite-amplitude periodic terrain or terrain-mode construction.

### Dependencies

Milestone 6.2.

### Deliverables

- Add `auditGlobalSmallTerrainPrimitive` without modifying the public wave–vortex coefficient layout.
- Use the compatible primitive spaces

```math
\hat u,\hat v\in P_N,
\qquad
\hat w\in P_{N+1}^{00},
\qquad
\hat\eta\in P_{N+1}^{0s},
\qquad
\hat p\in P_{N+1}.
```

- Assemble one globally coupled signed-Fourier saddle system with pressure retained through continuity and the bottom equation.
- Differentiate the complete mapped system analytically in the direction

```math
h=\delta\widetilde h
```

and compare it with centered differences of the unexpanded system.
- Construct the first-order physical-energy, exchange, APV, potential-enstrophy, and bottom maps independently.
- Test

```math
E_0L_1+E_1L_0=J_1,
```

```math
L_1^*E_0+E_0L_1+L_0^*E_1+E_1L_0=0,
```

```math
Q_0L_1+Q_1L_0=0,
\qquad
BL_1=R_1,
```

and the corresponding first-order quadratic-potential-enstrophy identity.
- Separate horizontally interior inputs, for which every first-order terrain sideband is retained, from Fourier-edge inputs whose sidebands leave the state space.
- Prohibit APV projection, empirical closure, energy symmetrization, or removal of failing states.

### Automated acceptance

- Analytic and centered first-order forms agree within $10^{-9}$.
- Raw energy and exchange tangents retain Hermitian and skew-Hermitian structure within $10^{-12}$.
- Weak evolution, physical energy, bottom evolution, and Fourier conjugacy close within $10^{-11}$.
- Pointwise APV and quadratic potential enstrophy close within $10^{-10}$ for the complete retained state.
- A single sinusoidal terrain component couples only the analytically selected neighboring Fourier blocks.
- Constant-$N$ flat frequencies converge spectrally to the analytic nonhydrostatic dispersion relation.
- The classification is unchanged under polynomial, quadrature, horizontal-oversampling, and smooth-stratification checks.
- Failure of the full APV identity after independent tangent agreement is recorded as a representation blocker; no corrected generator is constructed.

### Outcome

For constant $N^2$, resolution `[6 6 5]`, terrain $\widetilde h=20\cos(2\pi y/L_y)\ {\rm m}$, polynomial degree six, and horizontal oversampling factor two, the analytic and centered first-order operators agree to `2.97e-10`. The weak-evolution, physical-energy, bottom-evolution, and conjugacy defects are

```text
3.94e-15, 2.42e-15, 2.47e-13, 5.79e-16.
```

The energy and exchange Hermitian-structure defects remain below `4.6e-16`, Fourier coupling leakage is `1.03e-13`, and the flat mode-one frequency error is `1.64e-9`.

The pointwise-APV and quadratic-potential-enstrophy defects are instead

```text
2.34e-1, 3.28e-4.
```

The APV residual has two distinguishable pieces. For horizontally interior input modes, it decreases from `4.44e-3` at polynomial degree two to `5.38e-4` at degree six. Fourier-edge inputs remain near `4.14e-1` because multiplication by the terrain creates sidebands outside the retained prognostic state. That edge value is an unprojected external-sideband diagnostic; a finite Fourier support is not required to be closed under multiplication.

Milestone 6.3 historically takes the `representation-blocker` exit because that was the classification implemented by the audit. Retaining global terrain variation removes the analytic frozen-slope objection and produces the correct Fourier selection, but the audit does not apply a common restriction to its primitive and APV products. The later interpretation is therefore narrower: the trusted interior residual must be retested under refinement, while the edge residual must be reported as discarded support. Finite-amplitude terrain remains inactive until Milestone 6.7 performs that test.

## Milestone 6.4: Coupled volume–boundary PV frequency oracle

- [x] Complete — passed local dynamical gate; primitive periodic terrain remains blocked

### Purpose

Validate the coupled volume–boundary PV interpretation and the physical frequency eigenproblem in a controlled single-wavenumber QG problem before revisiting any primitive or periodic-terrain representation.

### Dependencies

Milestones 5–6.3 and the coupled dynamical formulation in `boundary-energy-enstrophy.tex` and `terrain-energy-galerkin.tex` at mathematical checkpoint `7527984`.

### Deliverables

- Add `auditCoupledPVFrequencyOracle` for one retained nonzero horizontal mode, a prescribed volume-PV gradient, and a prescribed bottom slope.
- Use the internal state

```math
\boldsymbol{\mathcal Q}_{\boldsymbol K}
=
\begin{pmatrix}
q_{\boldsymbol K}\\
r_{b,\boldsymbol K}
\end{pmatrix},
\qquad
r_b=-f\eta_b,
\qquad
\psi_{\boldsymbol K}
=
\mathcal G_{\boldsymbol K}[q_{\boldsymbol K},r_{b,\boldsymbol K}].
```

- Compute the projected gradients

```math
\beta_{\boldsymbol K}
=K_x\overline q_y-K_y\overline q_x,
\qquad
s_{\boldsymbol K}
=f(K_xh_y-K_yh_x),
```

and solve

```math
\partial_tq_{\boldsymbol K}
+i\beta_{\boldsymbol K}\psi_{\boldsymbol K}=0,
\qquad
\partial_tr_{b,\boldsymbol K}
+is_{\boldsymbol K}\psi_{b,\boldsymbol K}=0.
```

- Construct the direct volume–boundary inversion with a Legendre--Galerkin discretization and verify its physical-energy Green identity.
- Independently construct Yassin's eigenvalue-dependent endpoint problem with a Chebyshev--Lobatto tau discretization. When both projected gradients are nonzero, compare physical frequencies using

```math
\lambda_{\mathrm Y}
=-\kappa^2-\frac{\beta_{\boldsymbol K}}{\omega}.
```

- Construct the signed volume–boundary pseudoenstrophy only when both normalizing gradients are nonzero. Report its inertia and retain its sign.
- Reconstruct the QG modes as geostrophic velocity, displacement, and pressure, then project them into the existing `A0` plus bottom-coordinate subspace without changing the public coefficient layout.
- Classify the result as `passed`, `qg-formulation-blocker`, `endpoint-equivalence-blocker`, or `basis-representation-blocker`.
- Do not reuse the blocked frozen-slope primitive generator as the reference, alter the primitive generator, or introduce empirical closure.

### Automated acceptance

- The inversion and Green-identity defects are below $10^{-12}$.
- Physical-energy skew-adjointness is below $10^{-12}$.
- Signed-pseudoenstrophy conservation is below $10^{-12}$ whenever the form is nonsingular.
- Direct and endpoint physical frequencies agree below $10^{-10}$ for resolved modes, and matched eigenfunctions or invariant subspaces agree below $10^{-9}$ in physical energy.
- Frequencies are real to $10^{-11}$.
- When $\beta_{\boldsymbol K}=0$, every nonzero-frequency mode has negligible volume APV and the active boundary mode has nonzero $r_b$.
- When both projected gradients vanish, the complete PV state is stationary.
- Reversing both gradients reverses the frequencies without changing the modal subspaces.
- Opposite gradient signs retain the indefinite pseudoenstrophy without replacing it by a positive metric.
- The $\boldsymbol K$ and $-\boldsymbol K$ blocks satisfy Fourier conjugacy.
- Projection into the existing boundary-complete balanced subspace converges under vertical refinement and retains nonzero bottom participation for the boundary mode.
- The complete repository suite and static analysis pass.

### Outcome

The independent single-wavenumber oracle passes. At polynomial degree 20, the direct Legendre volume–boundary inversion has inversion, Green-identity, physical-energy, signed-pseudoenstrophy, and frequency-imaginary defects

```text
3.04e-16, 2.17e-16, 1.64e-16, 5.29e-17, 9.76e-18.
```

The four leading resolved frequencies agree with the independent Chebyshev--Lobatto Yassin endpoint problem to `4.82e-11`; their energy-normalized eigenfunction defect is at roundoff. Lower polynomial degrees converge spectrally toward this result. The weak Legendre inversion's strong bottom-flux residual decreases from `1.14e-1` at degree 8 to `1.30e-2` at degree 20, while its Green and dynamical identities remain at roundoff.

With zero volume-PV gradient, the only nonzero-frequency mode has zero volume APV to roundoff and dominant bottom participation. With zero bottom gradient, the bottom tendency row vanishes. With both gradients zero, the complete PV generator vanishes. Same-sign gradients give a definite pseudoenstrophy, opposite signs give the expected Pontryagin metric, and reversing the signed horizontal wavenumber reverses the frequencies.

Projection of the four leading QG modes into the existing `A0` plus bottom-coordinate subspace gives resolved energy-norm defects `1.64e-1`, `4.19e-2`, and `1.70e-2` at vertical resolutions 5, 7, and 9, respectively, while reproducing the bottom value exactly. Constant and exponential stratification and both antialias settings pass. The full repository suite passes 108 tests with zero failures, and `checkcode` reports no issues in the new source.

Milestone 6.4 therefore validates the coupled volume–boundary PV interpretation, confirms that the physical frequency comes from the dynamical eigenproblem, and shows that the present boundary-complete balanced basis can represent the known local QG modes. It does not repair the Milestone-6.3 Fourier-edge APV defect.

### Stopping condition

A passing local oracle validates the boundary-dynamical interpretation and identifies the current balanced basis as an admissible representation of the known QG modes. It does not determine whether the primitive APV law closes or converges under one common dealiased projection, and it does not activate Milestone 7. Any classified blocker is recorded without modifying the governing equations or fitting a correction.

## Milestone 6.5: Periodic coupled volume–boundary PV oracle

- [x] Complete — passed projected QG closure gate; primitive periodic terrain remains blocked

### Purpose

Determine whether the Fourier-edge failure in Milestone 6.3 is intrinsic to periodic terrain or belongs to the chosen primitive state and APV range. Evolve volume APV and bottom PV as explicit coordinates and use an orthogonal Fourier projection defined before discretization.

### Dependencies

Milestones 6.3–6.4 and the existing complete signed horizontal layout. This milestone uses the verified volume–boundary inversion from Milestone 6.4 but does not use or alter the blocked primitive generator.

### Deliverables

- Add `auditPeriodicCoupledPVOracle` on every retained nonzero signed horizontal mode.
- Use

```math
\boldsymbol{\mathcal Q}_{\boldsymbol K}
=
\begin{pmatrix}
q_{\boldsymbol K}\\
r_{b,\boldsymbol K}
\end{pmatrix},
\qquad
r_b=-f\eta_b,
\qquad
\psi_{\boldsymbol K}
=\mathcal G_{\boldsymbol K}[q_{\boldsymbol K},r_{b,\boldsymbol K}].
```

- Assemble the resting-$f$-plane dynamics

```math
\partial_tq=0,
\qquad
\partial_tr_b=-\mathcal P J(\psi_b,fh),
```

where $\mathcal P$ is the orthogonal projection onto the retained signed Fourier state.
- Construct the exact, non-aliased mode-number convolution

```math
\left[\partial_tr_b\right]_{\boldsymbol K}
=
f\sum_{\boldsymbol P+\boldsymbol Q=\boldsymbol K}
(P_xQ_y-P_yQ_x)
\psi_{b,\boldsymbol P}\widehat h_{\boldsymbol Q}.
```

- Compare it with an independent oversampled pseudospectral Jacobian and reject material terrain Nyquist coefficients.
- Assemble and test physical energy, stationary volume APV, physical potential enstrophy, projected bottom evolution, and Fourier conjugacy.
- Separate interior inputs from edge inputs. Report discarded external sidebands explicitly without aliasing them into the retained state or treating them as an internal APV defect.
- Reconstruct geostrophic velocity, displacement, and pressure and project the active modes into the existing `A0` plus bottom-coordinate subspace without changing the public coefficient layout.
- Classify the result as `passed`, `inversion-blocker`, `exact-convolution-blocker`, `pseudospectral-blocker`, or `basis-bridge-blocker`.

### Automated acceptance

- The volume–boundary inversion and Green identities close below $10^{-12}$ and the physical-energy form is positive.
- The exact projected generator satisfies

```math
L_{\rm PV}^*E_{\rm PV}+E_{\rm PV}L_{\rm PV}=0,
\qquad
Q_{\rm PV}L_{\rm PV}=0,
\qquad
B_{\rm PV}L_{\rm PV}=R_h,
```

and the physical-potential-enstrophy identity below $10^{-12}$.
- Exact convolution and oversampled pseudospectral application agree below $10^{-11}$.
- The mean bottom tendency, disallowed Fourier leakage, and conjugacy defects are below $10^{-12}$.
- A sinusoidal terrain component produces only its analytically selected neighboring Fourier blocks.
- At least one edge input has a nonzero discarded sideband while all projected identities remain at roundoff.
- Physical frequencies are real to $10^{-11}$ and every nonzero-frequency mode has negligible volume APV.
- Flat terrain produces a zero generator.
- The basis bridge reproduces the bottom value exactly and converges under simultaneous Legendre-oracle and WaveVortex vertical refinement.
- Constant and exponential stratification, both antialias settings, oversampling factors two and three, the complete repository suite, and static analysis pass.

### Outcome

The periodic QG oracle passes. For constant $N^2$, resolution `[6 6 9]`, terrain $20\cos(2\pi y/L_y)\ {\rm m}$, polynomial degree 12, and horizontal oversampling factor two, the maximum inversion and Green-identity defects are

```text
3.11e-16, 2.22e-16.
```

The physical-energy and frequency-imaginary defects are `2.05e-16` and `1.26e-16`. Stationary volume APV, physical potential enstrophy, projected bottom evolution, and Fourier conjugacy close to the reported numerical zero. The exact mode-number convolution and independent pseudospectral Jacobian agree to `7.25e-16`, and the mean bottom tendency is below `4e-17` on its normalized scale.

Fourteen horizontal inputs retain every sinusoidal sideband and ten are edge inputs with at least one discarded sideband. The largest discarded coupling coefficient is `1.70e-10`, while the projected APV defect remains zero. Thus the external sideband is removed by the declared Galerkin projection rather than cyclically aliased or included in the internal conservation audit.

The active modes have zero volume APV and are represented by the existing balanced-plus-bottom coordinates with a `4.0e-3` energy-norm defect at the reference settings. Under simultaneous vertical and Legendre refinement, the resolved defects decrease from approximately `5.5e-3` to `7e-4` to `1e-4`, with the bottom value reproduced exactly. Constant and exponential stratification and both antialias settings pass. The complete repository suite passes 115 tests with zero failures, and `checkcode` reports no issues in the new source.

Milestone 6.5 establishes that periodic terrain admits a closed finite QG volume–boundary PV projection even when terrain multiplication generates Fourier coefficients outside the retained state. Its external edge sidebands are discarded by the declared projection rather than treated as internal residuals. The next experiment at that point was a hybrid primitive/PV representation with a commuting reconstruction and tendency map. This result does not activate Milestone 7.

### Stopping condition

A passing result identifies the projected coupled-PV law as the horizontal closure oracle for a future hybrid representation. No full primitive periodic descriptor, finite-amplitude terrain, or modal construction is attempted in this milestone.

## Milestone 6.6: Hybrid primitive–PV commuting oracle

- [x] Complete — primitive-equivalence blocker established

### Purpose

Determine whether the Milestone-6.5 projected volume–boundary PV law can replace only the stationary rows of the first-order primitive system while retaining the ordinary flat wave rows. This is the smallest hybrid construction that could plausibly close spectral edges without changing the public coefficient layout.

### Dependencies

Milestones 6.3–6.5. The experiment uses the independently derived primitive polynomial oracle from Milestone 6.3 and the exact projected QG closure from Milestone 6.5. It remains fixed at one nonzero zonal Fourier block and first order in a zonally invariant sinusoidal terrain direction.

### Deliverables

- Add `auditHybridPrimitivePVOracle` without altering `horizontalLayout`, `stateLayout`, or the public Galerkin coefficients.
- Let \(W\) contain the energy-normalized nonzero-frequency flat primitive modes in the selected zonal block.
- Compress the sampled flat APV range without changing it:

```math
\overline Q_0=U_q^*Q_0,
\qquad
\overline Q_1=U_q^*Q_1,
```

where the columns of \(U_q\) are an orthonormal basis for \(\operatorname{range}Q_0\).
- Combine flat wave tests, projected APV, and the independent bottom coordinate into

```math
M_0=
\begin{pmatrix}
W^*E_0\\
\overline Q_0\\
B
\end{pmatrix},
\qquad
K_0=
\begin{pmatrix}
W^*J_0\\
0\\
0
\end{pmatrix},
```

and the first-order terrain rows

```math
M_1=
\begin{pmatrix}
W^*E_1\\
\overline Q_1\\
0
\end{pmatrix},
\qquad
K_1=
\begin{pmatrix}
W^*J_1\\
0\\
R_1
\end{pmatrix}.
```

- Form

```math
L_0=M_0^{-1}K_0,
\qquad
L_1=M_0^{-1}(K_1-M_1L_0)
```

only when \(M_0\) is square, full rank, and adequately conditioned.
- Test the defining wave, projected-APV, and bottom rows separately from the unmodified primitive weak equation

```math
E_0L_1+E_1L_0=J_1.
```

- Audit physical energy, full sampled APV, quadratic potential enstrophy, Fourier conjugacy, analytic-versus-centered tangents, the public flat wave spectrum, and the Milestone-6.5 projected QG closure.
- Classify the result without empirical symmetrization, fitted closure, APV projection of the final tendency, or alteration of the governing primitive forms.

### Automated acceptance

- The hybrid descriptor is square and full rank with reciprocal condition number above \(10^{-10}\).
- It recovers the independent flat primitive generator below \(10^{-10}\).
- Analytic and centered hybrid tangents agree below \(10^{-9}\).
- Its defining wave, projected-APV, and bottom rows close below \(10^{-12}\).
- The Milestone-6.5 physical-energy, APV, enstrophy, bottom, and conjugacy identities remain below \(10^{-12}\).
- Equivalence to the unmodified primitive weak equation, physical energy, full sampled APV, and potential enstrophy must each pass below \(10^{-10}\) before the milestone may activate periodic primitive terrain.
- The conjugate zonal block agrees below \(10^{-10}\).
- Edge and interior APV defects are reported independently under vertical refinement and arbitrary stationary stratification.
- The complete repository suite and static analysis pass.

### Outcome

The hybrid coordinate construction is algebraically complete but fails the primitive-equivalence gate. For constant \(N^2\), resolution `[6 6 7]`, polynomial degree 4, and terrain direction \(20\cos(2\pi y/L_y)\ {\rm m}\), the selected zonal block contains 70 coordinates: 40 wave rows, 25 independent projected-volume-APV rows, and 5 bottom rows. The descriptor reciprocal condition number is \(1.45\times10^{-8}\), and it recovers the independent flat generator to \(5.1\times10^{-15}\).

The analytic and independently centered hybrid tangents agree to \(1.6\times10^{-11}\). The wave weak rows, projected APV identity, bottom evolution, and Fourier conjugacy close at roundoff. The independently verified periodic QG volume–boundary PV law also retains its roundoff closure.

The complete primitive identities do not close:

| Audit | Relative defect |
|---|---:|
| full primitive weak equation | \(1.26\times10^{-3}\) |
| physical energy | \(1.49\times10^{-3}\) |
| full sampled volume APV | \(2.34\times10^{-1}\) |
| quadratic potential enstrophy | \(1.04\times10^{-14}\) |

The small enstrophy defect does not rescue the construction: it is an aggregate quadratic audit and does not imply the primitive weak equation. Per-horizontal-mode APV diagnostics are at roundoff for the three interior modes, while both edge modes contain an external sideband with defect near \(4.14\times10^{-1}\). That edge value is an unprojected truncation residual, not a failure of the projected QG law. The decisive hybrid defects are the complete primitive weak equation and physical energy, which remain nonzero under vertical refinement.

The result remains classified `primitive-equivalence-blocker`. Projected APV and bottom PV are valid coordinates, and they retain the QG closure, but replacing primitive stationary test equations by those coordinate rows changes the primitive dynamics. Milestone 6.7 therefore keeps the unmodified primitive rows and changes only the shared spectral evaluation and diagnostic projection.

The focused oracle tests and complete repository suite pass 120 tests with zero failures. `checkcode` reports no issues in the repository MATLAB source.

### Stopping condition

Milestone 7 remains inactive. Do not extend this hybrid descriptor to finite-amplitude periodic primitive terrain, construct terrain modes from it, or repair its residuals empirically. Proceed only through the dealiased projected-primitive tangent oracle in Milestone 6.7.

## Milestone 6.7: Dealiased projected primitive tangent oracle

- [x] Complete — `nonconvergent-projected` result resolved by Milestone 6.8

### Purpose

Test the ordinary Fourier–Galerkin construction that the earlier full-grid APV audit did not test. Retain the unmodified primitive weak forms while evaluating every terrain product with one zero-pad, multiply, and adjoint-restrict operation. Determine whether projected APV closes exactly or converges on a trusted physical band when the outer retained modes are treated as numerical support.

### Dependencies

Milestones 6.3, 6.5, and 6.6. Reuse the verified analytic primitive tangent, the exact projected QG convolution, and the hybrid row-equivalence diagnosis. Do not reuse the hybrid replacement rows.

### Deliverables

- Retain the primitive matrices \(H_0,H_1,J_0,J_1\) and construct

```math
L_0=H_0^{-1}J_0,
\qquad
L_1=H_0^{-1}(J_1-H_1L_0).
```

- Define a retained support space, an oversampled evaluation space, and a fixed trusted physical band.
- Implement one adjoint pair \(I_{N\to M}\), \(P_{M\to N}\) and evaluate every terrain multiplication as

```math
\mathcal M_{h,N}
=
P_{M\to N}\mathcal M_hI_{N\to M}.
```

- Use the same projection in the primitive weak forms, projected APV map, and bottom map.
- Compare exact mode-number convolution with the independently oversampled pseudospectral action.
- Report separately:
  - projected residuals on the complete retained support;
  - residuals restricted to the trusted physical band; and
  - discarded external sidebands before restriction.
- Preserve the existing public coefficient layout and Fourier conjugacy maps.
- Apply no APV-nullspace projection, replacement conservation rows, corrected generator, empirical symmetrization, or fitted closure.

### Automated acceptance

- Prolongation/restriction adjointness and exact-convolution agreement are below \(10^{-12}\).
- The analytic tangent and centered differences of the identically projected finite-terrain system agree below \(10^{-9}\).
- Primitive weak evolution, physical energy, projected bottom evolution, and Fourier conjugacy close below \(10^{-11}\).
- Padding factors two and three give trusted-band actions and residuals agreeing below \(10^{-10}\).
- Trusted-band APV and potential-enstrophy defects decrease across at least three independent horizontal and vertical refinements and reach \(10^{-8}\) or better.
- Discarded external sidebands are nonzero in an edge test, are not aliased into the retained state, and are excluded from the internal APV residual.
- If

```math
Q_{0,N}L_1+Q_{1,N}L_0=0
```

closes below \(10^{-10}\), classify the oracle as `exact-projected`.
- Otherwise classify it as `convergent-projected` only if the trusted-band APV and potential-enstrophy defects decrease by at least a factor of four per refinement and satisfy the final \(10^{-8}\) tolerance.
- A nonconvergent trusted-band result is a genuine blocker and keeps Milestone 7 inactive.
- The complete repository suite and `checkcode` pass.

### Implementation result

The public oracle is:

```matlab
audit = problem.auditDealiasedProjectedPrimitiveTangent( ...
    trustedModeBounds=[1 0], ...
    supportModeBounds=[1 1;1 2;1 3], ...
    polynomialDegrees=[4;8;12], ...
    paddingFactors=[2;3]);
```

It retains the unmodified primitive weak forms and applies one common Fourier prolongation, terrain multiplication, and adjoint restriction to the primitive forms, APV map, and bottom row. The deterministic divergence-free coordinate construction makes the retained coefficient basis independent of padding.

The structural gates pass:

| Test | Maximum defect |
|---|---:|
| Prolongation/restriction and exact convolution | `1.19e-15` |
| Analytic versus centered tangent | `6.10e-10` |
| Primitive weak evolution | `9.98e-15` |
| Physical energy | `4.57e-15` |
| Projected bottom evolution | `3.15e-12` |
| Fourier conjugacy | `2.37e-15` |
| Padding factors two versus three | `4.39e-12` |

The trusted-band refinement result is:

| Vertical degree and support | APV defect | Potential-enstrophy defect |
|---|---:|---:|
| `4`, `[1 1]` | `5.45e-4` | `3.34e-4` |
| `8`, `[1 2]` | `1.85e-4` | `1.12e-4` |
| `12`, `[1 3]` | `1.01e-4` | `6.08e-5` |

The successive APV reduction factors are `2.95` and `1.84`; the potential-enstrophy factors are `2.97` and `1.85`. Increasing horizontal support alone leaves the finest-degree trusted defects unchanged, showing that the first terrain sidebands already lie in the guard band. External edge sidebands remain nonzero, with maximum discarded terrain-multiplication norm `8.33e-3`, but exact convolution and the padding comparison show that they are discarded rather than aliased into the retained equations.

The result is therefore `nonconvergent-projected`. The ordinary dealiased projection removes the earlier false Fourier-edge obstruction, but it does not restore stationary trusted-band APV for this primitive discretization at the required rate or tolerance. No replacement rows, empirical correction, symmetrization, or APV-nullspace projection were applied.

The complete repository suite passes 125 tests with zero failures, and `checkcode` reports no issues in the repository MATLAB source.

### Stopping condition

This result stopped finite-amplitude work and motivated the continuous weak-eigenproblem derivation. Do not begin Milestone 7 from the Milestone-6.7 fixed-test construction alone.

## Milestone 6.8: Boundary-complete weak terrain eigenproblem

- [x] Complete — passing tangent oracle

### Purpose

Return to the continuous primitive weak equations and derive the terrain-dependent stationary geostrophic test states that make volume APV a consequence of those same equations. Test the resulting compatible weak sequence before any finite-amplitude terrain construction.

### Dependencies

Milestone 6.7 and the continuous derivation in `finite-terrain-projection-problem.tex`. Retain the unmodified primitive \(H_0,H_1,J_0,J_1\), the complete bottom coordinate, and the common dealiased terrain projection.

### Deliverables

- Construct the exact mapped geostrophic test state generated by a scalar \(\phi\), including its \(O(h)\) changes in \(\hat u,\hat v,\hat w,\hat\eta\).
- Use scalar tests satisfying

```math
\phi_b=0,
\qquad
\partial_\xi\phi(0)=0,
```

so the APV Green identity has no boundary term.
- Hold a trusted scalar vertical space fixed while independently enriching the primitive support space.
- Construct \(G_0\) and \(G_1\), the flat and first-terrain coefficients of the geostrophic inclusion into the primitive coefficient space.
- Verify the discrete Green identities

```math
G_0^*H_0=-M_0Q_0,
```

```math
G_1^*H_0+G_0^*H_1
=
-M_1Q_0-M_0Q_1,
```

and the stationary-row identities

```math
G_0^*J_0=0,
\qquad
G_1^*J_0+G_0^*J_1=0.
```

- Evaluate the APV moments independently from the strong mapped APV diagnostic and require agreement with the energy pairing.
- Solve the tangent generalized eigenproblem without deleting modes or projecting the generator. Classify nonzero-frequency modes only after solving for \(\omega\), and test energy orthogonality, weak APV, and the first-order bottom relation.

### Automated acceptance

- Scalar bottom values and upper derivatives close below \(10^{-12}\).
- The trusted terrain-dependent geostrophic test states are represented below \(10^{-10}\).
- Flat and first-terrain Green identities and stationary rows close below \(10^{-10}\).
- Weak APV moments, independently evaluated strong APV moments, and their evolution agree below \(10^{-10}\).
- Primitive physical energy, projected bottom evolution, and Fourier conjugacy retain their Milestone-6.7 tolerances.
- Tangent eigenproblem residuals and signed-frequency orthogonality identities close below \(10^{-10}\).
- Nonzero-frequency modes have weak volume APV below \(10^{-8}\), while modes with bottom participation are retained.
- Arbitrary stationary stratification converges when the primitive vertical support is enriched around a fixed trusted scalar space.
- The complete repository suite and `checkcode` pass.

### Implementation result

The public oracle is:

```matlab
audit = problem.auditBoundaryCompleteWeakEigenproblem( ...
    trustedModeBounds=[1 0], ...
    supportModeBounds=[1 2], ...
    polynomialDegrees=[4;8;12], ...
    paddingFactor=2);
```

For constant stratification and a single meridional terrain harmonic, the result is `compatible-boundary-complete-weak-oracle`. The finest degree-twelve diagnostics are:

| Test | Defect |
|---|---:|
| Scalar boundary conditions | `1.73e-16` |
| Trusted geostrophic-state representation | `4.51e-13` |
| Flat APV Green identity | `5.65e-15` |
| First-terrain APV Green identity | `8.01e-15` |
| Flat stationary geostrophic row | `2.02e-17` |
| First-terrain stationary geostrophic row | `6.63e-14` |
| Weak APV evolution | `7.03e-14` |
| Independent strong APV evolution | `1.11e-13` |
| Weak/strong APV agreement | `7.19e-14` |
| Physical energy | `4.28e-15` |
| Projected bottom evolution | `2.47e-12` |
| Fourier conjugacy | `2.01e-15` |
| Tangent eigenproblem | `9.49e-17` |
| Nonzero-frequency weak APV | `2.58e-11` |

The apparent slow APV convergence in Milestone 6.7 came from applying a fixed flat test space to a terrain-dependent Green identity. Once the geostrophic inclusion is differentiated consistently, the APV moments follow from the same primitive weak rows at roundoff. The construction does not replace primitive rows, modify the generator, symmetrize any operator, project into an APV nullspace, or delete modes.

For variable stratification, holding the scalar test degree fixed while increasing the primitive support from degrees two to four to eight reduces the geostrophic-state representation error spectrally and restores the Green and stationary identities. This confirms that vertical support enrichment, rather than a new conservation row, supplies the required product space.

The complete repository suite passes 130 tests with zero failures, and `checkcode` reports no issues in the repository MATLAB source.

### Stopping condition

Milestone 6.8 passes, but this goal stops here. Do not begin Milestone 7 or finite-amplitude terrain until that increment is separately planned and authorized.

## Milestone 7: Finite-amplitude projected primitive dense gate

- [x] Complete — passed finite-amplitude boundary-complete primitive gate

### Purpose

Extend the validated Milestone-6.7 projection and Milestone-6.8 boundary-complete weak sequence to finite-amplitude periodic terrain.

### Dependencies

Milestone 6.8 with status `compatible-boundary-complete-weak-oracle` or a documented convergent weak-sequence classification.

### Deliverables

- Extend the unmodified projected primitive construction to periodic sinusoidal terrain.
- Use the validated support, trusted-band, padding, and adjoint-restriction conventions from Milestone 6.7.
- Construct the finite-amplitude terrain-dependent geostrophic inclusion and its APV Green map from Milestone 6.8.
- Assemble the dense reference system without replacement APV rows, empirical correction, minimum-change closure, or post hoc skew-symmetrization.
- Require the exact finite-dimensional identities

```math
L^*E_\gamma+E_\gamma L=0,
\qquad
BL=R_{h,N},
```

- Compare the APV moments derived from the energy pairing with an independently evaluated strong mapped APV diagnostic.
- Require exact closure of the finite weak APV moments or documented convergence under independent scalar, primitive, support, padding, and vertical refinement.
- Report unprojected external sidebands separately.
- Recover pressure only as a strong-equation diagnostic after the projected primitive generator has been constructed.

### Automated acceptance

- Flat and uniform-depth cases pass exact applicable identities within \(10^{-12}\).
- Sinusoidal terrain passes primitive weak consistency, physical energy, projected bottom evolution, and Fourier conjugacy within \(10^{-10}\).
- Finite-amplitude APV Green, stationary-row, and weak/strong APV agreement defects close within \(10^{-10}\), or their trusted-band residuals decrease under independent refinement and reach \(10^{-8}\).
- Trusted-band results from padding factors two and three agree within \(10^{-10}\).
- Strong momentum, continuity, and physical bottom residuals decrease under refinement.
- Failure stops the roadmap before modal construction.

### Implementation result

The public finite-amplitude oracle is:

```matlab
audit = problem.auditFiniteAmplitudeBoundaryCompleteWeakSystem( ...
    trustedModeBounds=[1 0], ...
    supportModeBounds=[1 1;1 2;1 4], ...
    scalarPolynomialDegree=2, ...
    primitivePolynomialDegrees=[2;3;5], ...
    paddingFactors=[2;3], ...
    terrainScales=[0.25;0.5;1]);
```

The retained evolution is constructed only from the finite-terrain primitive energy and exchange forms,

```math
H_\gamma\dot{\boldsymbol A}
=
J_\gamma\boldsymbol A,
\qquad
L_\gamma
=
H_\gamma^{-1}J_\gamma.
```

The finite-amplitude geostrophic inclusion is evaluated independently at the same terrain amplitude. Its energy pairing is compared with the strong mapped APV diagnostic, while the physical bottom equation is evaluated as an independent residual. Pressure and the complete strong primitive solve are diagnostics; neither replaces an energy-weak evolution row.

For constant stratification, resolution `[6 10 5]`, terrain \(20\cos(2\pi y/L_y)\ {\rm m}\), trusted bounds `[1 0]`, and the settings above, the finest diagnostics are:

| Test | Defect |
|---|---:|
| Geostrophic-state representation | `2.84e-11` |
| APV Green identity | `4.33e-15` |
| Stationary geostrophic row | `9.86e-17` |
| Weak APV evolution | `1.32e-14` |
| Independent strong APV evolution | `4.47e-15` |
| Weak/strong APV agreement | `4.43e-15` |
| Primitive weak evolution | `2.62e-16` |
| Physical energy | `1.59e-15` |
| Projected bottom evolution | `2.01e-11` |
| Fourier conjugacy | `7.92e-15` |
| Quadratic potential enstrophy | `3.98e-14` |
| Padding factors two versus three | `2.18e-12` |

The trusted bottom residual decreases from `3.47e-5` to `2.89e-7` to `2.01e-11` under joint support and primitive-degree refinement. The complete strong primitive diagnostic also converges on the trusted band. Its unprojected outer-support residual is retained separately, as are terrain-generated external sidebands.

The exact energy, exchange, and geostrophic-inclusion forms were compared with their flat-plus-tangent approximations at terrain scales `0.25`, `0.5`, and `1`. Doubling the amplitude increases each remainder by a factor within `6.3e-4` of four, confirming the expected second-order finite-amplitude contribution.

Flat and uniform-depth controls pass. With \(N^2(z)=10^{-5}[1+\exp(z/1000)]\ {\rm s}^{-2}\), joint horizontal and vertical enrichment reaches geostrophic-state representation, strong APV, and bottom defects `7.18e-11`, `3.21e-15`, and `1.45e-11`. The result is `compatible-finite-amplitude-weak-oracle`.

No replacement APV row, empirical correction, post hoc symmetrization, APV-nullspace projection, or mode deletion is used.

The complete repository suite passes 139 tests with zero failures, and `checkcode` reports no issues in all 69 MATLAB files.

### Stopping condition

Milestone 7 passes, but this goal stops here. Terrain-mode construction, Milestone 8, and time integration require a separate approved increment.

## Milestone 8: Complete finite-terrain stationary balanced space

- [x] Complete — passed stationary-space gate

### Purpose

Extend the validated zero-bottom-value APV test family into the complete finite-terrain stationary balanced space $\mathcal G_{\gamma,N}$.

### Dependencies

Milestone 7.

### Deliverables

- Extend the finite-amplitude terrain-dependent geostrophic inclusion to nonzero bottom-value scalar states.
- Combine the zero-bottom-value APV test family, the complete balanced bottom inversion, bottom-buoyancy states, and compatible mean-density-anomaly states without changing the public coefficient layout.
- Treat every bottom inversion as a linked pressure–velocity–displacement state rather than a scalar displacement lift.
- Construct the global stationary space by enforcing bottom tangency:

```math
\boldsymbol u_{H,\mathrm g,b}\boldsymbol{\cdot}\nabla_Hh=0,
```

or equivalently

```math
\partial_x\phi_b\,\partial_yh
-
\partial_y\phi_b\,\partial_xh=0.
```

- Verify the complete APV Green identity, including

```math
c_b[\psi]
=
f\hat\eta_b+v_bh_x-u_bh_y.
```

- Determine the stationary-space rank, nullity, and physical subspace content under independent horizontal, vertical, support, and padding refinement.
- Recover the flat, uniform-depth, sinusoidal-symmetry, and variable-stratification limits without introducing an empirical metric or replacement evolution row.

### Automated acceptance

- Geostrophic-state representation and the complete APV Green identity close within $10^{-10}$.
- The stationary weak rows satisfy

```math
J_\gamma G_{\gamma,N}=0
```

within $10^{-10}$.
- Bottom tangency closes within $10^{-10}$ on the trusted band.
- The rank and nullity of $\mathcal G_{\gamma,N}$ stabilize under independent horizontal, vertical, support, and padding refinement.
- Flat and uniform-depth calculations reproduce the established stationary spaces within $10^{-10}$; sinusoidal-symmetry and variable-stratification results converge under refinement.
- Failure of the stationary dimension, bottom tangency, or boundary-complete Green identity blocks Milestone 9.

### Implementation result

The public stationary-space oracle is:

```matlab
audit = problem.auditCompleteStationaryBalancedSpace( ...
    trustedModeBounds=[1 1], ...
    supportModeBounds=[1 2;1 3;1 4], ...
    primitivePolynomialDegrees=[2;3;4], ...
    paddingFactors=[2;3]);
```

The first primitive degree defines the trusted scalar vertical space, while the later degrees provide vertical support. The finite-terrain geostrophic states include nonzero bottom values, enforce the commonly projected bottom-tangency condition, and retain the complete boundary term

```math
c_b[\psi]
=
f\hat\eta_b+v_bh_x-u_bh_y.
```

For the constant-stratification sinusoidal reference calculation, the finest trusted diagnostics are:

| Test | Defect |
|---|---:|
| Geostrophic-state representation | `1.29e-17` |
| Complete APV Green identity | `1.15e-14` |
| Stationary primitive row | `1.39e-13` |
| Projected bottom tangency | `1.44e-17` |
| Strong trusted-band bottom tangency | `2.77e-17` |
| Fourier conjugacy | `6.57e-14` |
| Padding factors two versus three | `9.37e-13` |

The trusted stationary space contains 11 states at this reduced resolution. Six independent bottom-streamfunction directions fail the tangency condition and remain in the full primitive state for subsequent dynamical classification; no mode is deleted or projected out. The complete balanced bottom inversion retained by the public coefficient layout continues to satisfy its flat APV and inversion residual gates, while the stationary polynomial oracle identifies the APV-bearing, zero-volume-APV candidate, bottom-participating, and mean-density sectors.

The raw energy-scaled exchange form also has a sharply separated numerical nullspace, but it is not used to define the physical stationary space. At the finest reference resolution the first nonzero boundary-mode singular value is only about `1e-14 s^-1`. Roundoff therefore rotates a hard SVD nullspace projector by about `9.33e-5` even though the physically derived stationary rows close to `1.39e-13`. This is the expected conditioning signature of active, very-low-frequency topographic boundary modes—not evidence that those modes are stationary or should be removed. Milestone 9 must use the Milestone-8 tangency construction to identify the stationary sector before classifying the nonzero-frequency complement.

Flat and uniform-depth controls pass with no non-tangent bottom complement. The sinusoidal calculation has the expected tangent bottom-streamfunction symmetry. Variable stratification passes after vertical enrichment, and the conjugacy and common-projection checks pass with both ordinary and antialiased originating layouts.

The result is `complete-stationary-space-oracle`. The public coefficient layout is unchanged, and no replacement APV rows, empirical corrections, symmetrization, APV-nullspace projection, or mode deletion is used.

The complete repository suite passes 147 tests with zero failures.

### Stopping condition

Milestone 8 passes, but this goal stops here. Milestone 9, nonzero-frequency terrain-mode classification, and time integration require a separate approved increment.

## Milestone 9: Dense physical-energy terrain modes

- [ ] Complete — blocking modal gate

### Purpose

Establish the complete dense terrain-mode oracle using the physical finite-terrain energy and identify the dynamical modes supported by the complete stationary space.

### Dependencies

Milestone 8.

### Deliverables

- Assemble the modal problem directly from the validated Milestone-7 primitive forms:

```math
iJ_\gamma\boldsymbol c_n
=
\omega_nH_\gamma\boldsymbol c_n.
```

- Use the positive physical-energy Gram matrix $H_\gamma$ and a Cholesky-scaled Hermitian eigensolve.
- Use the independently constructed Milestone-8 $\mathcal G_{\gamma,N}$ as the stationary sector rather than identifying it with a hard frequency cutoff.
- Construct its physical-energy orthogonal complement and solve the dynamical eigenproblem there; retain every nonzero eigenpair, including arbitrarily low-frequency topographic boundary modes.
- Define $\mathcal W_{\gamma,N}$ as the span of the resulting nonzero-frequency eigenvectors.
- Classify internal waves and topographic boundary waves using frequency, polarization, volume APV, and bottom participation; retain compatible mean-density and bottom-buoyancy states inside the stationary space.
- Evaluate

```math
-i\omega_n\eta_{n,b}
=
\boldsymbol u_{H,n,b}\boldsymbol{\cdot}\nabla_H h
```

or, equivalently,

```math
-i\omega_nB_N\boldsymbol c_n
=
R_{h,N}\boldsymbol c_n,
```

as an independent modal residual rather than a replacement eigenproblem row.
- Normalize and classify modes using physical energy, volume APV, and bottom participation.
- Recover pressure only after the eigensolve and report strong momentum, continuity, surface, pressure, and bottom residuals separately.

### Automated acceptance

- Frequencies are real and generalized eigen-residuals close within $10^{-10}$.
- Distinct-frequency modes are orthogonal under physical finite-terrain energy within $10^{-10}$.
- The Milestone-8 stationary basis has eigen-residual below $10^{-10}$, and the physical-energy complement is orthogonal to it within $10^{-10}$.
- Stationary and low-frequency boundary subspaces remain stable when the numerical frequency threshold is varied; no active mode may be reclassified or deleted solely because its frequency is close to zero.
- Nonzero-frequency modes have negligible projected volume APV at the validated exact or trusted-band convergence gate.
- Flat and uniform-depth frequencies reproduce the established nonhydrostatic solutions, and mode counts remain stable under refinement.
- Bottom participation and the internal-wave versus topographic-boundary-wave classification converge under refinement.
- Strong momentum, continuity, pressure, surface, and bottom residuals decrease under independent horizontal, vertical, support, and padding refinement.
- Failure of the modal classification or convergence tests blocks Milestone 10.

## Milestone 10: Residual-enriched terrain modes

- [ ] Complete

### Purpose

Construct selected terrain modes efficiently while preserving the complete stationary space and the physical-energy structure of the dense oracle.

### Dependencies

Milestone 9.

### Deliverables

- Use flat internal, geostrophic, bottom-buoyancy, and mean-density modes as block seeds for the physical-energy eigenproblem.
- Apply residual correction to complete resonant blocks while retaining every bottom coordinate and preserving the Milestone-8 stationary subspace.
- Orthogonalize accepted corrections using $H_\gamma$.
- Compare every enriched invariant subspace with the dense terrain oracle.

### Automated acceptance

- One correction reduces weak nonresonant $O(h)$ residuals to $O(h^2)$.
- Repeated enrichment reproduces dense eigenvalues and invariant subspaces within $10^{-9}$.
- Physical-energy structure and the complete stationary subspace remain invariant after every accepted enrichment; projected APV and bottom evolution retain their validated exact or convergent behavior.
- Resonant calculations converge only when the complete coupled stationary and wave block is retained.

## Milestone 11: Matrix-free terrain-energy operators

- [ ] Complete

### Purpose

Replace the dense physical-energy oracle with adjoint-consistent operator actions while preserving the complete stationary space, APV map, and bottom maps.

### Dependencies

Milestone 10.

### Deliverables

- Replace dense $H_\gamma$, $J_\gamma$, stationary-inclusion, APV, and bottom-map matrices with field reconstruction, oversampled terrain multiplication, and adjoint projection.
- Preserve the explicit bottom coefficients throughout packing, application, and reconstruction.
- Use the flat physical-energy operator and complete flat stationary space as preconditioning data.
- Keep diagnostic pressure recovery outside ordinary operator applications.

### Automated acceptance

- Matrix-free actions reproduce the dense physical-energy oracle within $10^{-10}$.
- Physical-energy adjoint identities close within $10^{-12}$.
- No diagnostic pressure solve occurs during an operator application.
- Bottom coefficients remain explicit throughout packing, application, and reconstruction.
- Matrix-free stationary-space, projected-APV, bottom-evolution, energy, and modal diagnostics reproduce the dense results and their trusted-band classification.

## Milestone 12: Energy-preserving evolution and scientific examples

- [ ] Complete

### Purpose

Evolve the terrain-energy system and demonstrate the physical role of its stationary and topographic boundary modes.

### Dependencies

Milestone 11.

### Deliverables

- Support exact phase evolution in a converged terrain-mode basis.
- Support a physical-energy midpoint or Cayley evolution of

```math
H_\gamma\dot{\boldsymbol A}
=
J_\gamma\boldsymbol A
```

for broad states.
- Add constant-slope, sinusoidal-terrain, and Gaussian-ridge examples showing wave scattering, stationary-space participation, topographic boundary-mode excitation, bottom displacement, APV, and physical energy.
- Compare dense exponentiation, modal phase evolution, and broad-state evolution at reference resolution.

### Automated acceptance

- Dense exponentiation, phase evolution, and Cayley evolution agree at reference resolution.
- Physical energy is conserved to solver tolerance.
- Projected APV either closes exactly or converges at the validated trusted-band rate; bottom evolution closes.
- The stationary and wave subspace projections agree with the dense Milestone-9 oracle.
- Example errors converge with time step and spatial resolution.

## Milestone 13: Research-production behavior

- [ ] Complete

### Purpose

Extend the validated terrain-energy formulation to repeatable research calculations without changing its scientific definition.

### Dependencies

Milestone 12.

### Deliverables

- Add arbitrary stationary stratification, broadband terrain, resolution rebuilding, restartable output, construction and evolution profiling, and symmetry or Bloch decomposition.
- Retain the pressure-free online path after the terrain-energy operator has been constructed.
- Persist the physical-energy basis, stationary and wave classification, bottom-coordinate convention, basis metadata, and construction version.
- Benchmark construction, memory, eigenanalysis, reconstruction, and online evolution at three resolutions.

### Automated acceptance

- Repeated construction is deterministic and resolution rebuilding preserves conjugacy, invariant normalization, APV classification, and bottom participation.
- Restart continuation matches uninterrupted evolution to $10^{-10}$ in invariant-normalized coefficients.
- Broadband and variable-stratification calculations retain every applicable Milestone-7 and Milestone-9 gate.
- Benchmarks report construction time, peak stored state, operator-application time, iteration counts, and reconstruction time.
- Ordinary reduced evolution performs no diagnostic pressure solve.

## Goal-oriented batch cadence

| Batch | Milestones | Goal-sized stopping condition |
|---|---:|---|
| **D1 — Boundary formulation** | 5–6 | Build the mixed descriptor, derive the Green identity, and record branch H, K, P, or incompatible. Do not begin periodic terrain. |
| **D1.1 — Primitive repair oracle** | 6.1 | Compare an independent eta-only polynomial discretization with the F–G descriptor. Stop without repair if the primitive physical identities do not all converge. |
| **D1.2 — Global first-order oracle** | 6.2–6.3 | Diagnose the frozen-slope APV source, restore global terrain coupling, and stop if the complete finite state remains incompatible with stationary APV. |
| **D1.3 — Coupled PV frequency oracle** | 6.4 | Validate the single-wavenumber volume–boundary PV dynamics, Yassin endpoint equivalence, and current-basis representation. Do not begin periodic terrain. |
| **D1.4 — Periodic coupled-PV closure oracle** | 6.5 | Test exact projected periodic QG closure, including Fourier-edge inputs, and stop before constructing a hybrid primitive/PV descriptor. |
| **D1.5 — Hybrid primitive–PV oracle** | 6.6 | Replace stationary primitive rows by projected APV and bottom rows in one fixed-zonal block; stop if the complete descriptor is not weakly equivalent to the primitive equations. |
| **D1.6 — Projected primitive spectral oracle** | 6.7 | Establish one common padded Galerkin projection and verify exact energy and bottom identities plus trusted-band APV convergence. Do not begin finite-amplitude terrain. |
| **D1.7 — Boundary-complete weak oracle** | 6.8 | Derive and discretize the terrain-dependent geostrophic test sequence, verify APV as a primitive weak consequence, and stop before finite-amplitude terrain. |
| **D2 — Periodic scientific gate** | 7 | Extend the validated projected primitive construction to finite-amplitude sinusoidal terrain. Stop if physical energy, bottom evolution, or trusted-band APV does not meet its branch gate. |
| **D3.1 — Complete stationary-space gate** | 8 | Construct the complete finite-terrain stationary balanced space and stop if its dimension, Green identity, or bottom tangency does not converge. Do not solve the complete terrain eigensystem. |
| **D3.2 — Dense terrain-mode gate** | 9 | Solve the physical-energy eigenproblem, classify the complete stationary and wave spaces, and stop before residual enrichment. |
| **E1 — Selective modes** | 10 | Implement residual enrichment and validate it exclusively against the dense oracle. |
| **E2 — Matrix-free production core** | 11 | Replace dense actions while preserving every Milestone-7 and Milestone-9 identity. |
| **F — Evolution and examples** | 12 | Add energy-preserving evolution and the three scientific examples. |
| **G — Research production** | 13 | Add persistence, rebuilding, broadband terrain, arbitrary stratification, and profiling. |

Each batch is suitable for a Codex goal of the form: “Implement Batch D1; do not proceed to Batch D2; continue until all acceptance criteria pass or a genuine scientific blocker is established.”

Use one focused commit per completed milestone and retain acceptance evidence in the automated tests and examples. Do not bypass a blocking gate by symmetrizing an incompatible operator, fitting a closure, or introducing an empirical correction.

## Definition of done

The scientific proof of concept is established when Milestones 5–9 pass: the boundary coordinate is represented as part of a complete linked state, the projected primitive system preserves physical energy, derives stationary volume APV, and converges to the strong bottom equation, the complete stationary balanced space is represented, and the physical-energy eigenproblem produces converged internal and topographic boundary modes.

The efficient research implementation is established when Milestones 10–12 pass: residual enrichment reproduces the dense physical-energy oracle, matrix-free actions preserve its stationary-space, energy, APV, and bottom identities, and energy-preserving evolution produces convergent constant-slope, sinusoidal-terrain, and Gaussian-ridge examples.

Milestone 13 completes research-production behavior through arbitrary stationary stratification, broadband terrain, resolution rebuilding, restartable output, symmetry decomposition, and profiling.

Nonlinear terrain dynamics, an additive terrain `WVForcing`, empirical closure, an MPM release, modifications to WaveVortexModel, independent surface buoyancy, and dynamic barotropic backreaction remain outside this roadmap.
