# Terrain-energy Galerkin milestones

## Objective

Develop `WVTerrainEnergyGalerkin`, a pressure-free Galerkin system for the linear rotating Boussinesq equations over stationary bottom topography. The formulation is linear in flow amplitude and exact in the resolved terrain. The goal is to discretize a boundary-complete state for which the unmodified finite-terrain weak forms conserve discrete energy and quadratic potential enstrophy while enforcing the resolved bottom evolution.

The semidiscrete evolution is

```math
E_\gamma\dot{\boldsymbol a}=J_\gamma\boldsymbol a,
\qquad
E_\gamma=E_\gamma^*>0,
\qquad
J_\gamma=-J_\gamma^*.
```

The mathematical specification begins with `terrain-energy-galerkin.tex` at commit `72c967c` in the `ape-apv-bottom-topography` literature repository. The boundary-energy and potential-enstrophy analysis in `boundary-energy-enstrophy.tex` at commit `716c216` supplies the subsequent correction to the numerical state: a nonzero bottom displacement cannot be represented as an isolated scalar lift if its dynamically associated balanced velocity, pressure, and APV structure are omitted.

The implementation will use ordinary hydrostatic modes as economical vertical coordinates. For every nonzero horizontal wavenumber, the bottom coefficient will multiply a complete zero-APV balanced bottom-inversion state, normalized by its bottom displacement. The flat weak diagonalization will then recover the complete nonhydrostatic waves and balanced nullspace in these coordinates. Finite-terrain modes will be constructed from the unmodified finite-terrain energy and exchange forms by generalized Hermitian diagonalization and residual-based enrichment.

The completed mean-depth generator and scattering implementation is retained as reusable engineering infrastructure. Its scientific roadmap is archived in [mean-depth-wave-generator-milestones.md](mean-depth-wave-generator-milestones.md). Neither existing forcing is used as the terrain-energy evolution operator.

The dense finite-terrain forms remain the scientific oracle. Milestones 4.5 and 4.6 are retained as completed negative results: neither pointwise-APV nor quadratic-enstrophy minimum-change closure can reconcile the strong bottom equation with the displacement-only bottom coordinate. Those audits diagnose an incomplete boundary state rather than motivate a modified evolution operator. Milestones 4.7 and 4.8 therefore replace that coordinate with a complete balanced bottom-inversion state and test the raw Galerkin generator directly. No constrained surrogate advances beyond those historical audits.

## Fixed conventions and boundaries

- The public target is `WVTransformBoussinesq`; WaveVortexModel remains an external dependency and is not modified.
- The initial balanced bottom-inversion construction requires a nonzero, spatially uniform Coriolis parameter $f$. Equatorial and variable-$f$ extensions are outside the present roadmap.
- `topographicHeight` is a real, finite, stationary, upward-positive, horizontally periodic field on the transform grid.
- The mapped geometry is

```math
\gamma=1-\frac{h}{D}>0.
```

- The prognostic weak state contains $(\hat u,\hat v,\hat w,\hat\eta)$. For each retained nonzero horizontal wavenumber, one additional balanced bottom-inversion coordinate supplies a linked velocity, pressure, displacement, and APV structure and is normalized so that its coefficient equals the bottom value of $\hat\eta$. Pressure is recovered only as a post-solve diagnostic.
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

- [ ] Complete — blocking scientific gate

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

## Milestone 5: Terrain modes and wave–balanced separation

- [ ] Complete — blocking scientific gate

### Purpose

Demonstrate that the finite-terrain weak system produces physically admissible oscillatory and balanced modes before production optimization.

### Dependencies

Milestone 4.8 with its compatible-and-convergent exit.

### Deliverables

- Solve the complete unmodified dense generalized Hermitian problem

```math
iJ_\gamma\boldsymbol c_n
=
\Omega_nE_\gamma\boldsymbol c_n.
```

- Classify nonzero-frequency modes by APV norm and projection onto the flat wave subspace. Verify that converged oscillatory modes lie in the nullspace of $Q_\gamma$.
- Retain the complete zero-frequency balanced subspace, including ordinary interior geostrophic flow, the bottom-inversion coordinates, and the mean-density anomaly.
- Within the zero-frequency subspace, solve

  ```math
  Z_\gamma\boldsymbol c_n
  =
  \lambda_nE_\gamma\boldsymbol c_n
  ```

  to obtain energy–enstrophy-orthogonal vortical modes and an energy-orthogonal zero-enstrophy remainder.
- Normalize modes in the finite-terrain energy and resolve degenerate eigenspaces by $E_\gamma$-orthogonalization.
- Use sinusoidal terrain to verify leading Fourier selection and terrain-induced modal coupling.
- Verify the oscillatory bottom relation

```math
\sigma i\Omega\hat\eta_b
=
\boldsymbol u_{H,b}\boldsymbol{\cdot}\nabla_Hh.
```

- Recover a zero-mean pressure after convergence by fitting the strong momentum equations. Treat this as a diagnostic only.
- Report strong horizontal and vertical momentum, continuity, surface, and physical-bottom residuals separately.

### Automated acceptance

- Terrain frequencies are real to $10^{-10}$ relative accuracy.
- Distinct-frequency modes are $E_\gamma$-orthogonal to $10^{-10}$.
- Total quadratic potential enstrophy is conserved to $10^{-12}$ and every converged nonzero-frequency mode has negligible APV.
- The balanced dimension, bottom-inversion freedom, and mean-density coordinate agree with the discrete state count.
- The secondary energy–enstrophy diagonalization spans the complete balanced nullspace and reproduces both quadratic forms to $10^{-10}$.
- Sinusoidal terrain produces the predicted Fourier couplings.
- Bottom and strong-equation residuals decrease under horizontal and vertical refinement.
- Failure of the energy, APV, quadratic potential-enstrophy, bottom, balanced-dimension, or strong-residual gates blocks Milestone 6.

## Milestone 6: Residual-enriched terrain dressing

- [ ] Complete

### Purpose

Construct selected terrain modes efficiently from flat nonhydrostatic seeds without diagonalizing the complete global dense problem.

### Dependencies

Milestone 5.

### Deliverables

- Use the Milestone-4.7 flat nonhydrostatic wave modes and complete balanced modes as initial vectors.
- Evaluate each seed's residual in the unmodified boundary-complete finite-terrain generalized eigenproblem approved by Milestone 4.8.
- Apply block residual corrections preconditioned by the flat signed-frequency operator.
- Remove the complete degenerate or near-resonant flat eigenspace before applying the complementary inverse.
- Include the complete bottom-inversion coordinate in the same correction, energy orthogonalization, and reduced Ritz solve.
- Rediagonalize the exact reduced $(iJ_\gamma,E_\gamma)$ pair after every enrichment step.
- Preserve the discrete APV nullspace during oscillatory-mode enrichment and retain the complete balanced block when a near-zero-frequency resonance is present.
- Compare selective dressed modes with the raw boundary-complete dense terrain oracle.

### Automated acceptance

- For weak nonresonant terrain, one correction reduces an $O(h)$ seed residual to $O(h^2)$, with observed terrain-amplitude order at least 1.8.
- Degenerate tests converge only after the complete coupled block is included.
- Residual norms decrease monotonically after accepted enrichment steps.
- Repeated enrichment reproduces targeted dense eigenvalues and energy-normalized eigenspaces within $10^{-9}$.
- Dressed oscillatory modes have negligible APV and satisfy the bottom relation at the dense-oracle tolerance.
- Dressed balanced blocks reproduce the dense oracle's energy–enstrophy eigenspaces.

## Milestone 7: Matrix-free operator actions

- [ ] Complete

### Purpose

Replace global dense terrain matrices with adjoint-consistent field reconstruction, terrain multiplication, and projection.

### Dependencies

Milestone 6.

### Deliverables

- Implement matrix-free applications of the raw boundary-complete $E_\gamma$, $J_\gamma$, and $Q_\gamma$.
- Unpack interior and balanced bottom-inversion coordinates, reconstruct the mapped and physical fields, multiply by terrain and stratification weights on an oversampled grid, and project with the exact adjoints.
- Use a default horizontal oversampling factor of two and permit larger factors for convergence studies.
- Implement flat signed-frequency preconditioning and block iterative eigensolves.
- Exploit preserved meridional wavenumber and Bloch classes when the terrain provides those symmetries.
- Keep diagnostic pressure recovery outside ordinary operator applications.

### Automated acceptance

- Matrix-free actions agree with the raw boundary-complete dense oracle within $10^{-10}$ at reference resolution.
- Matrix-free energy and exchange actions satisfy their adjoint identities within $10^{-12}$.
- Matrix-free actions satisfy the APV, quadratic-enstrophy, and bottom compatibility identities at the Milestone-4.8 tolerances.
- Matrix-free Ritz values, wave and balanced eigenspaces, quadratic potential enstrophy, modal APV diagnostics, and residual histories reproduce the dense results.
- Results converge independently with Fourier resolution, vertical modes, oversampling factor, and dressing iterations.
- Runtime operator applications allocate no global dense terrain matrix and perform no pressure solve.

## Milestone 8: Energy-preserving linear evolution and examples

- [ ] Complete

### Purpose

Evolve finite-terrain states without converting the weak system into an additive flat forcing.

### Dependencies

Milestone 7.

### Deliverables

- Evolve resolved terrain modes by exact phase multiplication:

```math
A_n(t)=A_n(0)e^{-i\Omega_nt}.
```

- Evolve broader mixed states with the midpoint/Cayley step

  ```math
  \left(E_\gamma-\frac{\Delta t}{2}J_\gamma\right)\boldsymbol a^{n+1}
  =
  \left(E_\gamma+\frac{\Delta t}{2}J_\gamma\right)\boldsymbol a^n.
  ```

- Use the flat energy and signed-frequency operator as the iterative-solve preconditioner.
- Convert between Galerkin states and WaveVortexModel field and coefficient conventions at initialization and output.
- Add uniform-depth, sinusoidal-terrain, and Gaussian-ridge examples showing modal scattering, physical vertical sections, bottom displacement, APV, and finite-terrain energy.
- Report operator, eigensolver, and time-integration errors separately.

### Automated acceptance

- Terrain-mode phase evolution agrees with direct matrix exponentiation at low resolution.
- Cayley evolution conserves $\boldsymbol a^*E_\gamma\boldsymbol a/2$ to the linear-solver tolerance.
- Cayley evolution conserves $\boldsymbol a^*Z_\gamma\boldsymbol a/2$ to the linear-solver tolerance.
- The uniform-depth evolution agrees with the independent exact solution at depth $H$.
- Sinusoidal and Gaussian examples converge with time step, horizontal resolution, vertical modes, and dressing iterations.
- Quadratic potential enstrophy remains constant and resolved APV is stationary to the matrix-free compatibility tolerance.
- Ordinary time stepping performs no diagnostic pressure solve.

## Milestone 9: Research-production behavior

- [ ] Complete

### Purpose

Add the resolution, persistence, performance, and broadband capabilities required for repeatable research calculations.

### Dependencies

Milestone 8.

### Deliverables

- Support arbitrary stationary stratification and broadband periodic terrain through the already validated generalized forms.
- Reuse the deterministic Goff terrain generator from the mean-depth implementation without reusing its forcing operator.
- Rebuild geometry, basis functions, and operator state when the target transform resolution changes.
- Add restartable NetCDF output for canonical inputs, Galerkin coefficients, bottom coefficients, time, basis metadata, and construction version.
- Add selective terrain-mode and broad-state restart continuation.
- Benchmark construction, memory, eigensolve, reconstruction, and one online evolution step at three resolutions.
- Preserve the pressure-free online path and distinguish diagonal terrain-mode evolution from matrix-free Cayley evolution in performance reports.

### Automated acceptance

- Repeated construction is deterministic and resolution rebuilding preserves conjugacy, energy normalization, and APV classification.
- Restart continuation matches uninterrupted evolution to $10^{-10}$ in energy-normalized coefficients.
- Broadband and variable-stratification calculations retain the structural, energy, APV, quadratic potential-enstrophy, and strong-residual gates.
- Benchmarks report construction time, peak stored state, operator-application time, iteration counts, and reconstruction time.
- Selective terrain-mode evolution reduces to phase multiplication and reconstruction.
- Broad-state evolution reduces to matrix-free $E_\gamma$ and raw boundary-complete $J_\gamma$ actions plus preconditioned Cayley solves.

## Planning and implementation cadence

- **Batch A — Milestones 1–4.6.** These completed milestones established the software foundation, the initial dense forms, and the two incompatible closure audits. Milestones 2–4 are retained as historical prototypes because their displacement-only bottom coordinate is superseded.
- **Batch B1 — Milestones 4.7–4.7.1.** Implement the balanced bottom-inversion basis, then resolve its flat stationary degeneracy with the common physical energy–enstrophy basis. Stop if either step fails to recover the complete flat nonhydrostatic and balanced problem with acceptable conditioning.
- **Batch B2 — Milestones 4.8–5.** Rebuild the raw finite-terrain forms in the revised basis and test their mutual compatibility before solving terrain modes. Milestone 4.8 is the principal blocking gate: only a compatible-and-convergent raw system proceeds to the modal classification and balanced energy–enstrophy diagonalization in Milestone 5.
- **Batch B3 — Milestone 6.** Implement residual enrichment only after the boundary-complete dense terrain eigensystem passes. Compare every selective result with the raw dense oracle.
- **Batch C — Milestones 7–9.** Plan only after the dense terrain modes and residual-enrichment tests pass. Implement matrix-free actions before online evolution, resolution rebuilding, or persistence.
- Use one focused commit per completed milestone and retain the acceptance evidence in the automated tests and examples.
- Do not bypass a blocking gate by symmetrizing a scientifically incorrect operator or by introducing empirical correction factors.

## Definition of done

The scientific proof of concept is established when Milestones 4.7–6 pass in addition to the completed foundation: the balanced bottom-inversion basis recovers the complete flat nonhydrostatic and balanced state, the unmodified finite-terrain forms satisfy energy, APV, quadratic-enstrophy, and bottom-evolution compatibility under refinement, terrain modes satisfy the physical bottom condition, and residual dressing converges to the raw boundary-complete dense oracle.

The research implementation is complete when Milestones 7–9 pass: matrix-free actions reproduce the oracle, linear evolution preserves finite-terrain energy and quadratic potential enstrophy, the uniform-depth and scattering examples converge, and restartable broadband calculations require no diagnostic pressure solve during ordinary evolution.

Nonlinear terrain dynamics, an additive terrain `WVForcing`, an MPM release, modifications to WaveVortexModel, independent surface buoyancy, and dynamic barotropic backreaction remain outside this roadmap.
