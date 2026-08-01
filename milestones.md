# Terrain-energy Galerkin milestones

> **Checkpoint before Milestone 10.2.5:** fixed-$\kappa$ internal waves remain the verified wave coordinates, but neither ordinary flat nor trained signed-Robin geostrophic coordinates economically represent the exact finite-terrain stationary inclusion. Global $G_1$ dressing of the bottom-tangent Robin combinations reduces stationary weak-row and bottom defects, yet the best stationary-projector defect remains $1.20\times10^{-2}$ and the finite-Robin Gram condition number reaches $1.33\times10^{15}$. The next experiment will enforce the global terrain-tangent trace constraint before compressing the homogeneous vertical interior. Milestone 10.3 and matrix-free work remain inactive unless Milestone 10.2.5 returns `global-tangent-acceleration`.

## Objective

Develop an economical, complete wave–vortex model for the linear rotating Boussinesq equations over stationary bottom topography. The formulation is linear in flow amplitude and exact in the resolved terrain. The immediate goal is to define a finite production state from explicit internal-wave, APV-bearing balanced, zero-APV bottom, and MDA coordinates; construct stationary, internal-wave, and topographic-wave projectors that exhaust that state; and evolve arbitrary production states while conserving physical finite-terrain energy and preserving the validated projected APV and bottom laws. The primitive polynomial representation remains the independent construction and verification oracle rather than the intended online state.

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

The verified dense oracle uses primitive polynomial coordinates. Milestone 9.2 compares them with a boundary-complete modal family containing fixed-\(\kappa\) nonhydrostatic waves, generalized-energy Robin APV modes, and an explicit zero-APV bottom inversion. The public coefficient ordering remains fixed, the bottom coordinate remains independent, and physical finite-terrain energy remains the evolution norm. The signed Robin length shapes only the candidate basis. At finite mode count it may strongly affect convergence; it is selected on a training oracle and frozen for independent validation. Only a stable complete limit is expected to be coordinate independent.

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

- [x] Complete — dense eigensystem established; physical-subspace classification continues in Milestone 9.1

### Purpose

Establish the complete dense physical-energy eigensystem, preserve the Milestone-8 stationary space and every non-tangent bottom direction, and expose the candidate physical and unresolved algebraic subspaces without requiring every auxiliary eigenvector to pass the continuum physical tests at one resolution.

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
- Construct its physical-energy orthogonal complement and solve the dynamical eigenproblem there; retain every eigenpair, including zero-frequency completion directions and arbitrarily low-frequency topographic boundary candidates.
- Record candidate internal-wave, topographic-boundary-wave, unresolved-zero, and unresolved-nonzero groups without treating a one-resolution label as physical validation.
- Retain compatible mean-density and bottom-buoyancy states inside the stationary space.
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
- Normalize modes using physical energy and diagnose frequency, polarization, volume APV, bottom participation, and backward uncertainty.
- Recover pressure only after the eigensolve and report strong momentum, continuity, surface, pressure, and bottom residuals separately.

### Automated acceptance

- Frequencies are real and generalized eigen-residuals close within $10^{-10}$.
- Distinct-frequency modes are orthogonal under physical finite-terrain energy within $10^{-10}$.
- The Milestone-8 stationary basis has eigen-residual below $10^{-10}$, and the physical-energy complement is orthogonal to it within $10^{-10}$.
- The stationary sector is defined exclusively by the Milestone-8 tangency construction. A mode is distinguishable from zero only when its absolute frequency exceeds its backward-error uncertainty; no direction is reclassified or deleted solely because its frequency is close to zero.
- Every candidate active bottom-mode frequency is reported with its numerical uncertainty and frequency-to-uncertainty ratio.
- Projected volume APV, bottom evolution, strong primitive residuals, and bottom participation are reported mode by mode for subsequent subspace-convergence testing.
- Flat and uniform-depth frequencies reproduce the established nonhydrostatic solutions, and mode counts remain stable under refinement.
- Strong momentum, continuity, pressure, surface, and bottom residuals are available separately rather than hidden in the energy eigen-residual.
- Failure to classify every auxiliary eigenvector does not invalidate a separately converged physical subspace; the required refinement and classification gate is Milestone 9.1.

### Numerical outcome

The dense physical-energy implementation is complete, but the blocking scientific gate does not pass.

The implementation:

- builds $H_\gamma$ and $J_\gamma$ directly from the unmodified primitive energy and Coriolis--buoyancy exchange forms;
- uses a Cholesky-scaled dense eigensolve without empirical symmetrization;
- defines the stationary sector only through the Milestone-8 tangency construction at each vertical refinement;
- retains every non-tangent bottom direction in the energy-orthogonal complement;
- diagnoses pressure, APV, bottom evolution, and the strong primitive equations only after the eigensolve; and
- assigns a backward-error uncertainty to every computed frequency.

For the constant-stratification sinusoidal reference, the generalized eigen-residual, physical-energy orthogonality, Fourier-conjugacy defect, and non-tangent-bottom retention defect are below approximately `2e-13`, `5e-15`, `7e-15`, and `3e-15`. Flat and uniform-depth first-mode frequencies converge below `1e-10` relative error. These results verify the dense energy eigensolver and the reference limits.

Two unresolved populations remain:

1. The physical-energy complement contains backward-error-indistinguishable zero-frequency directions outside the declared Milestone-8 stationary sector. They remain in the returned eigensystem as unresolved-zero directions; they are not added to the stationary sector or deleted.
2. Four bottom-dominated candidates have a minimum frequency-to-uncertainty ratio above $10^{12}$, comfortably exceeding $10^3$, but their reference trusted volume-APV, bottom-evolution, and strong-equation defects are approximately `1.60e-1`, `4.19e-6`, and `3.67e-1`. They remain unresolved pending subspace convergence.

A fixed-horizontal-support vertical study at polynomial degrees four, six, and eight shows that the dense result contains genuinely convergent physical branches. A tracked internal-wave branch gives:

| Degree | $\omega/f$ | APV defect | Bottom defect | Strong residual |
|---:|---:|---:|---:|---:|
| 4 | `4.428504` | `7.07e-6` | `3.49e-10` | `1.78e-3` |
| 6 | `4.428536` | `3.14e-7` | `4.88e-10` | `4.37e-5` |
| 8 | `4.428536` | `7.56e-9` | `5.22e-10` | `6.34e-6` |

The tracked subinertial bottom-dominated branch improves but is not yet physically classified:

| Degree | $\omega/f$ | APV defect | Bottom defect | Strong residual |
|---:|---:|---:|---:|---:|
| 4 | `0.168758` | `1.60e-1` | `4.19e-6` | `3.67e-1` |
| 6 | `0.169182` | `4.03e-2` | `4.74e-6` | `1.54e-1` |
| 8 | `0.168981` | `6.39e-3` | `4.81e-6` | `5.82e-2` |

The Milestone-8 geostrophic construction remains well resolved. Holding its trusted scalar degree fixed at four while enriching the primitive degree from four to eight keeps nine stationary states and two nonstationary bottom directions; representation, Green-identity, stationary-row, and strong bottom-tangency defects remain between approximately `1e-18` and `1e-14`. When the stationary polynomial degree grows with the primitive degree, the represented stationary dimension grows from 9 to 15 to 21 as expected for an increasingly rich APV space. The dense complement nevertheless contains 88, 120, and 152 additional unresolved-zero directions, respectively. Those directions are algebraic completion modes, not validated geostrophic states.

A variable-stratification control additionally fails the full-space eigen-residual gate at approximately `9.88e-7`, even though energy orthogonality, conjugacy, and bottom-direction retention remain accurate.

The existing public audit retains the historical status `dense-modal-classification-blocker` because its original all-mode gate does not pass. The revised interpretation is narrower: the dense energy eigensystem is valid, at least one internal-wave branch converges, and the unresolved-zero and subinertial candidate subspaces require the dedicated Milestone-9.1 convergence analysis. Milestone 10 remains inactive.

### Stopping condition

Proceed only to Milestone 9.1. Do not begin residual enrichment, matrix-free operators, or time integration until the physical spectral subspaces have been separated from the unresolved algebraic completion by nested-space convergence without changing the Milestone-8 stationary definition or deleting modes.

## Milestone 9.1: Converged physical-subspace classification

- [x] Complete — physical stationary and internal-wave subspaces validated; topographic-boundary candidates remain unresolved

### Purpose

Separate converged physical stationary, internal-wave, and topographic-boundary-wave subspaces from the unresolved algebraic completion of the dense finite basis.

### Dependencies

Milestone 9.

### Deliverables

- Extend the dense modal audit across nested horizontal supports, vertical polynomial degrees, padding factors, and terrain amplitudes.
- Compare degenerate and nearly degenerate eigenspaces through physical-energy spectral projectors and principal angles rather than individual eigenvector ordering.
- Split the Milestone-8 stationary space into APV-bearing and zero-volume-APV sectors using the APV map and balanced energy--enstrophy problem.
- Identify internal-wave subspaces through flat-limit continuation, wave polarization, frequency, and weak bottom participation.
- Identify topographic-boundary-wave subspaces through bottom participation, boundary-PV evolution, frequency, and continuation to stationary bottom states as $h\to0$.
- Return a conjugate-closed classification for every eigenvector: validated stationary family, validated internal-wave family, validated topographic-boundary-wave family, unresolved zero, unresolved nonzero, or demonstrated numerical.
- Retain every unresolved direction in the complete eigensystem and report the physical-energy fraction of any input state that lies in the unresolved remainder.
- Use constant-stratification sinusoidal terrain as the primary oracle, with flat, uniform-depth, variable-stratification, and small-terrain-amplitude controls.

The finite classification must satisfy

```math
\mathcal V_N
=
\mathcal G_N
\mathbin{\oplus_{H_N}}
\mathcal W_N
\mathbin{\oplus_{H_N}}
\mathcal R_N,
```

and

```math
\dim\mathcal V_N
=
\dim\mathcal G_N+\dim\mathcal W_N+\dim\mathcal R_N.
```

Here $\mathcal R_N$ is the retained unresolved algebraic completion. A direction is called numerical only after independent refinement demonstrates nonconvergence; otherwise it remains unresolved.

### Automated acceptance

- Fixed smooth interior-APV and bottom-buoyancy inversions converge in physical energy, reconstructed fields, APV, and bottom tangency.
- Stationary dimensions and energy projectors stabilize on a fixed trusted band while horizontal support and vertical resolution increase.
- Internal-wave spectral subspaces converge in principal angle and frequency; APV and bottom defects are below $10^{-8}$, and strong residuals decrease by at least a factor of four across successive refinements with a final value below $10^{-5}$.
- A topographic boundary branch is claimed only when:
  - its frequency exceeds its backward uncertainty by at least $10^3$;
  - its conjugate spectral subspace converges;
  - APV, bottom, and strong residuals all decrease across at least three refinements;
  - bottom participation approaches a nonzero limit; and
  - its frequency tends continuously to zero as $h\to0$.
- Accepted physical projectors are insensitive to padding factors two and three and to additional guard modes within $10^{-8}$.
- Physical and remainder projectors are $H_N$-orthogonal, Fourier-conjugate closed, and account for the complete finite dimension within $10^{-10}$.
- Modes that do not pass remain unresolved. The `demonstrated-numerical` label requires nonconvergence under independent horizontal, vertical, support, and padding refinement.
- Failure to establish a converged $\mathcal G_N$ and at least one converged dynamical subspace blocks Milestone 10.

### Numerical outcome

The public `auditConvergedPhysicalSubspaces` oracle preserves the complete dense eigensystem and compares conjugate-closed, nearly degenerate eigenspaces through physical-energy projectors. Spectral blocks are formed from complete positive- and negative-frequency eigenspaces rather than individual conjugate partners, which are not unique inside a degenerate block. Exact nested maps preserve the raw state, bottom value, and Fourier conjugacy below `2e-16`.

For constant stratification and weak sinusoidal terrain \(h=2.5\cos(2\pi y/L_y)\ {\rm m}\), the reference calculation uses trusted bounds `[1 0]`, nested supports `[1 2;1 3;1 4]`, primitive degrees `[4;6;8]`, padding factors two and three, and terrain scales `[0.125;0.25;0.5;1]`. The nine-dimensional stationary projector converges with successive principal sines `1.60e-9` and `1.78e-12`; its fixed-degree guard and padding defects are `1.67e-12` and `7.82e-13`. The APV split has stable rank nine, so this reference stationary family is entirely APV bearing.

One four-dimensional, conjugate-closed internal-wave subspace near \(4.4283f\) passes every physical gate:

| Degree | Principal sine | Relative frequency change | APV defect | Bottom defect | Strong residual |
|---:|---:|---:|---:|---:|---:|
| 4 | — | — | `1.10e-7` | `1.46e-11` | `1.78e-3` |
| 6 | `3.91e-3` | `6.91e-6` | `4.82e-9` | `9.73e-14` | `2.75e-5` |
| 8 | `5.91e-5` | `1.62e-9` | `1.07e-10` | `1.46e-13` | `2.42e-7` |

Its independent fixed-degree guard defect is `5.21e-9`, padding defect is `1.58e-13`, and minimum frequency-to-backward-uncertainty ratio exceeds \(3.6\times10^{13}\). The strong residual decreases by factors of approximately 65 and 113.

No topographic-boundary-wave subspace is claimed. Its leading candidate improves in APV, bottom evolution, and strong residual, but its final guard defect is `2.10e-2`; it therefore remains unresolved. No direction is labelled demonstrated numerical.

At the finest level,

```math
\dim\mathcal V_N=703,
\qquad
\dim\mathcal G_N=9,
\qquad
\dim\mathcal W_N=4,
\qquad
\dim\mathcal R_N=690.
```

The stationary, wave, and unresolved projectors are physical-energy orthogonal, Fourier-conjugate closed, and complete within `3.1e-13`. The unresolved space contains 168 zero-frequency and 522 nonzero-frequency algebraic directions; all remain in the returned eigensystem. The oracle returns `physical-subspace-classification-oracle`.

The focused Milestone-9.1 tests and complete repository suite pass 162 tests with zero failures. Static analysis reports no issues in all 79 MATLAB files.

### Stopping condition

Milestone 9.1 passes and its goal stops here. Milestone 9.2 is the next separately authorized compression experiment. Milestone 10, residual enrichment, matrix-free operators, and time integration remain inactive. A later topographic-boundary-wave claim must still satisfy the independent Milestone-9.1 convergence gates.

## Milestone 9.2: Boundary-complete vertical-mode compression oracle

- [x] Complete — `modal-incompatible`

### Purpose

Determine whether a vertical basis built from the known flat physical solution families reaches the Milestone-9.1 physical subspaces with substantially fewer degrees of freedom than the primitive polynomial oracle.

The candidate basis is

```math
\boxed{
\text{fixed-}\kappa\text{ nonhydrostatic waves}
+
\text{Robin APV-bearing geostrophic modes}
+
\text{an explicit zero-APV bottom inversion}
+
\text{the }\kappa=0\text{ MDA sector}.
}
```

The signed Robin length is a numerical basis parameter. It does not replace physical finite-terrain energy or introduce a new invariant.

### Dependencies

Milestone 9.1 and the boundary-complete vertical-basis derivation in `finite-terrain-projection-problem.tex` and `terrain-energy-galerkin.tex`.

### Deliverables

- Add a dense audit that constructs, for every retained nonzero horizontal wavenumber:
  - fixed-$\kappa$ nonhydrostatic wave modes;
  - generalized-energy Robin geostrophic modes, retaining every negative eigendepth and the requested positive eigendepths;
  - one independently normalized complete zero-APV bottom inversion.
- Retain the compatible MDA coordinates at $\kappa=0$.
- Represent every modal column first in the existing primitive polynomial oracle. Reuse the validated terrain quadrature, projection, $H_\gamma,J_\gamma,Q_\gamma,B_N,R_{h,N}$ forms, and Milestone-8 stationary construction without alteration.
- Test

```math
J_{\mathrm w}=J_{\mathrm g}\in\{1,2,4,8\}.
```

- Use $\ell_b=-D/4$ as the primary signed Robin length and sweep

```math
\frac{\ell_b}{D}
\in
\left\{
-\frac18,-\frac14,-\frac12,\infty
\right\}.
```

  The $\ell_b=\infty$ ordinary-boundary family is the convergence control.
- Use the Milestone-9.1 constant-$N$ sinusoidal terrain with trusted bounds `[1 0]`, support `[1 4]`, padding factors two and three, and a degree-12 primitive reference. Add flat, uniform-depth, variable-stratification, and small-terrain-amplitude controls.
- Construct the terrain-dependent stationary inclusion in the modal coordinates. Append the physical-energy-orthogonal representation residual whenever the truncated flat modal span does not contain a required stationary state.
- Solve the unmodified physical-energy eigenproblem in the resulting compatible modal space.
- Compare:
  - total state and vertical degrees of freedom;
  - physical-energy spectral-projector angles;
  - frequency errors;
  - APV, bottom-evolution, and strong primitive residuals;
  - stationary-space representation and Green identities;
  - Gram-matrix conditioning;
  - guard, padding, and Robin-length sensitivity;
  - convergence of the unresolved Milestone-9.1 subinertial bottom candidate.
- Preserve every direction in the primitive reference eigensystem. Do not use a frequency cutoff, replacement APV rows, empirical correction, symmetrization, APV-nullspace projection, or mode deletion.

The modal coordinate count at each nonzero wavenumber is

```math
N_{\mathrm{modal}}
=
2J_{\mathrm w}+J_{\mathrm g}+1,
```

to be compared with $3p+2$ admissible primitive coordinates at polynomial degree $p$.

### Automated acceptance

- Vertical EVP, endpoint, normalization, and flat reconstruction defects are below $10^{-11}$.
- The number of negative Robin eigendepths agrees with the inertia of the endpoint quadratic form.
- Every explicit bottom inversion has unit bottom displacement and volume-APV defect below $10^{-11}$.
- The modal coordinate Gram matrix is full rank without deleting or merging a physical state; its condition number and any near-dependence are reported.
- The assembled forms satisfy

```math
\frac{\lVert H_\gamma-H_\gamma^*\rVert}{\lVert H_\gamma\rVert}
\leq10^{-12},
\qquad
\frac{\lVert J_\gamma+J_\gamma^*\rVert}{\lVert J_\gamma\rVert}
\leq10^{-12}.
```

- The stationary Green identities and bottom evolution close below $10^{-10}$.
- Internal-wave frequencies and physical-energy projectors agree with the degree-12 primitive oracle within $10^{-8}$.
- Accepted physical projectors are insensitive to padding and guard support within $10^{-8}$. Dependence on $\ell_b$ is reported as a finite-order convergence diagnostic rather than used as a hard physical gate.
- A topographic boundary-wave claim must pass every Milestone-9.1 APV, bottom, strong-residual, backward-error, bottom-participation, guard, and $h\to0$ gate.
- Flat and uniform-depth calculations recover the established wave, APV-bearing geostrophic, zero-APV bottom, and MDA sectors.
- Variable-stratification results converge under independent modal-count and primitive-reference refinement.

### Outcome classification

- **`modal-acceleration`:** every physical-projector gate passes with at least a factor-two reduction in admissible vertical degrees of freedom.
- **`modal-equivalent`:** the modal family converges to the primitive physical subspaces but provides less than a factor-two reduction.
- **`modal-incompatible`:** the compatible Green, APV, or bottom identities do not converge.

The milestone completes after recording one of these outcomes. Only `modal-acceleration` becomes the preferred Milestone-10 seed. Otherwise Milestone 10 retains the primitive polynomial seed. A topographic boundary branch remains unresolved unless it independently passes the Milestone-9.1 classification gates.

### Outcome

The audit uses the unchanged Milestone-9.1 classification at primitive degrees four, six, and eight, embeds its validated internal-wave projector into one degree-12 primitive comparison space, and retains that comparison space's complete eigensystem. This avoids rebuilding a full degree-12 refinement sequence while preserving the primitive reference and public coefficient ordering.

The isolated `internal-modes-evp` checkout remains at `InternalModesEVP` commit `df86687e91faa31bf65941299062d125a96904b1`. Two provider qualifications were required without modifying that checkout:

- basis evaluation is performed in descending physical `z`, matching the provider's native Chebyshev ordering, and restored to the caller's ordering afterward;
- the Robin EVP is divided by its common factor \(f^2\), giving the equivalent but well-scaled coefficients \(p=1/N^2\) and \(r=1/g\).

Orders 96 and 192 return a spurious enormous negative Robin eigenvalue and are rejected. Orders 128 and 256 return the same physical Robin spectrum, including the single negative root for each negative \(\ell_b\). The remaining maximum order-to-order profile and endpoint defects are `1.02e-10` and `1.23e-10`. These narrowly miss the requested `1e-11` provider tolerance, but are far too small to explain the independent finite-terrain defects below.

For the primary \(\ell_b=-D/4\), padding factor two, and degree-12 primitive comparison, the count sweep is:

| \(J_{\mathrm w}=J_{\mathrm g}\) | Completed dimension | Compression factor | Internal projector defect | Frequency defect | Bottom-evolution defect | Internal APV defect | Strong residual |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 450 | `2.2822` | `2.50e-3` | `1.22e-6` | `6.14e-1` | `2.41e-9` | `5.80e-2` |
| 2 | 528 | `1.9451` | `2.90e-4` | `1.59e-8` | `8.97e-1` | `2.30e-7` | `6.29e-3` |
| 4 | 680 | `1.5103` | `2.76e-5` | `2.87e-10` | `8.93e-1` | `1.60e-6` | `6.82e-4` |
| 8 | 941 | `1.0914` | `7.86e-7` | `3.84e-13` | `3.10e-3` | `7.31e-11` | `4.36e-6` |

The primitive state dimension is 1026. The explicit stationary completion accounts for every missing stationary direction rather than deleting it; at count eight this leaves 941 coordinates and only a `1.0914` compression factor. The degree-8 to degree-12 nested primitive transfer closes to `2.51e-16`. Physical-energy and exchange structure remain at roundoff, and the trusted stationary representation is accurate to approximately `2.3e-12`.

Padding factors two and three change the internal projector by `1.58e-9`, which passes the padding gate. Changing the Robin length from \(-D/4\) produces projector changes between `1.39e-5` and `2.87e-5`. This is retained as evidence that \(\ell_b\) materially shapes a low-order trial space, not as a physical failure: future use must tune \(\ell_b\) on a training case and freeze it for validation. No topographic-boundary-wave claim is made.

The result is therefore **`modal-incompatible` for the tested compression strategy**. The reference modes themselves are useful coordinates—the frequency, APV, and strong residual improve substantially—but the flat Dirichlet wave family reaches the bottom identity too slowly and the required stationary completion removes essentially all compression before the internal projector reaches its tolerance. The cross-\(\ell_b\) variation is not part of this incompatibility classification. Milestone 9.3 tests a slope-compatible replacement for the wave coordinates before any Milestone-10 seed is selected. No mode was removed, merged, or classified by a frequency cutoff, and no APV row or generator was altered.

The six focused Milestone-9.2 tests and the complete repository suite pass 168 tests with zero failures. Static analysis reports no issues in all 83 MATLAB files.

### Stopping condition

Milestone 9.2 stops with the flat-Dirichlet boundary-complete modal seed rejected. Milestone 9.3 is the next separately authorized compression experiment. Residual enrichment, matrix-free operators, and time integration have not begun.

## Milestone 9.3: Slope-compatible wave-coordinate oracle

- [x] Complete — `per-wavenumber-slope-incompatible`

### Purpose

Determine whether carrying the leading active-bottom relation in the wave coordinates gives faster convergence than flat Dirichlet wave modes, and whether hydrostatic coordinates are as effective as nonhydrostatic coordinates once the final finite-terrain problem retains the unchanged primitive \(H_\gamma,J_\gamma\) forms.

The reference calculation leaves the interior flat and puts a constant reference slope \(\boldsymbol s\) only in the bottom equation. For hydrostatic pressure \(p/\rho_0=F(z)e^{i\boldsymbol K\cdot\boldsymbol x-i\omega t}\), the endpoint is

```math
\frac{i\omega}{N_b^2}F_z(-D)
=
\frac{
\omega\boldsymbol K\cdot\boldsymbol s
+if(\ell s_x-ks_y)}
{\omega^2-f^2}
F(-D).
```

Because this condition depends rationally on \(\omega\), the reference modes are constructed from a first-order primitive generalized eigenproblem with bottom displacement retained explicitly, not from a scalar real Robin problem.

### Dependencies

Milestones 9.1–9.2 and the slope-compatible derivation in `finite-terrain-projection-problem.tex` and `terrain-energy-galerkin.tex`. The isolated `internal-modes-evp` checkout may continue to supply the Robin geostrophic modes and zero-APV bottom inversion, but neither InternalModes checkout may be modified without approval.

### Deliverables

- Correct the Milestone-9.2 implementation so cross-\(\ell_b\) variation is diagnostic, not a hard pass condition.
- Add hydrostatic and nonhydrostatic boundary-only reference descriptors for each retained nonzero \(\boldsymbol K\):
  - the hydrostatic version makes vertical momentum diagnostic;
  - the nonhydrostatic version retains vertical acceleration;
  - both retain continuity, pressure gauge, surface normal flow, and the active bottom-displacement row;
  - slope appears only in the bottom row.
- Classify finite generalized eigenvalues with homogeneous generalized-Schur data, using \((\alpha,\beta)\) and normwise backward error. Do not use a frequency cutoff to remove diagnostic or constraint directions.
- Verify the reference descriptors independently before embedding their modes in the primitive oracle:
  - flat-limit frequencies and polarization;
  - hydrostatic or full vertical-momentum residual;
  - endpoint residual;
  - bottom displacement;
  - Fourier conjugacy.
- Retain the Milestone-9.2 Robin APV-bearing geostrophic family, explicit zero-APV bottom inversion, and MDA sector unchanged.
- Compare four wave-coordinate controls:
  - hydrostatic Dirichlet;
  - nonhydrostatic Dirichlet;
  - hydrostatic slope compatible;
  - nonhydrostatic slope compatible.
- Derive the reference-slope direction from the terrain-gradient covariance. For rank-one sinusoidal terrain, use its principal direction and symmetric positive and negative slope samples.
- Use the default parameter grids

```math
\frac{\ell_b}{D}
\in
\left\{
-\frac18,-\frac14,-\frac12,-1,\infty
\right\},
\qquad
\frac{s_{\rm ref}}{s_{\rm rms}}
\in
\left\{
0,\frac12,1,\sqrt2
\right\}.
```

- Select \(\ell_b\), reference-slope factor, and hydrostatic/nonhydrostatic family on the padding-two, largest-modal-count training calculation using

```math
\mathcal S
=
\max\left(
\frac{e_P}{10^{-8}},
\frac{e_b}{10^{-10}},
\frac{e_Q}{10^{-8}},
\frac{e_{\rm strong}}{10^{-5}}
\right).
```

  Candidates must first pass provider, endpoint, rank, and matrix-structure gates. Break scores within one percent by smaller completed dimension and then lower physical-energy Gram condition number.
- Freeze the selected parameters before validating padding factor three, other guard supports, smaller terrain amplitudes, and variable stratification. Do not retune on a validation case.
- Use the Milestone-9.1 constant-\(N\) sinusoidal terrain, trusted bounds `[1 0]`, supports `[1 2;1 3;1 4]`, stationary degree four, primitive degrees `[4;6;8]`, comparison degree 12, modal counts `[1;2;4;8]`, padding factors `[2;3]`, terrain scales `[0.125;0.25;0.5;1]`, and provider orders `[128;256]`.
- Assemble and solve only the unchanged primitive finite-terrain \(H_\gamma,J_\gamma\) system in the compressed coordinates. Preserve the primitive reference eigensystem, public coefficient ordering, stationary completion, and every unresolved direction.

### Automated acceptance

- Local reference descriptor, endpoint, and hydrostatic or vertical-momentum residuals are below \(10^{-11}\).
- The flat limit recovers the reference wave frequencies within \(10^{-10}\), and each accepted family is Fourier-conjugate closed.
- Every explicit zero-APV bottom inversion retains unit bottom displacement and volume-APV defect below \(10^{-11}\).
- The compressed terrain forms satisfy

```math
\frac{\lVert H_\gamma-H_\gamma^*\rVert}{\lVert H_\gamma\rVert}
\leq10^{-12},
\qquad
\frac{\lVert J_\gamma+J_\gamma^*\rVert}{\lVert J_\gamma\rVert}
\leq10^{-12}.
```

- Stationary representation and Green identities close below \(10^{-10}\).
- The accepted internal-wave projector and frequencies agree with the degree-12 primitive oracle within \(10^{-8}\).
- Accepted APV, bottom-evolution, and strong primitive residuals are below \(10^{-8}\), \(10^{-10}\), and \(10^{-5}\), respectively; strong residuals decrease by at least a factor of four across the final three modal counts.
- Frozen-parameter validation is insensitive to padding factors two and three and to additional guard support within \(10^{-8}\).
- The validation result remains convergent for smaller terrain amplitudes and variable stratification without retuning.
- A topographic boundary-wave claim remains subject to every Milestone-9.1 backward-error, projector, APV, bottom, strong-residual, bottom-participation, guard, padding, and \(h\to0\) gate.
- Every unresolved direction remains in the reference eigensystem. No frequency cutoff, replacement APV row, empirical correction, symmetrization, APV-nullspace projection, or mode deletion is permitted.

### Outcome classification

- **`hydrostatic-slope-acceleration`:** the frozen hydrostatic slope-compatible family passes all physical-projector gates with at least a factor-two reduction in admissible vertical degrees of freedom.
- **`nonhydrostatic-slope-acceleration`:** only the nonhydrostatic family passes with at least a factor-two reduction, or it provides more than a ten-percent compression advantage over a passing hydrostatic family.
- **`slope-modal-equivalent`:** a frozen slope-compatible family converges to the primitive physical subspaces but supplies less than a factor-two reduction.
- **`per-wavenumber-slope-incompatible`:** neither frozen family satisfies the compatible physical gates.

If both accelerated families differ in compression by no more than ten percent, prefer the hydrostatic coordinates. The milestone completes after recording one outcome and selecting the corresponding Milestone-10 seed. It does not implement Milestone 10.

### Outcome

The oracle evaluates the complete declared grid of five Robin lengths, seven signed or zero reference-slope samples for each of the hydrostatic and nonhydrostatic families, modal counts `[1;2;4;8]`, and padding factors two and three. The reference slope is obtained from the terrain-gradient covariance. Candidate parameters are scored only at padding two and the largest modal count; the selected parameters are then frozen while count and padding convergence are evaluated.

The project-local descriptor uses homogeneous generalized-Schur \((\alpha,\beta)\) data to distinguish finite dynamical modes from infinite constraint directions. Candidate wave subspaces are selected by physical-energy overlap with the flat wave projector rather than by a frequency threshold. Across the accepted local candidates, descriptor and momentum residuals are \(O(10^{-14})\), active-bottom residuals are \(O(10^{-13})\), and normwise backward errors are \(O(10^{-17})\). Neither InternalModes checkout is modified.

Both dynamical families select

```math
\frac{s_{\rm ref}}{s_{\rm rms}}=\frac12,
\qquad
\frac{\ell_b}{D}=\infty.
```

Their frozen validation results agree to the reported precision:

| Diagnostic | Hydrostatic | Nonhydrostatic | Gate |
|---|---:|---:|---:|
| Internal projector defect | `3.251e-7` | `3.251e-7` | `1e-8` |
| Internal frequency defect | `8.452e-14` | `8.333e-14` | `1e-8` |
| Internal APV defect | `3.139e-11` | `3.139e-11` | `1e-8` |
| Bottom-evolution defect | `4.597e-3` | `4.597e-3` | `1e-10` |
| Strong primitive residual | `2.073e-6` | `2.073e-6` | `1e-5` |
| Padding projector defect | `9.888e-10` | `9.763e-10` | `1e-8` |
| Local descriptor residual | `5.548e-14` | `5.067e-14` | `1e-11` |
| Local active-bottom residual | `8.121e-13` | `4.470e-13` | `1e-11` |
| Local backward error | `1.656e-17` | `1.481e-17` | `1e-11` |
| Compression factor | `1.075` | `1.075` | at least `2` |

The finite-terrain \(H_\gamma,J_\gamma\) forms are unchanged and retain their Hermitian/skew-Hermitian structure. The slope-compatible modes carry nonzero bottom displacement and satisfy their local active-bottom row, but the exact periodic-terrain bottom map is a global horizontal convolution. After the complete stationary space is appended, the compressed span remains nearly as large as the primitive reference and its global bottom identity does not converge. Vertical inertia at the coordinate-construction stage supplies no measurable advantage in this oracle.

The result is therefore **`per-wavenumber-slope-incompatible`** for both reference families. This does not invalidate the local endpoint derivation or imply that hydrostatic wave coordinates are physically wrong; it establishes that per-wavenumber slope tuning alone is not the missing compression mechanism. The primitive polynomial representation remains the Milestone-10 seed. Every unresolved direction is retained, and no frequency cutoff, replacement APV row, empirical correction, symmetrization, APV-nullspace projection, or mode deletion is used.

The five focused Milestone-9.3 tests and the complete repository suite pass 173 tests with zero failures. Static analysis reports no issues in all 86 MATLAB files.

### Stopping condition

Milestone 9.3 stops with `per-wavenumber-slope-incompatible`. Residual enrichment, matrix-free operators, and time integration have not begun.

## Milestone 9.4: Global first-order terrain-dressed block oracle

- [x] Complete — `global-dressing-seed`

### Purpose

Determine whether the verified global terrain tangent supplies the missing efficient coordinates. Unlike Milestones 9.2–9.3, the correction is a coupled horizontal block containing the complete \(O(h)\) terrain sidebands rather than a new vertical endpoint at one wavenumber.

For

```math
h=\delta\widetilde h,
\qquad
iJ_\gamma\boldsymbol c=\omega H_\gamma\boldsymbol c,
```

the first terrain coefficient satisfies

```math
\left(iJ_0-\omega_0H_0\right)\boldsymbol c_1
=-
\left(iJ_1-\omega_0H_1-\omega_1H_0\right)\boldsymbol c_0,
```

and

```math
-i\omega_0B\boldsymbol c_1-i\omega_1B\boldsymbol c_0
=
R_1\boldsymbol c_0.
```

The first-order construction selects coordinates only. Every finite-amplitude reduced solve continues to use the unchanged exact primitive \(H_\gamma,J_\gamma\) forms and the exact Milestone-8 stationary space.

### Dependencies

Milestones 6.8, 7, 8, 9.1, and 9.3. The primitive polynomial representation remains the ambient oracle. Neither InternalModes checkout participates in or may be modified by this milestone.

### Deliverables

- Add `auditGlobalFirstOrderTerrainDressing` with defaults:

  ```matlab
  audit = problem.auditGlobalFirstOrderTerrainDressing( ...
      trustedModeBounds=[1 0], ...
      supportModeBounds=[1 2;1 3;1 4], ...
      stationaryPolynomialDegree=4, ...
      primitivePolynomialDegrees=[4;6;8], ...
      comparisonPolynomialDegree=12, ...
      paddingFactors=[2;3], ...
      terrainScales=[1/16;1/8;1/4;1/2;1], ...
      tangentStep=1e-3);
  ```

- Reuse the analytic and independently centered \(H_1,J_1,R_1\) construction from the global primitive tangent oracle and the \(G_0,G_1\) stationary inclusion from the boundary-complete weak oracle.
- Dress every Milestone-9.1 validated internal-wave block and the complete trusted flat zero-frequency sector containing APV-bearing stationary states, zero-APV bottom directions, and compatible MDA states.
- Retain complete exactly degenerate or backward-error-inseparable blocks. A frequency gap smaller than \(10^3\) times the combined normwise uncertainty enlarges the block; it never deletes a direction.
- Inside a block \(X_D\), solve

  ```math
  X_D^*(iJ_1-\omega_0H_1)X_D\boldsymbol d
  =
  \omega_1X_D^*H_0X_D\boldsymbol d.
  ```

- Solve the complementary correction in the complete flat physical-energy eigenbasis. Retain every correction coefficient.
- For the zero-frequency block, solve

  ```math
  X_0^*iJ_1X_0\boldsymbol d
  =
  \omega_1X_0^*H_0X_0\boldsymbol d,
  ```

  while representing the tangent stationary nullspace directly by \(G_0+\delta G_1\).
- Verify the Fourier selection rule. For one sinusoidal terrain harmonic \(\boldsymbol q\), the first correction to a seed at \(\boldsymbol K\) may occupy only \(\boldsymbol K\pm\boldsymbol q\), apart from members already retained in the complete degenerate block.
- Form \(\boldsymbol c^{[1]}=\boldsymbol c_0+\delta\boldsymbol c_1\), insert the exact finite-terrain stationary space, project the dynamical columns into its \(H_\gamma\)-orthogonal complement, and assemble exact reduced \(H_r,J_r\).
- Preserve the complete primitive eigensystem and unresolved projector as the independent reference. No unrepresented primitive direction is relabelled or deleted.
- Report separately:
  - internal-wave block corrections;
  - the tangent stationary nullspace;
  - nonzero first-order zero-block pairs;
  - bottom participation;
  - \(\omega(\delta)/\delta\) convergence;
  - APV, bottom, and strong residuals;
  - dimension and compression.

### Automated acceptance

- Analytic and centered \(H_1,J_1,R_1\) agree below \(10^{-9}\).
- Flat block energy orthogonality, conjugacy, and sinusoidal Fourier selection close below \(10^{-11}\).
- Block solvability, complementary correction, and first-order bottom residuals close below \(10^{-10}\).
- The \(G_0,G_1\) Green and stationary-row identities close below \(10^{-10}\).
- Undressed exact weak and absolute bottom residuals scale as \(O(\delta)\).
- Dressed weak and absolute bottom residuals scale as \(O(\delta^2)\), with observed order at least `1.8` over the three smallest amplitudes.
- At the smallest amplitude, the dressed residual is at least four times smaller than its undressed counterpart.
- Exact reduced forms satisfy

  ```math
  \frac{\lVert H_r-H_r^*\rVert}{\lVert H_r\rVert}\leq10^{-12},
  \qquad
  \frac{\lVert J_r+J_r^*\rVert}{\lVert J_r\rVert}\leq10^{-12}.
  ```

- Internal frequency, APV, bottom, strong-equation, projector, padding, and guard-support diagnostics use the Milestone-9.1 tolerances.
- Internal projector and bottom errors decrease monotonically under support and vertical refinement.
- Stationary, dressed physical, and unresolved reference projectors are physical-energy orthogonal and account for the complete primitive dimension within \(10^{-10}\).
- A topographic boundary-wave claim requires every Milestone-9.1 frequency-uncertainty, projector, APV, bottom, strong, bottom-participation, guard, padding, and \(h\to0\) gate. Otherwise the branch remains unresolved.
- No replacement APV row, empirical correction, post hoc symmetrization, APV-nullspace projection, absolute frequency cutoff, or mode deletion is permitted.

### Outcome classification

- **`global-dressing-acceleration`:** every physical gate passes with at least factor-two dimension reduction.
- **`global-dressing-equivalent`:** every physical gate passes with less than factor-two dimension reduction.
- **`global-dressing-seed`:** coefficient-level \(O(h^2)\) weak and bottom gates pass and improve the physical projectors, but one correction does not pass every finite-amplitude gate.
- **`global-dressing-incompatible`:** the analytic correction fails the \(O(h^2)\) weak or bottom behavior, or violates the compatible stationary/APV structure.

The first three outcomes authorize Milestone 10 and use the primitive ambient representation with the dressed blocks as initial vectors. `global-dressing-incompatible` blocks Milestone 10.

### Stopping condition

Stop after one outcome is rigorously recorded. Do not begin repeated residual enrichment, matrix-free operators, or time integration.

### Completed result

The public `auditGlobalFirstOrderTerrainDressing` diagnostic constructs the complete signed-frequency corrections in the primitive polynomial ambient space. The nominal centered-difference step `1e-3` is the middle of the three-step validation sequence `[2e-3 1e-3 5e-4]`; the required analytic/centered \(H_1,J_1,R_1\) agreement is below `1e-9`. The broader \(L_1\) and reconstructed-tendency comparison is returned separately because it contains an additional ill-conditioned inverse and is not the coefficient gate.

For the degree-12, padding-two reference, the complete flat zero-frequency block has 377 directions. Backward uncertainty resolves 32 active \(O(h)\) directions and leaves 345 directions in its tangent stationary null sector. The latter are represented by \(G_0+\delta G_1\), rather than by an arbitrary basis of the singular zero-frequency inverse. No active pair is called a topographic boundary wave in this milestone.

The maximum first-order block-correction and active-bottom defects are `5.41e-14` and below `1.0e-11`. Fourier selection is `7.48e-13`; flat energy orthogonality is `1.27e-14`; and the compatible \(G_0,G_1\) Green and stationary-row defects are `1.34e-14` and `1.29e-13`. The undressed weak and absolute bottom residuals are first order. Their dressed counterparts have measured order `2.000` over the three smallest amplitudes and improve at the smallest amplitude by factors `2.09e3` and `1.87e3`.

Every exact finite-amplitude restriction retains the physical structure: the maximum \(H_r-H_r^*\) and \(J_r+J_r^*\) defects are `1.52e-14` and `4.50e-15`. One correction is not yet a converged finite-amplitude physical basis. At full terrain amplitude, the validated internal block has containment defect `1.58e-5`, APV defect `1.70e-11`, bottom defect `8.72e-8`, and strong primitive residual `2.46e-5`; the padding-two/three projector difference is `2.58e-8`. The selective stationary-plus-dressed span uses 45 of 1027 primitive coordinates, a nominal compression factor `22.82`, but it does not satisfy the Milestone-9.1 finite-amplitude projector, bottom, strong, and padding gates.

The dressed physical projector improves by at least a factor of two over the same undressed selective span. The outcome is therefore `global-dressing-seed`. The analytic global correction captures exactly the missing \(O(h)\) terrain convolution and supplies the approved initial vectors for Milestone 10, but repeated exact-residual enrichment is still required. Milestone 10 has not begun.

The eight focused Milestone-9.4 tests and the complete repository suite pass 181 tests with zero failures. Static analysis reports no issues in all 89 MATLAB files.

## Milestone 10: Residual-enriched terrain modes

- [x] Complete — `residual-enrichment-acceleration`

### Purpose

Construct selected terrain modes efficiently while preserving the complete stationary space and the physical-energy structure of the dense oracle.

### Dependencies

Milestone 9.4 with outcome `global-dressing-acceleration`, `global-dressing-equivalent`, or `global-dressing-seed`.

### Deliverables

- Retain the primitive ambient representation and initialize each selected block with the verified global first-order dressed coordinates from Milestone 9.4.
- Use the validated Milestone-9.1 physical stationary, internal-wave, and topographic-boundary-wave subspaces as the reference block projectors for the physical-energy eigenproblem.
- Apply repeated exact finite-amplitude residual correction to complete resonant blocks while retaining every bottom coordinate, the Milestone-8 stationary subspace, and the unresolved algebraic completion needed during construction.
- Recompute a complete block whenever its Ritz frequencies become backward-error inseparable; never continue an individual-vector correction through a newly detected resonance.
- Preserve the APV-bearing Robin sector, explicit zero-APV bottom sector, stationary sector, and unresolved reference completion throughout every accepted correction.
- Orthogonalize accepted corrections using $H_\gamma$.
- Compare every enriched physical invariant subspace with the dense Milestone-9.1 projectors.
- For an approximate block $X$ with Ritz vectors $XC$ and frequencies $\Theta$, evaluate the exact residual

  ```math
  R=iJ_\gamma XC-H_\gamma XC\Theta.
  ```

- Precondition internal-wave corrections with the flat signed-frequency eigensystem. Precondition directions emerging from the flat zero-frequency sector with the verified first-order zero-block splitting, because the unmodified flat zero block is singular.
- Use both preconditioners only to generate coordinates. Every residual, Ritz value, physical diagnostic, and accepted projector continues to use the exact finite-terrain forms.
- Compare an enriched block with an independently diagonalized dense block at the same resolution. Report cross-resolution and cross-padding drift separately so errors in the selective coordinates are not confused with errors already present in the dense discretization.

### Automated acceptance

- One correction reduces weak nonresonant $O(h)$ residuals to $O(h^2)$.
- Repeated enrichment reproduces validated dense eigenvalues and physical invariant subspaces within $10^{-9}$.
- Physical-energy structure, the complete stationary subspace, and the physical/remainder classification remain invariant after every accepted enrichment; projected APV and bottom evolution retain their validated exact or convergent behavior.
- Resonant calculations converge only when the complete coupled stationary and wave block is retained.

### Completed result

The public `auditResidualEnrichedTerrainModes` diagnostic reuses the complete Milestone-9.4 setup and preserves its result. The exact finite-amplitude iteration works in the primitive polynomial ambient space, keeps the nine-dimensional exact stationary space fixed, closes every correction under Fourier conjugacy, and expands any backward-error-inseparable Ritz cluster before applying a correction. Directions are removed from a proposed correction only when finite-terrain energy orthogonalization proves them algebraically dependent.

The degree-12, padding-two reference begins with 36 dynamical seed coordinates. Three exact-residual corrections increase the dynamical trial dimension to 156 inside the 1027-dimensional primitive ambient space. The physical stationary-plus-trial representation therefore has compression factor `6.224`.

The validated internal-wave block reaches:

| Diagnostic | Final defect |
|---|---:|
| Dense physical-energy projector | numerical zero |
| Relative frequency | below \(10^{-12}\) |
| Energy-scaled Ritz residual | \(1.37\times10^{-11}\) |
| Projected volume APV | below \(10^{-12}\) |
| Bottom evolution | \(1.55\times10^{-12}\) |
| Strong primitive equations | below \(2\times10^{-10}\) |
| Padding excess over the dense oracle | numerical zero |
| Nested-representation excess | \(1.08\times10^{-9}\) |

The internal Ritz residual decreases from `2.13e-6` for the globally dressed seed to `8.15e-9`, `8.37e-10`, and `1.37e-11` under successive exact corrections. Exact restrictions retain Hermitian physical energy and skew-Hermitian exchange at roundoff throughout.

The active directions emerging from the flat zero-frequency block remain explicitly retained as unresolved topographic-boundary candidates. Milestone 10 does not promote them to physical modes because it does not repeat the complete terrain-amplitude continuation required by Milestone 9.1. Their omission from the validated internal projector does not delete them from the primitive ambient eigensystem or unresolved complement.

The measured outcome is `residual-enrichment-acceleration` for the selected internal-wave block. This proves that global dressing followed by exact residual enrichment is an effective coordinate construction, but it does not define an economical complete state for arbitrary forward evolution. Milestone 10.1 must define that production state before matrix-free work begins.

The nine focused Milestone-10 tests and the complete repository suite pass 190 tests with zero failures. Static analysis reports no issues in all 92 MATLAB files.

## Milestone 10.1: Production physical-state contract

- [x] Complete — blocking gate passed

### Purpose

Define the complete finite-dimensional state that the forward wave–vortex model will evolve. The dense primitive polynomial system remains an overresolved independent oracle and is not the online state.

### Dependencies

Milestone 10.

### Deliverables

- For every retained horizontal wavenumber, construct fixed-$\kappa$ energy-orthogonal flat nonhydrostatic wave pairs, APV-bearing balanced inversion coordinates, one explicit zero-APV bottom inversion, and the compatible $\kappa=0$ MDA sector.
- Declare the production dimension

  ```math
  N_{\rm prod}
  =
  \sum_{\boldsymbol K\ne\boldsymbol 0}
  \left[
  2J_{\rm w}(\boldsymbol K)
  +
  J_{\rm q}(\boldsymbol K)
  +
  1
  \right]
  +
  2J_{\rm io}
  +
  N_{\rm MDA}
  +
  N_{\rm mean,b},
  ```

  where \(J_{\rm io}=1+J_{\rm w}(\boldsymbol 0)\) counts the retained zero-horizontal-wavenumber inertial vertical structures and \(N_{\rm mean,b}=1\) is the compatible mean-bottom coordinate.

- Separate trusted prognostic modes, outer guard modes used only for dealiased evaluation, and deliberately discarded primitive-oracle modes.
- Map every production coordinate into the primitive oracle and retain the public WaveVortexModel coefficient and Fourier-conjugacy conventions.
- Treat Robin length and local reference slope only as historical or optional coordinate-shaping parameters. They are not part of the production physical definition.
- Add a public audit that reports the degree-of-freedom budget, physical family, units, conjugate partner, native WaveVortexModel index, and primitive-oracle image of every production coordinate.

### Automated acceptance

- Coordinate counts, family counts, conjugacy maps, and the declared total dimension close exactly.
- Flat reconstruction/projection round trips, endpoint values, APV classification, and physical-energy normalization close within $10^{-12}$.
- The production flat wave, balanced, bottom, and MDA projectors converge toward the trusted primitive-oracle projectors under independent vertical and horizontal refinement.
- Guard modes do not appear in the prognostic layout and discarded oracle coordinates are reported only as truncation error.
- No coordinate inside the declared production state is labelled unresolved.

Failure blocks all finite-terrain production-basis work.

### Implementation result

The public `auditProductionPhysicalStateContract` diagnostic constructs the production layout without changing the existing public `horizontalLayout` or `stateLayout`. Native WaveVortexModel coefficients define the fixed-\(\kappa\) wave, inertial, balanced, and MDA coordinates. The existing complete balanced bottom inversion supplies the independent nonzero-wavenumber bottom coordinate, and the compatible linear mean-bottom profile supplies the \(\kappa=0\) bottom value. Their continuous profiles are energy normalized and mapped into the independent primitive polynomial oracle; neither Robin nor slope-compatible coordinates enter the contract.

For the reduced trusted band \(|k|\le1,\ell=0\), one retained wave mode, and balanced labels \(j\in\{0,1\}\), the complete accounting is:

| Flat physical family | Complex coordinates |
|---|---:|
| Internal-wave | 4 |
| APV-bearing balanced | 4 |
| Zero-APV bottom inversion | 2 |
| Inertial | 4 |
| MDA | 1 |
| Compatible mean-bottom | 1 |
| **Total** | **16** |

Fourier conjugacy makes these 16 complex signed-layout entries equivalent to 16 independent real degrees of freedom. No coordinate is labelled unresolved. Guard modes are returned separately, and the primitive directions not used by the production contract are reported only as discarded oracle coordinates.

Constant-stratification calculations with antialiasing disabled and enabled, together with an arbitrary-stratification control, all return `complete-production-state-contract`. Across those controls, reconstruction round trips and energy normalization close below \(5\times10^{-16}\), endpoint and zero-APV defects remain below \(4\times10^{-13}\), the native WaveVortexModel mode match is below \(6\times10^{-14}\), and the highest-degree primitive embedding defect is below \(2.5\times10^{-11}\). Independent vertical degrees \(8,16,24\), three nested horizontal supports, and padding factors two and three satisfy the convergence gates; support and padding Gram differences remain below \(5.3\times10^{-14}\).

The focused Milestone-10.1 suite passes eight tests for constant and variable stratification and both antialias conventions. The complete repository suite passes 198 tests with zero failures, and static analysis reports no issues in all 95 MATLAB files.

## Milestone 10.2: Complete internal-wave coverage

- [x] Complete — `internal-wave-numerical-blocker`; Milestone 10.3 remains blocked

### Purpose

Extend the successful Milestone-10 construction from one four-dimensional internal-wave block to every internal-wave degree of freedom declared by Milestone 10.1.

### Dependencies

Milestone 10.1.

### Deliverables

- Insert the exact finite-terrain stationary space before constructing dynamical coordinates.
- Globally dress every complete flat signed-frequency block using the verified $H_1,J_1,R_1,G_1$ derivatives.
- Apply repeated exact finite-amplitude residual enrichment with the unchanged $H_\gamma,J_\gamma$ forms.
- Expand every backward-error-inseparable or newly resonant block before applying a correction.
- Recycle shared correction directions across neighboring horizontal and vertical blocks so the total production span is not the direct sum of independently overbuilt trial spaces.
- Preserve Fourier conjugacy, the complete stationary space, every declared bottom coordinate, and the primitive oracle throughout construction.
- Construct the complete internal-wave projector $\mathcal W_{\rm int}$ and associate every declared flat wave coordinate with one converged finite-terrain invariant subspace.

### Automated acceptance

- The internal-wave dimension equals the declared $2J_{\rm w}$ total and remains stable under support, padding, and primitive-oracle refinement.
- Every accepted block satisfies the Milestone-9.1 frequency, physical-energy projector, APV, bottom, strong-equation, guard, and padding gates.
- The union of accepted blocks is $H_\gamma$-orthogonal to the exact stationary space, conjugate closed, and internally complete within $10^{-10}$.
- No retained internal-wave coordinate remains in an unresolved production remainder.
- Recycled enrichment reduces the total trial dimension relative to independently enriched blocks and reproduces the matched dense primitive projector within $10^{-9}$.

Failure blocks the bottom and topographic-wave construction.

### Implementation result

The public `auditCompleteInternalWaveCoverage` diagnostic maps every declared native WaveVortexModel wave into the unchanged primitive polynomial oracle, groups exact signed-frequency degeneracies, applies the verified global first-order dressing, and compares joint residual enrichment with independently enriched block controls. Secondary Fourier/branch columns are generated through the declared production conjugacy map. Bottom and unresolved primitive directions remain in the ambient reference space but are not used as selective Milestone-10.2 seeds.

The requested split controls retain their complete declared dimensions:

| Control | Declared internal-wave coordinates | Frequency blocks | Result |
|---|---:|---:|---|
| Zonal trusted band, \(j=1,2\) | 8 | 2 | The degree-12 physical projector passes, but vertical refinement is not stable and joint enrichment saves no directions |
| \(\lvert k\rvert,\lvert\ell\rvert\le1\), \(j=1\) | 16 | 3 | Two blocks converge; the meridional-axis block fails the residual/projector gates |

For the degree-12 zonal reference, the eight-dimensional dense target is recovered with energy-scaled Ritz residual \(9.81\times10^{-12}\), projected-APV defect \(1.16\times10^{-9}\), bottom defect \(7.12\times10^{-12}\), strong primitive residual \(3.23\times10^{-6}\), and numerical-zero dense-projector and frequency defects. The joint and independent constructions both require 48 trial directions; their physical projector difference is \(6.59\times10^{-10}\). Thus the required correction recycling produces no dimension reduction. At primitive degree eight, the multiblock iteration ceases to track the dense drift: its final projector defect is \(1.42\times10^{-2}\), giving nested excess \(1.42\times10^{-2}\).

The two-dimensional control separates sixteen coordinates into three native-frequency blocks. The first and third independent blocks converge, whereas the intermediate meridional-axis block reaches projector defect \(1\) and residual \(7.76\times10^{-2}\). The resulting recycled joint span is not a physical invariant subspace. This is a numerical limitation of the present flat-frequency residual correction and block-recycling strategy; it is not evidence that the corresponding continuum internal-wave family is absent.

The wider 2D contract also establishes two provider-level reproduction floors: \(1.92\times10^{-12}\) for the embedded zero-APV diagnostic and \(2.94\times10^{-10}\) for re-evaluating native modes on diagonal wavenumbers. The production-contract roundoff tolerances are therefore \(5\times10^{-12}\) and \(5\times10^{-10}\), respectively. No physical finite-terrain gate is relaxed.

The measured outcome is `internal-wave-numerical-blocker`. Milestone 10.3 must not begin. A future Milestone-10.2 reformulation needs a robust complete-block correction—most likely a projected Jacobi–Davidson or exact complementary solve—and an economical thick-restart/recycling rule that converges under vertical refinement and reduces the joint span below the independent sum without using dense eigenvectors as construction coordinates.

The focused Milestone-10.2 suite passes eight tests. The complete repository suite passes 206 tests with zero failures, and static analysis reports no issues in all 98 MATLAB files.

## Milestone 10.2.1: Coupled-block internal-wave correction

- [x] Complete — blocking gate

### Purpose

Determine whether the Milestone-10.2 failure came from generating separate flat-frequency residual corrections even though the reduced Ritz diagonalization used their joint span. Replace that frequency-local correction step by a genuinely coupled invariant-subspace correction without changing the production state, exact stationary space, primitive oracle, or finite-terrain forms.

### Dependencies

Milestone 10.2 with outcome `internal-wave-numerical-blocker`.

### Coupled formulation

Begin with every globally dressed production wave block in one finite-terrain-energy-orthonormal trial space \(X\), after projecting out the exact stationary space \(G_\gamma\). Solve the joint Ritz problem and define

```math
Y=XC,
\qquad
R=iJ_\gamma Y-H_\gamma Y\Theta,
\qquad
U=\begin{bmatrix}G_\gamma&Y\end{bmatrix}.
```

Treat \(Y\) as an invariant subspace. Its correction must therefore be unchanged by reordering or rotating the Ritz vectors inside that subspace. Replace the separate shifted-vector corrections by

```math
iJ_\gamma\Delta
-H_\gamma\Delta\Theta
+H_\gamma U\Lambda
=-R,
\qquad
U^*H_\gamma\Delta=0.
```

The Lagrange multiplier \(\Lambda\) enforces finite-terrain-energy orthogonality to both the exact stationary space and the current Ritz space. Every reduced solve, residual, and accepted projector continues to use the unchanged exact \(H_\gamma,J_\gamma\) forms.

### Two-stage oracle

1. **Exact complementary oracle.** Solve the coupled Sylvester or saddle system directly at reduced resolution. Determine whether the exact block correction repairs the failed meridional block and the nested vertical drift. Dense primitive eigenvectors remain validation data and must not become construction coordinates.
2. **Iterative production candidate.** Solve the same coupled equation with a block Krylov or Jacobi--Davidson method. Use the flat signed-frequency operator only as a preconditioner, share Krylov directions across all residual columns, and use thick restart while retaining the complete production Ritz subspace and every backward-error-inseparable cluster.

Compare the ordinary globally dressed production seed with the frozen \(s_{\rm ref}=s_{\rm rms}/2\) slope-compatible seed. The slope-compatible family is only an optional starting coordinate or preconditioner and does not alter the production contract or the final physical eigensystem.

### Controls

- Repeat the zonal trusted-band calculation with \(j=1,2\), giving eight declared wave coordinates.
- Repeat the two-dimensional \(\lvert k\rvert,\lvert\ell\rvert\le1\) calculation with \(j=1\), giving sixteen declared wave coordinates and retaining the previously failed meridional block.
- Use primitive degrees \(8,10,12\), padding factors two and three, nested guard supports, and at least three terrain amplitudes including the Milestone-10.2 reference amplitude.
- Retain the independent-block enrichment as the economy control and report

  ```math
  N_{\rm joint},
  \qquad
  N_{\rm independent},
  \qquad
  N_{\rm ambient}.
  ```

### Automated acceptance

- The dense physical-energy projector defect is below \(10^{-9}\).
- The energy-scaled Ritz residual is below \(10^{-10}\).
- Projected APV, bottom evolution, and strong primitive residuals are below \(10^{-8}\), \(10^{-10}\), and \(10^{-5}\), respectively.
- Stationary orthogonality and Fourier conjugacy close below \(10^{-10}\).
- Degree-eight to degree-twelve and padding-two to padding-three projector drift are below \(10^{-8}\).
- The formerly failed meridional block passes every physical, refinement, guard, and padding gate independently.
- The joint stationary-plus-wave representation achieves at least factor-two compression relative to the primitive ambient space.

### Outcome classification

- **`coupled-block-acceleration`:** every physical and refinement gate passes and \(N_{\rm joint}\le0.8N_{\rm independent}\).
- **`coupled-block-equivalent`:** every physical and refinement gate passes, but the twenty-percent correction-recycling target does not.
- **`coupled-block-iterative-blocker`:** the exact complementary oracle passes but the iterative block solver does not reproduce it robustly.
- **`coupled-block-formulation-blocker`:** the exact coupled correction fails to recover the dense physical subspace.
- **`coupled-block-isolation-blocker`:** the declared internal-wave sector cannot be separated consistently from another physical sector.

Only `coupled-block-acceleration` would have authorized Milestone 10.3 directly. The measured `coupled-block-isolation-blocker` outcome instead activates Milestone 10.2.2. No replacement APV rows, empirical corrections, post hoc symmetrization, APV-nullspace projection, frequency cutoff, mode deletion, matrix-free production work, or time integration is permitted.

### Implementation result

The public `auditCoupledBlockInternalWaveCorrection` diagnostic retains the production contract and complete stationary space, forms one physical-energy-orthonormal production wave block, and solves the constrained complementary Sylvester equation in Cholesky-scaled energy coordinates. The dense eigensystem is evaluated only afterward as an independent projector oracle. The exact correction updates the invariant subspace as \(Y+\Delta\); it does not append \(\Delta\) as an additional physical coordinate.

For the prescribed zonal \(j=1,2\) control at primitive degrees \(8,10,12\), one exact coupled correction reduces the energy-scaled Ritz residual to \(1.69\times10^{-12}\). At degree twelve, the dense-projector, projected-APV, bottom-evolution, correction-equation, and correction-orthogonality defects are \(4.03\times10^{-10}\), \(1.16\times10^{-9}\), \(3.50\times10^{-13}\), \(4.30\times10^{-15}\), and \(6.55\times10^{-17}\), respectively. Padding factors two and three agree to \(4.63\times10^{-13}\). The exact update retains eight joint physical coordinates inside a \(1027\)-dimensional ambient oracle; the Milestone-10.2 independent construction required \(48\) trial directions.

The blocking result is the nested degree-\(8,10,12\) physical-projector drift, \(5.96\times10^{-4}\), compared with the required \(10^{-8}\). Thus the coupled equation repairs the frequency-local correction algorithm at each fixed discretization, but the declared internal-wave sector is not isolated consistently from the changing primitive completion under vertical refinement. This is `coupled-block-isolation-blocker`, not `coupled-block-formulation-blocker`.

The milestone instructions require stopping as soon as the exact oracle establishes an isolation blocker. Consequently, the iterative coupled Krylov/Jacobi--Davidson realization, slope-seed comparison, and two-dimensional meridional control were not attempted. No numerical direction was removed or reclassified, and Milestone 10.3 remains blocked.

The focused Milestone-10.2.1 suite passes five tests. The complete repository suite passes 211 tests with zero failures, and `checkcode` reports no issues in all 101 MATLAB files.

## Milestone 10.2.2: Geometric-cascade and spectral-isolation audit

- [x] Complete — `cascade-slow`; blocking diagnostic gate

### Purpose

Determine whether the Milestone-10.2.1 projector drift is caused by an unresolved horizontal or vertical geometric-scattering tail, or because the declared eight-dimensional internal-wave sector becomes spectrally inseparable from a larger physical block.

The finite-terrain dynamics remain linear in the state. The cascade diagnosed here is produced by repeated multiplication by stationary geometric coefficients such as

```math
\gamma^{-1}
=
\frac{1}{1-h/D}
=
\sum_{n=0}^{\infty}
\left(\frac{h}{D}\right)^n.
```

It is therefore a representation and isolation audit, not a nonlinear wave--wave calculation and not a modification of the physical eigensystem.

### Dependencies

Milestone 10.2.1 with outcome `coupled-block-isolation-blocker`.

### Scattering and tail definitions

Let \(\mathcal H\) be the terrain Fourier support and \(\mathcal K_{\rm seed}\) the support of the flat production-wave block. Define

```math
\mathcal K_0=\mathcal K_{\rm seed},
\qquad
\mathcal K_{m+1}=\mathcal K_m+\mathcal H,
\qquad
\mathcal S_m=\mathcal K_m\setminus\mathcal K_{m-1}.
```

For sinusoidal terrain, \(\mathcal S_m\) contains the newly accessible \(\boldsymbol K\pm m\boldsymbol q\) sidebands. Let \(\mathcal V_{m,p}\) retain scattering order \(m\) and primitive vertical degree \(p\), and let \(\Pi_{m,p}\) be its physical-energy projector in one common comparison space. For an \(H_\gamma\)-orthonormal basis \(Y_\gamma\) of the continued physical block, measure

```math
\tau_{m,p}
=
\left\|
\left(I-\Pi_{m,p}\right)Y_\gamma
\right\|_{H_\gamma,2}.
```

This quantity is invariant under rotations or reorderings inside the physical block.

### Planned diagnostic

```matlab
audit = problem.auditGeometricCascadeIsolation( ...
    trustedModeBounds=[1 0], ...
    waveModeIndices=[1;2], ...
    scatteringOrders=[1;2;3;4;5], ...
    primitivePolynomialDegrees=[8;10;12;14], ...
    comparisonPolynomialDegree=16, ...
    paddingFactors=[2;3], ...
    terrainScales=[0;1/16;1/8;1/4;1/2;3/4;1]);
```

The oracle will:

1. Continue the flat eight-dimensional, symmetry-preserving, conjugate-closed wave projector through the prescribed terrain amplitudes.
2. Vary horizontal scattering order and vertical degree independently and embed all results into the degree-16 comparison space with the existing adjoint-consistent transfers.
3. Compare physical-energy spectral projectors rather than individual eigenvectors. Dense primitive eigenvectors are validation data only and never become construction coordinates.
4. Report energy by horizontal shell, vertical tail, cumulative \(\tau_{m,p}\), physical-energy principal angles, spectral gaps, normwise backward uncertainty, and padding sensitivity.
5. Decompose the worst drifting principal directions into exact stationary, declared internal-wave, bottom/MDA, higher-wave, and unresolved primitive sectors.
6. Enlarge the tracked block only when its separation from another direction is less than \(10^3\) times the combined backward uncertainty. The enlarged block must retain complete symmetry classes and Fourier conjugates.
7. Run the two-dimensional \(\lvert k\rvert,\lvert\ell\rvert\le1\), \(j=1\) control only after the zonal calculation reaches `cascade-resolved`.

The exact stationary space, primitive dense oracle, production coordinate contract, public coefficient layout, unresolved completion, and unmodified finite-terrain \(H_\gamma,J_\gamma\) forms remain fixed.

The implementation also provides a resumable content-addressed cache through `cacheDirectory`. Cache keys include the complete scientific configuration, geometry, stratification sample, coefficient layout, terrain field, and implementation signature. Only signature-validated audit artifacts are restored. The cache is opt-in, is restricted to the repository's ignored `output/` directory, and performs no disk writes when disabled. Each primitive configuration is checkpointed independently so an interrupted high-degree run can resume without invoking Milestone 10.2.1 or rebuilding completed forms.

### Automated acceptance

- Common-space transfer adjointness, physical-energy normalization, shell-energy accounting, and Fourier conjugacy close below \(10^{-11}\).
- Fixed-discretization Ritz, projected-APV, bottom-evolution, and strong primitive residuals retain the Milestone-10.2.1 tolerances \(10^{-10}\), \(10^{-8}\), \(10^{-10}\), and \(10^{-5}\).
- Padding factors two and three give physical projectors and tail measures agreeing below \(10^{-10}\).
- The final two independent horizontal-scattering refinements and the final two independent vertical refinements give projector drift and omitted-tail amplitude below \(10^{-8}\).
- An accepted block is separated from its complement by at least \(10^3\) times the combined normwise backward uncertainty.
- The accepted stationary-plus-wave representation retains at least factor-two compression relative to the primitive ambient space.
- Small-amplitude shell amplitude and shell energy are reported against their expected \(O(h^m)\) and \(O(h^{2m})\) onsets for sinusoidal terrain. These are diagnostics of the scattering path rather than gates imposed on resonant cases.
- If the zonal calculation passes, the two-dimensional control must satisfy the same physical, tail, isolation, padding, and compression gates.

### Outcome classification

- **`cascade-resolved`:** the original eight-dimensional projector is spectrally isolated, its horizontal and vertical tails converge, and every physical and economy gate passes.
- **`cascade-slow`:** the projector and tail measures decrease systematically, but the final tolerance or factor-two economy gate is not reached.
- **`resonant-block-required`:** the eight-dimensional projector loses isolation, but a smallest larger conjugate-closed physical block converges.
- **`cascade-nonconvergent`:** neither a stable physical projector nor systematic horizontal and vertical tail decay is established.

Only `cascade-resolved` authorizes Milestone 10.3. `resonant-block-required` requires an explicit revision of the production block before further work. The other outcomes stop for analysis. No iterative coupled solver, bottom-wave classification, matrix-free operator, time integration, frequency cutoff, replacement APV row, empirical correction, symmetrization, APV-nullspace projection, or mode deletion is included.

### Outcome

The constant-\(N\) zonal oracle used resolution `[6 14 5]`, trusted band `[1 0]`, the native \(j=1,2\) wave families, sinusoidal terrain of amplitude \(2.5\,\mathrm{m}\), scattering orders one through five, primitive degrees \(8,10,12,14\), degree-sixteen comparison, seven terrain amplitudes, and padding factors two and three. The production block retains eight directions, does not require backward-error enlargement, and is separated from the complementary spectrum by \(5.16\times10^{11}\) times the combined normwise uncertainty. With nine exact stationary directions in an ambient space of dimension \(1151\), the stationary-plus-wave compression factor is \(67.7\).

The final physical residuals are

```text
energy-scaled Ritz       2.74e-15
projected volume APV     1.18e-12
bottom evolution         1.01e-13
strong primitive         1.16e-07
padding projector        0
padding tail difference  2.46e-15
```

The independent horizontal tail amplitudes at scattering orders one through five are

```text
7.58e-4, 1.73e-5, 3.23e-7, 4.66e-9, 0,
```

so the geometric Fourier cascade is resolved through the declared support. The fitted shell-amplitude onsets for orders one through five are \(1.000,2.000,3.000,3.999,4.987\), matching the expected powers of terrain amplitude. The independent vertical tail amplitudes at degrees \(8,10,12,14,16\) are

```text
2.77e-4, 3.72e-5, 5.03e-6, 6.35e-7, 0.
```

They decrease systematically by factors near seven to eight, but degree fourteen remains above the \(10^{-8}\) gate. The corresponding projector defect is \(6.83\times10^{-7}\). The worst drifting principal direction contains energy fractions \(0.999999\) in the declared internal-wave sector, \(5.74\times10^{-7}\) in higher flat waves, and \(5.83\times10^{-10}\) in the unresolved primitive remainder; stationary and non-tangent-bottom fractions are negligible.

The outcome is therefore **`cascade-slow`**. The earlier \(5.96\times10^{-4}\) drift is not caused by loss of spectral isolation, a missing resonant block, padding, bottom physics, or an unresolved horizontal geometric cascade. It is a slowly converging vertical primitive tail. The two-dimensional control is not attempted because the approved cadence permits it only after `cascade-resolved`. Milestone 10.3 remains blocked. The next increment must decide whether a better vertically adapted coordinate family or a targeted higher-degree oracle can establish the same projector below \(10^{-8}\) without sacrificing economy.

The focused Milestone-10.2.2 suite passes seven tests. The complete repository suite passes 218 tests with zero failures, and static analysis reports no issues in all 104 MATLAB files.

## Milestone 10.2.3: Flat wave–vortex modal ambient oracle

- [x] Complete — `geostrophic-modal-blocker`; Milestone 10.3 remains blocked

### Purpose

Test whether the slow primitive vertical tail isolated by Milestone 10.2.2 can be replaced by a physically organized flat wave–vortex ambient space rather than by increasing generic polynomial degree. The candidate coordinates are

```math
\mathcal V_{\rm modal}
=
\mathcal W_{\kappa}^{\rm flat}
\oplus
\mathcal G_{\kappa}^{\rm flat}
\oplus
\mathcal B_{\kappa}^{0}
\oplus
\mathcal M_0,
```

where \(\mathcal W_{\kappa}^{\rm flat}\) contains fixed-\(\kappa\) nonhydrostatic internal-wave modes, \(\mathcal G_{\kappa}^{\rm flat}\) contains ordinary APV-bearing flat geostrophic modes, \(\mathcal B_{\kappa}^{0}\) is the explicit zero-APV bottom inversion, and \(\mathcal M_0\) is the compatible \(\kappa=0\) inertial, MDA, and mean-bottom sector. No Robin tuning or local-slope coordinate is used. The primitive polynomial system remains an independent validation and quadrature oracle rather than a candidate production basis.

### Dependencies

Milestone 10.2.2 with outcome `cascade-slow`, the Milestone-8 exact finite-terrain stationary construction, and the Milestone-10.1 production coordinate contract. The isolated InternalModesEVP checkout is pinned to commit `df86687` and is used only to generate and independently qualify the one-dimensional modal families; the WaveVortexModel InternalModes dependency is unchanged.

### Oracle

The public diagnostic is

```matlab
audit = problem.auditWaveVortexModalAmbient( ...
    trustedModeBounds=[1 0], ...
    targetWaveModeIndices=[1;2], ...
    waveGuardModeCounts=[2;4;6;8;12], ...
    geostrophicModeCounts=[2;4;6;8;12], ...
    scatteringOrders=[1;2;3;4;5], ...
    primitiveReferenceDegrees=[16;18;20], ...
    paddingFactors=[2;3], ...
    terrainScales=[0;1/8;1/4;1/2;3/4;1], ...
    internalModesEVPOrders=[128;256], ...
    shouldRunTwoDimensionalControl=true, ...
    cacheDirectory="output/milestone-10.2.3-cache");
```

For each horizontal coefficient the oracle constructs the modal families directly, embeds them into the primitive comparison space, removes the exact finite-terrain stationary space before identifying the dynamical wave projector, and restricts the unchanged exact finite-terrain \(H_\gamma,J_\gamma,Q_\gamma,B,R_h\) forms. It varies wave and geostrophic counts independently, compares padding and horizontal support, retains the complete primitive eigensystem as validation data, and never promotes a dense eigenvector to a construction coordinate. Every cache artifact is disposable, content addressed, signature checked, restricted to ignored `output/`, and written only when caching is enabled.

### Acceptance and outcome classification

The modal provider, endpoint conditions, flat reconstruction, conjugacy, energy normalization, and exact-form structure must close near roundoff. The finite-terrain internal-wave projector, frequencies, Ritz equation, APV, bottom evolution, and strong primitive equations must agree with the converged primitive oracle. The exact stationary space must be represented and stationary to the declared tolerances. Padding, guard support, polynomial reference degree, wave count, and geostrophic count must converge independently, and an accepted production ambient must retain at least factor-two compression.

The declared outcomes are:

- **`wave-vortex-modal-acceleration`:** all physical gates pass with at least factor-two compression;
- **`wave-vortex-modal-equivalent`:** all physical gates pass without factor-two compression;
- **`wave-modal-blocker`:** the fixed-\(\kappa\) wave family or its finite-terrain continuation does not converge;
- **`geostrophic-modal-blocker`:** the ordinary flat balanced/bottom family does not converge to the exact finite-terrain stationary space economically;
- **`primitive-reference-blocker`:** the independent polynomial validation sequence does not stabilize;
- **`provider-blocker`:** the isolated one-dimensional EVP provider fails its qualification.

The measured `geostrophic-modal-blocker` outcome routes the roadmap to Milestone 10.2.4 before Milestone 10.3. No frequency cutoff, replacement APV row, empirical correction, post hoc symmetrization, APV-nullspace projection, mode deletion, matrix-free operator, or time integration is included.

### Outcome

The approved constant-\(N\) zonal oracle used resolution `[6 14 5]`, trusted band `[1 0]`, wave and ordinary geostrophic counts \(2,4,6,8,12\), primitive reference degrees \(16,18,20\), scattering orders one through five, terrain amplitude \(2.5\,\mathrm{m}\), and padding factors two and three. The isolated provider orders 128 and 256 agree in physical-energy projector to \(1.51\times10^{-11}\), with maximum endpoint defect \(3.32\times10^{-12}\). The flat modal control is correspondingly clean: its wave-projector defect is zero, and its Ritz, APV, bottom, and strong residuals are \(1.68\times10^{-12}\), \(6.04\times10^{-14}\), \(1.00\times10^{-13}\), and \(1.47\times10^{-13}\).

The finite-terrain count sweep is

| Wave/geostrophic count | Modal dimension | Primitive/modal compression | Wave-projector defect | Stationary representation defect |
|---:|---:|---:|---:|---:|
| 2 | 160 | 8.92 | \(2.97\times10^{-2}\) | \(5.82\times10^{-1}\) |
| 4 | 298 | 4.79 | \(1.68\times10^{-3}\) | \(6.41\times10^{-2}\) |
| 6 | 436 | 3.27 | \(6.92\times10^{-4}\) | \(2.36\times10^{-2}\) |
| 8 | 574 | 2.49 | \(4.32\times10^{-4}\) | \(1.17\times10^{-2}\) |
| 12 | 850 | 1.68 | \(2.33\times10^{-4}\) | \(4.26\times10^{-3}\) |

At the largest count the frequency, Ritz, projected-APV, bottom, and strong defects are \(1.11\times10^{-8}\), \(1.23\times10^{-5}\), \(2.53\times10^{-8}\), \(9.70\times10^{-6}\), and \(2.83\times10^{-6}\). The exact stationary weak-row and bottom-tangency defects are \(2.07\times10^{-7}\) and \(4.07\times10^{-8}\). Padding factors two and three give the same \(2.33\times10^{-4}\) wave-projector discrepancy, and nested horizontal support stabilizes at that value. Thus neither dealiasing nor missing horizontal guard modes explains the stop.

The outcome is **`geostrophic-modal-blocker`**. The fixed-\(\kappa\) internal-wave coordinates and isolated InternalModesEVP provider are not the failure. The ordinary flat APV-bearing geostrophic modes plus the explicit bottom inversion converge systematically, but they approximate the exact terrain-dependent stationary inclusion too slowly: the factor-two economy gate is lost before the stationary or wave projector reaches tolerance. The two-dimensional control is not attempted because the zonal acceleration gate does not pass. Milestone 10.3 remains inactive; a future increment must construct a more terrain-adapted balanced/stationary coordinate family while retaining the successful fixed-\(\kappa\) wave sector and the primitive oracle.

The focused Milestone-10.2.3 suite passes six tests. The complete repository suite passes 224 tests with zero failures, and static analysis reports no issues in all 107 MATLAB files.

## Milestone 10.2.4: Terrain-dressed Robin stationary-coordinate oracle

- [x] Complete — blocking stationary-compression gate

### Purpose

Determine whether bottom-localized Robin geostrophic modes become an economical representation of the exact finite-terrain stationary space when combined with the global first terrain correction verified in Milestones 6.8 and 9.4. This experiment retains the successful fixed-$\kappa$ internal-wave sector from Milestone 10.2.3 and isolates the failed balanced-coordinate sector.

For scalar Robin columns $\Phi_{\ell_b}$, define the globally dressed geostrophic coordinates

```math
X_{\rm g}^{[1]}(\ell_b,\delta)
=
\left(G_0+\delta G_1\right)\Phi_{\ell_b},
```

and the candidate ambient space

```math
\mathcal V_{\ell_b}^{[1]}
=
\mathcal W_{\kappa}^{\rm flat}
\oplus
X_{\rm g}^{[1]}(\ell_b,\delta)
\oplus
\mathcal B_{\kappa}^{0}
\oplus
\mathcal M_0.
```

Here $\mathcal W_{\kappa}^{\rm flat}$ is the qualified fixed-$\kappa$ wave family, $G_1$ retains every terrain-generated horizontal sideband, $\mathcal B_{\kappa}^{0}$ remains the independent zero-APV bottom inversion, and $\mathcal M_0$ is the compatible mean sector. The Robin condition shapes vertical convergence only. The exact finite-terrain $H_\gamma,J_\gamma,Q_\gamma,B,R_h$ forms and the exact Milestone-8 stationary space remain the physical evolution and validation oracles.

The implementation refines this schematic expression to preserve admissibility. If $T_{\ell_b}$ spans the bottom-tangent scalar combinations and $C_{\ell_b}$ is its complementary scalar subspace,

```math
R_1G_0\Phi_{\ell_b}T_{\ell_b}=0,
```

then the equal-dimension balanced ambient is

```math
X_{\rm g}^{[1]}
=
\left[
(G_0+\delta G_1)\Phi_{\ell_b}T_{\ell_b}
\quad
G_0\Phi_{\ell_b}C_{\ell_b}
\right].
```

Only the tangent combinations may receive the stationary $G_1$ correction: applying $G_0+\delta G_1$ to a non-tangent scalar column creates nonzero mapped normal velocity at the bottom. The complementary Robin directions remain admissible flat coordinates and are retained for later dynamical classification; none is removed.

### Dependencies

Milestone 10.2.3 with outcome `geostrophic-modal-blocker`, the Milestone-6.8 $G_0,G_1$ construction, the Milestone-8 exact finite-terrain stationary space, the Milestone-9.2 qualified signed-Robin provider, and the Milestone-10.1 production coordinate contract. The isolated InternalModesEVP checkout remains pinned to commit `df86687` and neither InternalModes checkout may be modified.

### Training and ablations

Use the existing constant-$N$, zonal sinusoidal-terrain case with padding factor two, primitive reference degree 18, geostrophic count eight, and

```math
\frac{\ell_b}{D}
\in
\left\{
-\frac18,-\frac14,-\frac12,-1,\infty
\right\}.
```

For every $\ell_b$, compare at equal retained dimension:

1. ordinary flat geostrophic modes $G_0\Phi_\infty$;
2. undressed Robin modes $G_0\Phi_{\ell_b}$;
3. globally dressed Robin modes $(G_0+\delta G_1)\Phi_{\ell_b}$.

Select $\ell_b$ using the smallest exact stationary-projector representation defect, subject to provider, endpoint, conditioning, conjugacy, Green-identity, and bottom-tangency gates. Break numerical ties using the smaller internal-wave projector defect and then the smaller physical-energy Gram condition number. Freeze the selected $\ell_b$ before every validation calculation. Cross-$\ell_b$ variation is a finite-order convergence diagnostic rather than a physical invariance requirement.

Retain the explicit zero-APV bottom inversion separately for every $\ell_b$. A Robin eigenfunction must not be identified with, substituted for, or allowed to duplicate that independent boundary coordinate.

### Oracle

The planned public diagnostic is

```matlab
audit = problem.auditTerrainDressedRobinStationaryAmbient( ...
    trustedModeBounds=[1 0], ...
    supportModeBounds=[1 2;1 3;1 4], ...
    targetWaveModeIndices=[1;2], ...
    geostrophicModeCounts=[2;4;6;8;12], ...
    robinLengthRatios=[-1/8;-1/4;-1/2;-1;Inf], ...
    trainingGeostrophicModeCount=8, ...
    primitiveReferenceDegrees=[16;18;20], ...
    paddingFactors=[2;3], ...
    terrainScales=[0;1/8;1/4;1/2;3/4;1], ...
    internalModesEVPOrders=[128;256], ...
    shouldRunTwoDimensionalControl=true, ...
    cacheDirectory="output/milestone-10.2.4-cache");
```

Training uses only the declared padding-two, degree-18, count-eight case. Validation uses the frozen $\ell_b$ while varying modal count, primitive reference degree, horizontal support, padding, and terrain amplitude. The two-dimensional trusted-band control runs only after the zonal acceleration gate passes.

Apply $G_1$ globally with the common dealiased terrain projection and retain every generated sideband inside the declared support. Compare the dressed span with the exact Milestone-8 stationary projector without appending missing stationary directions. Attach the unchanged fixed-$\kappa$ wave family and measure the resulting internal-wave projector, frequency, APV, bottom-evolution, and strong primitive residuals. Reuse signature-compatible Milestone-10.2.3 primitive-oracle cache artifacts when available; every new artifact remains disposable, content addressed, signature validated, opt-in, and restricted to ignored `output/`.

No frequency cutoff, replacement APV row, empirical correction, post hoc symmetrization, APV-nullspace projection, mode deletion, dense-oracle construction coordinate, matrix-free operator, or time integration is included.

### Automated acceptance

- Robin provider and endpoint defects are below $10^{-10}$.
- Analytic and centered $G_1$ actions agree below $10^{-9}$.
- Fourier selection, conjugacy, and the Hermitian/skew-Hermitian finite-terrain structure close below $10^{-11}$.
- Tangent Green identities and stationary weak rows close below $10^{-10}$.
- The frozen dressed basis represents the exact stationary projector within $10^{-8}$.
- Stationary bottom tangency and stationary weak-row residuals are below $10^{-10}$.
- Internal-wave projector and frequency defects are below $10^{-8}$.
- Projected APV, bottom-evolution, and strong primitive residuals are below $10^{-8}$, $10^{-10}$, and $10^{-5}$.
- Padding, support, and primitive-reference-degree projector drifts are below $10^{-8}$.
- The accepted ambient space retains at least factor-two compression relative to the primitive reference.
- Report separately the improvement from replacing the ordinary flat modes by Robin modes and the additional improvement from applying $G_1$.

Classify the result as:

- **`dressed-robin-acceleration`:** a finite negative $\ell_b$ passes every physical and economy gate;
- **`dressed-flat-acceleration`:** $\ell_b=\infty$ passes after global dressing, showing that Robin localization is unnecessary for the tested production basis;
- **`dressed-stationary-equivalent`:** every physical gate passes but factor-two compression does not;
- **`dressed-stationary-seed`:** global dressing materially improves stationary and wave projectors but does not pass the finite-amplitude physical gates;
- **`dressed-stationary-blocker`:** neither Robin localization nor global first-order dressing gives a convergent stationary representation.

Under the original Milestone-10.2.4 gate, only `dressed-robin-acceleration` or `dressed-flat-acceleration` would have authorized Milestone 10.3. Neither occurred; the measured blocker routes the roadmap to Milestone 10.2.5.

### Outcome

The declared constant-$N$ zonal training oracle used resolution `[6 14 5]`, support `[1 4]`, primitive reference degree 18, geostrophic count eight, padding factor two, and isolated provider orders 128 and 256. No Robin length passed the prerequisite provider, centered-inclusion, conjugacy, APV Green-identity, and stationary-compression gates, so no value of $\ell_b$ was frozen and the independent finite-amplitude validation sweep was correctly skipped.

The dressed stationary results are:

| $\ell_b/D$ | Stationary representation | Stationary weak row | Bottom tangency | APV Green identity | Gram condition number |
|---:|---:|---:|---:|---:|---:|
| $-1/8$ | $1.66\times10^{-2}$ | $7.48\times10^{-5}$ | $6.54\times10^{-8}$ | $4.48\times10^{-3}$ | $2.89\times10^{12}$ |
| $-1/4$ | $1.39\times10^{-2}$ | $1.36\times10^{-4}$ | $1.23\times10^{-7}$ | $2.12\times10^{-3}$ | $5.09\times10^{12}$ |
| $-1/2$ | $1.26\times10^{-2}$ | $1.54\times10^{-4}$ | $1.41\times10^{-7}$ | $1.04\times10^{-3}$ | $3.54\times10^{13}$ |
| $-1$ | $1.20\times10^{-2}$ | $1.61\times10^{-4}$ | $1.48\times10^{-7}$ | $5.25\times10^{-4}$ | $1.33\times10^{15}$ |
| $\infty$ | $1.26\times10^{-2}$ | $1.68\times10^{-4}$ | $1.59\times10^{-7}$ | $5.35\times10^{-5}$ | $1.99\times10^{7}$ |

The tangent scalar combinations themselves satisfy bottom tangency below $3\times10^{-20}$. Thus the failure is not caused by dressing inadmissible columns. Relative to the undressed coordinates, $G_1$ lowers the stationary weak-row and bottom defects by factors between approximately $2.5$ and $5.6$, but it does not materially lower the exact stationary-projector error. Finite negative Robin lengths give at most a small projector improvement over the ordinary value $1.26\times10^{-2}$ while worsening the physical-energy Gram conditioning by five to eight orders of magnitude. Their provider and endpoint comparisons are $1.05\times10^{-10}$--$2.75\times10^{-10}$, and analytic-versus-centered $G_1$ defects are approximately $2.0\times10^{-9}$; these narrow numerical misses do not explain the much larger $10^{-2}$--$10^{-3}$ physical defects.

The outcome is **`dressed-stationary-blocker`**. Bottom localization and the first stationary terrain derivative improve selected residuals, but a count-eight Robin family still does not economically span the exact finite-terrain stationary inclusion. Milestone 10.3 remains inactive. The next formulation, if pursued, must construct the stationary coordinates from the global tangent constraint itself rather than dress a truncated collection of independently generated per-wavenumber Robin functions.

The focused Milestone-10.2.4 suite passes six tests. The complete repository suite passes 230 tests with zero failures; static analysis covers all 109 MATLAB files.

## Milestone 10.2.5: Global tangent-scalar stationary compression oracle

- [ ] Complete — blocking stationary-compression gate

### Purpose

Reverse the unsuccessful per-wavenumber construction order:

```math
\boxed{
\text{enforce the global terrain constraint first}
\quad\longrightarrow\quad
\text{compress its vertical dependence second}.
}
```

The exact stationary scalar space will be represented as

```math
\mathcal S_{\gamma,N}
=
\mathcal S_{0,N}
\oplus
\mathcal L_{\gamma,N}\mathcal C_{h,N},
```

where the zero-bottom-trace sector carries the interior APV structure and the second sector contains energy-minimal lifts of globally terrain-tangent bottom traces. The successful fixed-$\kappa$ waves, explicit zero-APV bottom inversion, compatible mean sector, public coefficient layout, and exact finite-terrain forms remain unchanged.

### Dependencies

Milestone 10.2.4 with outcome `dressed-stationary-blocker`, the Milestone-7 exact finite-terrain primitive forms, the Milestone-8 polynomial stationary oracle, the Milestone-10.1 production contract, and the qualified isolated InternalModesEVP checkout pinned to commit `df86687`. The polynomial oracle remains the independent validation space and may not be inserted as a construction coordinate.

### Global trace factorization

Construct the projected bottom-trace operator with the common dealiased product:

```math
\mathcal T_{h,N}b
\equiv
P_{M\to N}
\left[
(\partial_xb)h_y-(\partial_yb)h_x
\right].
```

Use a conjugate-consistent singular-value decomposition to obtain the gauge-free tangent kernel and its retained complement,

```math
\mathcal C_{h,N}
\equiv
\ker\mathcal T_{h,N},
\qquad
\mathcal C_{h,N}^{\perp}.
```

For each $b\in\mathcal C_{h,N}$, construct the lift $\mathcal L_{\gamma,N}b$ with trace $b$ and finite-terrain-energy orthogonality to every zero-bottom-trace geostrophic state. By the APV Green identity, this is the energy-minimal projected-zero-volume-APV representative of the trace. Retain the lifted non-tangent traces in $\mathcal C_{h,N}^{\perp}$ as complementary coordinates for later dynamical classification; do not delete them or call them stationary.

### Two-stage oracle

1. **Polynomial factorization control**
   - Factor the existing global polynomial stationary scalar space into zero-bottom-trace and lifted tangent-trace sectors.
   - Reconstruct the Milestone-8 stationary projector without adding or deleting a direction.
   - Preserve every non-tangent trace direction in a separately identified complementary sector.
   - Verify the factorization, trace, stationary Green, weak-row, bottom-tangency, conjugacy, and dimension identities before any modal compression is attempted.

2. **Independent modal compression**
   - Generate mixed Dirichlet--Neumann interior modes with the isolated InternalModesEVP checkout pinned to `df86687`:

     ```math
     \phi_b=0,
     \qquad
     \partial_\xi\phi(0)=0.
     ```

   - Use the public explicit zero-APV bottom inversion as the boundary pivot. Introduce no additional bottom degree of freedom.
   - Build the tangent and non-tangent lifts in the retained zero-trace modal space, map them through the exact terrain geostrophic-state construction, and physical-energy orthogonalize the result.
   - Append the verified fixed-$\kappa$ wave family and compatible $\kappa=0$ sector.
   - Train only the smallest acceptable interior geostrophic-mode count on the declared constant-$N$ zonal case. Freeze that count before independent changes in support, padding, primitive reference degree, terrain amplitude, stratification, and terrain dimensionality.
   - Compare with the polynomial stationary and dynamical projectors without inserting dense polynomial stationary vectors into the modal candidate.

### Public audit

```matlab
audit = problem.auditGlobalTangentStationaryCompression( ...
    trustedModeBounds=[1 0], ...
    supportModeBounds=[1 2;1 3;1 4], ...
    targetWaveModeIndices=[1;2], ...
    interiorGeostrophicModeCounts=[2;4;6;8;12], ...
    stationaryReferenceDegree=4, ...
    primitiveReferenceDegrees=[16;18;20], ...
    paddingFactors=[2;3], ...
    terrainScales=[0;1/8;1/4;1/2;3/4;1], ...
    internalModesEVPOrders=[128;256], ...
    shouldRunTwoDimensionalControl=true, ...
    cacheDirectory="output/milestone-10.2.5-cache");
```

The audit must make the training record, frozen modal count, independent validation cases, polynomial-factorization result, stationary and complementary trace dimensions, and full physical degree-of-freedom accounting explicit. Cache artifacts are disposable, content addressed, signature validated, opt-in, and confined to ignored `output/`; disabling caching performs no hidden disk writes.

### Automated acceptance

- Polynomial factorization, projector reconstruction, and complete dimension accounting close below $10^{-12}$.
- The trace-nullspace rank is stable under the declared support and padding changes, and Fourier conjugacy closes below $10^{-11}$.
- The analytic zonal-ridge tangent-trace kernel is recovered below $10^{-12}$.
- The isolated mixed-boundary provider closes its endpoint, physical-energy, and potential-enstrophy checks below $10^{-11}$ at qualified orders 128 and 256.
- Lift trace and finite-terrain energy-minimization defects are below $10^{-10}$.
- Every tangent lift has projected volume APV below $10^{-10}$.
- The compressed stationary projector agrees with the independent polynomial stationary projector below $10^{-8}$.
- Stationary APV Green, weak-row, and bottom-tangency defects are below $10^{-10}$.
- The internal-wave physical-energy projector and frequency defects are below $10^{-8}$.
- Internal-wave APV, bottom-evolution, and strong primitive residuals are below $10^{-8}$, $10^{-10}$, and $10^{-5}$.
- Padding, support, primitive-reference-degree, and terrain-amplitude projector drifts are below $10^{-8}$ after the modal count is frozen.
- The flat explicit bottom-inversion projector is recovered without duplication.
- Stationary and non-tangent trace sectors are $H_\gamma$-orthogonal and account for the declared trace dimension within $10^{-10}$.
- The accepted modal space retains at least factor-two compression relative to the primitive ambient reference.

### Outcome classification

- **`global-tangent-acceleration`:** every physical gate passes with at least factor-two compression;
- **`global-tangent-equivalent`:** every physical gate passes without factor-two compression;
- **`global-tangent-seed`:** the polynomial factorization passes and the modal construction materially improves convergence, but the final physical gates do not all pass;
- **`global-tangent-modal-blocker`:** the exact polynomial factorization passes, but the mixed-boundary modes or their lifts do not converge;
- **`global-tangent-factorization-blocker`:** the polynomial stationary space does not admit the proposed finite factorization within tolerance;
- **`global-tangent-provider-blocker`:** the isolated one-dimensional provider cannot represent or qualify the required mixed endpoint problem.

Only `global-tangent-acceleration` authorizes Milestone 10.3. Every other outcome stops for analysis or repair. No bottom- or topographic-wave classification, frequency cutoff, replacement APV row, empirical correction, post hoc symmetrization, APV-nullspace projection, mode deletion, matrix-free work, or time integration is included.

## Milestone 10.3: Complete bottom and topographic-wave sector

- [ ] Complete — blocking gate

### Purpose

Partition every declared flat bottom and zero-frequency coordinate into the exact finite-terrain stationary space or a converged dynamical topographic-wave space.

### Dependencies

Milestone 10.2.5 with outcome `global-tangent-acceleration`. Milestone 10.2.3 established the fixed-$\kappa$ wave sector, while Milestone 10.2.4 established that independently truncated Robin families do not economically represent the stationary space; neither result by itself satisfies this dependency.

### Deliverables

- Begin with the complete physical flat zero-frequency block: APV-bearing balanced states, explicit zero-APV bottom inversions, and MDA states.
- Construct the exact finite-terrain bottom-tangent stationary inclusion $G_{\gamma,N}$.
- Define the $H_\gamma$-orthogonal non-tangent bottom complement without using a frequency threshold.
- Use the verified first-order zero-block splitting only as the preconditioner for global dressing and exact residual enrichment of the complete non-tangent block.
- Continue the resulting invariant subspaces in terrain amplitude and report backward uncertainty, bottom participation, projected APV, bottom evolution, strong equations, padding, guard support, and the $h\to0$ limit.
- Retain complete degenerate pairs and Fourier conjugates throughout. Do not promote or delete individual small-frequency eigenvectors.

### Automated acceptance

- Every declared bottom and zero-frequency coordinate belongs to either the exact stationary space or the converged topographic-wave space.
- Stationary and dynamical-bottom projectors are $H_\gamma$-orthogonal, conjugate closed, and dimensionally complete within $10^{-10}$.
- Every claimed topographic-wave frequency exceeds its backward uncertainty by at least $10^3$.
- Accepted topographic subspaces satisfy the Milestone-9.1 APV, bottom, strong-equation, guard, padding, bottom-participation, and terrain-amplitude gates.
- The topographic-wave projector approaches the appropriate flat stationary bottom block continuously as $h\to0$.
- No residual small-frequency or unclassified bottom direction remains in the production state.

Failure blocks assembly of the complete forward basis.

## Milestone 10.4: Complete production-basis gate

- [ ] Complete — principal blocking gate

### Purpose

Combine the stationary, internal-wave, and topographic-wave sectors into one complete and economical state for arbitrary linear forward evolution.

### Dependencies

Milestone 10.3.

### Deliverables

- Assemble

  ```math
  \mathcal V_{\rm prod}
  =
  \mathcal G_\gamma
  \mathbin{\oplus_{H_\gamma}}
  \mathcal W_{\rm int}
  \mathbin{\oplus_{H_\gamma}}
  \mathcal W_{\rm topo}.
  ```

- Require the production projectors to satisfy

  ```math
  P_{\rm g}
  +
  P_{\rm int}
  +
  P_{\rm topo}
  =
  I_{\rm prod}.
  ```

- Compare the complete production projector with the trusted physical part of the overresolved primitive oracle:

  ```math
  \left\|
  (I-P_{\rm prod})
  P_{\rm oracle,trusted}
  \right\|_{H_\gamma}
  \longrightarrow0.
  ```

- Test arbitrary conjugate-symmetric production states, not only individual eigenmodes.
- Report deliberately omitted oracle energy as spectral truncation error rather than unresolved production energy.
- Record the selected horizontal scattering order, vertical support, physical-energy tail estimate, spectral-isolation margin, and unresolved-energy fraction with the production basis.
- Diagonalize the complete reduced physical-energy pencil and record the stationary, internal, and topographic mode counts.

### Automated acceptance and exits

- The three physical projectors are pairwise $H_\gamma$-orthogonal, conjugate closed, and complete within $10^{-10}$.
- Arbitrary production states reconstruct in the primitive oracle and return to production coordinates within $10^{-10}$.
- The complete terrain eigensystem has real frequencies, $H_\gamma$-orthogonality, projected APV, bottom evolution, and strong residuals at the existing physical tolerances.
- Production observables and projectors converge independently with wave count, APV count, horizontal support, guard width, padding, and primitive-oracle degree.
- The production state uses no unresolved coordinate and achieves at least factor-two dimension reduction relative to the primitive oracle at matched trusted-band accuracy.

Classify the result as:

- `complete-basis-acceleration` — every gate passes with at least factor-two reduction;
- `complete-basis-equivalent` — physical completeness passes but the economy target fails;
- `complete-basis-blocker` — dimensional or physical completeness fails.

Only `complete-basis-acceleration` authorizes Milestone 11.

## Milestone 11: Matrix-free production basis

- [ ] Complete

### Purpose

Replace dense production-basis construction and operator application with adjoint-consistent matrix-free algorithms without returning to the primitive ambient state online.

### Dependencies

Milestone 10.4 with outcome `complete-basis-acceleration`.

### Public interface

Introduce an immutable basis object:

```matlab
basis = WVTerrainEnergyBasis.fromGalerkin(problem, ...
    trustedModeBounds=trustedModeBounds, ...
    waveModeIndices=waveModeIndices, ...
    apvModeIndices=apvModeIndices, ...
    guardModeBounds=guardModeBounds, ...
    paddingFactor=paddingFactor);
```

`WVTerrainEnergyGalerkin` remains the scientific construction and oracle object. `WVTerrainEnergyBasis` stores the production layout, exact reduced forms, terrain modes, family projectors, reconstruction maps, selected scattering and vertical supports, physical-energy tail estimate, spectral-isolation margin, unresolved-energy fraction, and construction version.

### Deliverables

- Apply $H_\gamma,J_\gamma,Q_\gamma,B,R_h$, and the stationary inclusion through field reconstruction, common dealiased terrain multiplication, and adjoint restriction.
- Reuse global dressing, recycled residual-enrichment spaces, and complete resonant blocks during construction.
- Use dense reduced diagonalization at oracle sizes and symmetry/Bloch-block iterative eigensolves at production sizes.
- Preserve explicit bottom coordinates throughout packing, application, projection, and reconstruction.
- Keep diagnostic pressure recovery outside ordinary operator applications.

### Automated acceptance

- Matrix-free actions reproduce the complete dense production oracle within $10^{-10}$.
- Hermitian and skew-Hermitian identities close within $10^{-12}$.
- Stationary, internal-wave, topographic-wave, and complete-basis projectors agree with Milestone 10.4 within $10^{-9}$.
- Matrix-free truncation metadata reproduces the Milestone-10.4 scattering support, vertical support, tail estimate, and isolation margin.
- Matrix-free construction retains the Milestone-10.4 dimension reduction and performs no global primitive-matrix assembly online.
- No diagnostic pressure solve occurs during an operator application.

## Milestone 12: Forward wave–vortex evolution

- [ ] Complete

### Purpose

Evolve arbitrary states in the complete production basis.

### Dependencies

Milestone 11.

### Public interface

Add a separate mutable model:

```matlab
model = WVTerrainEnergyModel.fromBasis(basis, ...
    Ap=Ap, ...
    Am=Am, ...
    A0=A0, ...
    bottomDisplacement=etaB, ...
    time=t0);
```

The basis and Galerkin construction remain immutable. The model owns time, production coefficients, integration configuration, diagnostics, and output callbacks.

### Deliverables

- Support complete terrain-mode phase evolution as the primary path.
- Support the physical-energy Cayley step

  ```math
  \left(
  H_\gamma-\frac{\Delta t}{2}J_\gamma
  \right)
  \boldsymbol A^{n+1}
  =
  \left(
  H_\gamma+\frac{\Delta t}{2}J_\gamma
  \right)
  \boldsymbol A^n
  ```

  as the validation and fallback path.
- Convert to and from WaveVortexModel `Ap`, `Am`, `A0`, and bottom-displacement conventions.
- Reconstruct physical fields and expose stationary, internal-wave, topographic-wave, physical-energy, APV, bottom-displacement, and truncation diagnostics.
- Keep pressure recovery diagnostic and absent from ordinary evolution.

### Automated acceptance

- Phase evolution, Cayley evolution, and dense matrix exponentiation agree within $10^{-10}$ at reference resolution.
- Arbitrary production states round-trip through native WaveVortexModel coefficients and reconstructed fields within $10^{-10}$.
- Physical energy is conserved to solver tolerance, projected APV remains stationary, and bottom evolution closes.
- Stationary coefficients remain constant and internal/topographic mode energies agree between evolution paths.
- Ordinary phase and Cayley evolution perform no pressure solve.

## Milestone 13: Forward scientific benchmarks and examples

- [ ] Complete

### Purpose

Demonstrate that the complete forward model reproduces the expected stationary, internal-wave, topographic-wave, and scattering dynamics.

### Dependencies

Milestone 12.

### Deliverables

- Add flat and uniform-depth recovery benchmarks.
- Add stationary APV-bearing and bottom-buoyancy examples.
- Add an initialized topographic boundary-wave example.
- Add sinusoidal-terrain internal-wave scattering and Gaussian-ridge packet examples.
- Show physical fields, modal-family energy transfers, bottom displacement, projected APV, physical energy, and truncation error.
- Compare the production model with the dense primitive oracle at reduced resolution and compare phase with Cayley evolution.

### Automated acceptance

- Flat and uniform-depth forward solutions recover the corresponding WaveVortexModel solutions within $10^{-10}$.
- Stationary examples remain stationary; topographic-wave examples retain their converged invariant subspaces.
- Sinusoidal and Gaussian calculations converge with production-mode count, support, guard width, padding, and time step.
- Physical energy, projected APV, and bottom evolution retain their Milestone-12 gates.
- Every example reports convergent truncation error and uses no unresolved production subspace.

## Milestone 14: Research-production behavior

- [ ] Complete

### Purpose

Extend the validated linear forward model to repeatable research calculations without changing its scientific definition.

### Dependencies

Milestone 13.

### Deliverables

- Add arbitrary stationary stratification, broadband terrain, resolution rebuilding, deterministic basis reconstruction, restartable NetCDF output, and symmetry or Bloch decomposition.
- Persist the physical truncation, terrain basis, family classification, bottom-coordinate convention, construction version, and evolution method.
- Add construction, eigensolve, reconstruction, phase-evolution, and Cayley-evolution profiling at three resolutions.
- Retain the pressure-free online path and the common dealiased terrain-product convention.

### Automated acceptance

- Repeated construction is deterministic and resolution rebuilding preserves conjugacy, physical normalization, family dimensions, projectors, APV classification, and bottom participation.
- Restart continuation matches uninterrupted evolution within $10^{-10}$ in physical-energy-normalized coefficients.
- Broadband and variable-stratification calculations retain every Milestone-10.4 through Milestone-13 gate.
- Benchmarks report construction time, peak stored state, operator-application time, eigensolver iterations, reconstruction time, and online evolution time.
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
| **D3.2 — Dense terrain eigensystem** | 9 | Solve the complete physical-energy eigenproblem, preserve every direction, and expose candidate physical and unresolved subspaces. |
| **D3.3 — Physical-subspace classification** | 9.1 | Establish converged stationary and dynamical physical projectors, quantify the unresolved remainder, and stop before vertical compression or residual enrichment. |
| **D3.4 — Boundary-complete modal compression** | 9.2 | Compare the Robin, zero-APV-bottom, and fixed-$\kappa$ modal coordinates with the primitive oracle, select the Milestone-10 seed, and stop before residual enrichment. |
| **D3.5 — Slope-compatible wave coordinates** | 9.3 | Train and freeze hydrostatic and nonhydrostatic active-bottom reference families, validate them independently against the primitive oracle, select the Milestone-10 seed, and stop before residual enrichment. |
| **D3.6 — Global terrain dressing** | 9.4 | Construct the complete global $O(h)$ internal and zero-frequency block corrections, verify $O(h^2)$ weak and bottom residuals, classify the seed, and stop before repeated enrichment. |
| **E1 — Selective modes** | 10 | Implement residual enrichment and validate it exclusively against the dense oracle. |
| **E1.1 — Production-state contract** | 10.1 | Establish a complete physical coordinate set with exact degree-of-freedom accounting. |
| **E1.2 — Complete internal waves** | 10.2 | Continue and enrich every declared internal-wave block. |
| **E1.2.1 — Coupled-block internal waves** | 10.2.1 | Exact correction passed at fixed resolution but exposed a nested vertical-isolation blocker; the iterative stage was not attempted. |
| **E1.2.2 — Geometric-cascade isolation** | 10.2.2 | Resolve horizontal and vertical terrain-scattering tails independently, test spectral isolation of the complete internal-wave block, and stop before bottom-wave classification. |
| **E1.2.3 — Flat wave–vortex modal ambient** | 10.2.3 | Replace generic polynomial candidate coordinates by fixed-\(\kappa\) waves, ordinary flat geostrophic modes, explicit bottom inversions, and the compatible mean sector; stop if the exact terrain stationary space is not represented economically. |
| **E1.2.4 — Terrain-dressed stationary coordinates** | 10.2.4 | Train and freeze the Robin length, apply the global $G_0+\delta G_1$ stationary dressing, and stop unless the exact stationary and internal-wave projectors pass with factor-two compression. |
| **E1.2.5 — Global tangent stationary compression** | 10.2.5 | Factor the global terrain-tangent trace space before vertical compression, validate the modal lifts against the polynomial stationary oracle, and stop unless every physical gate passes with factor-two compression. |
| **E1.3 — Bottom/topographic sector** | 10.3 | Partition every bottom coordinate into the exact stationary or converged topographic-wave space. |
| **E1.4 — Complete basis gate** | 10.4 | Prove complete, economical coverage of the declared production state. |
| **E2 — Matrix-free production basis** | 11 | Reproduce the complete dense production oracle without global primitive matrices. |
| **F1 — Forward evolution** | 12 | Evolve arbitrary production states with phase and Cayley methods. |
| **F2 — Scientific examples** | 13 | Validate the complete forward model in physical benchmarks. |
| **G — Research production** | 14 | Add persistence, rebuilding, broadband terrain, arbitrary stratification, and profiling. |

Each batch is suitable for a separate Codex goal of the form: “Implement Batch E1.1; do not proceed to Batch E1.2; continue until all acceptance criteria pass or a genuine scientific or numerical blocker is established.” Substitute the batch being pursued and its immediate successor.

Use one focused commit per completed milestone and retain acceptance evidence in the automated tests and examples. Do not bypass a blocking gate by symmetrizing an incompatible operator, fitting a closure, or introducing an empirical correction.

## Definition of done

The continuum and dense-oracle scientific proof of concept is established by Milestones 5–10: the boundary coordinate is represented as part of a complete linked state, the projected primitive system preserves physical energy, derives stationary volume APV, converges to the strong bottom equation, constructs the complete stationary balanced space, and demonstrates that global dressing followed by exact residual enrichment efficiently recovers a selected internal-wave subspace.

Scientific completeness of the production basis requires Milestones 10.1–10.4, including the intervening Milestone-10.2.1 coupled-block, Milestone-10.2.2 geometric-cascade, Milestone-10.2.3 flat wave–vortex ambient, Milestone-10.2.4 terrain-dressed stationary-coordinate, and Milestone-10.2.5 global tangent-scalar gates. The exact stationary space must be represented economically before bottom/topographic-wave classification begins. Every declared wave, APV, bottom, and MDA coordinate must belong to the stationary, internal-wave, or topographic-wave projector, with no unresolved coordinate inside the production state. The selected horizontal scattering support, vertical support, modal family and count, tangent-trace rank, frozen interior-mode count, omitted-tail estimate, and unresolved-energy fraction must be recorded with the production basis. Guard modes and omitted primitive-oracle directions are numerical support and truncation diagnostics rather than prognostic degrees of freedom.

A complete forward wave–vortex model requires Milestones 11–13: matrix-free construction of the complete production basis, arbitrary-state phase and Cayley evolution, and convergent stationary, topographic-wave, sinusoidal-terrain, and Gaussian-ridge benchmarks.

Milestone 14 completes research-production behavior through arbitrary stationary stratification, broadband terrain, resolution rebuilding, restartable output, symmetry decomposition, and profiling.

Nonlinear wave–wave interactions, nonlinear terrain dynamics, an additive terrain `WVForcing`, empirical closure, an MPM release, modifications to WaveVortexModel, independent surface buoyancy, and dynamic barotropic backreaction remain outside this roadmap.
