# Terrain-energy Galerkin milestones

## Objective

Develop `WVTerrainEnergyGalerkin`, a pressure-free Galerkin system for the linear rotating Boussinesq equations over stationary bottom topography. The formulation is linear in flow amplitude and exact in the resolved terrain. Its raw Galerkin form conserves the discrete finite-terrain energy; the approved closed evolution must also conserve the discrete quadratic potential enstrophy and enforce the resolved bottom evolution.

The semidiscrete evolution is

```math
E_\gamma\dot{\boldsymbol a}=J_\gamma\boldsymbol a,
\qquad
E_\gamma=E_\gamma^*>0,
\qquad
J_\gamma=-J_\gamma^*.
```

The mathematical specification is `terrain-energy-galerkin.tex` at commit `72c967c` in the `ape-apv-bottom-topography` literature repository. The implementation will use ordinary hydrostatic modes as economical vertical coordinates, retain one explicit bottom-displacement coefficient per horizontal wavenumber, recover the complete flat nonhydrostatic modes through weak diagonalization, and construct finite-terrain modes through a generalized Hermitian eigensolve and residual-based enrichment.

The completed mean-depth generator and scattering implementation is retained as reusable engineering infrastructure. Its scientific roadmap is archived in [mean-depth-wave-generator-milestones.md](mean-depth-wave-generator-milestones.md). Neither existing forcing is used as the terrain-energy evolution operator.

The dense finite-terrain forms are the unmodified Galerkin oracle. Milestone 4.5 records that requiring every resolved APV sample to remain fixed is incompatible with the strong bottom evolution in the present finite-dimensional state space. Milestone 4.6 therefore tests the physically weaker quadratic requirement: exact conservation of total discrete quadratic potential enstrophy together with energy and strong bottom evolution. Any adopted closure is a conservative discrete surrogate for the raw finite-dimensional Galerkin operator; its correction magnitude and convergence must remain explicit diagnostics.

## Fixed conventions and boundaries

- The public target is `WVTransformBoussinesq`; WaveVortexModel remains an external dependency and is not modified.
- `topographicHeight` is a real, finite, stationary, upward-positive, horizontally periodic field on the transform grid.
- The mapped geometry is

```math
\gamma=1-\frac{h}{D}>0.
```

- The prognostic weak state contains $(\hat u,\hat v,\hat w,\hat\eta)$ and one bottom-displacement coefficient per retained horizontal wavenumber. Pressure is recovered only as a post-solve diagnostic.
- The velocity satisfies mapped continuity and homogeneous normal-flow conditions,

```math
\nabla_\xi\boldsymbol{\cdot}\hat{\boldsymbol u}=0,
\qquad
\hat w(-D)=\hat w(0)=0.
```

- Surface displacement vanishes. The bottom value of $\hat\eta$ is unrestricted and evolves through the displacement equation.
- Independent surface-buoyancy anomalies are excluded. Independent bottom buoyancy is retained in the zero-frequency balanced state space.
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

## Milestone 2: Mixed hydrostatic and bottom-displacement basis

- [x] Complete

### Purpose

Construct the complete coordinate space required by the finite-terrain weak equations before forming any dynamical matrices.

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

## Milestone 3: Flat nonhydrostatic dense oracle

- [x] Complete — blocking scientific gate

### Purpose

Prove that the mixed hydrostatic coordinates recover the complete flat nonhydrostatic wave–vortex problem before introducing terrain.

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
- Nonzero-frequency modes have negligible discrete APV and the zero-frequency dimension matches the expected balanced, MDA, and bottom-buoyancy coordinates.
- Failure to recover the flat problem or acceptable conditioning blocks Milestone 4.

## Milestone 4: Exact finite-terrain dense forms

- [x] Complete

### Purpose

Establish the finite-terrain weak system, exact in the resolved terrain height $h$, as a low-resolution dense scientific oracle.

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

- [x] Complete — incompatible exit; Milestone 5 remains blocked

### Purpose

Determine whether the current complete mixed state space admits an energy-preserving evolution that also conserves discrete finite-terrain APV and enforces the strong bottom-displacement evolution exactly. Do not assume that these constraints are jointly feasible.

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
- Test the constraints on the complete state space, including the APV-bearing balanced and bottom-buoyancy coordinates. Do not obtain feasibility by silently deleting balanced columns or restricting the audit to a selected wave subspace.

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
  - **Feasible:** the exact constraints close at the stated tolerances, the correction is reproducible, and Milestone 5 proceeds using $J_\gamma^c$ while retaining the raw operator as a diagnostic reference.
  - **Incompatible:** the rank-revealing solve establishes a nonzero minimum constraint residual, identifies the failed compatibility condition and implicated state subspace, and Milestone 5 remains blocked pending a revised state or constraint formulation.
- Approximate feasibility obtained by relaxed tolerances, empirical correction factors, or removal of physical state coordinates does not pass this gate.

### Outcome

The dense audit preserves the complete mixed state and constructs the resolved bottom row by projecting the oversampled physical product into the retained bottom Fourier space. Flat and uniform-depth cases are admissible: the APV nullity equals the 56-dimensional flat wave space in the automated reference problem, all three constraints close at roundoff, and the minimum-change closure is the raw generator.

The 20 m sinusoidal-terrain reference is incompatible. Under the documented numerical-rank criterion, its APV map has rank 65 and nullity 35, leaving 21 fewer resolved zero-APV directions than the flat wave count. Its minimum bottom-constraint residual is approximately $2.46\times10^{-5}$, compared with the $10^{-12}$ gate; both the APV–bottom right compatibility and skew compatibility conditions fail materially. The audit therefore returns no constrained terrain generator. This establishes the incompatible exit without deleting balanced coordinates or relaxing the APV rank. Milestone 5 remains blocked under the pointwise-APV formulation; Milestone 4.6 tests the revised quadratic-invariant constraint.

## Milestone 4.6: Quadratic energy–enstrophy closure audit

- [ ] Complete — blocking scientific gate

### Purpose

Determine whether the current complete mixed state space admits a minimum-change evolution that preserves the finite-terrain energy and total discrete quadratic potential enstrophy while enforcing the resolved strong bottom-displacement evolution. This replaces the pointwise APV constraint that Milestone 4.5 proved incompatible; it does not alter the energy or potential-enstrophy norms.

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
  {\max(1,\lVert K_c\rVert_F)}
  \leq10^{-13},
  ```

  ```math
  \frac{\lVert G_\gamma K_c-K_cG_\gamma\rVert_F}
  {\max(1,\lVert G_\gamma\rVert_F\lVert K_c\rVert_F)}
  \leq10^{-12},
  ```

  ```math
  \frac{\lVert\overline B K_c-\overline R\rVert_F}
  {\max(1,\lVert\overline R\rVert_F)}
  \leq10^{-12}.
  ```

- Random complex and conjugate-symmetric real states must have zero instantaneous finite-terrain energy and quadratic potential-enstrophy tendencies at the stated tolerances.
- The constrained generator must preserve real-field conjugacy, have real frequencies to eigensolver tolerance, and satisfy the resolved bottom relation for random states and eigenvectors.
- Report $Q_\gamma K_c$ and the APV content of every eigenvector without requiring either to vanish. Verify instead that APV redistribution leaves the total quadratic potential enstrophy unchanged.
- The audit has three scientifically valid exits:
  - **Feasible and convergent:** all exact constraints close, the correction tends to zero with terrain amplitude and decreases or stabilizes under resolution refinement, and the corrected dynamics remain close to the raw oracle. Milestone 5 may proceed with $J_\gamma^c$.
  - **Feasible but dynamically poor:** the algebraic constraints close, but the correction remains order one, fails to converge, or produces modal structure inconsistent with the raw oracle. Retain the result as a diagnostic and keep Milestone 5 blocked.
  - **Incompatible:** the rank-revealing solve establishes a nonzero minimum bottom or conjugacy residual within the exact commutant of $G_\gamma$. Keep Milestone 5 blocked and next investigate revised boundary forms or state coordinates.
- Approximate feasibility obtained by merging spectrally distinct eigenvalues, relaxing invariant tolerances, deleting physical coordinates, or changing the energy or potential-enstrophy norm does not pass this gate.

## Milestone 5: Terrain modes and wave–balanced separation

- [ ] Complete — blocking scientific gate

### Purpose

Demonstrate that the finite-terrain weak system produces physically admissible oscillatory and balanced modes before production optimization.

### Dependencies

Milestone 4.6 with its feasible-and-convergent exit.

### Deliverables

- Solve the complete quadratic-invariant-preserving dense generalized Hermitian problem

```math
iJ_\gamma^c\boldsymbol c_n
=
\Omega_nE_\gamma\boldsymbol c_n.
```

- Retain the raw $J_\gamma$ eigensystem and correction norm as diagnostic comparisons; do not describe $J_\gamma^c$ as the unmodified Galerkin operator.
- Classify modes by frequency, quadratic potential-enstrophy content, and projection onto the flat wave and balanced subspaces. Retain the zero-frequency complement, including bottom buoyancy.
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
- Total quadratic potential enstrophy is conserved to $10^{-12}$, and modal APV content and APV redistribution converge under refinement.
- The balanced dimension and bottom-buoyancy freedom agree with the discrete state count.
- Sinusoidal terrain produces the predicted Fourier couplings.
- Bottom and strong-equation residuals decrease under horizontal and vertical refinement.
- Failure of the quadratic potential-enstrophy, energy, bottom, or strong-residual gates blocks Milestone 6. Nonzero modal APV alone is reported but is not a failure when the exact quadratic invariant and other gates pass.

## Milestone 6: Residual-enriched terrain dressing

- [ ] Complete

### Purpose

Construct selected terrain modes efficiently from flat nonhydrostatic seeds without diagonalizing the complete global dense problem.

### Dependencies

Milestone 5.

### Deliverables

- Use the Milestone-3 flat nonhydrostatic Ritz modes as initial vectors.
- Evaluate each seed's residual in the quadratic-invariant-preserving finite-terrain generalized eigenproblem approved by Milestone 4.6.
- Apply block residual corrections preconditioned by the flat signed-frequency operator.
- Remove the complete degenerate or near-resonant flat eigenspace before applying the complementary inverse.
- Include the bottom coefficient in the same correction, energy orthogonalization, and reduced Ritz solve.
- Rediagonalize the exact reduced $(iJ_\gamma^c,E_\gamma)$ pair after every enrichment step.
- Compare selective dressed modes with the constrained dense terrain oracle while retaining raw-operator residuals as diagnostics.

### Automated acceptance

- For weak nonresonant terrain, one correction reduces an $O(h)$ seed residual to $O(h^2)$, with observed terrain-amplitude order at least 1.8.
- Degenerate tests converge only after the complete coupled block is included.
- Residual norms decrease monotonically after accepted enrichment steps.
- Repeated enrichment reproduces targeted dense eigenvalues and energy-normalized eigenspaces within $10^{-9}$.
- Dressed modes reproduce the dense oracle's quadratic potential-enstrophy content and satisfy the bottom relation at the dense-oracle tolerance.

## Milestone 7: Matrix-free operator actions

- [ ] Complete

### Purpose

Replace global dense terrain matrices with adjoint-consistent field reconstruction, terrain multiplication, and projection.

### Dependencies

Milestone 6.

### Deliverables

- Implement matrix-free applications of $E_\gamma$, raw $J_\gamma$, $Q_\gamma$, and the feasible quadratic-invariant-preserving closure selected in Milestone 4.6.
- Unpack interior and bottom coordinates, reconstruct the mapped and physical fields, multiply by terrain and stratification weights on an oversampled grid, and project with the exact adjoints.
- Use a default horizontal oversampling factor of two and permit larger factors for convergence studies.
- Implement flat signed-frequency preconditioning and block iterative eigensolves.
- Exploit preserved meridional wavenumber and Bloch classes when the terrain provides those symmetries.
- Keep diagnostic pressure recovery outside ordinary operator applications.

### Automated acceptance

- Matrix-free raw and constrained actions agree with their dense oracles within $10^{-10}$ at reference resolution.
- Matrix-free energy and exchange actions satisfy their adjoint identities within $10^{-12}$.
- Matrix-free Ritz values, eigenspaces, quadratic potential enstrophy, modal APV diagnostics, and residual histories reproduce the dense results.
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
\left(E_\gamma-\frac{\Delta t}{2}J_\gamma^c\right)\boldsymbol a^{n+1}
=
\left(E_\gamma+\frac{\Delta t}{2}J_\gamma^c\right)\boldsymbol a^n.
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
- Quadratic potential enstrophy remains constant, while resolved APV redistribution agrees with the approved dense generator.
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
- Broadband and variable-stratification calculations retain the structural, energy, quadratic potential-enstrophy, and strong-residual gates; modal APV leakage is reported separately.
- Benchmarks report construction time, peak stored state, operator-application time, iteration counts, and reconstruction time.
- Selective terrain-mode evolution reduces to phase multiplication and reconstruction.
- Broad-state evolution reduces to matrix-free $E_\gamma$ and constrained $J_\gamma^c$ actions plus preconditioned Cayley solves.

## Planning and implementation cadence

- **Batch A — Milestones 1–3.** Plan these milestones together and implement them sequentially. Milestone 3 is a blocking scientific gate: stop if the mixed hydrostatic basis does not accurately and stably recover the flat nonhydrostatic transform.
- **Batch B — Milestones 4–6.** Plan only after the flat gate passes. Establish the complete dense finite-terrain oracle. Preserve Milestone 4.5 as the documented incompatible pointwise-APV audit, then perform the Milestone-4.6 quadratic energy–enstrophy closure audit. Only a feasible-and-convergent Milestone-4.6 exit permits Milestone 5, which remains the second blocking modal gate.
- **Batch C — Milestones 7–9.** Plan only after the dense terrain modes and residual-enrichment tests pass. Implement matrix-free actions before online evolution, resolution rebuilding, or persistence.
- Use one focused commit per completed milestone and retain the acceptance evidence in the automated tests and examples.
- Do not bypass a blocking gate by symmetrizing a scientifically incorrect operator or by introducing empirical correction factors.

## Definition of done

The scientific proof of concept is established when Milestones 1–6 pass: the mixed basis recovers the flat nonhydrostatic modes, the dense finite-terrain forms have the correct energy structure, the quadratic energy–enstrophy closure has a feasible-and-convergent exit, terrain modes conserve both quadratic invariants and satisfy the physical bottom condition, modal APV diagnostics converge under refinement, and residual dressing converges to the constrained dense oracle.

The research implementation is complete when Milestones 7–9 pass: matrix-free actions reproduce the oracle, linear evolution preserves finite-terrain energy and quadratic potential enstrophy, the uniform-depth and scattering examples converge, and restartable broadband calculations require no diagnostic pressure solve during ordinary evolution.

Nonlinear terrain dynamics, an additive terrain `WVForcing`, an MPM release, modifications to WaveVortexModel, independent surface buoyancy, and dynamic barotropic backreaction remain outside this roadmap.
