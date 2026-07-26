# Terrain-energy Galerkin milestones

## Objective

Develop `WVTerrainEnergyGalerkin`, a pressure-free Galerkin system for the linear rotating Boussinesq equations over stationary bottom topography. The formulation is linear in flow amplitude, exact in the resolved terrain, and conserves the discrete finite-terrain energy.

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
- The formulation is not an additive `WVForcing`; ordinary evolution acts through $E_\gamma^{-1}J_\gamma$, terrain-mode phases, or an energy-preserving implicit step.

## Milestone 1: Repository and scientific contract

- [ ] Complete

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

- [ ] Complete

### Purpose

Construct the complete coordinate space required by the finite-terrain weak equations before forming any dynamical matrices.

### Dependencies

Milestone 1.

### Deliverables

- Construct a matched `WVTransformHydrostatic` using the target domain, grid, stratification, Coriolis parameter, density, gravity, and antialias convention.
- Use its ordinary $F$–$G$ modes for the homogeneous interior fields.
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

- [ ] Complete — blocking scientific gate

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

- [ ] Complete

### Purpose

Establish the exact-in-resolved-$h$ finite-terrain weak system as a low-resolution dense scientific oracle.

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
\operatorname{Re}
\left(\boldsymbol a^*J_\gamma\boldsymbol a\right)=0.
```

- Use constant $\gamma$ as an exact mapped-coordinate oracle and compare against an independent flat transform of physical depth $H=\gamma D$.
- Separate quadrature, truncation, and eigensolver errors in the diagnostics.

### Automated acceptance

- Raw finite-terrain structural residuals are below $10^{-12}$.
- $E_\gamma$ is positive definite and remains well conditioned at the documented reference terrains.
- Random-state finite-terrain energy tendencies close to roundoff.
- For $h=0$, all finite-terrain matrices reproduce the Milestone-3 flat matrices within $10^{-12}$.
- Uniform-depth frequencies, mapped eigenfunctions, and energy normalization converge to the independent depth-$H$ solution within $10^{-9}$.
- Repeated dense construction is deterministic to $10^{-13}$.

## Milestone 5: Terrain modes and wave–balanced separation

- [ ] Complete — blocking scientific gate

### Purpose

Demonstrate that the finite-terrain weak system produces physically admissible oscillatory and balanced modes before production optimization.

### Dependencies

Milestone 4.

### Deliverables

- Solve the complete dense generalized Hermitian problem

```math
iJ_\gamma\boldsymbol c_n
=
\Omega_nE_\gamma\boldsymbol c_n.
```

- Classify nonzero-frequency modes with $Q_\gamma$ and retain the zero-frequency balanced complement, including bottom buoyancy.
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
- Nonzero-frequency modes have negligible APV at the spatial-discretization tolerance.
- The balanced dimension and bottom-buoyancy freedom agree with the discrete state count.
- Sinusoidal terrain produces the predicted Fourier couplings.
- Bottom and strong-equation residuals decrease under horizontal and vertical refinement.
- Failure of the APV, energy, bottom, or strong-residual gates blocks Milestone 6.

## Milestone 6: Residual-enriched terrain dressing

- [ ] Complete

### Purpose

Construct selected terrain modes efficiently from flat nonhydrostatic seeds without diagonalizing the complete global dense problem.

### Dependencies

Milestone 5.

### Deliverables

- Use the Milestone-3 flat nonhydrostatic Ritz modes as initial vectors.
- Evaluate each seed's residual in the exact finite-terrain generalized eigenproblem.
- Apply block residual corrections preconditioned by the flat signed-frequency operator.
- Remove the complete degenerate or near-resonant flat eigenspace before applying the complementary inverse.
- Include the bottom coefficient in the same correction, energy orthogonalization, and reduced Ritz solve.
- Rediagonalize the exact reduced $(iJ_\gamma,E_\gamma)$ pair after every enrichment step.
- Compare selective dressed modes with the complete dense terrain oracle.

### Automated acceptance

- For weak nonresonant terrain, one correction reduces an $O(h)$ seed residual to $O(h^2)$, with observed terrain-amplitude order at least 1.8.
- Degenerate tests converge only after the complete coupled block is included.
- Residual norms decrease monotonically after accepted enrichment steps.
- Repeated enrichment reproduces targeted dense eigenvalues and energy-normalized eigenspaces within $10^{-9}$.
- Dressed nonzero-frequency modes retain negligible APV and satisfy the bottom relation at the dense-oracle tolerance.

## Milestone 7: Matrix-free operator actions

- [ ] Complete

### Purpose

Replace global dense terrain matrices with adjoint-consistent field reconstruction, terrain multiplication, and projection.

### Dependencies

Milestone 6.

### Deliverables

- Implement matrix-free applications of $E_\gamma$, $J_\gamma$, and $Q_\gamma$.
- Unpack interior and bottom coordinates, reconstruct the mapped and physical fields, multiply by terrain and stratification weights on an oversampled grid, and project with the exact adjoints.
- Use a default horizontal oversampling factor of two and permit larger factors for convergence studies.
- Implement flat signed-frequency preconditioning and block iterative eigensolves.
- Exploit preserved meridional wavenumber and Bloch classes when the terrain provides those symmetries.
- Keep diagnostic pressure recovery outside ordinary operator applications.

### Automated acceptance

- Matrix-free actions agree with the dense oracle within $10^{-10}$ at reference resolution.
- Matrix-free energy and exchange actions satisfy their adjoint identities within $10^{-12}$.
- Matrix-free Ritz values, eigenspaces, APV, and residual histories reproduce the dense results.
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
- Add uniform-depth, sinusoidal-terrain, and Gaussian-ridge examples showing modal scattering, physical $x$–$z$ fields, bottom displacement, APV, and finite-terrain energy.
- Report operator, eigensolver, and time-integration errors separately.

### Automated acceptance

- Terrain-mode phase evolution agrees with direct matrix exponentiation at low resolution.
- Cayley evolution conserves $\boldsymbol a^*E_\gamma\boldsymbol a/2$ to the linear-solver tolerance.
- The uniform-depth evolution agrees with the independent exact depth-$H$ solution.
- Sinusoidal and Gaussian examples converge with time step, horizontal resolution, vertical modes, and dressing iterations.
- Nonzero-frequency evolution remains in the discrete zero-APV wave space.
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
- Broadband and variable-stratification calculations retain the structural, energy, APV, and strong-residual gates.
- Benchmarks report construction time, peak stored state, operator-application time, iteration counts, and reconstruction time.
- Selective terrain-mode evolution reduces to phase multiplication and reconstruction.
- Broad-state evolution reduces to matrix-free $E_\gamma$ and $J_\gamma$ actions plus preconditioned Cayley solves.

## Planning and implementation cadence

- **Batch A — Milestones 1–3.** Plan these milestones together and implement them sequentially. Milestone 3 is a blocking scientific gate: stop if the mixed hydrostatic basis does not accurately and stably recover the flat nonhydrostatic transform.
- **Batch B — Milestones 4–6.** Plan only after the flat gate passes. Establish the complete dense finite-terrain oracle before implementing selective residual dressing. Milestone 5 is the second blocking gate.
- **Batch C — Milestones 7–9.** Plan only after the dense terrain modes and residual-enrichment tests pass. Implement matrix-free actions before online evolution, resolution rebuilding, or persistence.
- Use one focused commit per completed milestone and retain the acceptance evidence in the automated tests and examples.
- Do not bypass a blocking gate by symmetrizing a scientifically incorrect operator or by introducing empirical correction factors.

## Definition of done

The scientific proof of concept is established when Milestones 1–6 pass: the mixed basis recovers the flat nonhydrostatic modes, the dense finite-terrain forms have the correct energy structure, oscillatory modes have negligible APV and satisfy the physical bottom condition, and residual dressing converges to the dense oracle.

The research implementation is complete when Milestones 7–9 pass: matrix-free actions reproduce the oracle, linear evolution preserves finite-terrain energy, the uniform-depth and scattering examples converge, and restartable broadband calculations require no diagnostic pressure solve during ordinary evolution.

Nonlinear terrain dynamics, an additive terrain `WVForcing`, an MPM release, modifications to WaveVortexModel, independent surface buoyancy, and dynamic barotropic backreaction remain outside this roadmap.
