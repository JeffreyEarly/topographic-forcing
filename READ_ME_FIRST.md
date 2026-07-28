# Read me first: terrain-energy Galerkin status

## Current status

**This research direction is paused at a scientific representation blocker.** The project is not abandoned, but Milestone 7 must not begin from the current finite primitive representation.

The authoritative global numerical checkpoint is commit [`060560a`](https://github.com/JeffreyEarly/topographic-forcing/tree/060560a) on the [`terrain-energy-galerkin`](https://github.com/JeffreyEarly/topographic-forcing/tree/terrain-energy-galerkin) branch. The matching mathematical checkpoint is commit `7527984` in the `ape-apv-bottom-topography` literature repository.

The objective was to construct a flow-linear system, exact in the resolved stationary terrain, that satisfies the continuum requirements

```math
\frac{d\mathcal E_\gamma^{(2)}}{dt}=0,
\qquad
\partial_t q_\gamma=0,
\qquad
\partial_t\hat\eta_b
=
\boldsymbol u_{H,b}\boldsymbol{\cdot}\nabla_H h.
```

For a discrete generator `L_gamma`, energy matrix `E_gamma`, APV map `Q_gamma`, bottom extraction `B`, and bottom-kinematic map `R_h`, the corresponding tests are

```math
L_\gamma^*E_\gamma+E_\gamma L_\gamma=0,
\qquad
Q_\gamma L_\gamma=0,
\qquad
BL_\gamma=R_h.
```

Pointwise APV conservation implies the quadratic potential-enstrophy identity

```math
L_\gamma^*Z_\gamma+Z_\gamma L_\gamma=0,
\qquad
Z_\gamma=Q_\gamma^*Q_\gamma.
```

The continuous equations satisfy these requirements simultaneously. The blocker described here belongs to the current finite representation, not to the physical conservation laws.

## What is established

The Milestone-6.3 global small-terrain oracle differentiates a fully coupled signed-Fourier primitive saddle system about flat topography. It retains pressure until continuity and the strong bottom equation have been imposed, constructs the first terrain coefficient analytically, and verifies it independently with centered differences of the unexpanded mapped equations.

For constant stratification, resolution `[6 6 5]`, terrain `20*cos(2*pi*y/Ly) m`, polynomial degree six, and horizontal oversampling factor two, the verified normalized defects are:

| Test | Defect |
|---|---:|
| Analytic versus centered terrain tangent | `2.97e-10` |
| Weak evolution | `3.94e-15` |
| Physical energy | `2.42e-15` |
| Strong bottom evolution | `2.47e-13` |
| Fourier conjugacy | `5.79e-16` |
| Pointwise APV | `2.34e-1` |
| Quadratic potential enstrophy | `3.28e-4` |

The energy and exchange forms retain their Hermitian structure, the flat nonhydrostatic dispersion converges, and a single terrain Fourier component produces only the expected neighboring horizontal couplings.

The last complete repository verification at checkpoint `060560a` passed 101 tests with zero failures.

Milestone 6.4 subsequently validates the coupled volume–boundary PV frequency problem in a local single-wavenumber QG oracle. The direct $L^2\oplus\mathbb C$ system and an independent discretization of Yassin's endpoint problem agree to `4.82e-11` in the leading resolved frequencies, while physical energy and signed pseudoenstrophy close at roundoff. The $f$-plane boundary mode has zero volume APV and active bottom PV, and its representation in the existing balanced-plus-bottom basis converges with vertical resolution. The complete suite now contains 108 passing tests.

This is positive evidence for the boundary-dynamical formulation, but it does not alter the global Milestone-6.3 blocker: the periodic finite Fourier state remains unclosed under the terrain convolution required by the volume-APV cancellation.

## The blocker

The APV defect separates into two parts:

- For horizontally interior input modes, whose first-order sidebands remain in the state, the normalized defect decreases from `4.44e-3` at polynomial degree two to `5.38e-4` at degree six.
- For Fourier-edge inputs, the defect remains near `4.14e-1` because terrain multiplication creates sidebands outside the retained prognostic state.

Thus vertical refinement improves the interior approximation but cannot make a fixed finite horizontal Fourier state closed under terrain convolution. The cancellation required by

```math
Q_0L_1+Q_1L_0=0
```

cannot be represented at the spectral edge by the current state and APV range. The oracle classifies this result as `representation-blocker`.

This conclusion is finite-dimensional. It is not evidence that stationary APV is incompatible with topography in the continuous equations.

## What the investigation ruled out

The sequence of independent audits removed the following explanations:

1. A displacement-only bottom function was incomplete. It was replaced by a complete balanced bottom inversion carrying its linked pressure, velocity, displacement, and APV structure.
2. The complete bottom inversion repaired the flat state space but did not make the finite-terrain volume projection compatible with APV and the strong bottom equation.
3. Retaining pressure and the independent bottom row in a primitive descriptor showed that premature pressure elimination was not the only problem.
4. An independent polynomial primitive oracle showed that hydrostatic `F`–`G` coordinates and the public modal reconstruction were not the sole causes.
5. The frozen constant-slope approximation introduced an analytic cross-slope APV source. Restoring the globally varying terrain coefficients removed that approximation and recovered the correct Fourier selection.
6. Increasing vertical polynomial degree reduced the interior error but left the Fourier-edge failure unchanged.

Changing the bottom basis alone, changing a quadratic norm without a derived invariant, or tuning the existing truncation therefore does not address the remaining issue.

## Do not resume by

Do not proceed by:

- starting Milestone 7 with the current representation;
- projecting the generator into an APV nullspace;
- empirically symmetrizing the generator;
- fitting a minimum-change closure;
- adding an endpoint energy or enstrophy term without an exact closure identity;
- deleting failing balanced or spectral-edge states; or
- interpreting energy conservation alone as validation of the state tendency.

Those operations either answer a different dynamical question or hide the diagnosed incompatibility.

## Conditions for resuming

Resume this direction only after deriving one of the following:

- a finite horizontal state closed under the terrain products required by both the tendency and APV maps;
- a rigorously projected discrete APV law whose state, range, and quadratic invariant use the same truncation; or
- a different compatible discretization, such as an enriched mixed or exact-sequence formulation, that preserves the continuum identities without empirical correction.

Before finite-amplitude terrain or modal construction resumes, the replacement must:

1. reproduce the analytic first terrain coefficient independently;
2. satisfy weak evolution, physical energy, pointwise APV, strong bottom evolution, and conjugacy on the complete retained state;
3. reduce both interior and Fourier-edge APV defects below `1e-10`;
4. converge under independent horizontal, vertical, quadrature, and oversampling refinement; and
5. pass without post hoc projection or symmetrization.

## Reproducing the numerical result

The public entry point is [`auditGlobalSmallTerrainPrimitive`](<@WVTerrainEnergyGalerkin/auditGlobalSmallTerrainPrimitive.m>):

```matlab
audit = problem.auditGlobalSmallTerrainPrimitive( ...
    polynomialDegree=6, ...
    tangentStep=1e-3);
```

The acceptance tests are in [`TestWVTerrainEnergyGlobalSmallTerrainPrimitive`](UnitTests/TestWVTerrainEnergyGlobalSmallTerrainPrimitive.m), and the full numerical record is in [Milestone 6.3](milestones.md#milestone-63-global-small-terrain-primitive-compatibility-oracle).

The local coupled-PV oracle is exposed by [`auditCoupledPVFrequencyOracle`](<@WVTerrainEnergyGalerkin/auditCoupledPVFrequencyOracle.m>), tested by [`TestWVTerrainEnergyCoupledPVFrequencyOracle`](UnitTests/TestWVTerrainEnergyCoupledPVFrequencyOracle.m), and recorded in [Milestone 6.4](milestones.md#milestone-64-coupled-volumeboundary-pv-frequency-oracle).

## Mathematical sources

The mathematical repository is `ape-apv-bottom-topography` at checkpoint `7527984`. Its relevant sources and equation labels are:

- `main.tex`: `eq:finite-terrain-linear-displacement-boundaries`, `eq:finite-terrain-linear-total-energy`, `eq:finite-terrain-linear-apv`, `eq:finite-terrain-linear-potential-enstrophy`, and `eq:finite-terrain-pressure-free-weak-equation`;
- `terrain-energy-galerkin.tex`: `eq:galerkin-linear-apv`, `eq:galerkin-discrete-compatibility-identities`, `eq:galerkin-boundary-descriptor-pencil`, and `eq:galerkin-first-order-conservation-identities`;
- `boundary-energy-enstrophy.tex`: the exact energy, potential-enstrophy, bottom-Casimir, and complete zero-APV state-space derivations.

The literature repository is hosted privately, so this public repository identifies it by repository name, commit, source file, and equation label rather than publishing its Overleaf project URL.

## Practical alternative

The [`mean-depth-wave-generator`](https://github.com/JeffreyEarly/topographic-forcing/tree/mean-depth-wave-generator) branch, currently at checkpoint [`753432d`](https://github.com/JeffreyEarly/topographic-forcing/tree/753432d), remains the useful near-term implementation. It provides validated first-order bottom wave generation and scattering without direct linear PV generation.

That implementation is a controlled first-order approximation. It does not resolve the exact terrain-energy Galerkin representation problem documented here.
