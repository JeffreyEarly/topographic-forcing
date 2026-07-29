# Read me first: terrain-energy Galerkin status

## Current status

**This branch is at a passing Milestone-6.8 tangent checkpoint.** Milestone 6.7 showed that a fixed flat test space does not recover the trusted-band APV cancellation even under a common padded projection. Milestone 6.8 returns to the continuous Green identity and includes the terrain derivative of the stationary geostrophic test states. The resulting primitive weak APV moments, physical energy, projected bottom evolution, and tangent eigenproblem pass without replacement APV rows or corrections. Milestone 7 and finite-amplitude terrain have not begun.

The pre-Milestone-6.7 checkpoint is commit [`d5ab18f`](https://github.com/JeffreyEarly/topographic-forcing/tree/d5ab18f) on the [`terrain-energy-galerkin`](https://github.com/JeffreyEarly/topographic-forcing/tree/terrain-energy-galerkin) branch. The current branch contains the completed projected and boundary-complete weak oracles. The matching mathematics is maintained in `finite-terrain-projection-problem.tex` and `terrain-energy-galerkin.tex` in the `ape-apv-bottom-topography` literature repository.

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

For a retained state, the standard terrain multiplication is

```math
\mathcal M_{h,N}
=
P_{M\to N}\mathcal M_hI_{N\to M},
```

where `I_N->M` zero-pads into an oversampled evaluation space and `P_M->N` is its quadrature adjoint. The mandatory finite identities are

```math
L_\gamma^*E_\gamma+E_\gamma L_\gamma=0,
\qquad
BL_\gamma=R_{h,N}.
```

Projected APV is tested with `Q_gamma,N = P_q,N Q_gamma I_N->M`. Exact closure `Q_gamma,N L_gamma = 0` is the strongest result. Otherwise its residual and the corresponding potential-enstrophy residual must converge on a fixed trusted physical band under independent support, padding, and vertical refinement.

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

Milestone 6.5 then constructs the periodic QG state directly from volume APV and bottom PV. Its exact projected terrain convolution and an independent oversampled pseudospectral Jacobian agree below `8e-16`. Physical energy, stationary volume APV, physical potential enstrophy, projected bottom evolution, and Fourier conjugacy close at roundoff, including for inputs whose unprojected terrain sidebands leave the retained Fourier state. The complete suite now contains 115 passing tests.

Milestone 6.6 combines flat primitive wave rows with independent projected-volume-APV and bottom rows. The resulting descriptor is square, full rank, and well enough conditioned; it recovers the flat primitive generator and satisfies all three selected row families at roundoff. It nevertheless fails the complete primitive weak equation by approximately `1.26e-3` and physical energy by `1.49e-3`. Those are genuine row-equivalence failures. Its `4.14e-1` full-grid APV value at the meridional edges is instead the norm of an external unprojected sideband.

The complete repository suite now contains 120 passing tests, including the Milestone-6.6 blocker audit.

Milestone 6.7 performs that next experiment without replacing any primitive row. For trusted bounds `[1 0]`, supports `[1 1;1 2;1 3]`, vertical degrees `[4;8;12]`, and padding factors two and three, its structural defects are:

| Test | Maximum defect |
|---|---:|
| Prolongation/restriction and exact convolution | `1.19e-15` |
| Analytic versus centered terrain tangent | `6.10e-10` |
| Primitive weak evolution | `9.98e-15` |
| Physical energy | `4.57e-15` |
| Projected bottom evolution | `3.15e-12` |
| Fourier conjugacy | `2.37e-15` |
| Padding factors two versus three | `4.39e-12` |

The trusted-band APV defects are `5.45e-4`, `1.85e-4`, and `1.01e-4`; the potential-enstrophy defects are `3.34e-4`, `1.12e-4`, and `6.08e-5`. Their successive reduction factors are approximately `(2.95,1.84)` and `(2.97,1.85)`, rather than the required factor four, and their finest values remain far above `1e-8`. Increasing horizontal support alone does not change the finest-degree trusted result.

This establishes that the earlier `4.14e-1` Fourier-edge value was not the decisive obstruction: a common dealiased projection handles and separately reports those external sidebands. Milestone 6.8 then establishes that the remaining trusted-band defect was caused by omitting the terrain derivative of the geostrophic test inclusion.

For every scalar test with \(\phi_b=0\) and \(\partial_\xi\phi(0)=0\), the exact terrain-dependent geostrophic state belongs to the primitive test space and satisfies

```math
2\rho_0
\langle\Psi_{\mathrm g}[\phi],\psi\rangle_{H_\gamma}
=
-\frac{\rho_0}{A}
\int\gamma\phi^*q_\gamma\,dV.
```

Its first terrain coefficient supplies \(G_1\), which was absent from the fixed-basis Milestone-6.7 audit. The verified tangent identities are

```math
G_1^*H_0+G_0^*H_1
=
-M_1Q_0-M_0Q_1,
```

```math
G_1^*J_0+G_0^*J_1=0.
```

At vertical degree twelve, the APV Green-identity defect is `8.01e-15`, the stationary-row defect is `6.63e-14`, independently evaluated APV evolution is `1.11e-13`, physical energy is `4.28e-15`, and projected bottom evolution is `2.47e-12`. The tangent eigenproblem has residual `9.49e-17`, and its nonzero-frequency modes have weak volume APV `2.58e-11`.

This result changes the interpretation of Milestone 6.7: ordinary padding was necessary but not sufficient because the correct weak scalar-to-primitive inclusion is itself terrain dependent. Once that inclusion and its vertical support are represented, APV follows from the original primitive weak equations rather than being imposed as a replacement row.

The complete repository suite now contains 130 passing tests, including the Milestone-6.8 constant- and variable-stratification weak-sequence and eigenproblem audits.

## What the spectral-edge audit means

The APV defect separates into two parts:

- For horizontally interior input modes, whose first-order sidebands remain in the state, the normalized defect decreases from `4.44e-3` at polynomial degree two to `5.38e-4` at degree six.
- For Fourier-edge inputs, the full-grid defect remains near `4.14e-1` because terrain multiplication creates sidebands outside the retained prognostic state.

The edge value is an expected truncation residual when an unprojected oversampled product is compared with a finite state. A finite Fourier space is not required to be closed under multiplication. Zero-padding computes the external sideband without aliasing, and adjoint restriction removes it from the retained equations.

The fixed-test projected coefficient identity

```math
Q_{0,N}L_1+Q_{1,N}L_0=0
```

does not close exactly for the unmodified primitive Galerkin operator and does not satisfy the approved trusted-band convergence gate. The historical `representation-blocker` status belongs to the unprojected audit; the common projected experiment is classified `nonconvergent-projected`.

## What the investigation ruled out

The sequence of independent audits removed the following explanations:

1. A displacement-only bottom function was incomplete. It was replaced by a complete balanced bottom inversion carrying its linked pressure, velocity, displacement, and APV structure.
2. The complete bottom inversion repaired the flat state space but did not make the finite-terrain volume projection compatible with APV and the strong bottom equation.
3. Retaining pressure and the independent bottom row in a primitive descriptor showed that premature pressure elimination was not the only problem.
4. An independent polynomial primitive oracle showed that hydrostatic `F`–`G` coordinates and the public modal reconstruction were not the sole causes.
5. The frozen constant-slope approximation introduced an analytic cross-slope APV source. Restoring the globally varying terrain coefficients removed that approximation and recovered the correct Fourier selection.
6. Increasing vertical polynomial degree reduced the interior error; the unchanged Fourier-edge value was subsequently recognized as discarded external support.
7. Retaining the terrain derivative of the geostrophic test inclusion restores the APV Green identity and stationary primitive rows at roundoff.

Changing the bottom basis alone, changing a quadratic norm without a derived invariant, or merely enlarging a fixed flat test space does not reproduce the derived weak sequence.

## Do not resume by

Do not proceed by:

- starting Milestone 7 without carrying the Milestone-6.8 terrain-dependent geostrophic test sequence into the finite-amplitude construction;
- projecting the generator into an APV nullspace;
- empirically symmetrizing the generator;
- fitting a minimum-change closure;
- adding an endpoint energy or enstrophy term without an exact closure identity;
- deleting failing balanced states or hiding external sidebands; or
- interpreting energy conservation alone as validation of the state tendency.

Those operations either answer a different dynamical question or hide the diagnosed incompatibility.

## Conditions for continuing

Milestone 6.8 supplies the derived primitive representation that Milestone 6.7 lacked. The next increment may test finite-amplitude terrain only by constructing the terrain-dependent geostrophic inclusion and APV Green map at the same amplitude as the primitive \(H_\gamma,J_\gamma\) forms. It must retain exact physical energy and projected bottom evolution, compare the energy-derived APV moments with an independent strong APV diagnostic, and refine the trusted scalar and primitive support spaces independently.

External sidebands must remain separately reported. No corrected generator, fitted closure, empirical symmetrization, replacement APV rows, APV-nullspace projection, or mode deletion is permitted.

## Reproducing the numerical result

The public entry point is [`auditGlobalSmallTerrainPrimitive`](<@WVTerrainEnergyGalerkin/auditGlobalSmallTerrainPrimitive.m>):

```matlab
audit = problem.auditGlobalSmallTerrainPrimitive( ...
    polynomialDegree=6, ...
    tangentStep=1e-3);
```

The acceptance tests are in [`TestWVTerrainEnergyGlobalSmallTerrainPrimitive`](UnitTests/TestWVTerrainEnergyGlobalSmallTerrainPrimitive.m), and the full numerical record is in [Milestone 6.3](milestones.md#milestone-63-global-small-terrain-primitive-compatibility-oracle).

The local coupled-PV oracle is exposed by [`auditCoupledPVFrequencyOracle`](<@WVTerrainEnergyGalerkin/auditCoupledPVFrequencyOracle.m>), tested by [`TestWVTerrainEnergyCoupledPVFrequencyOracle`](UnitTests/TestWVTerrainEnergyCoupledPVFrequencyOracle.m), and recorded in [Milestone 6.4](milestones.md#milestone-64-coupled-volumeboundary-pv-frequency-oracle).

The periodic projected oracle is exposed by [`auditPeriodicCoupledPVOracle`](<@WVTerrainEnergyGalerkin/auditPeriodicCoupledPVOracle.m>), tested by [`TestWVTerrainEnergyPeriodicCoupledPVOracle`](UnitTests/TestWVTerrainEnergyPeriodicCoupledPVOracle.m), and recorded in [Milestone 6.5](milestones.md#milestone-65-periodic-coupled-volumeboundary-pv-oracle).

The hybrid primitive–PV oracle is exposed by [`auditHybridPrimitivePVOracle`](<@WVTerrainEnergyGalerkin/auditHybridPrimitivePVOracle.m>), tested by [`TestWVTerrainEnergyHybridPrimitivePVOracle`](UnitTests/TestWVTerrainEnergyHybridPrimitivePVOracle.m), and recorded in [Milestone 6.6](milestones.md#milestone-66-hybrid-primitivepv-commuting-oracle).

The common projected primitive oracle is exposed by [`auditDealiasedProjectedPrimitiveTangent`](<@WVTerrainEnergyGalerkin/auditDealiasedProjectedPrimitiveTangent.m>), tested by [`TestWVTerrainEnergyDealiasedProjectedPrimitiveTangent`](UnitTests/TestWVTerrainEnergyDealiasedProjectedPrimitiveTangent.m), and recorded in [Milestone 6.7](milestones.md#milestone-67-dealiased-projected-primitive-tangent-oracle).

The boundary-complete weak oracle is exposed by [`auditBoundaryCompleteWeakEigenproblem`](<@WVTerrainEnergyGalerkin/auditBoundaryCompleteWeakEigenproblem.m>), tested by [`TestWVTerrainEnergyBoundaryCompleteWeakEigenproblem`](UnitTests/TestWVTerrainEnergyBoundaryCompleteWeakEigenproblem.m), and recorded in [Milestone 6.8](milestones.md#milestone-68-boundary-complete-weak-terrain-eigenproblem).

## Mathematical sources

The mathematical repository is `ape-apv-bottom-topography`. Its relevant sources and equation labels are:

- `main.tex`: `eq:finite-terrain-linear-displacement-boundaries`, `eq:finite-terrain-linear-total-energy`, `eq:finite-terrain-linear-apv`, `eq:finite-terrain-linear-potential-enstrophy`, and `eq:finite-terrain-pressure-free-weak-equation`;
- `finite-terrain-projection-problem.tex`: `eq:terrain-projection-apv-green-identity`, `eq:terrain-projection-volume-apv-weak-conservation`, and `eq:terrain-projection-discrete-green-commutation`;
- `terrain-energy-galerkin.tex`: `eq:galerkin-linear-apv`, `eq:galerkin-discrete-compatibility-identities`, `eq:galerkin-boundary-descriptor-pencil`, and `eq:galerkin-first-order-conservation-identities`;
- `boundary-energy-enstrophy.tex`: the exact energy, potential-enstrophy, bottom-Casimir, and complete zero-APV state-space derivations.

The literature repository is hosted privately, so this public repository identifies it by repository name, commit, source file, and equation label rather than publishing its Overleaf project URL.

## Practical alternative

The [`mean-depth-wave-generator`](https://github.com/JeffreyEarly/topographic-forcing/tree/mean-depth-wave-generator) branch, currently at checkpoint [`753432d`](https://github.com/JeffreyEarly/topographic-forcing/tree/753432d), remains the useful near-term implementation. It provides validated first-order bottom wave generation and scattering without direct linear PV generation.

That implementation is a controlled first-order approximation. It does not resolve the exact terrain-energy Galerkin representation problem documented here.
