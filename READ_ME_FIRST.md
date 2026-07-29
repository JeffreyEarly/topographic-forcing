# Read me first: terrain-energy Galerkin status

## Current status

**This research direction is paused before Milestone 6.7, a defined dealiased projected-primitive tangent experiment.** Milestone 6.5 shows that a rigorously projected QG volume–boundary PV state handles discarded Fourier sidebands correctly. Milestone 6.6 shows that replacing primitive stationary rows by those projected conservation rows is not equivalent to the unmodified primitive weak equations. The next experiment will therefore keep the primitive energy and exchange rows and apply one common padded Galerkin projection to the primitive forms, APV diagnostic, and bottom equation.

The authoritative completed numerical checkpoint is commit [`d5ab18f`](https://github.com/JeffreyEarly/topographic-forcing/tree/d5ab18f) on the [`terrain-energy-galerkin`](https://github.com/JeffreyEarly/topographic-forcing/tree/terrain-energy-galerkin) branch. The matching mathematics is maintained in `finite-terrain-projection-problem.tex` and `terrain-energy-galerkin.tex` in the `ape-apv-bottom-topography` literature repository.

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

This establishes that the Fourier edge is not an obstruction to a rigorously projected coupled-PV law. It also shows why the next experiment must retain the unmodified primitive weak evolution instead of substituting projected PV rows.

## What the spectral-edge audit means

The APV defect separates into two parts:

- For horizontally interior input modes, whose first-order sidebands remain in the state, the normalized defect decreases from `4.44e-3` at polynomial degree two to `5.38e-4` at degree six.
- For Fourier-edge inputs, the full-grid defect remains near `4.14e-1` because terrain multiplication creates sidebands outside the retained prognostic state.

The edge value is an expected truncation residual when an unprojected oversampled product is compared with a finite state. A finite Fourier space is not required to be closed under multiplication. Zero-padding computes the external sideband without aliasing, and adjoint restriction removes it from the retained equations.

The still-unresolved question is whether the projected coefficient identity

```math
Q_{0,N}L_1+Q_{1,N}L_0=0
```

closes for the unmodified primitive Galerkin operator or converges on a trusted band. The existing `representation-blocker` status records the historical unprojected audit; it is not the classification of the new projected experiment.

## What the investigation ruled out

The sequence of independent audits removed the following explanations:

1. A displacement-only bottom function was incomplete. It was replaced by a complete balanced bottom inversion carrying its linked pressure, velocity, displacement, and APV structure.
2. The complete bottom inversion repaired the flat state space but did not make the finite-terrain volume projection compatible with APV and the strong bottom equation.
3. Retaining pressure and the independent bottom row in a primitive descriptor showed that premature pressure elimination was not the only problem.
4. An independent polynomial primitive oracle showed that hydrostatic `F`–`G` coordinates and the public modal reconstruction were not the sole causes.
5. The frozen constant-slope approximation introduced an analytic cross-slope APV source. Restoring the globally varying terrain coefficients removed that approximation and recovered the correct Fourier selection.
6. Increasing vertical polynomial degree reduced the interior error; the unchanged Fourier-edge value was subsequently recognized as discarded external support.

Changing the bottom basis alone, changing a quadratic norm without a derived invariant, or tuning the existing truncation therefore does not address the remaining issue.

## Do not resume by

Do not proceed by:

- starting Milestone 7 before the Milestone-6.7 projected tangent gate;
- projecting the generator into an APV nullspace;
- empirically symmetrizing the generator;
- fitting a minimum-change closure;
- adding an endpoint energy or enstrophy term without an exact closure identity;
- deleting failing balanced states or hiding external sidebands; or
- interpreting energy conservation alone as validation of the state tendency.

Those operations either answer a different dynamical question or hide the diagnosed incompatibility.

## Conditions for resuming

Resume with Milestone 6.7, which must:

1. retain the unmodified primitive `E_0,E_1,J_0,J_1` construction;
2. use one zero-pad, multiply, and adjoint-restrict operation for every terrain product;
3. separate a trusted physical band from an evolved outer support band;
4. reproduce exact Fourier convolution with the oversampled pseudospectral calculation below `1e-12`;
5. retain weak evolution, physical energy, bottom evolution, and conjugacy below `1e-11`;
6. make padding factors two and three agree on the trusted band within `1e-10`; and
7. either close projected APV below `1e-10` or reduce its trusted-band residual by at least a factor of four per refinement to below `1e-8`.

External sidebands must be reported separately. No corrected generator, fitted closure, empirical symmetrization, or APV-nullspace projection is permitted.

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

## Mathematical sources

The mathematical repository is `ape-apv-bottom-topography`. Its relevant sources and equation labels are:

- `main.tex`: `eq:finite-terrain-linear-displacement-boundaries`, `eq:finite-terrain-linear-total-energy`, `eq:finite-terrain-linear-apv`, `eq:finite-terrain-linear-potential-enstrophy`, and `eq:finite-terrain-pressure-free-weak-equation`;
- `terrain-energy-galerkin.tex`: `eq:galerkin-linear-apv`, `eq:galerkin-discrete-compatibility-identities`, `eq:galerkin-boundary-descriptor-pencil`, and `eq:galerkin-first-order-conservation-identities`;
- `boundary-energy-enstrophy.tex`: the exact energy, potential-enstrophy, bottom-Casimir, and complete zero-APV state-space derivations.

The literature repository is hosted privately, so this public repository identifies it by repository name, commit, source file, and equation label rather than publishing its Overleaf project URL.

## Practical alternative

The [`mean-depth-wave-generator`](https://github.com/JeffreyEarly/topographic-forcing/tree/mean-depth-wave-generator) branch, currently at checkpoint [`753432d`](https://github.com/JeffreyEarly/topographic-forcing/tree/753432d), remains the useful near-term implementation. It provides validated first-order bottom wave generation and scattering without direct linear PV generation.

That implementation is a controlled first-order approximation. It does not resolve the exact terrain-energy Galerkin representation problem documented here.
