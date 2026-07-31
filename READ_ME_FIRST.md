# Read me first: terrain-energy Galerkin status

## Current status

**This branch has completed Milestone 10.2 with a numerical blocker and is paused before Milestone 10.2.1.** Milestones 7 and 8 provide a passing finite-amplitude primitive weak oracle and a tangency-defined stationary balanced sector. Milestones 9 and 9.1 preserve that sector, construct the complete physical-energy eigensystem, and identify a converged stationary projector and one conjugate-closed internal-wave projector while retaining every unvalidated direction in an explicit unresolved projector. Milestones 9.2 and 9.3 show why modifying separate vertical problems or local bottom endpoints is insufficient: the leading terrain correction is a global horizontal convolution. Milestone 9.4 constructs that correction from the verified \(H_1,J_1,R_1,G_1\) derivatives. Milestone 10 applies repeated exact finite-amplitude residual corrections and reproduces one independent dense internal-wave block with compression factor `6.22`. Milestone 10.1 declares the smaller state intended for production evolution. Milestone 10.2 then tests complete wave coverage and finds that the current multiblock correction is not robust or economical enough for that state. Its reduced Ritz solve is joint, but its correction columns are still generated through separate flat-frequency shifted inverses.

**The production contract contains only physically declared coordinates.** It uses native WaveVortexModel fixed-\(\kappa\) wave and inertial modes, APV-bearing balanced modes, one existing complete zero-APV bottom inversion per retained nonzero horizontal wavenumber, and compatible MDA and mean-bottom coordinates at \(\kappa=0\). Every coordinate has a physical family, units, native index where applicable, Fourier partner, and converged primitive-oracle image. Guard modes remain evaluation infrastructure, omitted primitive directions are truncation diagnostics, and no unresolved coordinate enters the production layout. Milestones 10.2.1–10.4 must still assign every production coordinate to the finite-terrain stationary, internal-wave, or topographic-wave projector before matrix-free construction or time integration begins.

The pre-Milestone-6.7 checkpoint is commit [`d5ab18f`](https://github.com/JeffreyEarly/topographic-forcing/tree/d5ab18f) on the [`terrain-energy-galerkin`](https://github.com/JeffreyEarly/topographic-forcing/tree/terrain-energy-galerkin) branch. The current branch contains the completed projected, tangent boundary-complete, and finite-amplitude weak oracles. The matching mathematics is maintained in `finite-terrain-projection-problem.tex` and `terrain-energy-galerkin.tex` in the `ape-apv-bottom-topography` literature repository.

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

Milestone 7 removes the small-terrain expansion. It evaluates \(H_\gamma\), \(J_\gamma\), the mapped APV operator, and the terrain-dependent geostrophic inclusion at the same finite amplitude, then advances the retained coefficients with

```math
L_\gamma
=
H_\gamma^{-1}J_\gamma.
```

For the constant-stratification reference calculation with terrain \(20\cos(2\pi y/L_y)\ {\rm m}\), trusted bounds `[1 0]`, supports `[1 1;1 2;1 4]`, primitive degrees `[2;3;5]`, and padding factors two and three, the finest diagnostics are:

| Test | Defect |
|---|---:|
| Primitive weak evolution | `2.62e-16` |
| Physical energy | `1.59e-15` |
| APV Green identity | `4.33e-15` |
| Stationary geostrophic row | `9.86e-17` |
| Weak APV evolution | `1.32e-14` |
| Independent strong APV evolution | `4.47e-15` |
| Weak/strong APV agreement | `4.43e-15` |
| Projected bottom evolution | `2.01e-11` |
| Fourier conjugacy | `7.92e-15` |
| Quadratic potential enstrophy | `3.98e-14` |
| Padding factors two versus three | `2.18e-12` |

The projected bottom defect decreases from `3.47e-5` to `2.89e-7` to `2.01e-11`. The finite forms' remainders relative to their flat-plus-tangent approximations decrease by a factor of approximately four whenever terrain amplitude is halved, confirming the expected \(O(h^2)\) term rather than a first-order inconsistency. An exponential-stratification calculation passes after joint horizontal and vertical enrichment, with geostrophic-state representation, strong APV, and bottom defects of `7.18e-11`, `3.21e-15`, and `1.45e-11`.

Milestone 7 is therefore classified `compatible-finite-amplitude-weak-oracle`. This validates the dense finite-amplitude weak structure but does not yet establish a terrain-mode construction or online evolution method.

Milestone 8 constructs the trusted stationary scalar family with nonzero bottom values and the projected tangency condition

```math
\partial_x\phi_b\,\partial_yh
-
\partial_y\phi_b\,\partial_xh
=0.
```

The resulting primitive states satisfy the complete APV Green identity with

```math
c_b[\psi]
=
f\hat\eta_b+v_bh_x-u_bh_y.
```

For the constant-stratification sinusoidal reference calculation, the geostrophic representation, complete Green identity, stationary row, projected bottom tangency, strong trusted bottom tangency, conjugacy, and padding defects are respectively `1.29e-17`, `1.15e-14`, `1.39e-13`, `1.44e-17`, `2.77e-17`, `6.57e-14`, and `9.37e-13`. Six non-tangent bottom-streamfunction directions remain in the full primitive state.

The raw energy-scaled exchange form contains active topographic boundary directions with singular values near `1e-14 s^-1`. A hard SVD nullspace is therefore ill-conditioned even when the physically derived stationary rows close. The Milestone-8 tangency construction, not a numerical frequency threshold, therefore defines the stationary sector in the Milestone-9 audit.

Milestone 8 is classified `complete-stationary-space-oracle`. Flat and uniform-depth controls pass, variable stratification converges, and ordinary and antialiased layouts preserve the common projection and Fourier conjugacy.

Milestone 9 solves the Cholesky-scaled physical-energy pencil without modifying \(H_\gamma\) or \(J_\gamma\). For constant stratification, the eigen-residual, energy orthogonality, Fourier conjugacy, and retention of all non-tangent bottom directions close below `2e-13`, `5e-15`, `7e-15`, and `3e-15`. Flat and uniform-depth wave frequencies converge below `1e-10` relative error.

That numerical eigensolve is not sufficient by itself to establish physical terrain modes. The backward-error zero-frequency projector contains directions outside the Milestone-8 tangency sector. Four bottom-dominated candidates are spectrally resolved—the smallest ratio of frequency to backward uncertainty exceeds $10^{12}$—but their reference trusted APV, bottom, and strong-equation defects are approximately `1.60e-1`, `4.19e-6`, and `3.67e-1`. The public audit therefore retains the historical status `dense-modal-classification-blocker`; a variable-stratification control also fails the full-space eigen-residual gate.

A fixed-horizontal-support study at vertical degrees four, six, and eight changes the interpretation of that all-mode result. A tracked internal branch near \(4.43f\) converges in frequency, APV, bottom evolution, and the strong primitive residual:

| Degree | $\omega/f$ | APV defect | Bottom defect | Strong residual |
|---:|---:|---:|---:|---:|
| 4 | `4.428504` | `7.07e-6` | `3.49e-10` | `1.78e-3` |
| 6 | `4.428536` | `3.14e-7` | `4.88e-10` | `4.37e-5` |
| 8 | `4.428536` | `7.56e-9` | `5.22e-10` | `6.34e-6` |

The tracked subinertial bottom-dominated branch near \(0.169f\) improves in APV and strong residual, but its bottom defect remains near `5e-6`. It is therefore unresolved, not demonstrated numerical. The constructed geostrophic space continues to satisfy representation, Green-identity, stationary-row, and bottom-tangency tests near roundoff. Additional zero directions outside that construction grow from 88 to 120 to 152 as the vertical space is enriched; they remain unresolved algebraic completion directions rather than validated geostrophic states.

The finite classification target is

```math
\mathcal V_N
=
\mathcal G_N
\mathbin{\oplus_{H_N}}
\mathcal W_N
\mathbin{\oplus_{H_N}}
\mathcal R_N,
```

where `G_N` is the validated stationary balanced space, `W_N` contains converged internal and topographic boundary spectral subspaces, and `R_N` is the retained unresolved completion. All candidates, unresolved zero directions, and non-tangent bottom coordinates remain in the returned basis. They are not reassigned by a frequency cutoff, projected into an APV nullspace, or deleted.

The complete repository suite contains 154 passing tests, including seven Milestone-9 structure, retention, backward-error, reference-limit, and blocker tests. Static analysis reports no issues in all 76 MATLAB files.

Milestone 9.1 implements the required projector-level classification. For weak sinusoidal terrain \(h=2.5\cos(2\pi y/L_y)\ {\rm m}\), nested supports `[1 2;1 3;1 4]`, vertical degrees `[4;6;8]`, and padding factors two and three, the nine-dimensional stationary projector has final nested, fixed-degree guard, and padding defects `1.78e-12`, `1.67e-12`, and `7.82e-13`.

One four-dimensional internal-wave subspace near \(4.4283f\) passes. Its final principal sine, relative frequency change, APV defect, bottom defect, and strong residual are `5.91e-5`, `1.62e-9`, `1.07e-10`, `1.46e-13`, and `2.42e-7`. The principal sine and strong residual decrease by factors greater than four at both refinement steps. Independent guard and padding defects are `5.21e-9` and `1.58e-13`.

The leading topographic-boundary candidate is retained but unresolved because its fixed-degree guard defect is `2.10e-2`. The finest physical-energy dimensions are

```math
\dim\mathcal V_N=703,
\qquad
\dim\mathcal G_N=9,
\qquad
\dim\mathcal W_N=4,
\qquad
\dim\mathcal R_N=690.
```

The physical and remainder projectors are mutually energy orthogonal, Fourier-conjugate closed, and complete within `3.1e-13`. The remainder contains 168 unresolved-zero and 522 unresolved-nonzero directions. No direction is labelled demonstrated numerical.

The complete repository suite now contains 162 passing tests. Static analysis reports no issues in all 79 MATLAB files.

Milestone 9.2 embeds an independently generated boundary-complete vertical family into a degree-12 primitive reference space. The family contains fixed-\(\kappa\) nonhydrostatic waves, signed-Robin APV-bearing geostrophic modes, one explicit zero-APV bottom inversion per retained nonzero horizontal coefficient, and the complete mean sector. The Robin and wave eigenproblems are evaluated with the isolated `internal-modes-evp` checkout at commit `df86687e91faa31bf65941299062d125a96904b1`; the WaveVortexModel dependency remains the main InternalModes checkout.

At the primary Robin length \(\ell_b=-D/4\), increasing the retained wave and geostrophic counts from one to eight improves the internal-wave projector defect from `2.50e-3` to `7.86e-7`, the frequency defect from `1.22e-6` to `3.84e-13`, and the strong primitive residual from `5.80e-2` to `4.36e-6`. The bottom-evolution defect reaches `3.10e-3`, however, and the completed modal dimension grows from 450 to 941 compared with the 1026-dimensional primitive reference. The finest compression factor is therefore only `1.09`.

Padding factors two and three agree to `1.58e-9`. Changing \(\ell_b/D\) over `[-1/8,-1/4,-1/2,Inf]` changes the accepted projector by as much as `2.87e-5`; this is expected for distinct low-order trial spaces and is retained as a convergence-tuning diagnostic, not a physical gate. The isolated provider reproduces its physical Robin spectrum at qualified orders 128 and 256; orders 96 and 192 are rejected because that provider revision returns spurious very-large negative eigenvalues at those orders. The provider and endpoint defects near `1e-10` are far too small to explain the failed bottom and projector gates.

The result is therefore `modal-incompatible` for the tested flat-Dirichlet compression strategy, not a failure of the underlying primitive weak oracle or of the individual reference-mode ideas. The classification follows from slow bottom-evolution and physical-projector convergence together with ineffective compression, not from cross-\(\ell_b\) variation. No eigenvector was removed, merged, or classified by a frequency cutoff, and no evolution row or generator was altered.

The complete repository suite now contains 168 passing tests. Static analysis reports no issues in all 83 MATLAB files.

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

## Do not continue by

Do not proceed by:

- starting Milestone 11 before Milestones 10.1–10.4 define and validate a complete production state whose stationary, internal-wave, and topographic-wave projectors exhaust its coordinates;
- projecting the generator into an APV nullspace;
- empirically symmetrizing the generator;
- fitting a minimum-change closure;
- adding an endpoint energy or enstrophy term without an exact closure identity;
- deleting failing balanced states or hiding external sidebands; or
- interpreting energy conservation alone as validation of the state tendency.

Those operations either answer a different dynamical question or hide the diagnosed incompatibility.

## Conditions for continuing

Milestones 7–10.1 supply the finite-amplitude dense weak oracle, trusted stationary balanced space, complete dense eigensystem, one validated internal-wave projector, controlled rejections of both flat-Dirichlet and slope-compatible per-wavenumber compression, a verified global first-order dressed seed, a converged exact-residual enrichment of the internal-wave block, and an explicit production coordinate contract. The primitive polynomial representation remains the ambient oracle during Milestones 10.2–10.4, including the intervening Milestone 10.2.1, so that omitted directions and truncation errors can be measured independently. It is not the intended online state, and its unresolved algebraic completion will not be carried into matrix-free evolution.

Milestone 10.2.1 specifies that reformulation. With the complete Ritz basis \(Y\), residual \(R=iJ_\gamma Y-H_\gamma Y\Theta\), and \(U=[G_\gamma\;Y]\), it first solves the exact constrained coupled correction

```math
iJ_\gamma\Delta-H_\gamma\Delta\Theta+H_\gamma U\Lambda=-R,
\qquad
U^*H_\gamma\Delta=0.
```

Only after this complementary oracle passes will the same problem be approximated by a block Krylov or Jacobi--Davidson method with shared directions and thick restart. Milestone 10.3 must not begin unless the result is `coupled-block-acceleration` on both split controls. Later increments resolve the bottom/topographic-wave sector and then verify that those dynamical subspaces and the stationary space exhaust every production coordinate. Guard modes remain construction infrastructure, while primitive directions omitted from the declared state remain visible through oracle reconstruction errors and unresolved-energy diagnostics. The unresolved topographic-boundary candidate must not be presented as a physical mode unless a later independent refinement study satisfies every Milestone-9.1 frequency, projector, APV, bottom, strong-equation, guard, padding, and bottom-participation gate. A direction is called numerical only after independent refinement demonstrates nonconvergence. No corrected generator, fitted closure, empirical symmetrization, replacement APV rows, APV-nullspace projection, frequency cutoff, or mode deletion is permitted. Milestone 11 and all later work require separate authorization.

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

The finite-amplitude weak oracle is exposed by [`auditFiniteAmplitudeBoundaryCompleteWeakSystem`](<@WVTerrainEnergyGalerkin/auditFiniteAmplitudeBoundaryCompleteWeakSystem.m>), tested by [`TestWVTerrainEnergyFiniteAmplitudeBoundaryCompleteWeakSystem`](UnitTests/TestWVTerrainEnergyFiniteAmplitudeBoundaryCompleteWeakSystem.m), and recorded in [Milestone 7](milestones.md#milestone-7-finite-amplitude-projected-primitive-dense-gate).

The complete stationary-space oracle is exposed by [`auditCompleteStationaryBalancedSpace`](<@WVTerrainEnergyGalerkin/auditCompleteStationaryBalancedSpace.m>), tested by [`TestWVTerrainEnergyCompleteStationaryBalancedSpace`](UnitTests/TestWVTerrainEnergyCompleteStationaryBalancedSpace.m), and recorded in [Milestone 8](milestones.md#milestone-8-complete-finite-terrain-stationary-balanced-space).

The dense modal gate is exposed by [`auditDensePhysicalEnergyTerrainModes`](<@WVTerrainEnergyGalerkin/auditDensePhysicalEnergyTerrainModes.m>), tested by [`TestWVTerrainEnergyDensePhysicalEnergyTerrainModes`](UnitTests/TestWVTerrainEnergyDensePhysicalEnergyTerrainModes.m), and recorded in [Milestone 9](milestones.md#milestone-9-dense-physical-energy-terrain-modes). Its historical all-mode blocker status motivates [Milestone 9.1](milestones.md#milestone-91-converged-physical-subspace-classification); it does not authorize Milestone 10.

The physical-subspace classifier is exposed by [`auditConvergedPhysicalSubspaces`](<@WVTerrainEnergyGalerkin/auditConvergedPhysicalSubspaces.m>), tested by [`TestWVTerrainEnergyConvergedPhysicalSubspaces`](UnitTests/TestWVTerrainEnergyConvergedPhysicalSubspaces.m), and recorded in [Milestone 9.1](milestones.md#milestone-91-converged-physical-subspace-classification). Its passing result completes the dense classification gate but does not itself implement Milestone 10.

The boundary-complete vertical compression experiment is exposed by [`auditBoundaryCompleteVerticalModeCompression`](<@WVTerrainEnergyGalerkin/auditBoundaryCompleteVerticalModeCompression.m>), tested by [`TestWVTerrainEnergyBoundaryCompleteVerticalModeCompression`](UnitTests/TestWVTerrainEnergyBoundaryCompleteVerticalModeCompression.m), and recorded in [Milestone 9.2](milestones.md#milestone-92-boundary-complete-vertical-mode-compression-oracle). Its `modal-incompatible` result motivates the slope-compatible wave-coordinate experiment in [Milestone 9.3](milestones.md#milestone-93-slope-compatible-wave-coordinate-oracle); it no longer selects the Milestone-10 seed by itself.

The slope-compatible experiment is exposed by [`auditSlopeCompatibleWaveModeCompression`](<@WVTerrainEnergyGalerkin/auditSlopeCompatibleWaveModeCompression.m>), tested by [`TestWVTerrainEnergySlopeCompatibleWaveModeCompression`](UnitTests/TestWVTerrainEnergySlopeCompatibleWaveModeCompression.m), and recorded in [Milestone 9.3](milestones.md#milestone-93-slope-compatible-wave-coordinate-oracle). Its `per-wavenumber-slope-incompatible` result retains the primitive polynomial representation as the ambient oracle for Milestone 9.4. Neither InternalModes checkout was modified.

The global terrain-dressing oracle is exposed by [`auditGlobalFirstOrderTerrainDressing`](<@WVTerrainEnergyGalerkin/auditGlobalFirstOrderTerrainDressing.m>), tested by [`TestWVTerrainEnergyGlobalFirstOrderTerrainDressing`](UnitTests/TestWVTerrainEnergyGlobalFirstOrderTerrainDressing.m), and recorded in [Milestone 9.4](milestones.md#milestone-94-global-first-order-terrain-dressed-block-oracle). Its `global-dressing-seed` result verifies the analytic \(O(h)\) coordinates and exact reduced energy structure while reserving repeated exact-residual enrichment for Milestone 10.

The residual-enrichment oracle is exposed by [`auditResidualEnrichedTerrainModes`](<@WVTerrainEnergyGalerkin/auditResidualEnrichedTerrainModes.m>), tested by [`TestWVTerrainEnergyResidualEnrichedTerrainModes`](UnitTests/TestWVTerrainEnergyResidualEnrichedTerrainModes.m), and recorded in [Milestone 10](milestones.md#milestone-10-residual-enriched-terrain-modes). Its `residual-enrichment-acceleration` result validates the compressed internal-wave block while leaving active bottom candidates in the explicit unresolved complement.

The production-state contract is exposed by [`auditProductionPhysicalStateContract`](<@WVTerrainEnergyGalerkin/auditProductionPhysicalStateContract.m>), tested by [`TestWVTerrainEnergyProductionPhysicalStateContract`](UnitTests/TestWVTerrainEnergyProductionPhysicalStateContract.m), and recorded in [Milestone 10.1](milestones.md#milestone-101-production-physical-state-contract). It preserves the public coefficient layouts, uses the native WaveVortexModel flat families and existing complete bottom inversion, and returns `complete-production-state-contract` for constant and variable stratification and both antialias conventions.

The focused Milestone-10.1 suite passes eight tests, the complete repository suite passes 198 tests, and `checkcode` reports no issues in all 95 MATLAB files.

The focused Milestone-10.2 suite passes eight tests. The latest complete repository suite passes 206 tests with zero failures, and `checkcode` reports no issues in all 98 MATLAB files.

Milestone 10.2.1 is a planned oracle and has no implementation entry point yet. Its exact complementary solve must precede any iterative block implementation; the completed Milestone-10.2 audit remains the independent failure control.

## Mathematical sources

The mathematical repository is `ape-apv-bottom-topography`. Its relevant sources and equation labels are:

- `main.tex`: `eq:finite-terrain-linear-displacement-boundaries`, `eq:finite-terrain-linear-total-energy`, `eq:finite-terrain-linear-apv`, `eq:finite-terrain-linear-potential-enstrophy`, and `eq:finite-terrain-pressure-free-weak-equation`;
- `finite-terrain-projection-problem.tex`: `eq:terrain-projection-apv-green-identity`, `eq:terrain-projection-volume-apv-weak-conservation`, `eq:terrain-projection-discrete-green-commutation`, and `eq:terrain-projection-finite-physical-remainder-decomposition`;
- `terrain-energy-galerkin.tex`: `eq:galerkin-linear-apv`, `eq:galerkin-discrete-compatibility-identities`, `eq:galerkin-boundary-descriptor-pencil`, and `eq:galerkin-first-order-conservation-identities`;
- `boundary-energy-enstrophy.tex`: the exact energy, potential-enstrophy, bottom-Casimir, and complete zero-APV state-space derivations.

The literature repository is hosted privately, so this public repository identifies it by repository name, commit, source file, and equation label rather than publishing its Overleaf project URL.

## Practical alternative

The [`mean-depth-wave-generator`](https://github.com/JeffreyEarly/topographic-forcing/tree/mean-depth-wave-generator) branch, currently at checkpoint [`753432d`](https://github.com/JeffreyEarly/topographic-forcing/tree/753432d), remains the useful near-term implementation. It provides validated first-order bottom wave generation and scattering without direct linear PV generation.

That implementation is a controlled first-order approximation. It does not resolve the exact terrain-energy Galerkin representation problem documented here.
