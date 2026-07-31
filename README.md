# Topographic forcing research implementations

> **Milestone-10 scientific checkpoint:** repeated exact-residual enrichment is classified `residual-enrichment-acceleration`. Starting from the globally dressed Milestone-9.4 seed, three corrections reproduce the independent dense internal-wave projector at numerical precision while reducing the full-scale Ritz, bottom, and strong residuals to `1.37e-11`, `1.55e-12`, and below `2e-10`. The exact stationary-plus-dynamical representation uses 165 directions inside the 1027-dimensional primitive ambient space, for compression factor `6.22`. The primitive eigensystem and unresolved complement remain intact, and active bottom directions remain unresolved rather than being deleted or promoted without the Milestone-9.1 continuation gates. Milestone 11 has not begun. Read [Read me first: terrain-energy Galerkin status](READ_ME_FIRST.md) before continuing. The implemented mean-depth generator and scattering classes remain a validated first-order baseline.

Potential upstream Fourier and modal-layout additions are prioritized in [Missing WaveVortexModel Infrastructure](MISSING_WAVEVORTEXMODEL_INFRASTRUCTURE.md).

## Branch-P discrete oracle

The local constant-slope audit can be reproduced with:

```matlab
audit = problem.auditBranchPDiscreteOracle( ...
    bottomSlope=[0.01 0], ...
    polynomialDegree=6);
```

The oracle is independent of the hydrostatic F–G reconstruction. It uses an eta-only bottom coordinate, retains pressure through a constrained saddle solve, and differentiates the saddle system analytically with respect to slope. The coupled physical-energy weak form closes energy and bottom kinematics at roundoff and recovers the analytic flat nonhydrostatic dispersion relation spectrally. Pointwise APV and quadratic potential enstrophy do not close or converge under polynomial refinement, so `audit.status` is `"scientific-blocker"` and `audit.repairCandidate` is `"none"`.

The audit also reduces the existing F–G descriptor directly in native coordinates and maps the same reduced system to the unchanged public coefficient ordering. Those comparisons reproduce the prior defects and show that the eigenvector-based public reconstruction was not their sole cause. The full acceptance evidence and numerical values are recorded in [Milestone 6.1](milestones.md#milestone-61-branch-p-primitive-oracle-and-fg-repair-audit).

## APV-compatible local primitive audit

The follow-up audit compares a compatible mixed polynomial descriptor with vorticity--divergence and explicit-APV formulations:

```matlab
audit = problem.auditAPVCompatiblePrimitive( ...
    bottomSlope=[0 0.01], ...
    horizontalMode=[1 0], ...
    polynomialDegree=6);
```

For the frozen Branch-P equations, direct differentiation gives

```math
\partial_t q_{\boldsymbol s}
=
\frac{i(k s_y-\ell s_x)}{\rho_0D}\,p.
```

The source is nonzero whenever the retained horizontal wavenumber and local slope are not parallel. The energy-weak candidate uses `u,v` in `P_N`, homogeneous mapped vertical velocity in the endpoint-zero part of `P_{N+1}`, displacement in the surface-zero part of `P_{N+1}`, pressure in `P_{N+1}`, and the complete surface-zero space to test vertical momentum. It closes physical energy, the weak equations, continuity, and strong bottom evolution at roundoff and reproduces the analytic APV source. The vorticity--divergence candidate and its explicit-APV coordinate form instead make `q_t=0` at roundoff, but they have nonzero physical-energy and weak-equation defects.

Consequently, `audit.status` is `"mathematical-blocker"` for a cross-slope mode and `"compatible-aligned-special-case"` when `k*s_y-l*s_x=0`. This result is specific to the frozen local approximation; it does not contradict stationary APV in the exact spatially varying terrain equations. The complete evidence is recorded in [Milestone 6.2](milestones.md#milestone-62-apv-compatible-local-primitive-investigation).

## Global small-terrain primitive audit

The global follow-up differentiates a fully coupled signed-Fourier primitive saddle system in the direction supplied by `topographicHeight`:

```matlab
audit = problem.auditGlobalSmallTerrainPrimitive( ...
    polynomialDegree=6, ...
    tangentStep=1e-3);
```

The audit retains pressure until continuity and the bottom equation have been imposed. It constructs the first terrain coefficient analytically and compares it with centered differences of the unexpanded mapped equations. For sinusoidal terrain, the two constructions agree below `3e-10`, the weak equation, physical energy, bottom evolution, and conjugacy close near roundoff, and coupling is confined to the expected neighboring Fourier blocks.

The historical audit gives `audit.status="representation-blocker"`. At the reference resolution its normalized full-grid APV defect is approximately `2.34e-1`. For horizontally interior inputs the defect decreases from `4.44e-3` to `5.38e-4` as the polynomial degree increases from two to six, while Fourier-edge inputs remain near `4.14e-1`. The edge value is an external sideband on the oversampled grid and does not test a consistently restricted APV law. No APV projection, energy symmetrization, or fitted correction is applied. The complete evidence is recorded in [Milestone 6.3](milestones.md#milestone-63-global-small-terrain-primitive-compatibility-oracle).

## Coupled volume–boundary PV frequency oracle

Milestone 6.4 validates the dynamical interpretation in a controlled single-wavenumber QG problem:

```matlab
audit = problem.auditCoupledPVFrequencyOracle( ...
    horizontalMode=[1 0], ...
    volumePVGradient=[0 2e-11], ...
    bottomSlope=[0 0.01], ...
    polynomialDegree=20);
```

The direct $L^2\oplus\mathbb C$ volume–boundary inversion conserves physical energy and the applicable signed pseudoenstrophy at roundoff. Its leading physical frequencies agree with an independent Chebyshev discretization of Yassin's eigenvalue-dependent endpoint problem to `4.82e-11`. The $f$-plane limit contains an active bottom mode with zero volume APV, and the existing balanced-plus-bottom Galerkin basis represents the leading modes with a residual that decreases under vertical refinement.

This passing local oracle does not establish the projected primitive APV law. Milestone 6.5 isolates the horizontal projection by evolving volume and boundary PV directly.

## Periodic coupled volume–boundary PV oracle

Milestone 6.5 applies the coupled-PV formulation to periodic terrain:

```matlab
audit = problem.auditPeriodicCoupledPVOracle( ...
    polynomialDegree=12);
```

For the resting $f$-plane problem it evolves

```math
\partial_tq=0,
\qquad
\partial_tr_b=-\mathcal P J(\psi_b,fh),
\qquad
\psi=\mathcal G[q,r_b],
```

where $\mathcal P$ is the orthogonal projection onto the retained signed Fourier layout. The exact mode-number convolution agrees with an independently oversampled pseudospectral Jacobian below `8e-16`. Physical energy, stationary volume APV, physical potential enstrophy, projected bottom evolution, and Fourier conjugacy close at roundoff for both interior and Fourier-edge inputs.

Edge inputs still produce nonzero sidebands outside the retained state. Those sidebands are explicitly measured and discarded by the Galerkin projection rather than aliased back into the state or counted as an internal APV defect. The existing balanced-plus-bottom basis represents the active topographic modes increasingly accurately under simultaneous oracle and vertical refinement.

This establishes a closed periodic QG oracle and identifies an explicit PV-coordinate formulation as the promising horizontal representation. Milestone 6.6 tests that idea in a complete hybrid descriptor; the result below explains why Milestone 7 was inactive at that stage.

Milestone 6.6 tests the smallest such hybrid descriptor:

```matlab
audit = problem.auditHybridPrimitivePVOracle( ...
    zonalMode=1,polynomialDegree=4);
```

The wave, projected-volume-APV, and bottom rows form a complete coordinate system and close at roundoff. Nevertheless, the resulting first-order generator fails the unmodified primitive weak equation and physical energy. The full sampled APV audit also contains external meridional sidebands. The oracle therefore classifies the result as `primitive-equivalence-blocker`: using projected APV as a replacement coordinate is not by itself an equivalent discretization of the primitive equations.

## Dealiased projected primitive tangent oracle

Milestone 6.7 returns to the unmodified primitive energy and exchange rows:

```matlab
audit = problem.auditDealiasedProjectedPrimitiveTangent( ...
    trustedModeBounds=[1 0], ...
    supportModeBounds=[1 1;1 2;1 3], ...
    polynomialDegrees=[4;8;12], ...
    paddingFactors=[2;3]);
```

All terrain products use one zero-pad, multiply, and adjoint-restrict operation. Exact mode-number convolution and the pseudospectral action agree to `1.19e-15`; the primitive weak, physical-energy, projected-bottom, conjugacy, and padding defects are all below their acceptance tolerances. External edge sidebands are measured separately and are not aliased into the retained equations.

The trusted-band APV defect decreases from `5.45e-4` to `1.85e-4` to `1.01e-4` at polynomial degrees `4`, `8`, and `12`. The corresponding potential-enstrophy defect decreases from `3.34e-4` to `1.12e-4` to `6.08e-5`. These rates and final values fail the required factor-four and `1e-8` gates, while horizontal support refinement alone leaves the result unchanged. The oracle therefore returns `audit.status="nonconvergent-projected"`. No projected replacement rows or corrective closure are used. This result kept finite-amplitude work inactive until Milestone 6.8 supplied the terrain-dependent geostrophic inclusion.

The complete repository suite passes 125 tests with zero failures.

## Boundary-complete weak terrain eigenproblem

Milestone 6.8 uses the terrain-dependent geostrophic test states derived from the primitive APV Green identity:

```matlab
audit = problem.auditBoundaryCompleteWeakEigenproblem( ...
    trustedModeBounds=[1 0], ...
    supportModeBounds=[1 2], ...
    polynomialDegrees=[4;8;12], ...
    paddingFactor=2);
```

The oracle constructs both the flat and first-terrain geostrophic inclusions, evaluates APV moments independently from the strong mapped diagnostic, and solves the unmodified tangent generalized eigenproblem. The first-terrain Green-identity and stationary-row defects are `8.01e-15` and `6.63e-14`; independently evaluated APV evolution is `1.11e-13`. Physical energy and projected bottom evolution retain defects `4.28e-15` and `2.47e-12`.

The result is `compatible-boundary-complete-weak-oracle`. APV follows from the original primitive weak equations once the terrain derivative of the test space is retained. No APV replacement rows, empirical corrections, symmetrization, APV-nullspace projection, or mode deletion are used. Milestone 7 extends this construction without making a small-terrain approximation.

The complete repository suite passes 130 tests with zero failures.

## Finite-amplitude boundary-complete primitive gate

Milestone 7 constructs the finite-terrain energy and exchange forms and the terrain-dependent geostrophic inclusion at the same amplitude:

```matlab
audit = problem.auditFiniteAmplitudeBoundaryCompleteWeakSystem( ...
    trustedModeBounds=[1 0], ...
    supportModeBounds=[1 1;1 2;1 4], ...
    scalarPolynomialDegree=2, ...
    primitivePolynomialDegrees=[2;3;5], ...
    paddingFactors=[2;3], ...
    terrainScales=[0.25;0.5;1]);
```

The coefficient evolution is the unmodified primitive weak system

```math
H_\gamma\dot{\boldsymbol A}
=
J_\gamma\boldsymbol A,
\qquad
L_\gamma
=
H_\gamma^{-1}J_\gamma.
```

Pressure recovery and the strong bottom equation are independent diagnostics; neither replaces an evolution row. At the finest constant-stratification reference calculation, the weak equation and physical-energy defects are `2.62e-16` and `1.59e-15`. The APV Green identity, stationary geostrophic row, weak APV evolution, independently evaluated strong APV evolution, and their agreement are respectively `4.33e-15`, `9.86e-17`, `1.32e-14`, `4.47e-15`, and `4.43e-15`.

The trusted-band bottom defect decreases from `3.47e-5` to `2.89e-7` to `2.01e-11` as the primitive vertical and horizontal support spaces are enriched. Padding factors two and three agree to `2.18e-12`; external sidebands remain reported separately. Variable stratification reaches a geostrophic-state representation defect of `7.18e-11`, APV evolution of `3.21e-15`, and bottom defect of `1.45e-11` under joint horizontal and vertical enrichment.

The exact finite-amplitude forms differ from their flat-plus-tangent approximations by \(O(h^2)\): halving terrain amplitude reduces the energy, exchange, and geostrophic-inclusion remainders by a factor of approximately four. The oracle returns `compatible-finite-amplitude-weak-oracle`. It does not construct terrain modes or advance a state in time.

The complete repository suite passes 139 tests with zero failures, and `checkcode` reports no issues in all 69 MATLAB files.

## Complete finite-terrain stationary balanced space

Milestone 8 extends the geostrophic inclusion to nonzero bottom values and enforces the finite-terrain tangency condition:

```matlab
audit = problem.auditCompleteStationaryBalancedSpace( ...
    trustedModeBounds=[1 1], ...
    supportModeBounds=[1 2;1 3;1 4], ...
    primitivePolynomialDegrees=[2;3;4], ...
    paddingFactors=[2;3]);
```

The oracle verifies the full APV Green identity, including

```math
c_b[\psi]
=
f\hat\eta_b+v_bh_x-u_bh_y,
```

and returns `complete-stationary-space-oracle`. At the reduced constant-stratification reference resolution, geostrophic representation, the complete Green identity, stationary rows, trusted bottom tangency, conjugacy, and padding agreement are all below `1.4e-13`. Six non-tangent bottom-streamfunction directions remain in the full state for dynamical classification.

The raw exchange matrix contains extremely low-frequency active boundary directions next to its exact stationary kernel. Consequently, the Milestone-9 audit starts from the tangency-defined Milestone-8 space rather than a hard SVD or frequency cutoff, and retains the complete physical-energy complement.

The complete repository suite passes 147 tests with zero failures.

## Dense physical-energy terrain-mode audit

Milestone 9 uses the Milestone-8 tangency construction at every vertical refinement, retains every non-tangent bottom direction, and solves

```math
iJ_\gamma\boldsymbol c_n
=
\omega_nH_\gamma\boldsymbol c_n
```

by a Cholesky-scaled dense eigensolve:

```matlab
audit = problem.auditDensePhysicalEnergyTerrainModes( ...
    trustedModeBounds=[1 0], ...
    supportModeBounds=[1 2;1 3;1 4], ...
    stationaryPolynomialDegree=2, ...
    primitivePolynomialDegrees=[2;3;4], ...
    paddingFactors=[2;3], ...
    terrainScales=1);
```

For the constant-stratification sinusoidal reference, the generalized eigen-residual, physical-energy orthogonality, Fourier-conjugacy defect, and non-tangent-bottom retention defect are respectively below `2e-13`, `5e-15`, `7e-15`, and `3e-15`. Flat and uniform-depth first-mode frequencies converge below `1e-10` relative error.

The original all-mode classification gate does not pass. The physical-energy complement contains backward-error-indistinguishable zero-frequency directions outside the declared Milestone-8 stationary sector. Four bottom-dominated candidates have frequencies resolved by more than $10^{12}$ times their backward uncertainty, but their reference trusted APV, bottom, and strong-equation defects are approximately `1.60e-1`, `4.19e-6`, and `3.67e-1`.

A fixed-support degree-\(4,6,8\) study shows that a tracked internal branch near \(4.43f\) converges to APV, bottom, and strong defects `7.56e-9`, `5.22e-10`, and `6.34e-6`. The subinertial branch near \(0.169f\) improves to APV and strong defects `6.39e-3` and `5.82e-2`, but its bottom defect remains near `5e-6`; it therefore remains unresolved rather than demonstrated numerical. Every candidate and unclassified direction remains in the returned eigensystem.

The audit retains the historical status `dense-modal-classification-blocker`. Milestone 9.1 replaces the all-mode interpretation with nested physical-energy spectral-projector convergence and an explicit unresolved remainder. Milestone 10 remains inactive until that gate passes.

The complete repository suite passes 154 tests with zero failures, and `checkcode` reports no issues in all 76 MATLAB files.

## Converged physical-subspace classification

Milestone 9.1 compares complete conjugate-closed spectral subspaces rather than individual eigenvectors:

```matlab
audit = problem.auditConvergedPhysicalSubspaces( ...
    trustedModeBounds=[1 0], ...
    supportModeBounds=[1 2;1 3;1 4], ...
    stationaryPolynomialDegree=4, ...
    primitivePolynomialDegrees=[4;6;8], ...
    paddingFactors=[2;3], ...
    terrainScales=[0.125;0.25;0.5;1]);
```

For the \(2.5\ {\rm m}\) constant-stratification sinusoidal terrain control, the nine-dimensional stationary projector has a final nested principal sine of `1.78e-12`; its independent fixed-degree guard and padding defects are `1.67e-12` and `7.82e-13`.

One four-dimensional internal-wave subspace near \(4.4283f\) is validated. At degrees four, six, and eight its strong residual decreases from `1.78e-3` to `2.75e-5` to `2.42e-7`; its final APV and bottom defects are `1.07e-10` and `1.46e-13`. Its independent guard and padding defects are `5.21e-9` and `1.58e-13`.

The topographic-boundary candidate is not promoted to a physical mode: its final guard defect remains `2.10e-2`. It stays in the unresolved projector along with every other unvalidated direction. The finest decomposition is

```math
\dim\mathcal V_N=703,
\qquad
\dim\mathcal G_N=9,
\qquad
\dim\mathcal W_N=4,
\qquad
\dim\mathcal R_N=690.
```

The three projectors are physical-energy orthogonal, Fourier-conjugate closed, and complete within `3.1e-13`. The oracle returns `physical-subspace-classification-oracle`; Milestone 10 remains a separate, unimplemented residual-enrichment step.

The complete repository suite passes 162 tests with zero failures, and `checkcode` reports no issues in all 79 MATLAB files.

## Boundary-complete vertical-mode compression

Milestone 9.2 tests a smaller vertical representation built from fixed-\(\kappa\) nonhydrostatic waves, signed-Robin APV-bearing geostrophic modes, explicit zero-APV bottom inversions, and the complete mean sector:

```matlab
audit = problem.auditBoundaryCompleteVerticalModeCompression( ...
    trustedModeBounds=[1 0], ...
    supportModeBounds=[1 2;1 3;1 4], ...
    stationaryPolynomialDegree=4, ...
    primitivePolynomialDegrees=[4;6;8], ...
    comparisonPolynomialDegree=12, ...
    modalCounts=[1;2;4;8], ...
    robinLengthRatios=[-1/8;-1/4;-1/2;Inf], ...
    paddingFactors=[2;3], ...
    terrainScales=[0.125;0.25;0.5;1], ...
    internalModesEVPOrders=[128;256]);
```

The one-dimensional eigenproblems come from the isolated `internal-modes-evp` checkout pinned to commit `df86687`; WaveVortexModel continues to use the main InternalModes dependency. Every candidate column is first embedded in the primitive oracle, missing stationary directions are appended by physical-energy orthogonal completion, and the reduced system uses the unchanged finite-terrain \(H_\gamma\) and \(J_\gamma\).

The family converges toward the validated internal-wave projector, but at eight wave and Robin modes its projector defect is `7.86e-7`, bottom-evolution defect is `3.10e-3`, and compression factor is only `1.0914`. Robin-length changes produce projector differences up to `2.87e-5`; this is retained as finite-order convergence evidence and is no longer treated as a hard physical gate. The audit remains `modal-incompatible` because the flat Dirichlet wave coordinates do not close bottom evolution before stationary completion removes useful compression. Milestone 9.3 subsequently tests slope-compatible hydrostatic and nonhydrostatic replacements while retaining every unresolved direction and the unchanged primitive \(H_\gamma,J_\gamma\) evolution. Milestone 10 has not begun.

The complete repository suite passes 168 tests with zero failures, and `checkcode` reports no issues in all 83 MATLAB files.

## Slope-compatible wave-coordinate compression

Milestone 9.3 replaces only the wave coordinates in the Milestone-9.2 experiment. For each horizontal wavenumber, it compares hydrostatic and nonhydrostatic first-order descriptors that leave the interior flat and retain a constant reference slope only in the active bottom-displacement row:

```matlab
audit = problem.auditSlopeCompatibleWaveModeCompression( ...
    trustedModeBounds=[1 0], ...
    supportModeBounds=[1 2;1 3;1 4], ...
    stationaryPolynomialDegree=4, ...
    primitivePolynomialDegrees=[4;6;8], ...
    comparisonPolynomialDegree=12, ...
    modalCounts=[1;2;4;8], ...
    robinLengthRatios=[-1/8;-1/4;-1/2;-1;Inf], ...
    referenceSlopeFactors=[0;0.5;1;sqrt(2)], ...
    dynamics=["hydrostatic";"nonhydrostatic"], ...
    paddingFactors=[2;3], ...
    terrainScales=[0.125;0.25;0.5;1], ...
    internalModesEVPOrders=[128;256]);
```

The Robin length and reference slope are trained jointly at padding two and then frozen. Finite descriptor modes are classified using homogeneous generalized-Schur data; the wave subspace is tracked by physical-energy projector overlap rather than a frequency cutoff. The final periodic-terrain calculation continues to use the unchanged primitive \(H_\gamma,J_\gamma\) forms.

Both families select \(s_{\rm ref}/s_{\rm rms}=0.5\) and \(\ell_b/D=\infty\). Their local descriptor residuals are about `5e-14`, active-bottom residuals are below `9e-13`, and backward errors are below `2e-17`. Globally, both give an internal-projector defect of `3.251e-7`, bottom-evolution defect of `4.597e-3`, strong residual of `2.073e-6`, and compression factor `1.075`. The result is `per-wavenumber-slope-incompatible`. Hydrostatic and nonhydrostatic reference dynamics are indistinguishable at this resolution; the unresolved bottleneck is the global periodic-terrain stationary/bottom coupling, not vertical inertia in the local wave coordinate.

The complete repository suite passes 173 tests with zero failures, and `checkcode` reports no issues in all 86 MATLAB files.

## Global first-order terrain dressing

Milestone 9.4 uses the verified analytic terrain derivatives to construct the horizontal sidebands and vertical corrections of complete flat signed-frequency blocks:

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

The first-order block calculation uses \(H_1,J_1,R_1,G_1\) only to choose coordinates. At every terrain amplitude, the reduced eigensystem is assembled from the unchanged exact primitive \(H_\gamma,J_\gamma\) forms and the exact stationary space. At the degree-12 comparison resolution, the complete flat zero-frequency block has 377 directions: 345 remain in its tangent stationary sector and 32 form backward-error-resolved nonzero first-order pairs. These pairs remain topographic-boundary candidates rather than claimed physical modes.

The analytic and centered \(H_1,J_1,R_1\) calculations agree below `1e-9`; the complementary correction, Fourier selection, active first-order bottom equation, and stationary Green identities close below `1e-10`. Dressed weak and bottom residuals have observed order `2.000`, and the dressed physical projector improves by at least a factor of two over the same undressed selective span. The full-amplitude internal APV defect is `1.70e-11`, but the physical projector and strong residual still exceed the Milestone-9.1 gates. The result is therefore `global-dressing-seed`: it authorizes repeated exact-residual enrichment as a separate Milestone 10, not a finite-amplitude mode claim.

The complete repository suite passes 181 tests with zero failures, and `checkcode` reports no issues in all 89 MATLAB files.

## Terrain-energy Galerkin flat oracle

`WVTerrainEnergyGalerkin` is a standalone linear scientific system, not a `WVForcing`. It uses ordinary hydrostatic wave-vortex modes as coordinates, adds one complete balanced bottom-inversion state per retained nonzero horizontal Fourier coefficient, and projects the complete flat nonhydrostatic weak equations into that mixed basis. The public bottom coefficient remains the bottom displacement.

```matlab
problem = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
    topographicHeight=h, ...
    verticalModeIndices=wvt.j, ...
    horizontalOversamplingFactor=2);
```

The state uses a full-complex Fourier representation internally. Explicit maps connect it to WaveVortexModel's nonredundant real-field layout without forcing arbitrary complex eigenvectors through a symmetric inverse transform. The source transform's `shouldAntialias` convention is inherited and may be either true or false.

For every retained horizontal wavenumber, the flat oracle constructs the pressure-free forms

```math
E_0\dot{\boldsymbol a}=J_0\boldsymbol a,
\qquad
iJ_0\boldsymbol c=\omega E_0\boldsymbol c.
```

The raw matrices satisfy the required Hermitian identities at roundoff. Constant-stratification frequencies recover the analytic nonhydrostatic dispersion relation, and arbitrary-stratification frequencies and eigenfunctions converge to the directly computed wavenumber-dependent modes as hydrostatic vertical coordinates are added. Nonzero-frequency modes have negligible flat QGPV and bottom displacement.

The dynamical solve does not uniquely select vectors inside its degenerate stationary subspace. The flat oracle therefore also forms

```math
Z_0=Q_0^*W_\xi Q_0
```

and diagonalizes physical potential enstrophy within that subspace. The resulting common basis separates APV-bearing geostrophic modes from one stationary zero-APV bottom inversion for every nonzero horizontal wavenumber. Each flat block exposes the physical-energy projector onto that bottom state. No generalized bottom energy or additional boundary invariant is introduced.

The dense finite-terrain oracle evaluates the mapped geometry and stratification on a common horizontally oversampled grid and stores

```math
E_\gamma,\qquad J_\gamma,\qquad Q_\gamma.
```

The raw energy and exchange forms retain their Hermitian and skew-Hermitian structure before roundoff cleanup. Their flat limit reproduces the per-wavenumber oracle, and the constant-$\gamma$ problem agrees with an independent flat transform of physical depth $H=\gamma D$.

## Raw finite-terrain compatibility audit

Audit the unmodified generator without constructing a closure:

```matlab
audit = problem.auditFiniteTerrainCompatibility();
```

The audit forms

```math
L_\gamma=E_\gamma^{-1}J_\gamma,
\qquad
Z_\gamma=Q_\gamma^*Q_\gamma
```

and directly tests energy, pointwise APV, quadratic potential enstrophy, strong bottom evolution, and Fourier conjugacy. It reports both pre-restoration and operational matrices, deterministic random-state tendencies, APV rank, flat-common-subspace residuals, vertical APV rows, horizontal bottom rows, and the unresolved part of the oversampled bottom product. It never re-skews, projects, or corrects the generator.

Flat and uniform-depth references pass within $10^{-12}$. For the documented 20 m sinusoidal terrain at resolution $[4,4,5]$ and oversampling factor two, the normalized APV, bottom, and enstrophy defects are respectively

```math
1.89\times10^{-7},
\qquad
2.99\times10^{-2},
\qquad
4.78\times10^{-8}.
```

Energy and conjugacy remain at roundoff. Increasing horizontal oversampling from two to three does not change the incompatible defects, and the finest horizontal refinement pair changes them by less than $13\%$. The result therefore satisfies the roadmap's incompatible exit rather than its compatible-and-convergent exit. Milestone 5 is blocked; no corrected closure or relaxed invariant is substituted.

## Constrained-closure audit

The dense finite-terrain forms conserve their discrete energy, but energy skew-Hermiticity alone does not guarantee exact discrete APV conservation or the resolved strong bottom equation. Audit those additional requirements with:

```matlab
audit = problem.auditConstrainedClosure();
```

The audit transforms the raw generator into energy coordinates and asks whether a skew-Hermitian operator can satisfy

```math
\overline QK_c=0,
\qquad
\overline BK_c=\overline R.
```

It preserves the raw Galerkin matrices, reports the numerical APV rank and bottom-constraint compatibility, and computes the minimum-change closure only when the complete constraint system is feasible. The bottom equation is enforced in the retained Fourier space; the unrepresented part of the oversampled terrain product is reported separately.

Flat and uniform-depth reference problems pass: the APV nullity equals the 56-dimensional flat wave space in the automated reference case, and the minimum-change closure is the raw generator to roundoff. The corresponding 20 m sinusoidal-terrain problem does not pass. Under the documented numerical-rank criterion, its discrete APV nullity is 35, and its minimum bottom-constraint residual is approximately $2.46\times10^{-5}$. The audit therefore returns `status="incompatible"` and no constrained terrain generator.

The weaker quadratic-invariant question is available separately:

```matlab
audit = problem.auditQuadraticInvariantClosure();
```

This audit uses independent real physical coordinates and seeks a skew-symmetric energy-coordinate generator satisfying

```math
[G_\gamma,K_c]=0,
\qquad
\overline B K_c=\overline R,
\qquad
G_\gamma=S^{-*}Q_\gamma^*Q_\gamma S^{-1}.
```

The commutator conserves total discrete quadratic potential enstrophy without requiring every APV sample to remain fixed. Flat and uniform-depth problems again return the raw generator with zero correction. The sinusoidal reference remains incompatible: at 20 m amplitude its minimum bottom residual is $2.46\times10^{-5}$, or $88.5\%$ of the bottom target, and 23 of 44 resolved enstrophy eigenspaces fail the bottom constraint. The residual scales linearly with terrain amplitude and remains between approximately $81\%$ and $89\%$ under the tested horizontal, vertical, and oversampling refinements. These results describe the superseded displacement-only state and are retained as historical diagnostics; the raw Milestone-4.8 audit above is the authoritative result for the boundary-complete basis.

`topographic-forcing` provides a fast, first-order bottom wave generator for WaveVortexModel. The formulation retains the ordinary rigid-lid wave--vortex basis and represents weak topography through the mean-depth bottom condition.

For a prescribed, horizontally uniform barotropic current, the bottom velocity is

```math
g_b(\boldsymbol{x},t)
=
\boldsymbol{U}_{\mathrm{bt}}(t)\boldsymbol{\cdot}\nabla_H h(\boldsymbol{x}).
```

The wave forcing is obtained directly from the bottom pressure of each complete phase-inclusive mode:

```math
\dot A_\alpha
=
\frac{1}{A E_\alpha}
\int_A p_{\alpha,d}^*(\boldsymbol{x},t)g_b(\boldsymbol{x},t)\,dA,
\qquad
\alpha\in\{+,-\}.
```

This construction adds only to the wave coefficients $A_+$ and $A_-$. It therefore produces no direct linear interior QGPV tendency and leaves the balanced forcing coefficient $F_0$ unchanged. It requires no finite-terrain pressure solve, modal terrain matrix, or artificial bottom-localized vertical envelope.

The prescribed barotropic current is an external energy reservoir. Wave energy is not conserved by itself; its rate of increase must equal the bottom pressure work supplied by the prescribed current.

## Quick start

Construct a `WVTransformBoussinesq`, prescribe the terrain and complex barotropic velocity amplitude, and register the spectral forcing:

```matlab
forcing = WVBottomWaveGenerationForcing(wvt, ...
    topographicHeight=h, ...
    barotropicVelocityAmplitude=[0.05; 0]);
wvt.removeAllForcing();
wvt.addForcing(forcing);
```

The velocity amplitude is the complex two-component vector $\widehat{\boldsymbol U}_{\mathrm{bt}}$ in meters per second. The default frequency is M2; `frequency`, `rampDuration`, `startTime`, and `name` are optional constructor arguments.

The constructor precomputes the bottom-pressure projection on the transform's native spectral layout. It supports any stationary stratification represented by `WVTransformBoussinesq`. Each subsequent forcing call evaluates the prescribed current, combines two response arrays per wave branch, and applies WaveVortexModel's interaction phases. There is no runtime pressure solve, FFT, spatial projection, or modal coupling matrix. Transforms with either value of `shouldAntialias` are supported.

`forcingWithResolutionOfTransform` spectrally transfers the terrain and rebuilds all modal responses for the new transform. The physical forcing configuration is also included when its parent transform or model is written to NetCDF; transform-derived response arrays are rebuilt after restoration.

Milestones 1--6, 8, and 9 of the archived [mean-depth development roadmap](mean-depth-wave-generator-milestones.md) are implemented on the `mean-depth-wave-generator` branch. The optional comparison in Milestone 7 was intentionally skipped.

A separate [second-order roadmap](second-order-milestones.md) records a possible extension of the mean-depth scattering approach. That hierarchy is planned work and is not part of the current forcing classes or the terrain-energy Galerkin roadmap.

The scientific and computational status of the earlier repositories and branches is summarized in [PRIOR_APPROACHES.md](PRIOR_APPROACHES.md).

## Exact sinusoidal-ridge example

Run the known-solution benchmark and create its two diagnostic figures with:

```matlab
benchmark = SinusoidalRidgeWaveGenerationBenchmark;
```

The example uses a uniform M2 current over one sinusoidal terrain component. It compares three adaptive `ode78` integrations with the analytically integrated interaction coefficients, verifies that wave-energy growth equals the work done by bottom pressure, and shows the generated vertical velocity and displacement fields. Pass `shouldMakeFigures=false` for diagnostics without graphics or set `resolution` and `relativeTolerances` explicitly.

At the automated reference resolution, the coefficient errors for tolerances $10^{-6}$, $10^{-8}$, and $10^{-10}$ are approximately $6.4\times10^{-8}$, $5.6\times10^{-10}$, and $5.8\times10^{-12}$. The tight-run energy/work error is approximately $2.5\times10^{-12}$. The balanced tendency is exactly zero, while an independent physical-space calculation confirms the linear-QGPV source vanishes to the accuracy of the discrete derivative transforms.

## Goff abyssal-hill example

Generate deterministic periodic Goff topography with:

```matlab
[virtualDepth,h,diagnostics] = ...
    WVBottomWaveGenerationForcing.goffAbyssalHillTopography( ...
        wvt,minimumWavelength=40e3);
```

The generator uses a local random stream, returns upward-positive zero-mean terrain with the requested post-filter RMS, and reports realized slope and radial-spectrum diagnostics. It does not change MATLAB's global random state.

Run the variable-stratification broadband example with:

```matlab
result = GoffAbyssalHillWaveGenerationExample;
```

The example drives a seed-2023, $100$ m RMS Goff field with a uniform M2 current and shows the terrain spectrum, wave-energy distribution, and bottom-work budget. Use `BottomWaveGenerationPerformanceBenchmark` to report construction, storage, and forcing-application costs at several resolutions.

## Autonomous first-order wave scattering

`WVBottomWaveScatteringForcing` is a separate forcing for scattering an existing wave field. It evaluates the wave-only mean-depth bottom velocity

```math
g_b
=
\boldsymbol u_{H,d}\boldsymbol{\cdot}\nabla_Hh
-
h\,\partial_z w_d
```

and applies the same bottom-pressure projection to the two wave branches:

```matlab
scattering = WVBottomWaveScatteringForcing( ...
    wvt,topographicHeight=h);
wvt.removeAllForcing();
wvt.addForcing(scattering);
```

The ordinary forcing call reconstructs only the three required bottom fields, performs horizontal pseudospectral products, and projects the resulting bottom velocity. It does not evaluate full three-dimensional fields, solve for pressure, assemble a modal terrain matrix, or modify the balanced tendency.

Bottom displacement is a postprocessed diagnostic rather than a model state. Save ordinary WaveVortexModel output and integrate the saved bottom velocity with:

```matlab
diagnostics = ...
    WVBottomWaveScatteringForcing.bottomDisplacementFromFile( ...
        "scattering-output.nc");
```

The method reads `t`, `Ap`, and `Am` from the standard `wave-vortex` group, reconstructs $g_b$ at the selected accepted output times, and applies cumulative trapezoidal quadrature to $\partial_t\eta_d=g_b$. Use `iTime` to select output records, `forcingName` when the file contains more than one scattering forcing, and `initialBottomDisplacement` for a nonzero initial condition.

Two examples exercise the autonomous formulation:

```matlab
uniform = UniformDepthWaveScatteringBenchmark;
ridge = SinusoidalRidgeWaveScatteringExample;
gaussian = GaussianRidgeWaveScatteringExample;
```

The uniform-depth benchmark verifies the first-order frequency correction against the exact depth-$D-h_0$ dispersion relation. The sinusoidal-ridge example writes standard output, displays the scattered sidebands and energy exchange, and reconstructs the bottom displacement from that file.

The Gaussian-ridge example follows a localized rightward mode-one M2 wave packet as it crosses a gentle subcritical ridge. Its default terrain has $h_0/D=0.05$ and $\max|h_x|/\mu=0.05$, keeping the calculation in the controlled first-order regime. It uses adaptive `ode78`, returns the sampled `Ap`, `Am`, and `A0` trajectory and diagnostics in memory, and creates no output file unless one is requested:

```matlab
result = GaussianRidgeWaveScatteringExample( ...
    outputPath="gaussian-ridge.nc");
```

The output is standard restartable WaveVortexModel NetCDF. A separate renderer reads only that saved file and creates an $x$--$z$ movie and PNG poster:

```matlab
movie = GaussianRidgeWaveScatteringMovie( ...
    "gaussian-ridge.nc",field="u");
```

The movie masks the reconstructed field beneath the physical bottom, uses wet points only for its fixed color scale, and shows the first-order depth-integrated wave-energy profile alongside the evolving rightward, reflected, and higher-mode energy. Its default section uses a labelled 20-times vertical exaggeration, configurable with `verticalExaggeration`. The renderer also supports `field="w"` and `field="eta"`, a subset of saved records through `iTime`, and explicit `videoPath` and `posterPath` values. Existing files are never replaced unless `shouldOverwriteExisting=true`.

## Scientific scope

Both forcings are accurate through first order in terrain height. The prescribed generator supports broadband terrain with arbitrary stationary stratification, transform-resolution rebuilding, and restart persistence. The autonomous forcing adds wave--wave scattering with the same production behavior. Its strict perturbation hierarchy conserves the first-order physical energy; when the first-order operator is iterated autonomously, the unresolved energy tendency is $O(h^2)$.

The planned strict second Born model will retain the zeroth-, first-, and second-order coefficients separately and stop before uncontrolled higher-order feedback. See the [second-order milestones](second-order-milestones.md) for its scientific gates and implementation sequence.

The following remain outside the initial proof of concept:

- exact finite-amplitude terrain dynamics;
- nonlinear terrain-aware advection;
- independent bottom-buoyancy dynamics;
- backreaction on a dynamically evolving barotropic tide;
- exact finite-amplitude energy conservation under autonomous first-order scattering.

WaveVortexModel and InternalModes remain external dependencies and are not modified or vendored by this repository. The terrain-energy Galerkin branch additionally uses the boundary-mode solver from `InternalModesEVP`, currently pinned to commit `df86687e91faa31bf65941299062d125a96904b1` in an isolated sibling checkout named `internal-modes-evp`. This keeps WaveVortexModel on the main InternalModes implementation while exposing the uniquely named `IMSurfaceGeostrophicModes` and `IMSolverSpectral` classes.

Create the isolated checkout without changing the main InternalModes working tree:

```bash
git -C ../internal-modes worktree add --detach ../internal-modes-evp df86687e91faa31bf65941299062d125a96904b1
```

Run the automated suite with:

```matlab
results = runTests;
```

The runner resolves WaveVortexModel, main InternalModes, and InternalModesEVP from explicit `waveVortexModelRoot`, `internalModesRoot`, and `internalModesEVPRoot` options; the corresponding environment variables; or the sibling repositories.
