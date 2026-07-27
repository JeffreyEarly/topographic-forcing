# Topographic forcing research implementations

> **Development branch:** `terrain-energy-galerkin` contains the pressure-free finite-terrain-energy Galerkin prototypes, the completed closure audits, and the boundary-complete flat basis. Milestone 4.8 found that the unmodified finite-terrain forms remain incompatible with simultaneous APV, quadratic-enstrophy, and strong bottom-evolution conservation, so development stops before Milestone 5. The implemented mean-depth generator and scattering classes documented below are retained as a validated baseline; their completed roadmap is archived in [mean-depth-wave-generator-milestones.md](mean-depth-wave-generator-milestones.md).

Potential upstream Fourier and modal-layout additions are prioritized in [Missing WaveVortexModel Infrastructure](MISSING_WAVEVORTEXMODEL_INFRASTRUCTURE.md).

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
