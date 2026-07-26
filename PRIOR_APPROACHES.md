# Prior topographic forcing approaches

This note records the scientific and computational status of the terrain-forcing approaches that preceded the mean-depth bottom wave generator. The reviewed snapshots are fixed so that later changes in the other repositories do not obscure what informed this roadmap.

| Approach | Reviewed snapshot | Scientific target | Role going forward |
| --- | --- | --- | --- |
| Pseudo-topography | [`Pseudo-topography` at `c9b65c7`](https://github.com/JeffreyEarly/Pseudo-topography/tree/c9b65c7) | Represent omitted bottom motion as a prescribed displacement tendency projected by WaveVortexModel | External baseline and source of terrain, tide, restart, and example patterns |
| Conservative terrain forcing | [`conservative-terrain-forcing` at `c543ed8`](https://github.com/JeffreyEarly/conservative-terrain-forcing/tree/c543ed8) | Replace the raw linear terrain response by a modal surrogate preserving flat-reference energy and potential enstrophy | Reference for spectral forcing, phase handling, persistence, and operator diagnostics |
| Mapped strong-form forcing | [`topographic-forcing` at `a2c2e91`](https://github.com/JeffreyEarly/topographic-forcing/tree/a2c2e91) | Apply exact-in-height mapped terrain terms through the ordinary flat spatial projection | Preserved scientific checkpoint; not the active formulation |
| Terrain-weighted weak forcing | [`terrain-weighted-weak-form` at `3c75ba8`](https://github.com/JeffreyEarly/topographic-forcing/tree/3c75ba8) | Evolve the mapped system using the finite-terrain energy inner product | Preserved alternative weak formulation; not the active implementation |

## Pseudo-topography

`PseudoTopographyForcing` supplies a spatial displacement tendency and lets the ordinary WaveVortexModel transform project it onto wave and balanced modes. The repository contains mature rough-topography generation, harmonic barotropic forcing, restart support, eddy--tide examples, and extensive energy and PV diagnostics.

That projection retains a nonzero \(F_0\); the repository's unit tests explicitly verify that the balanced contribution is present. A vertically localized displacement tendency also has a nonzero linear QGPV source before wave-only projection. This makes the existing implementation a useful physical and numerical comparison, but not the desired no-direct-PV wave generator.

The new implementation should reuse or adapt the following patterns only after its scientific kernel is validated:

- deterministic periodic terrain generation;
- barotropic-tide phase and startup ramps;
- restart configuration and external-path test setup;
- rough-topography examples and diagnostic reporting.

It should not reuse the artificial vertical envelope or the unfiltered balanced projection. The new boundary projection obtains the vertical modal distribution directly from output-mode bottom pressure.

## Conservative terrain forcing

`conservative-terrain-forcing` constructs an exact-in-height raw linear terrain response, energy-scales its real modal matrix, removes unequal enstrophy-to-energy couplings, and takes the skew-symmetric part. The resulting surrogate preserves the flat-reference quadratic energy and potential enstrophy in continuous time. The repository includes dense and matrix-free pressure construction, sparse runtime storage, persistence, resolution rebuilding, and linear scattering examples.

This is a coherent solution to a different requirement. Its construction is expensive, its invariant-preserving closure deliberately changes the raw terrain response, and its modal matrix is unnecessary for prescribed bottom generation.

Useful implementation precedents include:

- `WVForcingType("Spectral")` accumulation without disturbing other forcings;
- componentwise interaction-phase handling;
- resolution-specific derived state and restart persistence;
- conjugacy, invariant-tendency, and matrix-exponential tests.

The new generator should not carry forward the finite-terrain pressure solve, dense or sparse terrain matrix, equal-\(\lambda\) restriction, or flat-invariant surrogate.

## Earlier topographic-forcing branches

The mapped strong-form checkpoint evaluates terrain metric terms using pressure reconstructed from the ordinary flat modal state and then applies the ordinary flat spatial projection. The work established that integrated pressure cancellation does not justify substituting the flat diagnostic pressure into finite-terrain strong-form metric accelerations: a pressure mismatch can be converted by the terrain metric into a non-gradient modal tendency.

The `terrain-weighted-weak-form` branch instead applies the finite-terrain energy inner product, under which pressure work vanishes. That formulation exposes the terrain-dependent Gram problem and the missing bottom displacement value, but it does not reduce to the inexpensive ordinary flat projection sought for the wave generator.

Both branches remain valuable derivational checkpoints. Neither should be incrementally modified into the new generator.

## Active mean-depth direction

The active formulation keeps the flat interior equations and expands the physical lower boundary about the mean depth. For prescribed depth-independent barotropic flow,

\[
g_b
=
\boldsymbol U_{\mathrm{bt}}\boldsymbol{\cdot}\nabla_Hh.
\]

A Green identity projects this boundary velocity directly onto each flat wave mode through its bottom pressure:

\[
\dot A_\alpha
=
\frac{1}{A E_\alpha}
\int_Ap_{\alpha,d}^*g_b\,dA,
\qquad
\alpha\in\{+,-\}.
\]

This is a first-order boundary formulation rather than an exact finite-height mapped formulation. It is deliberately restricted to the wave subspace, so it supplies no direct balanced tendency or linear interior QGPV. For a prescribed harmonic barotropic current, the spectral response can be precomputed and ordinary runtime calls require neither a pressure solve nor a modal matrix multiplication.

The repositories remain independent. Comparisons may load another checkout explicitly, but no terrain-forcing repository becomes a runtime dependency of another.
