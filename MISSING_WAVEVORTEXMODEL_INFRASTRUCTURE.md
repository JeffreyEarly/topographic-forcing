# Missing WaveVortexModel Infrastructure

This document records Fourier and modal-layout capabilities that would simplify the terrain-energy Galerkin implementation and may be broadly useful within WaveVortexModel. These are proposed future additions to the WaveVortexModel authoring repository; they are not requirements for the completed Batch A flat oracle.

## Priority definitions

- **P0 — scientific prerequisite:** implement before relying on the associated operation in finite-terrain forms.
- **P1 — production prerequisite:** the current Galerkin code has a correct local implementation, but WaveVortexModel support should exist before production matrix-free work.
- **P2 — general hardening:** useful for clarity, validation, or broader reuse, but not currently blocking the Galerkin roadmap.

## P0: Adjoint-consistent horizontal Fourier oversampling

### Proposed functions

```matlab
targetDFT = prolongDFTGridToGeometry(sourceGeometry,targetGeometry,sourceDFT);
sourceDFT = restrictDFTGridFromGeometryAdjoint(sourceGeometry,targetGeometry,targetDFT);
```

These functions would transfer DFT-grid variables between compatible doubly periodic geometries by matching integer horizontal mode numbers. Prolongation would zero-pad unresolved modes, while restriction would be the exact discrete adjoint under the documented Fourier normalization.

### Why it is needed

Finite-terrain products require oversampled horizontal quadrature. The same prolongation and restriction must be used for reconstruction and projection so the discrete terrain energy form remains Hermitian and the exchange form remains skew-Hermitian. Interpolation or unmatched truncation would compromise those identities.

`spectralVariableWithResolution` is not sufficient for this purpose: it operates on WV modal arrays, copies a prefix of the stored ordering, and does not expose an adjoint full-DFT transfer. `waveVortexTransformWithDoubleResolution` also changes the vertical resolution when only horizontal oversampling is required.

Milestone 4 now contains a validated project-local implementation for the dense oracle. It reconstructs retained signed modes on the oversampled periodic grid, uses the exact horizontal mean as the adjoint projection, and treats terrain interpolation separately from the prognostic truncation. The upstream API remains P0 because Milestone 7 will need the same operation as a reusable matrix-free prolongation/restriction pair rather than as an explicitly materialized dense reconstruction.

### Acceptance criteria

- Preserve every horizontal Fourier coefficient common to the two geometries.
- Handle rectangular domains, antialiasing enabled or disabled, and all retained signed modes.
- Satisfy the discrete adjoint identity to roundoff for random complex arrays.
- Preserve real-field Hermitian conjugacy.
- Make normalization independent of the source and target grid sizes.

## P1: Complex DFT-grid inverse transformation

### Proposed function

```matlab
u = transformToSpatialDomainFromComplexDFTGrid(geometry,uBar);
```

This function would apply the normalized inverse two-dimensional DFT without the `symmetric` assumption and return a complex spatial field. Its forward partner can remain `transformFromSpatialDomainToDFTGrid`, provided that method formally documents support for complex input.

### Why it is needed

`transformToSpatialDomainFromDFTGrid` intentionally reconstructs a real field and discards non-Hermitian imaginary content. Individual Galerkin trial vectors, eigenvectors, and residual vectors are complex and need not be paired with their conjugates during operator construction. The current Galerkin class therefore contains a small private complex inverse transform.

An explicit complex API would remove that fallback and prevent callers from accidentally using a real-field transform for complex modal calculations.

### Acceptance criteria

- Round-trip arbitrary complex DFT arrays to numerical precision.
- Retain the existing WaveVortexModel Fourier normalization.
- Support two- and three-dimensional arrays with arbitrary trailing vertical length.
- Agree with the real-field transform for Hermitian input.
- Work for MATLAB and FFTW transform backends.

## P1: Complex horizontal spectral derivatives

### Proposed functions

```matlab
uX = diffXComplex(geometry,u,n=1);
uY = diffYComplex(geometry,u,n=1);
```

These functions would differentiate complex spatial fields without invoking a symmetric inverse transform. They would use the same wavenumber arrays, Fourier normalization, and Nyquist convention as the real-valued `diffX` and `diffY` methods.

### Why it is needed

The existing derivative methods return real fields and therefore remove the imaginary component of an unpaired complex Fourier mode. Complex derivatives are required when evaluating a single Galerkin basis vector, an eigen-residual, or a matrix-free operator action.

The terrain-energy implementation currently supplies a private fallback built from `transformFromSpatialDomainToDFTGrid`, the public DFT wavenumbers, and a complex inverse transform.

### Acceptance criteria

- Differentiate analytic complex Fourier modes to roundoff.
- Support arbitrary derivative order consistently with `diffX` and `diffY`.
- Apply a documented Nyquist convention.
- Reduce to the existing real derivative for real input without changing normalization.

## P1: Complete retained full-signed Fourier layout

### Proposed function

```matlab
layout = fullComplexHorizontalLayout(geometry);
```

The returned layout would describe every retained signed horizontal Fourier coefficient after applying the geometry’s antialiasing and Nyquist policies. It should include:

- DFT linear index;
- zonal and meridional mode numbers;
- dimensional wavenumbers;
- whether the coefficient belongs to the ordinary nonredundant WV layout;
- its corresponding primary WV index;
- its conjugate partner within the full layout.

### Why it is needed

The terrain Galerkin system couples positive and negative horizontal wavenumbers and therefore needs an explicit full-signed layout. The required information can currently be assembled from `dftPrimaryIndices2D`, `dftConjugateIndices2D`, and `indicesOfFourierConjugates`, but every client must repeat the indexing logic.

Constructing `WVGeometryDoublyPeriodic` with `shouldExludeConjugates=false` provides the signed modes, but `wvConjugateIndex` is a fill map for the ordinary real transform rather than a complete conjugate involution in that configuration.

### Acceptance criteria

- Return a complete conjugate involution for every retained mode.
- Respect antialiasing, Nyquist exclusion, and `conjugateDimension`.
- Work for rectangular even grids.
- Reproduce the ordinary WV ordering and indices without mode-number searches.

## P2: Full-signed wave–vortex coefficient conversion

### Proposed functions

```matlab
state = transformWaveVortexCoefficientsToFullSignedLayout(wvt,Ap,Am,A0);
[Ap,Am,A0] = transformFullSignedLayoutToWaveVortexCoefficients(wvt,state);
```

These functions would convert the nonredundant WaveVortexModel coefficient matrices into a complete signed modal representation and back. The conversion must apply the wave-branch conjugacy

```text
Aplus(-K)  = conjugate(Aminus(K))
Aminus(-K) = conjugate(Aplus(K))
A0(-K)     = conjugate(A0(K))
```

along with the special zero-wavenumber inertial and mean-density relations.

### Why it is needed

Scalar Fourier conjugacy is insufficient for wave–vortex coefficients because conjugation exchanges the two wave branches. WaveVortexModel’s flow-component masks identify where the physical degrees of freedom live, but there is no general pack/unpack API for a full-signed modal state.

The current implementation is concise and terrain-specific, so this should move upstream only if another subsystem needs the same representation.

### Acceptance criteria

- Round-trip every valid `Ap`, `Am`, and `A0` state.
- Preserve inertial conjugacy and real zero-wavenumber balanced coefficients.
- Respect inactive coefficient masks and antialiasing.
- Reconstruct identical physical fields before and after conversion.

## P2: Explicit horizontal-grid and Nyquist contract

### Proposed change

Either validate even `Nx` and `Ny` in `WVGeometryDoublyPeriodic`, or generalize the conjugacy masks, degrees-of-freedom counts, and Nyquist utilities to odd grid sizes.

### Why it is needed

Several geometry utilities index `Nx/2+1` and `Ny/2+1`, which assumes even horizontal dimensions, while the constructor currently validates only positivity. A clear contract would turn an implicit implementation assumption into an actionable public error or add genuine odd-grid support.

The terrain-energy implementation currently inherits the standard even-grid WaveVortexModel configuration and verifies that no Nyquist mode enters its retained signed layout.

### Acceptance criteria

- Give a structured error for unsupported sizes, or pass complete odd-grid Fourier tests.
- Keep `maskForNyquistModes`, conjugate masks, and degrees-of-freedom counts mutually consistent.
- Document the behavior in the geometry constructor.

## P2: Fourier-geometry regression suite

### Proposed tests

Add focused WaveVortexModel tests for:

- WV-to-DFT-to-WV round trips;
- primary and conjugate index maps;
- complete signed-layout conjugacy;
- rectangular domains;
- antialiasing enabled and disabled;
- both conjugate dimensions;
- Nyquist exclusion;
- real and complex Fourier transforms;
- real and complex spectral derivatives.

### Why it is needed

The Fourier geometry is foundational infrastructure used by transforms, forcings, diagnostics, and the terrain Galerkin system. Existing tests cover ordinary real Fourier reconstruction but do not exercise the complete combination of layouts, conjugacy choices, complex trial fields, and Nyquist behavior needed by operator-based methods.

Each new infrastructure function above should land with its corresponding tests rather than relying only on downstream terrain tests.

## Recommended implementation order

1. Implement adjoint-consistent oversampling before assembling finite-terrain forms.
2. Add complex inverse transforms and derivatives before production matrix-free operator actions.
3. Add the full-signed layout descriptor and migrate the Galerkin indexing adapter to it.
4. Promote wave–vortex coefficient conversion only if a second use case appears.
5. Clarify the even-grid contract and expand the geometry regression suite alongside the preceding changes.
