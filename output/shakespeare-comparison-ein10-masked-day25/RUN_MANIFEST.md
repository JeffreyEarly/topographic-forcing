# Masked Ein10 day-25 validation pair

## Purpose

Rerun the Shakespeare-comparison Ein10 eddy and no-initial-eddy cases with bottom-wave generation excluded from the exact nonzero support of `WVAdaptiveDamping`. Compare against the first 25 days of the completed unmasked Ein10 pair.

## Software

- MATLAB: 25.2.0.3150157 (R2025b) Update 4
- `topographic-forcing`: `7f43573f83673cf73ff3ef3f30c58674758433c2`, branch `mean-depth-wave-generator`
- `wave-vortex-model`: `52bdcdccaa65fa4c388adb39bfc23a5cbd7e7d2f`, branch `main`
- `wave-vortex-model-diagnostics`: `f1dcedeeba32b1474eb1cdcf9ea0140072e9aa69`, branch `main`

The WaveVortexModel and diagnostics authoring repositories were clean at launch. The only untracked content in `topographic-forcing` was `.DS_Store` and the existing `output/` tree.

## Configuration

- Cases: shallow Gaussian initial eddy and matched no-initial-eddy control
- Domain: 500 km by 500 km by 2 km
- Grid: `Nxy=256`, `Nz=45`, 29 wave modes
- Stratification: constant \(N^2=2\times10^{-5}\ \mathrm{s^{-2}}\)
- Latitude: 45 degrees north
- Integration: 25 days
- Output interval: 6 hours; 101 records including day 0
- Barotropic forcing: zonal M2 current, \(0.05/\sqrt{10}\ \mathrm{m\,s^{-1}}\)
- Ramp: one M2 period, 12.420602 hours
- Terrain: seed-2023 Goff realization, 100 m RMS, \(10^{-4}\ \mathrm{m^{-1}}\) corner wavenumber, 6 km minimum wavelength
- Dynamics: nonlinear advection, adaptive damping, transform-level antialiasing
- Topographic forcing: generation only; scattering disabled
- Generation mask: `shouldAvoidAdaptiveDamping=true`
- Manual bounds: `maximumForcedHorizontalWavenumber=Inf`, `maximumForcedVerticalMode=Inf`

## Output files

- `eddy-tide-topographic-eddy-generation-only-Lxy500km-Nxy256-Ein10-masked-day25-hrms100m-lmin6km-seed2023.nc`
- `eddy-tide-topographic-no-eddy-generation-only-Lxy500km-Nxy256-Ein10-masked-day25-hrms100m-lmin6km-seed2023.nc`
- Companion diagnostics use the standard `-diagnostics.nc` suffix.
- Figures and summaries are retained under `checkpoints/day-0025`.
- Case-specific logs and success/failure markers are retained in this directory.

## References

- Eddy: `../eddy-tide-topographic-eddy-generation-only-Lxy500km-Nxy256-Ein10-hrms100m-lmin6km-seed2023.nc`
- Control: `../eddy-tide-topographic-no-eddy-generation-only-Lxy500km-Nxy256-Ein10-hrms100m-lmin6km-seed2023.nc`

The reference files and their completed diagnostics are read-only inputs and must not be modified.

## Storage

- Free space at staging: 99,113,340 KiB, approximately 94.5 GiB.
- Projected pair plus diagnostics: approximately 3.4 GiB.
- Required reserve: 20 GiB.

Each runner checks the reserve before constructing its production output.

## Preflight and execution record

The focused forcing, production/persistence, and eddy-tide workflow test suites
all passed before production:

- `TestWVBottomWaveGenerationForcing`
- `TestWVBottomWaveGenerationProduction`
- `TestEddyTideTopographicForcingSimulation`

Both accepted simulations passed exact terrain and initial-state comparison
against their corresponding original Ein10 files. Their effective masks were
the wave-valid region intersected with `adaptiveDamping.damp == 0`, and both
wave-source tendencies were exactly zero on the damping support.

The first control launch was interrupted after an orchestration probe attempted
to open its live writer. That partial 532 MiB file, combined log, preflight, and
structured `NC_EHDFERR` record are preserved under
`attempts/control-attempt1/`. It is not part of the accepted pair. The control
was then rerun from day 0 with a clean isolated log.

The accepted eddy and control runs completed on 2026-07-30 with exactly 101
records at 0, 6, ..., 600 hours. Both final transforms contained finite
coefficients, restored the automatic mask, and reopened successfully after
their writers closed. The accepted case and analysis logs contain no
`NC_EHDFERR`, close failure, or open-NetCDF warning.

## Day-25 results

Full-cadence diagnostics and all requested figures were completed under
`checkpoints/day-0025`. Spatially integrated adaptive-damping work agreed with
the exact diagnostics to maximum relative errors of
`3.47e-16` (original) and `9.82e-16` (masked).

For the no-initial-eddy control:

- Final masked geostrophic energy: `0.01311724096 m^3 s^-2`
- Final masked wave energy: `0.4172671911 m^3 s^-2`
- `E_g/E_w`: `0.0314361`
- `E_g/(E_g+E_w)`: `0.0304780`
- Final geostrophic energy change from unmasked: `-61.46%`
- Final wave energy change from unmasked: `-6.91%`
- Maximum horizontal wave-speed change from unmasked: `-10.40%`
- Cumulative generation change from unmasked: `-20.79%`
- Cumulative adaptive-damping magnitude change from unmasked: `-99.16%`

The days-20--25 mean exact damping work changed from
`-4.00519e-8 m^3 s^-3` to `-1.21858e-9 m^3 s^-3`.

The complete raw series, final metrics, validation state, diagnostics paths,
and figure paths are stored in
`checkpoints/day-0025/eddy-tide-Lxy500km-Nxy256-Ein10-masked-day25-validation-summary.mat`.
