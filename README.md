# Exact topographic forcing

`topographic-forcing` is a research add-on for WaveVortexModel that evaluates the exact flow-linear bottom-topography terms in mapped coordinates. It reuses the rigid-lid modes and the existing spatial-forcing projection of `WVTransformBoussinesq`.

The current proof of concept supports stationary periodic terrain, discretely constant stratification, and `shouldAntialias=false`. It is exact in prescribed terrain height and linear in flow amplitude.

```matlab
wvt = WVTransformBoussinesq([20e3 20e3 2e3],[8 4 5], ...
    N2=@(z) 2e-5*ones(size(z)),latitude=45,shouldAntialias=false);
x = repmat(reshape(wvt.x,[],1),1,wvt.Ny);
h = 50*cos(2*pi*x/wvt.Lx);

forcing = WVExactTopographicForcing(wvt,topographicHeight=h);
wvt.removeAllForcing();
wvt.addForcing(forcing);
model = WVModel(wvt);
```

The transform coordinate `wvt.z` is interpreted as the mapped coordinate \(\xi\), and `wvt.u`, `v`, `w`, `eta`, and `p` are interpreted as the hatted projection-ready fields.

## Tests

Run the suite with WaveVortexModel in the sibling repository:

```matlab
results = runTests();
```

Alternatively, pass `waveVortexModelRoot` explicitly or set `WAVE_VORTEX_MODEL_ROOT`. The test runner restores the original MATLAB path when it exits.

See [milestones.md](milestones.md) for the scientific validation roadmap. Uniform-depth and sloping-terrain solution oracles follow in Milestones 5–7.
