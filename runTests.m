function results = runTests(options)
%RUNTESTS Run the topographic-forcing automated test suite.
%
% WaveVortexModel and WaveVortexModelDiagnostics are resolved, in order,
% from explicit options, their corresponding environment variables, or
% sibling authoring repositories. The MATLAB path is restored when the
% test run finishes.
arguments
    options.waveVortexModelRoot (1,1) string = ""
    options.waveVortexModelDiagnosticsRoot (1,1) string = ""
end

repositoryRoot = string(fileparts(mfilename("fullpath")));
waveVortexModelRoot = options.waveVortexModelRoot;
if strlength(waveVortexModelRoot) == 0
    waveVortexModelRoot = string(getenv("WAVE_VORTEX_MODEL_ROOT"));
end
if strlength(waveVortexModelRoot) == 0
    waveVortexModelRoot = fullfile(fileparts(repositoryRoot),"wave-vortex-model");
end

if ~isfolder(waveVortexModelRoot) || ~isfile(fullfile(waveVortexModelRoot,"WVForcing.m"))
    error("BottomWaveGenerationForcing:WaveVortexModelNotFound", ...
        "WaveVortexModel was not found at '%s'. Pass waveVortexModelRoot or set WAVE_VORTEX_MODEL_ROOT.", waveVortexModelRoot)
end

waveVortexModelDiagnosticsRoot = options.waveVortexModelDiagnosticsRoot;
if strlength(waveVortexModelDiagnosticsRoot) == 0
    waveVortexModelDiagnosticsRoot = string(getenv("WAVE_VORTEX_MODEL_DIAGNOSTICS_ROOT"));
end
if strlength(waveVortexModelDiagnosticsRoot) == 0
    waveVortexModelDiagnosticsRoot = fullfile(fileparts(repositoryRoot),"wave-vortex-model-diagnostics");
end

if ~isfolder(waveVortexModelDiagnosticsRoot) || ~isfile(fullfile(waveVortexModelDiagnosticsRoot,"@WVDiagnostics","WVDiagnostics.m"))
    error("BottomWaveGenerationForcing:WaveVortexModelDiagnosticsNotFound", ...
        "WaveVortexModelDiagnostics was not found at '%s'. Pass waveVortexModelDiagnosticsRoot or set WAVE_VORTEX_MODEL_DIAGNOSTICS_ROOT.", waveVortexModelDiagnosticsRoot)
end

originalPath = path;
pathCleanup = onCleanup(@()path(originalPath));
addpath(repositoryRoot)
addpath(fullfile(repositoryRoot,"Examples"))
addpath(genpath(waveVortexModelRoot))
addpath(genpath(waveVortexModelDiagnosticsRoot))

suite = matlab.unittest.TestSuite.fromFolder(fullfile(repositoryRoot,"UnitTests"));
runner = matlab.unittest.TestRunner.withTextOutput;
results = runner.run(suite);
if any([results.Failed])
    error("BottomWaveGenerationForcing:TestsFailed", "%d of %d tests failed.", nnz([results.Failed]), numel(results))
end
end
