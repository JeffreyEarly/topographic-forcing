function results = runTests(options)
%RUNTESTS Run the topographic-forcing automated test suite.
%
% WaveVortexModel is resolved, in order, from the explicit
% `waveVortexModelRoot` option, the `WAVE_VORTEX_MODEL_ROOT` environment
% variable, or the sibling `wave-vortex-model` repository. The MATLAB path
% is restored when the test run finishes.
arguments
    options.waveVortexModelRoot (1,1) string = ""
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
    error("ExactTopographicForcing:WaveVortexModelNotFound", ...
        "WaveVortexModel was not found at '%s'. Pass waveVortexModelRoot or set WAVE_VORTEX_MODEL_ROOT.", waveVortexModelRoot)
end

originalPath = path;
pathCleanup = onCleanup(@()path(originalPath));
addpath(repositoryRoot)
addpath(fullfile(repositoryRoot,"Examples"))
addpath(genpath(waveVortexModelRoot))

suite = matlab.unittest.TestSuite.fromFolder(fullfile(repositoryRoot,"UnitTests"));
runner = matlab.unittest.TestRunner.withTextOutput;
results = runner.run(suite);
if any([results.Failed])
    error("ExactTopographicForcing:TestsFailed", "%d of %d tests failed.", nnz([results.Failed]), numel(results))
end
end
