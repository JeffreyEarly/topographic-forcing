function results = runTests(options)
%RUNTESTS Run the topographic-forcing automated test suite.
%
% WaveVortexModel, main InternalModes, and the isolated InternalModesEVP
% checkout are resolved from explicit options, environment variables, or
% sibling repositories. The MATLAB path is restored when the run finishes.
arguments
    options.waveVortexModelRoot (1,1) string = ""
    options.internalModesRoot (1,1) string = ""
    options.internalModesEVPRoot (1,1) string = ""
end

repositoryRoot = string(fileparts(mfilename("fullpath")));
waveVortexModelRoot = options.waveVortexModelRoot;
if strlength(waveVortexModelRoot) == 0
    waveVortexModelRoot = string(getenv("WAVE_VORTEX_MODEL_ROOT"));
end
if strlength(waveVortexModelRoot) == 0
    waveVortexModelRoot = fullfile(fileparts(repositoryRoot),"wave-vortex-model");
end
internalModesRoot = resolveRoot(options.internalModesRoot,"INTERNAL_MODES_ROOT", ...
    fullfile(fileparts(repositoryRoot),"internal-modes"));
internalModesEVPRoot = resolveRoot(options.internalModesEVPRoot,"INTERNAL_MODES_EVP_ROOT", ...
    fullfile(fileparts(repositoryRoot),"internal-modes-evp"));

if ~isfolder(waveVortexModelRoot) || ~isfile(fullfile(waveVortexModelRoot,"WVForcing.m"))
    error("BottomWaveGenerationForcing:WaveVortexModelNotFound", ...
        "WaveVortexModel was not found at '%s'. Pass waveVortexModelRoot or set WAVE_VORTEX_MODEL_ROOT.", waveVortexModelRoot)
end
if ~isfile(fullfile(internalModesRoot,"@InternalModesSpectral","InternalModesSpectral.m"))
    error("WVTerrainEnergyGalerkin:InternalModesNotFound", ...
        "Main InternalModes was not found at '%s'. Pass internalModesRoot or set INTERNAL_MODES_ROOT.",internalModesRoot)
end
if ~isfile(fullfile(internalModesEVPRoot,"@IMSolverSpectral","IMSolverSpectral.m")) || ...
        ~isfile(fullfile(internalModesEVPRoot,"@IMSurfaceGeostrophicModes","IMSurfaceGeostrophicModes.m"))
    error("WVTerrainEnergyGalerkin:InternalModesEVPNotFound", ...
        "InternalModesEVP was not found at '%s'. Pass internalModesEVPRoot or set INTERNAL_MODES_EVP_ROOT.",internalModesEVPRoot)
end

originalPath = path;
pathCleanup = onCleanup(@()path(originalPath));
addpath(repositoryRoot)
addpath(fullfile(repositoryRoot,"Examples"))
addpath(genpath(waveVortexModelRoot))
addpath(internalModesRoot)
addpath(internalModesEVPRoot,"-end")
if ~startsWith(string(which("InternalModesSpectral")),internalModesRoot)
    error("WVTerrainEnergyGalerkin:InternalModesPathConflict", ...
        "InternalModesSpectral must resolve from the main InternalModes checkout.")
end
if ~startsWith(string(which("IMSolverSpectral")),internalModesEVPRoot)
    error("WVTerrainEnergyGalerkin:InternalModesEVPPathConflict", ...
        "IMSolverSpectral must resolve from the isolated InternalModesEVP checkout.")
end

suite = matlab.unittest.TestSuite.fromFolder(fullfile(repositoryRoot,"UnitTests"));
runner = matlab.unittest.TestRunner.withTextOutput;
results = runner.run(suite);
if any([results.Failed])
    error("BottomWaveGenerationForcing:TestsFailed", "%d of %d tests failed.", nnz([results.Failed]), numel(results))
end

function root = resolveRoot(explicitRoot,environmentVariable,siblingRoot)
root = explicitRoot;
if strlength(root) == 0
    root = string(getenv(environmentVariable));
end
if strlength(root) == 0
    root = siblingRoot;
end
end
end
