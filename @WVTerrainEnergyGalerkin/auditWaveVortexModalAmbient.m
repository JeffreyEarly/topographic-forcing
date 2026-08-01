function audit = auditWaveVortexModalAmbient(self,options)
% Audit a boundary-complete flat wave-vortex modal ambient space.
%
% This Milestone-10.2.3 diagnostic evaluates fixed-wavenumber internal-wave
% modes, ordinary APV-bearing geostrophic modes, the complete zero-APV
% bottom inversion, and the compatible mean sector directly on the common
% finite-terrain quadrature. The primitive polynomial system remains an
% independent validation oracle rather than the candidate evolution basis.
%
% ```matlab
% audit = problem.auditWaveVortexModalAmbient( ...
%     trustedModeBounds=[1 0], ...
%     targetWaveModeIndices=[1;2], ...
%     waveGuardModeCounts=[2;4;6;8;12], ...
%     geostrophicModeCounts=[2;4;6;8;12], ...
%     scatteringOrders=[1;2;3;4;5], ...
%     primitiveReferenceDegrees=[16;18;20], ...
%     paddingFactors=[2;3], ...
%     terrainScales=[0;1/8;1/4;1/2;3/4;1], ...
%     internalModesEVPOrders=[128;256], ...
%     cacheDirectory="output/milestone-10.2.3-cache");
% ```
%
% - Topic: Audit a flat wave-vortex modal ambient
% - Declaration: audit = auditWaveVortexModalAmbient(options)
% - Parameter trustedModeBounds: nonnegative integer `[kMax lMax]`
% - Parameter targetWaveModeIndices: positive production wave labels
% - Parameter waveGuardModeCounts: increasing fixed-wavenumber wave counts
% - Parameter geostrophicModeCounts: increasing APV-bearing mode counts
% - Parameter scatteringOrders: consecutive positive terrain-scattering orders
% - Parameter primitiveReferenceDegrees: increasing polynomial oracle degrees
% - Parameter paddingFactors: distinct horizontal oversampling factors
% - Parameter terrainScales: increasing terrain multipliers from zero to one
% - Parameter internalModesEVPOrders: increasing isolated-provider orders
% - Parameter shouldRunTwoDimensionalControl: run the two-dimensional control after a zonal pass
% - Parameter cacheDirectory: optional untracked directory for resumable artifacts
% - Parameter quadratureOrder: optional common vertical quadrature order
% - Returns audit: modal tails, physical projectors, and outcome classification
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.trustedModeBounds (1,2) double {mustBeInteger,mustBeNonnegative} = [1 0]
    options.targetWaveModeIndices (:,1) double {mustBeInteger,mustBePositive} = [1;2]
    options.waveGuardModeCounts (:,1) double {mustBeInteger,mustBePositive} = [2;4;6;8;12]
    options.geostrophicModeCounts (:,1) double {mustBeInteger,mustBePositive} = [2;4;6;8;12]
    options.scatteringOrders (:,1) double {mustBeInteger,mustBePositive} = [1;2;3;4;5]
    options.primitiveReferenceDegrees (:,1) double {mustBeInteger,mustBePositive} = [16;18;20]
    options.paddingFactors (:,1) double {mustBeInteger,mustBePositive} = [2;3]
    options.terrainScales (:,1) double {mustBeNonnegative} = [0;1/8;1/4;1/2;3/4;1]
    options.internalModesEVPOrders (:,1) double {mustBeInteger,mustBePositive} = [128;256]
    options.shouldRunTwoDimensionalControl (1,1) logical = true
    options.cacheDirectory (1,1) string = ""
    options.quadratureOrder double {mustBeInteger} = []
end

if isempty(options.targetWaveModeIndices) ...
        || ~issorted(options.targetWaveModeIndices) ...
        || numel(unique(options.targetWaveModeIndices)) ...
        ~=numel(options.targetWaveModeIndices) ...
        || any(~ismember(options.targetWaveModeIndices, ...
        self.verticalModeIndices))
    error("WVTerrainEnergyGalerkin:InvalidModalAmbientTargetWaves", ...
        "targetWaveModeIndices must be sorted unique positive labels retained by the originating transform.")
end
if numel(options.waveGuardModeCounts) < 3 || any(diff(options.waveGuardModeCounts) <= 0) || options.waveGuardModeCounts(1) < max(options.targetWaveModeIndices)
    error("WVTerrainEnergyGalerkin:InvalidModalAmbientWaveCounts", ...
        "waveGuardModeCounts must contain at least three increasing counts and retain every target wave label.")
end
if numel(options.geostrophicModeCounts) < 3 || any(diff(options.geostrophicModeCounts) <= 0)
    error("WVTerrainEnergyGalerkin:InvalidModalAmbientGeostrophicCounts", ...
        "geostrophicModeCounts must contain at least three strictly increasing counts.")
end
if numel(options.scatteringOrders) < 3 || ~isequal(options.scatteringOrders,(1:max(options.scatteringOrders)).')
    error("WVTerrainEnergyGalerkin:InvalidModalAmbientScatteringOrders", ...
        "scatteringOrders must contain the consecutive positive orders one through its maximum value.")
end
if numel(options.primitiveReferenceDegrees) < 3 || any(diff(options.primitiveReferenceDegrees) <= 0)
    error("WVTerrainEnergyGalerkin:InvalidModalAmbientReferenceDegrees", ...
        "primitiveReferenceDegrees must contain at least three strictly increasing degrees.")
end
if numel(options.paddingFactors) < 2 || any(options.paddingFactors < 2) || numel(unique(options.paddingFactors)) ~= numel(options.paddingFactors)
    error("WVTerrainEnergyGalerkin:InvalidModalAmbientPadding", ...
        "paddingFactors must contain at least two unique values of two or greater.")
end
if numel(options.terrainScales) < 5 || options.terrainScales(1) ~= 0 || options.terrainScales(end) ~= 1 || any(diff(options.terrainScales) <= 0)
    error("WVTerrainEnergyGalerkin:InvalidModalAmbientTerrainScales", ...
        "terrainScales must contain at least five increasing values beginning at zero and ending at one.")
end
if numel(options.internalModesEVPOrders) < 2 || any(options.internalModesEVPOrders < 16) || any(diff(options.internalModesEVPOrders) <= 0)
    error("WVTerrainEnergyGalerkin:InvalidModalAmbientEVPOrders", ...
        "internalModesEVPOrders must contain at least two increasing values of sixteen or greater.")
end
if any(ismember(options.internalModesEVPOrders,[96;192]))
    error("WVTerrainEnergyGalerkin:UnqualifiedModalAmbientEVPOrder", ...
        "InternalModesEVP orders 96 and 192 are excluded by the existing provider qualification.")
end
if isempty(options.quadratureOrder)
    quadratureOrder = max(2*max(options.primitiveReferenceDegrees)+7,2*max([options.waveGuardModeCounts;options.geostrophicModeCounts])+15);
elseif ~isscalar(options.quadratureOrder) || options.quadratureOrder < 2*max(options.primitiveReferenceDegrees)+3
    error("WVTerrainEnergyGalerkin:InsufficientModalAmbientQuadrature", ...
        "quadratureOrder must be empty or at least 2*max(primitiveReferenceDegrees)+3.")
else
    quadratureOrder = options.quadratureOrder;
end
if max(abs(self.topographicHeight),[],"all") == 0
    error("WVTerrainEnergyGalerkin:ModalAmbientRequiresTerrain", ...
        "topographicHeight must provide a nonzero terrain direction.")
end
if any(1-options.terrainScales*max(self.topographicHeight,[],"all")/self.originatingTransform.Lz <= 0)
    error("WVTerrainEnergyGalerkin:NonpositiveModalAmbientGamma", ...
        "Every requested terrain scale must satisfy gamma>0.")
end
cacheDirectory = validatedCacheDirectory(options.cacheDirectory);

audit = buildWaveVortexModalAmbientAudit(self,options.trustedModeBounds, ...
    options.targetWaveModeIndices,options.waveGuardModeCounts, ...
    options.geostrophicModeCounts,options.scatteringOrders, ...
    options.primitiveReferenceDegrees,options.paddingFactors, ...
    options.terrainScales,options.internalModesEVPOrders, ...
    options.shouldRunTwoDimensionalControl,cacheDirectory,quadratureOrder);
end

function directory = validatedCacheDirectory(directory)
if directory == ""
    return
end
repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
outputRoot = string(java.io.File(char(fullfile(repositoryRoot,"output"))).getCanonicalPath());
candidate = java.io.File(char(directory));
if ~candidate.isAbsolute()
    candidate = java.io.File(char(fullfile(repositoryRoot,directory)));
end
directory = string(candidate.getCanonicalPath());
if directory ~= outputRoot && ~startsWith(directory,outputRoot+filesep)
    error("WVTerrainEnergyGalerkin:InvalidModalAmbientCacheLocation", ...
        "cacheDirectory must lie inside the repository's untracked output directory.")
end
if isfile(directory)
    error("WVTerrainEnergyGalerkin:InvalidModalAmbientCacheDirectory", ...
        "cacheDirectory must be empty or name a directory rather than an existing file.")
end
end
