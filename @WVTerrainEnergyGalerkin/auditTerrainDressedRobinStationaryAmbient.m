function audit = auditTerrainDressedRobinStationaryAmbient(self,options)
% Audit globally dressed Robin stationary coordinates.
%
% This Milestone-10.2.4 diagnostic trains a signed Robin length on the
% stationary sector alone, freezes that coordinate choice, and validates
% the globally dressed geostrophic inclusion $$G_0+\delta G_1$$ against the
% exact finite-terrain stationary space. Fixed-wavenumber waves, the
% explicit zero-APV bottom inversion, and the compatible mean sector remain
% separate coordinates. The primitive polynomial system is an independent
% validation oracle.
%
% ```matlab
% audit = problem.auditTerrainDressedRobinStationaryAmbient( ...
%     trustedModeBounds=[1 0], ...
%     supportModeBounds=[1 2;1 3;1 4], ...
%     targetWaveModeIndices=[1;2], ...
%     waveGuardModeCounts=[2;4;6;8;12], ...
%     geostrophicModeCounts=[2;4;6;8;12], ...
%     robinLengthRatios=[-1/8;-1/4;-1/2;-1;Inf], ...
%     trainingGeostrophicModeCount=8, ...
%     primitiveReferenceDegrees=[16;18;20], ...
%     paddingFactors=[2;3], ...
%     terrainScales=[0;1/8;1/4;1/2;3/4;1], ...
%     internalModesEVPOrders=[128;256], ...
%     tangentStep=1e-3, ...
%     cacheDirectory="output/milestone-10.2.4-cache");
% ```
%
% - Topic: Audit terrain-dressed Robin stationary coordinates
% - Declaration: audit = auditTerrainDressedRobinStationaryAmbient(options)
% - Parameter trustedModeBounds: nonnegative integer `[kMax lMax]`
% - Parameter supportModeBounds: increasing retained support bounds
% - Parameter targetWaveModeIndices: positive production wave labels
% - Parameter waveGuardModeCounts: increasing fixed-wavenumber wave counts
% - Parameter geostrophicModeCounts: increasing APV-bearing mode counts
% - Parameter robinLengthRatios: signed Robin lengths divided by depth
% - Parameter trainingGeostrophicModeCount: count used only to select the Robin length
% - Parameter primitiveReferenceDegrees: increasing polynomial oracle degrees whose middle entry is used for training
% - Parameter paddingFactors: distinct horizontal oversampling factors including two and three
% - Parameter terrainScales: increasing terrain multipliers from zero to one
% - Parameter internalModesEVPOrders: increasing isolated-provider orders
% - Parameter tangentStep: positive centered-difference step for the geostrophic inclusion
% - Parameter shouldRunTwoDimensionalControl: run the two-dimensional control after a zonal pass
% - Parameter cacheDirectory: optional untracked directory for resumable artifacts
% - Parameter quadratureOrder: optional common vertical quadrature order
% - Returns audit: training ablations, validation sweeps, and outcome classification
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.trustedModeBounds (1,2) double {mustBeInteger,mustBeNonnegative} = [1 0]
    options.supportModeBounds (:,2) double {mustBeInteger,mustBeNonnegative} = [1 2;1 3;1 4]
    options.targetWaveModeIndices (:,1) double {mustBeInteger,mustBePositive} = [1;2]
    options.waveGuardModeCounts (:,1) double {mustBeInteger,mustBePositive} = [2;4;6;8;12]
    options.geostrophicModeCounts (:,1) double {mustBeInteger,mustBePositive} = [2;4;6;8;12]
    options.robinLengthRatios (:,1) double {mustBeReal} = [-1/8;-1/4;-1/2;-1;Inf]
    options.trainingGeostrophicModeCount (1,1) double {mustBeInteger,mustBePositive} = 8
    options.primitiveReferenceDegrees (:,1) double {mustBeInteger,mustBePositive} = [16;18;20]
    options.paddingFactors (:,1) double {mustBeInteger,mustBePositive} = [2;3]
    options.terrainScales (:,1) double {mustBeNonnegative} = [0;1/8;1/4;1/2;3/4;1]
    options.internalModesEVPOrders (:,1) double {mustBeInteger,mustBePositive} = [128;256]
    options.tangentStep (1,1) double {mustBePositive} = 1e-3
    options.shouldRunTwoDimensionalControl (1,1) logical = true
    options.cacheDirectory (1,1) string = ""
    options.quadratureOrder double {mustBeInteger} = []
end

if size(options.supportModeBounds,1) < 3 || any(any(diff(options.supportModeBounds,1,1) < 0)) || any(all(diff(options.supportModeBounds,1,1) == 0,2))
    error("WVTerrainEnergyGalerkin:InvalidDressedRobinSupport", ...
        "supportModeBounds must contain at least three distinct componentwise-nondecreasing rows.")
end
if any(options.supportModeBounds(1,:) < options.trustedModeBounds)
    error("WVTerrainEnergyGalerkin:InsufficientDressedRobinSupport", ...
        "Every support bound must contain trustedModeBounds.")
end
if isempty(options.targetWaveModeIndices) || ~issorted(options.targetWaveModeIndices) || numel(unique(options.targetWaveModeIndices)) ~= numel(options.targetWaveModeIndices) || any(~ismember(options.targetWaveModeIndices,self.verticalModeIndices))
    error("WVTerrainEnergyGalerkin:InvalidDressedRobinTargetWaves", ...
        "targetWaveModeIndices must be sorted unique labels retained by the originating transform.")
end
if numel(options.waveGuardModeCounts) < 3 || any(diff(options.waveGuardModeCounts) <= 0) || options.waveGuardModeCounts(1) < max(options.targetWaveModeIndices)
    error("WVTerrainEnergyGalerkin:InvalidDressedRobinWaveCounts", ...
        "waveGuardModeCounts must contain at least three increasing counts that retain every target wave label.")
end
if numel(options.geostrophicModeCounts) < 3 || any(diff(options.geostrophicModeCounts) <= 0) || ~ismember(options.trainingGeostrophicModeCount,options.geostrophicModeCounts)
    error("WVTerrainEnergyGalerkin:InvalidDressedRobinGeostrophicCounts", ...
        "geostrophicModeCounts must contain at least three increasing counts including trainingGeostrophicModeCount.")
end
if numel(options.robinLengthRatios) < 3 || any(isnan(options.robinLengthRatios)) || any(options.robinLengthRatios == 0) || any(isinf(options.robinLengthRatios) & options.robinLengthRatios < 0) || numel(unique(options.robinLengthRatios)) ~= numel(options.robinLengthRatios) || ~any(isinf(options.robinLengthRatios)) || nnz(isfinite(options.robinLengthRatios) & options.robinLengthRatios < 0) < 2
    error("WVTerrainEnergyGalerkin:InvalidDressedRobinLengths", ...
        "robinLengthRatios must contain unique nonzero finite values or Inf, including Inf and at least two negative values.")
end
if numel(options.primitiveReferenceDegrees) < 3 || any(diff(options.primitiveReferenceDegrees) <= 0)
    error("WVTerrainEnergyGalerkin:InvalidDressedRobinReferenceDegrees", ...
        "primitiveReferenceDegrees must contain at least three increasing degrees.")
end
if numel(options.paddingFactors) < 2 || any(options.paddingFactors < 2) || numel(unique(options.paddingFactors)) ~= numel(options.paddingFactors) || ~all(ismember([2;3],options.paddingFactors))
    error("WVTerrainEnergyGalerkin:InvalidDressedRobinPadding", ...
        "paddingFactors must contain unique values of two or greater, including two and three.")
end
if numel(options.terrainScales) < 5 || options.terrainScales(1) ~= 0 || options.terrainScales(end) ~= 1 || any(diff(options.terrainScales) <= 0)
    error("WVTerrainEnergyGalerkin:InvalidDressedRobinTerrainScales", ...
        "terrainScales must contain at least five increasing values beginning at zero and ending at one.")
end
if numel(options.internalModesEVPOrders) < 2 || any(options.internalModesEVPOrders < 16) || any(diff(options.internalModesEVPOrders) <= 0) || any(ismember(options.internalModesEVPOrders,[96;192]))
    error("WVTerrainEnergyGalerkin:InvalidDressedRobinEVPOrders", ...
        "internalModesEVPOrders must contain at least two increasing qualified values of sixteen or greater; orders 96 and 192 are excluded.")
end
if options.tangentStep >= 0.1
    error("WVTerrainEnergyGalerkin:InvalidDressedRobinTangentStep", ...
        "tangentStep must be positive and less than 0.1.")
end
if isempty(options.quadratureOrder)
    quadratureOrder = max(2*max(options.primitiveReferenceDegrees)+7,2*max([options.waveGuardModeCounts;options.geostrophicModeCounts])+15);
elseif ~isscalar(options.quadratureOrder) || options.quadratureOrder < 2*max(options.primitiveReferenceDegrees)+3
    error("WVTerrainEnergyGalerkin:InsufficientDressedRobinQuadrature", ...
        "quadratureOrder must be empty or at least 2*max(primitiveReferenceDegrees)+3.")
else
    quadratureOrder = options.quadratureOrder;
end
if max(abs(self.topographicHeight),[],"all") == 0
    error("WVTerrainEnergyGalerkin:DressedRobinRequiresTerrain", ...
        "topographicHeight must provide a nonzero terrain direction.")
end
if any(1-options.terrainScales*max(self.topographicHeight,[],"all")/self.originatingTransform.Lz <= 0)
    error("WVTerrainEnergyGalerkin:NonpositiveDressedRobinGamma", ...
        "Every requested terrain scale must satisfy gamma>0.")
end
cacheDirectory = validatedCacheDirectory(options.cacheDirectory);

experiment = struct("supportModeBounds",options.supportModeBounds, ...
    "robinLengthRatios",options.robinLengthRatios, ...
    "trainingGeostrophicModeCount",options.trainingGeostrophicModeCount, ...
    "tangentStep",options.tangentStep);
orders = (1:size(options.supportModeBounds,1)).';
audit = buildWaveVortexModalAmbientAudit(self,options.trustedModeBounds, ...
    options.targetWaveModeIndices,options.waveGuardModeCounts, ...
    options.geostrophicModeCounts,orders,options.primitiveReferenceDegrees, ...
    options.paddingFactors,options.terrainScales, ...
    options.internalModesEVPOrders,options.shouldRunTwoDimensionalControl, ...
    cacheDirectory,quadratureOrder,experiment);
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
    error("WVTerrainEnergyGalerkin:InvalidDressedRobinCacheLocation", ...
        "cacheDirectory must lie inside the repository's untracked output directory.")
end
if isfile(directory)
    error("WVTerrainEnergyGalerkin:InvalidDressedRobinCacheDirectory", ...
        "cacheDirectory must be empty or name a directory rather than an existing file.")
end
end
