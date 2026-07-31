function audit = auditGeometricCascadeIsolation(self,options)
% Audit geometric-scattering tails and internal-wave spectral isolation.
%
% This Milestone-10.2.2 diagnostic continues the complete production-wave
% projector from flat terrain, measures its horizontal and vertical tails
% in one common primitive comparison space, and tests whether the declared
% internal-wave block remains isolated under independent refinement.
%
% ```matlab
% audit = problem.auditGeometricCascadeIsolation( ...
%     trustedModeBounds=[1 0], ...
%     waveModeIndices=[1;2], ...
%     scatteringOrders=[1;2;3;4;5], ...
%     primitivePolynomialDegrees=[8;10;12;14], ...
%     comparisonPolynomialDegree=16, ...
%     paddingFactors=[2;3], ...
%     terrainScales=[0;1/16;1/8;1/4;1/2;3/4;1], ...
%     cacheDirectory="output/milestone-10.2.2-cache");
% ```
%
% - Topic: Audit geometric-cascade isolation
% - Declaration: audit = auditGeometricCascadeIsolation(options)
% - Parameter trustedModeBounds: nonnegative integer `[kMax lMax]`
% - Parameter waveModeIndices: retained positive native wave labels
% - Parameter apvModeIndices: retained balanced labels including zero
% - Parameter stationaryPolynomialDegree: stationary-space polynomial degree
% - Parameter scatteringOrders: consecutive positive terrain-scattering orders
% - Parameter primitivePolynomialDegrees: increasing primitive vertical degrees
% - Parameter comparisonPolynomialDegree: common primitive comparison degree
% - Parameter paddingFactors: horizontal oversampling factors
% - Parameter terrainScales: increasing terrain multipliers from zero to one
% - Parameter shouldRunTwoDimensionalControl: run the two-dimensional control after a zonal pass
% - Parameter cacheDirectory: optional untracked directory for resumable audit artifacts
% - Parameter quadratureOrder: optional common vertical quadrature order
% - Returns audit: cascade tails, isolation diagnostics, and outcome classification
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.trustedModeBounds (1,2) double {mustBeInteger,mustBeNonnegative} = [1 0]
    options.waveModeIndices (:,1) double {mustBeInteger} = self.verticalModeIndices(self.verticalModeIndices > 0 & self.verticalModeIndices <= 2)
    options.apvModeIndices (:,1) double {mustBeInteger} = self.verticalModeIndices
    options.stationaryPolynomialDegree (1,1) double {mustBeInteger,mustBePositive} = 4
    options.scatteringOrders (:,1) double {mustBeInteger,mustBePositive} = [1;2;3;4;5]
    options.primitivePolynomialDegrees (:,1) double {mustBeInteger,mustBePositive} = [8;10;12;14]
    options.comparisonPolynomialDegree (1,1) double {mustBeInteger,mustBePositive} = 16
    options.paddingFactors (:,1) double {mustBeInteger,mustBePositive} = [2;3]
    options.terrainScales (:,1) double {mustBeNonnegative} = [0;1/16;1/8;1/4;1/2;3/4;1]
    options.shouldRunTwoDimensionalControl (1,1) logical = true
    options.cacheDirectory (1,1) string = ""
    options.quadratureOrder double {mustBeInteger} = []
end

orders = options.scatteringOrders;
if numel(orders) < 3 || ~isequal(orders,(1:max(orders)).')
    error("WVTerrainEnergyGalerkin:InvalidGeometricCascadeOrders", ...
        "scatteringOrders must contain the consecutive positive orders one through its maximum value.")
end
degrees = options.primitivePolynomialDegrees;
if numel(degrees) < 3 || any(diff(degrees) <= 0) || any(degrees < options.stationaryPolynomialDegree)
    error("WVTerrainEnergyGalerkin:InvalidGeometricCascadeDegrees", ...
        "primitivePolynomialDegrees must contain at least three strictly increasing values no smaller than stationaryPolynomialDegree.")
end
if options.comparisonPolynomialDegree <= max(degrees)
    error("WVTerrainEnergyGalerkin:InvalidGeometricCascadeComparisonDegree", ...
        "comparisonPolynomialDegree must exceed every primitivePolynomialDegrees entry.")
end
if isempty(options.waveModeIndices) || any(options.waveModeIndices <= 0) || ~issorted(options.waveModeIndices) || numel(unique(options.waveModeIndices)) ~= numel(options.waveModeIndices) || any(~ismember(options.waveModeIndices,self.verticalModeIndices))
    error("WVTerrainEnergyGalerkin:InvalidGeometricCascadeWaveIndices", ...
        "waveModeIndices must be sorted unique positive labels retained by the originating transform.")
end
if isempty(options.apvModeIndices) || options.apvModeIndices(1) ~= 0 || ~issorted(options.apvModeIndices) || numel(unique(options.apvModeIndices)) ~= numel(options.apvModeIndices) || any(~ismember(options.apvModeIndices,self.verticalModeIndices))
    error("WVTerrainEnergyGalerkin:InvalidGeometricCascadeAPVIndices", ...
        "apvModeIndices must be sorted unique retained labels and include zero.")
end
if numel(options.paddingFactors) < 2 || any(options.paddingFactors < 2) || numel(unique(options.paddingFactors)) ~= numel(options.paddingFactors)
    error("WVTerrainEnergyGalerkin:InvalidGeometricCascadePadding", ...
        "paddingFactors must contain at least two unique values of two or greater.")
end
scales = options.terrainScales;
if numel(scales) < 5 || any(~isfinite(scales)) || scales(1) ~= 0 || scales(end) ~= 1 || any(diff(scales) <= 0)
    error("WVTerrainEnergyGalerkin:InvalidGeometricCascadeScales", ...
        "terrainScales must contain at least five finite increasing values beginning at zero and ending at one.")
end
if any(1-scales*max(self.topographicHeight,[],"all")/self.originatingTransform.Lz <= 0)
    error("WVTerrainEnergyGalerkin:NonpositiveGeometricCascadeGamma", ...
        "Every requested terrain scale must satisfy gamma>0.")
end
if options.shouldRunTwoDimensionalControl && ~ismember(1,self.verticalModeIndices)
    error("WVTerrainEnergyGalerkin:MissingGeometricCascadeTwoDimensionalMode", ...
        "The two-dimensional control requires native wave mode one.")
end
if isempty(options.quadratureOrder)
    quadratureOrder = max(2*options.comparisonPolynomialDegree+7,18);
elseif ~isscalar(options.quadratureOrder) || options.quadratureOrder < 2*options.comparisonPolynomialDegree+3
    error("WVTerrainEnergyGalerkin:InsufficientGeometricCascadeQuadrature", ...
        "quadratureOrder must be empty or at least 2*comparisonPolynomialDegree+3.")
else
    quadratureOrder = options.quadratureOrder;
end
cacheDirectory = validatedCacheDirectory(options.cacheDirectory);
if max(abs(self.topographicHeight),[],"all") == 0
    error("WVTerrainEnergyGalerkin:GeometricCascadeRequiresTerrain", ...
        "topographicHeight must provide a nonzero terrain direction.")
end

audit = buildGeometricCascadeIsolationAudit(self, ...
    options.trustedModeBounds,options.waveModeIndices, ...
    options.apvModeIndices,options.stationaryPolynomialDegree,orders, ...
    degrees,options.comparisonPolynomialDegree,options.paddingFactors, ...
    scales,options.shouldRunTwoDimensionalControl, ...
    cacheDirectory,quadratureOrder);
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
    error("WVTerrainEnergyGalerkin:InvalidGeometricCascadeCacheLocation", ...
        "cacheDirectory must lie inside the repository's untracked output directory.")
end
if isfile(directory)
    error("WVTerrainEnergyGalerkin:InvalidGeometricCascadeCacheDirectory", ...
        "cacheDirectory must be empty or name a directory rather than an existing file.")
end
end
