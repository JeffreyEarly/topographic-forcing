function audit = auditGlobalFirstOrderTerrainDressing(self,options)
% Audit globally coupled first-order terrain-dressed spectral blocks.
%
% This Milestone-9.4 diagnostic uses the independently verified terrain
% derivatives of the primitive physical-energy forms to dress every
% validated internal-wave block and the complete flat zero-frequency
% block. The first-order calculation constructs coordinates only. Every
% finite-amplitude reduced eigensystem uses the unchanged exact primitive
% energy and exchange matrices and the exact stationary balanced space.
%
% ```matlab
% audit = problem.auditGlobalFirstOrderTerrainDressing( ...
%     trustedModeBounds=[1 0], ...
%     supportModeBounds=[1 2;1 3;1 4], ...
%     stationaryPolynomialDegree=4, ...
%     primitivePolynomialDegrees=[4;6;8], ...
%     comparisonPolynomialDegree=12, ...
%     paddingFactors=[2;3], ...
%     terrainScales=[1/16;1/8;1/4;1/2;1], ...
%     tangentStep=1e-3);
% ```
%
% - Topic: Audit global first-order terrain dressing
% - Declaration: audit = auditGlobalFirstOrderTerrainDressing(options)
% - Parameter trustedModeBounds: nonnegative integer `[kMax lMax]`
% - Parameter supportModeBounds: nested signed Fourier bounds
% - Parameter stationaryPolynomialDegree: trusted stationary-space degree
% - Parameter primitivePolynomialDegrees: primitive classification degrees
% - Parameter comparisonPolynomialDegree: independent primitive comparison degree
% - Parameter paddingFactors: horizontal oversampling factors
% - Parameter terrainScales: positive increasing terrain multipliers ending at one
% - Parameter tangentStep: positive nominal centered-difference terrain step
% - Parameter quadratureOrder: optional common primitive quadrature order
% - Returns audit: block corrections, finite-amplitude audits, and outcome classification
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.trustedModeBounds (1,2) double {mustBeInteger,mustBeNonnegative} = [1 0]
    options.supportModeBounds (:,2) double {mustBeInteger,mustBeNonnegative} = [1 2;1 3;1 4]
    options.stationaryPolynomialDegree (1,1) double {mustBeInteger,mustBePositive} = 4
    options.primitivePolynomialDegrees (:,1) double {mustBeInteger,mustBePositive} = [4;6;8]
    options.comparisonPolynomialDegree (1,1) double {mustBeInteger,mustBePositive} = 12
    options.paddingFactors (:,1) double {mustBeInteger,mustBePositive} = [2;3]
    options.terrainScales (:,1) double {mustBePositive} = [1/16;1/8;1/4;1/2;1]
    options.tangentStep (1,1) double {mustBePositive} = 1e-3
    options.quadratureOrder double {mustBeInteger} = []
end

degrees = options.primitivePolynomialDegrees;
supports = options.supportModeBounds;
if numel(degrees) < 3 || any(diff(degrees) <= 0) || any(degrees < options.stationaryPolynomialDegree)
    error("WVTerrainEnergyGalerkin:InvalidGlobalDressingDegrees", ...
        "primitivePolynomialDegrees must contain at least three strictly increasing values no smaller than stationaryPolynomialDegree.")
end
if size(supports,1) ~= numel(degrees) || any(diff(supports,1,1) < 0,"all") || any(~any(diff(supports,1,1) > 0,2))
    error("WVTerrainEnergyGalerkin:InvalidGlobalDressingSupports", ...
        "supportModeBounds must contain one strictly expanding row per primitivePolynomialDegrees entry.")
end
if any(options.trustedModeBounds > supports(1,:))
    error("WVTerrainEnergyGalerkin:GlobalDressingTrustedBandOutsideSupport", ...
        "trustedModeBounds must lie inside every supportModeBounds row.")
end
if options.comparisonPolynomialDegree <= max(degrees)
    error("WVTerrainEnergyGalerkin:InvalidGlobalDressingComparisonDegree", ...
        "comparisonPolynomialDegree must exceed every primitivePolynomialDegrees entry.")
end
if numel(options.paddingFactors) < 2 || any(options.paddingFactors < 2) || numel(unique(options.paddingFactors)) ~= numel(options.paddingFactors)
    error("WVTerrainEnergyGalerkin:InvalidGlobalDressingPadding", ...
        "paddingFactors must contain at least two unique values of two or greater.")
end
scales = options.terrainScales;
if numel(scales) < 4 || any(~isfinite(scales)) || any(diff(scales) <= 0) || scales(end) ~= 1
    error("WVTerrainEnergyGalerkin:InvalidGlobalDressingTerrainScales", ...
        "terrainScales must contain at least four finite increasing values ending at one.")
end
if any(1-scales*max(self.topographicHeight,[],"all")/self.originatingTransform.Lz <= 0)
    error("WVTerrainEnergyGalerkin:NonpositiveGlobalDressingGamma", ...
        "Every requested terrain scale must satisfy gamma>0.")
end
if ~isfinite(options.tangentStep) || options.tangentStep >= 0.1
    error("WVTerrainEnergyGalerkin:InvalidGlobalDressingTangentStep", ...
        "tangentStep must be finite and less than 0.1.")
end
if isempty(options.quadratureOrder)
    quadratureOrder = [];
elseif ~isscalar(options.quadratureOrder) || options.quadratureOrder < 2*options.comparisonPolynomialDegree+3
    error("WVTerrainEnergyGalerkin:InsufficientGlobalDressingQuadrature", ...
        "quadratureOrder must be empty or at least 2*comparisonPolynomialDegree+3.")
else
    quadratureOrder = options.quadratureOrder;
end
if max(abs(self.topographicHeight),[],"all") == 0
    error("WVTerrainEnergyGalerkin:GlobalDressingRequiresTerrainDirection", ...
        "topographicHeight must provide a nonzero terrain direction.")
end

audit = buildGlobalFirstOrderTerrainDressingAudit(self, ...
    options.trustedModeBounds,supports, ...
    options.stationaryPolynomialDegree,degrees, ...
    options.comparisonPolynomialDegree,options.paddingFactors, ...
    scales,options.tangentStep,quadratureOrder);
end
