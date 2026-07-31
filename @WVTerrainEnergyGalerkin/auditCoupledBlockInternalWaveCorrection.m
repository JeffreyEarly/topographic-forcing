function audit = auditCoupledBlockInternalWaveCorrection(self,options)
% Audit coupled invariant-subspace correction of all production waves.
%
% This Milestone-10.2.1 diagnostic preserves the complete production state,
% exact stationary space, and finite-terrain physical-energy forms. It first
% solves the complementary Sylvester correction directly. The iterative
% realization is attempted only after that exact oracle passes.
%
% - Topic: Audit coupled-block internal-wave correction
% - Declaration: audit = auditCoupledBlockInternalWaveCorrection(options)
% - Parameter trustedModeBounds: nonnegative integer `[kMax lMax]`
% - Parameter supportModeBounds: nested signed-Fourier primitive bounds
% - Parameter waveModeIndices: retained positive native wave labels
% - Parameter apvModeIndices: retained balanced labels including zero
% - Parameter stationaryPolynomialDegree: stationary-space polynomial degree
% - Parameter primitivePolynomialDegrees: increasing primitive audit degrees
% - Parameter paddingFactors: horizontal oversampling factors
% - Parameter terrainScales: positive increasing terrain multipliers ending at one
% - Parameter tangentStep: positive centered terrain-derivative step
% - Parameter maximumOuterIterations: maximum coupled Ritz corrections
% - Parameter maximumKrylovRestarts: maximum iterative correction restarts
% - Parameter shouldCompareSlopeSeed: compare the frozen Milestone-9.3 seed
% - Parameter quadratureOrder: optional common primitive quadrature order
% - Returns audit: exact, iterative, economy, and physical diagnostics
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.trustedModeBounds (1,2) double {mustBeInteger,mustBeNonnegative} = [1 0]
    options.supportModeBounds (:,2) double {mustBeInteger,mustBeNonnegative} = [1 2;1 3;1 4]
    options.waveModeIndices (:,1) double {mustBeInteger} = self.verticalModeIndices(self.verticalModeIndices > 0)
    options.apvModeIndices (:,1) double {mustBeInteger} = self.verticalModeIndices
    options.stationaryPolynomialDegree (1,1) double {mustBeInteger,mustBePositive} = 4
    options.primitivePolynomialDegrees (:,1) double {mustBeInteger,mustBePositive} = [8;10;12]
    options.paddingFactors (:,1) double {mustBeInteger,mustBePositive} = [2;3]
    options.terrainScales (:,1) double {mustBePositive} = [1/16;1/8;1/4;1/2;1]
    options.tangentStep (1,1) double {mustBePositive} = 1e-3
    options.maximumOuterIterations (1,1) double {mustBeInteger,mustBePositive} = 6
    options.maximumKrylovRestarts (1,1) double {mustBeInteger,mustBePositive} = 4
    options.shouldCompareSlopeSeed (1,1) logical = true
    options.quadratureOrder double {mustBeInteger} = []
end

degrees = options.primitivePolynomialDegrees;
supports = options.supportModeBounds;
if numel(degrees) < 3 || any(diff(degrees) <= 0) || any(degrees < options.stationaryPolynomialDegree)
    error("WVTerrainEnergyGalerkin:InvalidCoupledBlockDegrees", ...
        "primitivePolynomialDegrees must contain at least three strictly increasing values no smaller than stationaryPolynomialDegree.")
end
if size(supports,1) ~= numel(degrees) || any(diff(supports,1,1) < 0,"all") || any(~any(diff(supports,1,1) > 0,2))
    error("WVTerrainEnergyGalerkin:InvalidCoupledBlockSupports", ...
        "supportModeBounds must contain one strictly expanding row per primitivePolynomialDegrees entry.")
end
if any(options.trustedModeBounds > supports(1,:))
    error("WVTerrainEnergyGalerkin:CoupledBlockTrustedBandOutsideSupport", ...
        "trustedModeBounds must lie inside every supportModeBounds row.")
end
if isempty(options.waveModeIndices) || any(options.waveModeIndices <= 0) || ~issorted(options.waveModeIndices) || numel(unique(options.waveModeIndices)) ~= numel(options.waveModeIndices) || any(~ismember(options.waveModeIndices,self.verticalModeIndices))
    error("WVTerrainEnergyGalerkin:InvalidCoupledBlockWaveIndices", ...
        "waveModeIndices must be sorted unique positive labels retained by the originating transform.")
end
if isempty(options.apvModeIndices) || options.apvModeIndices(1) ~= 0 || ~issorted(options.apvModeIndices) || numel(unique(options.apvModeIndices)) ~= numel(options.apvModeIndices) || any(~ismember(options.apvModeIndices,self.verticalModeIndices))
    error("WVTerrainEnergyGalerkin:InvalidCoupledBlockAPVIndices", ...
        "apvModeIndices must be sorted unique retained labels and include zero.")
end
if numel(options.paddingFactors) < 2 || any(options.paddingFactors < 2) || numel(unique(options.paddingFactors)) ~= numel(options.paddingFactors)
    error("WVTerrainEnergyGalerkin:InvalidCoupledBlockPadding", ...
        "paddingFactors must contain at least two unique values of two or greater.")
end
scales = options.terrainScales;
if numel(scales) < 4 || any(~isfinite(scales)) || any(diff(scales) <= 0) || scales(end) ~= 1
    error("WVTerrainEnergyGalerkin:InvalidCoupledBlockTerrainScales", ...
        "terrainScales must contain at least four finite increasing values ending at one.")
end
if any(1-scales*max(self.topographicHeight,[],"all")/self.originatingTransform.Lz <= 0)
    error("WVTerrainEnergyGalerkin:NonpositiveCoupledBlockGamma", ...
        "Every requested terrain scale must satisfy gamma>0.")
end
if ~isfinite(options.tangentStep) || options.tangentStep >= 0.1
    error("WVTerrainEnergyGalerkin:InvalidCoupledBlockTangentStep", ...
        "tangentStep must be finite and less than 0.1.")
end
if isempty(options.quadratureOrder)
    quadratureOrder = [];
elseif ~isscalar(options.quadratureOrder) || options.quadratureOrder < 2*max(degrees)+3
    error("WVTerrainEnergyGalerkin:InsufficientCoupledBlockQuadrature", ...
        "quadratureOrder must be empty or at least 2*max(primitivePolynomialDegrees)+3.")
else
    quadratureOrder = options.quadratureOrder;
end

contractSupports = options.trustedModeBounds+[0 0;0 1;0 2];
contract = self.auditProductionPhysicalStateContract( ...
    trustedModeBounds=options.trustedModeBounds, ...
    supportModeBounds=contractSupports, ...
    waveModeIndices=options.waveModeIndices, ...
    apvModeIndices=options.apvModeIndices, ...
    primitivePolynomialDegrees=[8;16;24], ...
    paddingFactors=options.paddingFactors);
if ~contract.isCompatible
    error("WVTerrainEnergyGalerkin:IncompleteCoupledBlockProductionContract", ...
        "Milestone 10.2.1 requires a passing production physical-state contract.")
end

audit = buildCoupledBlockInternalWaveCorrectionAudit(self, ...
    options.trustedModeBounds,supports, ...
    options.stationaryPolynomialDegree,degrees, ...
    options.paddingFactors,scales,options.tangentStep, ...
    options.maximumOuterIterations,options.maximumKrylovRestarts, ...
    options.shouldCompareSlopeSeed,quadratureOrder,contract);
end
