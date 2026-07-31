function audit = auditProductionPhysicalStateContract(self,options)
% Audit the flat physical-state contract for future forward evolution.
%
% The production state contains native fixed-wavenumber Boussinesq wave
% pairs, APV-bearing balanced coordinates, one complete zero-APV bottom
% inversion at every nonzero horizontal wavenumber, and the compatible
% inertial and mean-density sector at zero horizontal wavenumber. The
% primitive polynomial representation remains an independent oracle.
%
% ```matlab
% audit = problem.auditProductionPhysicalStateContract( ...
%     trustedModeBounds=[1 0], ...
%     supportModeBounds=[1 0;1 1;1 2], ...
%     primitivePolynomialDegrees=[8;16;24], ...
%     paddingFactors=[2;3]);
% ```
%
% - Topic: Audit the production physical-state contract
% - Declaration: audit = auditProductionPhysicalStateContract(options)
% - Parameter trustedModeBounds: nonnegative integer `[kMax lMax]`
% - Parameter supportModeBounds: nested signed-Fourier oracle bounds
% - Parameter waveModeIndices: positive retained fixed-wavenumber mode labels
% - Parameter apvModeIndices: retained balanced labels, including mode zero
% - Parameter primitivePolynomialDegrees: increasing primitive oracle degrees
% - Parameter paddingFactors: distinct horizontal oversampling factors
% - Parameter quadratureOrder: optional common vertical quadrature order
% - Returns audit: production layout, primitive maps, and contract diagnostics
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.trustedModeBounds (1,2) double {mustBeInteger,mustBeNonnegative}
    options.supportModeBounds (:,2) double {mustBeInteger,mustBeNonnegative}
    options.waveModeIndices (:,1) double {mustBeInteger} = self.verticalModeIndices(self.verticalModeIndices > 0)
    options.apvModeIndices (:,1) double {mustBeInteger} = self.verticalModeIndices
    options.primitivePolynomialDegrees (:,1) double {mustBeInteger,mustBePositive}
    options.paddingFactors (:,1) double {mustBeInteger,mustBePositive} = [2;3]
    options.quadratureOrder double {mustBeInteger} = []
end

if numel(options.primitivePolynomialDegrees) < 3 || any(diff(options.primitivePolynomialDegrees) <= 0)
    error("WVTerrainEnergyGalerkin:InvalidProductionPrimitiveDegrees", ...
        "primitivePolynomialDegrees must contain at least three strictly increasing values.")
end
if size(options.supportModeBounds,1) < 3 || any(diff(options.supportModeBounds,1,1) < 0,"all") || any(~any(diff(options.supportModeBounds,1,1) > 0,2))
    error("WVTerrainEnergyGalerkin:InvalidProductionSupportBounds", ...
        "supportModeBounds must contain at least three strictly expanding nested rows.")
end
if any(options.trustedModeBounds > options.supportModeBounds(1,:))
    error("WVTerrainEnergyGalerkin:ProductionTrustedBandOutsideSupport", ...
        "trustedModeBounds must lie inside every supportModeBounds row.")
end
if isempty(options.waveModeIndices) || any(options.waveModeIndices <= 0) || ~issorted(options.waveModeIndices) || numel(unique(options.waveModeIndices)) ~= numel(options.waveModeIndices) || any(~ismember(options.waveModeIndices,self.verticalModeIndices))
    error("WVTerrainEnergyGalerkin:InvalidProductionWaveModeIndices", ...
        "waveModeIndices must be sorted, unique positive labels retained by the problem.")
end
if isempty(options.apvModeIndices) || options.apvModeIndices(1) ~= 0 || ~issorted(options.apvModeIndices) || numel(unique(options.apvModeIndices)) ~= numel(options.apvModeIndices) || any(~ismember(options.apvModeIndices,self.verticalModeIndices))
    error("WVTerrainEnergyGalerkin:InvalidProductionAPVModeIndices", ...
        "apvModeIndices must be sorted, unique retained labels and include mode zero.")
end
if numel(options.paddingFactors) < 2 || any(options.paddingFactors < 2) || numel(unique(options.paddingFactors)) ~= numel(options.paddingFactors)
    error("WVTerrainEnergyGalerkin:InvalidProductionPaddingFactors", ...
        "paddingFactors must contain at least two unique values of two or greater.")
end
if ~isempty(options.quadratureOrder) && (~isscalar(options.quadratureOrder) || options.quadratureOrder < 2*max(options.primitivePolynomialDegrees)+3)
    error("WVTerrainEnergyGalerkin:InsufficientProductionQuadrature", ...
        "quadratureOrder must be empty or at least 2*max(primitivePolynomialDegrees)+3.")
end

audit = buildProductionPhysicalStateContractAudit(self, ...
    options.trustedModeBounds,options.supportModeBounds, ...
    options.waveModeIndices,options.apvModeIndices, ...
    options.primitivePolynomialDegrees,options.paddingFactors, ...
    options.quadratureOrder);
end
