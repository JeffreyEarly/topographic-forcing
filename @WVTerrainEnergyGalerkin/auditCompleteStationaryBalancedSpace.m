function audit = auditCompleteStationaryBalancedSpace(self,options)
% Audit the complete finite-terrain stationary balanced space.
%
% This diagnostic constructs the terrain-dependent geostrophic states with
% unrestricted bottom displacement, enforces bottom tangency, and compares
% the resulting physical-energy subspace with the stationary nullspace of
% the unmodified finite-amplitude primitive exchange form. Non-tangent
% bottom directions remain in the full primitive state.
%
% ```matlab
% audit = problem.auditCompleteStationaryBalancedSpace( ...
%     trustedModeBounds=[1 0], ...
%     supportModeBounds=[1 1;1 2;1 4], ...
%     primitivePolynomialDegrees=[2;3;5], ...
%     paddingFactors=[2;3]);
% ```
%
% - Topic: Audit the complete stationary balanced space
% - Declaration: audit = auditCompleteStationaryBalancedSpace(options)
% - Parameter trustedModeBounds: nonnegative integer `[kMax lMax]`
% - Parameter supportModeBounds: one retained signed Fourier bound per refinement
% - Parameter primitivePolynomialDegrees: increasing primitive support degrees; the first is the trusted scalar degree
% - Parameter paddingFactors: horizontal oversampling factors
% - Parameter quadratureOrder: optional common Gauss--Legendre order
% - Returns audit: stationary-space, Green-identity, and refinement diagnostics
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.trustedModeBounds (1,2) double {mustBeInteger,mustBeNonnegative}
    options.supportModeBounds (:,2) double {mustBeInteger,mustBeNonnegative}
    options.primitivePolynomialDegrees (:,1) double {mustBeInteger,mustBePositive}
    options.paddingFactors (:,1) double {mustBeInteger,mustBePositive} = [2;3]
    options.quadratureOrder double {mustBeInteger} = []
end

degrees = options.primitivePolynomialDegrees;
supports = options.supportModeBounds;
if numel(degrees) < 3 || any(diff(degrees) <= 0) || any(degrees < 2)
    error("WVTerrainEnergyGalerkin:InvalidStationarySpaceDegrees", ...
        "primitivePolynomialDegrees must contain at least three strictly increasing values of two or greater.")
end
if size(supports,1) ~= numel(degrees) || any(diff(supports,1,1) < 0,"all")
    error("WVTerrainEnergyGalerkin:InvalidStationarySpaceSupports", ...
        "supportModeBounds must have one nondecreasing row per primitivePolynomialDegrees entry.")
end
if any(options.trustedModeBounds > supports(1,:))
    error("WVTerrainEnergyGalerkin:StationaryTrustedBandOutsideSupport", ...
        "trustedModeBounds must lie inside every supportModeBounds row.")
end
if numel(options.paddingFactors) < 2 || any(options.paddingFactors < 2) ...
        || numel(unique(options.paddingFactors)) ~= numel(options.paddingFactors)
    error("WVTerrainEnergyGalerkin:InvalidStationarySpacePadding", ...
        "paddingFactors must contain at least two unique values of two or greater.")
end
if ~isempty(options.quadratureOrder) ...
        && (~isscalar(options.quadratureOrder) || options.quadratureOrder < 2*max(degrees)+3)
    error("WVTerrainEnergyGalerkin:InsufficientQuadratureOrder", ...
        "quadratureOrder must be empty or at least 2*max(primitivePolynomialDegrees)+3.")
end

audit = buildCompleteStationaryBalancedSpaceAudit(self, ...
    options.trustedModeBounds,supports,degrees,options.paddingFactors, ...
    options.quadratureOrder);
end
