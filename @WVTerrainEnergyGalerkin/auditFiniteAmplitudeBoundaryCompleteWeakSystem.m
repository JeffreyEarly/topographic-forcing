function audit = auditFiniteAmplitudeBoundaryCompleteWeakSystem(self,options)
% Audit the finite-amplitude boundary-complete primitive weak system.
%
% This diagnostic constructs the exact terrain-dependent geostrophic test
% states at each requested terrain scale. It compares the APV moments
% derived from the primitive energy rows with an independent strong APV
% evaluation, and audits physical energy and projected bottom evolution
% without replacing any primitive evolution row.
%
% ```matlab
% audit = problem.auditFiniteAmplitudeBoundaryCompleteWeakSystem( ...
%     trustedModeBounds=[1 0], ...
%     supportModeBounds=[1 2;1 3;1 4], ...
%     scalarPolynomialDegree=4, ...
%     primitivePolynomialDegrees=[4;8;12], ...
%     paddingFactors=[2;3], ...
%     terrainScales=[0.25;0.5;1]);
% ```
%
% - Topic: Audit a finite-amplitude boundary-complete weak system
% - Declaration: audit = auditFiniteAmplitudeBoundaryCompleteWeakSystem(options)
% - Parameter trustedModeBounds: nonnegative integer `[kMax lMax]`
% - Parameter supportModeBounds: one retained signed Fourier bound per refinement
% - Parameter scalarPolynomialDegree: fixed APV-test polynomial degree
% - Parameter primitivePolynomialDegrees: increasing primitive support degrees
% - Parameter paddingFactors: horizontal oversampling factors
% - Parameter terrainScales: nonnegative multipliers of `topographicHeight`
% - Parameter quadratureOrder: optional common Gauss--Legendre order
% - Returns audit: finite-amplitude weak-sequence and refinement diagnostics
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.trustedModeBounds (1,2) double {mustBeInteger,mustBeNonnegative}
    options.supportModeBounds (:,2) double {mustBeInteger,mustBeNonnegative}
    options.scalarPolynomialDegree (1,1) double {mustBeInteger,mustBePositive}
    options.primitivePolynomialDegrees (:,1) double {mustBeInteger,mustBePositive}
    options.paddingFactors (:,1) double {mustBeInteger,mustBePositive} = [2;3]
    options.terrainScales (:,1) double {mustBeNonnegative} = [0.25;0.5;1]
    options.quadratureOrder double {mustBeInteger} = []
end

degrees = options.primitivePolynomialDegrees;
supports = options.supportModeBounds;
if numel(degrees) < 3 || any(diff(degrees) <= 0) || any(degrees < options.scalarPolynomialDegree)
    error("WVTerrainEnergyGalerkin:InvalidFiniteAmplitudeDegrees", ...
        "primitivePolynomialDegrees must contain at least three strictly increasing values no smaller than scalarPolynomialDegree.")
end
if size(supports,1) ~= numel(degrees) || any(diff(supports,1,1) < 0,"all")
    error("WVTerrainEnergyGalerkin:InvalidFiniteAmplitudeSupports", ...
        "supportModeBounds must have one nondecreasing row per primitivePolynomialDegrees entry.")
end
if any(options.trustedModeBounds > supports(1,:))
    error("WVTerrainEnergyGalerkin:FiniteAmplitudeTrustedBandOutsideSupport", ...
        "trustedModeBounds must lie inside every supportModeBounds row.")
end
if numel(options.paddingFactors) < 2 || any(options.paddingFactors < 2) || numel(unique(options.paddingFactors)) ~= numel(options.paddingFactors)
    error("WVTerrainEnergyGalerkin:InvalidFiniteAmplitudePadding", ...
        "paddingFactors must contain at least two unique values of two or greater.")
end
if isempty(options.terrainScales) || any(~isfinite(options.terrainScales)) || ~any(options.terrainScales == 1)
    error("WVTerrainEnergyGalerkin:InvalidFiniteAmplitudeScales", ...
        "terrainScales must contain finite nonnegative values and include scale one.")
end
if any(options.terrainScales ~= sort(options.terrainScales)) || numel(unique(options.terrainScales)) ~= numel(options.terrainScales)
    error("WVTerrainEnergyGalerkin:InvalidFiniteAmplitudeScales", ...
        "terrainScales must be strictly increasing.")
end
if any(1-options.terrainScales*max(self.topographicHeight,[],"all")/self.originatingTransform.Lz <= 0)
    error("WVTerrainEnergyGalerkin:NonpositiveFiniteAmplitudeGamma", ...
        "Every requested terrain scale must satisfy gamma>0.")
end
if ~isempty(options.quadratureOrder) ...
        && (~isscalar(options.quadratureOrder) || options.quadratureOrder < 2*max(degrees)+3)
    error("WVTerrainEnergyGalerkin:InsufficientQuadratureOrder", ...
        "quadratureOrder must be empty or at least 2*max(primitivePolynomialDegrees)+3.")
end

audit = buildFiniteAmplitudeBoundaryCompleteWeakAudit(self, ...
    options.trustedModeBounds,supports,options.scalarPolynomialDegree,degrees, ...
    options.paddingFactors,options.terrainScales,options.quadratureOrder);
end
