function audit = auditBoundaryCompleteWeakEigenproblem(self,options)
% Audit the boundary-complete weak terrain eigenproblem.
%
% This diagnostic constructs the terrain-dependent geostrophic test states
% derived from the finite-terrain APV Green identity. It tests whether those
% states lie in the primitive trial space, whether their weak rows reproduce
% the independently evaluated APV moments, and whether the resulting tangent
% eigenproblem preserves energy, APV moments, and bottom evolution.
%
% ```matlab
% audit = problem.auditBoundaryCompleteWeakEigenproblem( ...
%     trustedModeBounds=[1 0], ...
%     supportModeBounds=[1 2], ...
%     polynomialDegrees=[4;8;12], ...
%     paddingFactor=2);
% ```
%
% - Topic: Audit a boundary-complete weak eigenproblem
% - Declaration: audit = auditBoundaryCompleteWeakEigenproblem(options)
% - Parameter trustedModeBounds: nonnegative integer `[kMax lMax]`
% - Parameter supportModeBounds: retained signed Fourier support `[kMax lMax]`
% - Parameter polynomialDegrees: increasing support degrees; the first is the trusted scalar degree
% - Parameter paddingFactor: positive horizontal oversampling factor
% - Parameter quadratureOrder: optional common Gauss--Legendre order
% - Returns audit: weak-sequence, tangent-mode, and refinement diagnostics
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.trustedModeBounds (1,2) double {mustBeInteger,mustBeNonnegative}
    options.supportModeBounds (1,2) double {mustBeInteger,mustBeNonnegative}
    options.polynomialDegrees (:,1) double {mustBeInteger,mustBePositive} = [4;8;12]
    options.paddingFactor (1,1) double {mustBeInteger,mustBePositive} = 2
    options.quadratureOrder double {mustBeInteger} = []
end

degrees = options.polynomialDegrees;
if numel(degrees) < 3 || any(diff(degrees) <= 0) || any(degrees < 2)
    error("WVTerrainEnergyGalerkin:InvalidWeakEigenproblemDegrees", ...
        "polynomialDegrees must contain at least three strictly increasing values of two or greater.")
end
if any(options.trustedModeBounds > options.supportModeBounds)
    error("WVTerrainEnergyGalerkin:WeakEigenproblemTrustedBandOutsideSupport", ...
        "trustedModeBounds must lie inside supportModeBounds.")
end
if options.paddingFactor < 2
    error("WVTerrainEnergyGalerkin:InsufficientWeakEigenproblemPadding", ...
        "paddingFactor must be at least two.")
end
if ~isempty(options.quadratureOrder) ...
        && (~isscalar(options.quadratureOrder) || options.quadratureOrder < 2*max(degrees)+3)
    error("WVTerrainEnergyGalerkin:InsufficientQuadratureOrder", ...
        "quadratureOrder must be empty or at least 2*max(polynomialDegrees)+3.")
end
if max(abs(self.topographicHeight),[],"all") == 0
    error("WVTerrainEnergyGalerkin:WeakEigenproblemRequiresTerrainDirection", ...
        "topographicHeight must provide a nonzero terrain direction.")
end

audit = buildBoundaryCompleteWeakEigenproblemAudit(self, ...
    options.trustedModeBounds,options.supportModeBounds,degrees, ...
    options.paddingFactor,options.quadratureOrder);
end
