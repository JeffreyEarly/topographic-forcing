function audit = auditDealiasedProjectedPrimitiveTangent(self,options)
% Audit the dealiased projected primitive terrain tangent.
%
% This diagnostic retains the unmodified primitive weak equations and
% evaluates every first-order terrain product with one zero-pad,
% multiply, and adjoint-restrict Fourier projection. Nested horizontal
% support, vertical polynomial, and padding refinements distinguish exact
% projected APV closure from trusted-band convergence and from a genuine
% nonconvergent result.
%
% ```matlab
% audit = problem.auditDealiasedProjectedPrimitiveTangent( ...
%     trustedModeBounds=[1 0], ...
%     supportModeBounds=[1 1;1 2;1 3]);
% ```
%
% - Topic: Audit a dealiased projected primitive tangent
% - Declaration: audit = auditDealiasedProjectedPrimitiveTangent(options)
% - Parameter trustedModeBounds: nonnegative integer `[kMax lMax]`
% - Parameter supportModeBounds: nested integer `[kMax lMax]` rows
% - Parameter polynomialDegrees: increasing vertical polynomial degrees
% - Parameter paddingFactors: horizontal padding factors including 2 and 3
% - Parameter quadratureOrder: optional common Gauss--Legendre order
% - Parameter tangentStep: positive dimensionless centered-difference step
% - Returns audit: projected primitive refinement results and classification
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.trustedModeBounds (1,2) double {mustBeInteger,mustBeNonnegative}
    options.supportModeBounds (:,2) double {mustBeInteger,mustBeNonnegative}
    options.polynomialDegrees (:,1) double {mustBeInteger,mustBePositive} = [4;8;12]
    options.paddingFactors (:,1) double {mustBeInteger,mustBePositive} = [2;3]
    options.quadratureOrder double {mustBeInteger} = []
    options.tangentStep (1,1) double {mustBePositive} = 1e-3
end

support = options.supportModeBounds;
degrees = options.polynomialDegrees;
padding = sort(unique(options.paddingFactors(:)));
if size(support,1) < 3 || numel(degrees) ~= size(support,1)
    error("WVTerrainEnergyGalerkin:InvalidProjectedPrimitiveRefinement", ...
        "supportModeBounds and polynomialDegrees must define the same number of at least three refinement levels.")
end
if any(any(diff(support,1,1) < 0)) || any(all(diff(support,1,1) == 0,2))
    error("WVTerrainEnergyGalerkin:InvalidProjectedPrimitiveSupport", ...
        "supportModeBounds must be nested and each row must enlarge at least one bound.")
end
if any(diff(degrees) <= 0) || any(degrees < 2)
    error("WVTerrainEnergyGalerkin:InvalidProjectedPrimitiveDegrees", ...
        "polynomialDegrees must be strictly increasing and at least two.")
end
if ~all(ismember([2;3],padding))
    error("WVTerrainEnergyGalerkin:InvalidProjectedPrimitivePadding", ...
        "paddingFactors must include factors two and three.")
end
if ~isempty(options.quadratureOrder) ...
        && (~isscalar(options.quadratureOrder) || options.quadratureOrder < max(degrees)+3)
    error("WVTerrainEnergyGalerkin:InsufficientQuadratureOrder", ...
        "quadratureOrder must be empty or at least max(polynomialDegrees)+3.")
end
if ~isfinite(options.tangentStep) || options.tangentStep >= 0.1
    error("WVTerrainEnergyGalerkin:InvalidTangentStep", ...
        "tangentStep must be finite and less than 0.1.")
end
if any(options.trustedModeBounds > support(1,:))
    error("WVTerrainEnergyGalerkin:TrustedBandOutsideSupport", ...
        "trustedModeBounds must lie inside the smallest supportModeBounds row.")
end
if max(abs(self.topographicHeight),[],"all") == 0
    error("WVTerrainEnergyGalerkin:ProjectedPrimitiveAuditRequiresTerrainDirection", ...
        "topographicHeight must provide a nonzero terrain direction.")
end

audit = buildDealiasedProjectedPrimitiveTangentAudit(self, ...
    options.trustedModeBounds,support,degrees,padding, ...
    options.quadratureOrder,options.tangentStep);
end
