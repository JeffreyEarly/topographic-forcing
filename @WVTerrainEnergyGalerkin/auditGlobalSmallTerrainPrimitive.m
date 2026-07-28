function audit = auditGlobalSmallTerrainPrimitive(self,options)
% Audit the global first-order primitive terrain equations.
%
% This diagnostic differentiates a globally coupled primitive-variable
% saddle system about flat topography in the direction supplied by
% `topographicHeight`. It compares an analytic terrain derivative with
% centered differences of the unexpanded mapped equations and audits
% physical energy, APV, potential enstrophy, and bottom evolution before
% any terrain modes are constructed.
%
% ```matlab
% audit = problem.auditGlobalSmallTerrainPrimitive( ...
%     polynomialDegree=6,tangentStep=1e-4);
% ```
%
% - Topic: Audit the Branch-P discrete oracle
% - Declaration: audit = auditGlobalSmallTerrainPrimitive(options)
% - Parameter polynomialDegree: positive Legendre polynomial degree
% - Parameter quadratureOrder: optional Gauss--Legendre quadrature order
% - Parameter tangentStep: positive dimensionless centered-difference step
% - Returns audit: global first-order forms, diagnostics, and classification
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.polynomialDegree (1,1) double {mustBeInteger,mustBePositive} = 6
    options.quadratureOrder double {mustBeInteger} = []
    options.tangentStep (1,1) double {mustBePositive} = 1e-3
end

if options.polynomialDegree < 2
    error("WVTerrainEnergyGalerkin:InvalidPolynomialDegree", ...
        "polynomialDegree must be at least two.")
end
if ~isempty(options.quadratureOrder) && (~isscalar(options.quadratureOrder) || options.quadratureOrder < 1)
    error("WVTerrainEnergyGalerkin:InvalidQuadratureOrder", ...
        "quadratureOrder must be empty or a positive integer.")
end
if isempty(options.quadratureOrder)
    quadratureOrder = max(2*options.polynomialDegree+7,18);
else
    quadratureOrder = options.quadratureOrder;
end
if quadratureOrder < options.polynomialDegree+3
    error("WVTerrainEnergyGalerkin:InsufficientQuadratureOrder", ...
        "quadratureOrder must be at least polynomialDegree+3.")
end
if ~isfinite(options.tangentStep) || options.tangentStep >= 0.1
    error("WVTerrainEnergyGalerkin:InvalidTangentStep", ...
        "tangentStep must be finite and less than 0.1.")
end
if max(abs(self.topographicHeight),[],"all") == 0
    error("WVTerrainEnergyGalerkin:GlobalPrimitiveAuditRequiresTerrainDirection", ...
        "topographicHeight must provide a nonzero terrain direction.")
end

audit = buildGlobalSmallTerrainPrimitiveAudit(self,options.polynomialDegree, ...
    quadratureOrder,options.tangentStep);
end
