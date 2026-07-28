function audit = auditBranchPDiscreteOracle(self,options)
% Audit an independent primitive Branch-P discretization.
%
% The oracle uses Legendre polynomial trial spaces, an eta-only bottom
% coordinate, and a pressure-retaining constrained saddle solve. It is
% independent of the hydrostatic F--G reconstruction used by the public
% Galerkin state.
%
% ```matlab
% audit = problem.auditBranchPDiscreteOracle( ...
%     bottomSlope=[0.01 0],polynomialDegree=6);
% ```
%
% - Topic: Audit boundary-dynamical evolution
% - Declaration: audit = auditBranchPDiscreteOracle(options)
% - Parameter bottomSlope: real two-component vector `[h_x h_y]`
% - Parameter polynomialDegree: positive Legendre polynomial degree
% - Parameter quadratureOrder: optional Gauss--Legendre quadrature order
% - Returns audit: primitive oracle, F--G comparison, and repair diagnosis
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.bottomSlope (1,2) double = [0.01 0]
    options.polynomialDegree (1,1) double {mustBeInteger,mustBePositive} = 6
    options.quadratureOrder double {mustBeInteger} = []
end

if ~isreal(options.bottomSlope) || any(~isfinite(options.bottomSlope))
    error("WVTerrainEnergyGalerkin:InvalidBottomSlope", ...
        "bottomSlope must contain two real finite values.")
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
    quadratureOrder = max(2*options.polynomialDegree+5,16);
else
    quadratureOrder = options.quadratureOrder;
end
if quadratureOrder < options.polynomialDegree+2
    error("WVTerrainEnergyGalerkin:InsufficientQuadratureOrder", ...
        "quadratureOrder must be at least polynomialDegree+2.")
end

audit = buildBranchPDiscreteOracleAudit(self,options.bottomSlope, ...
    options.polynomialDegree,quadratureOrder);
end
