function audit = auditPeriodicCoupledPVOracle(self,options)
% Audit the projected periodic volume--boundary PV dynamics.
%
% This diagnostic evolves volume APV and active bottom PV as explicit
% coordinates. It compares an exact, non-aliased Fourier convolution with
% an independently oversampled pseudospectral bottom Jacobian and audits
% physical energy, stationary volume APV, potential enstrophy, projected
% bottom evolution, Fourier conjugacy, and representation in the existing
% balanced-plus-bottom coordinates.
%
% ```matlab
% audit = problem.auditPeriodicCoupledPVOracle( ...
%     polynomialDegree=12);
% ```
%
% - Topic: Audit periodic coupled PV closure
% - Declaration: audit = auditPeriodicCoupledPVOracle(options)
% - Parameter polynomialDegree: positive vertical polynomial degree
% - Parameter quadratureOrder: optional Gauss--Legendre quadrature order
% - Returns audit: periodic coupled-PV forms and closure diagnostics
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.polynomialDegree (1,1) double {mustBeInteger,mustBePositive} = 12
    options.quadratureOrder double {mustBeInteger} = []
end

if options.polynomialDegree < 3
    error("WVTerrainEnergyGalerkin:InvalidPolynomialDegree", ...
        "polynomialDegree must be at least three.")
end
if ~isempty(options.quadratureOrder) && (~isscalar(options.quadratureOrder) || options.quadratureOrder < options.polynomialDegree+3)
    error("WVTerrainEnergyGalerkin:InsufficientQuadratureOrder", ...
        "quadratureOrder must be empty or at least polynomialDegree+3.")
end
if isempty(options.quadratureOrder)
    quadratureOrder = max(2*options.polynomialDegree+7,32);
else
    quadratureOrder = options.quadratureOrder;
end

audit = buildPeriodicCoupledPVOracleAudit(self,options.polynomialDegree,quadratureOrder);
end
