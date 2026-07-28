function audit = auditCoupledPVFrequencyOracle(self,options)
% Audit the coupled volume--boundary PV frequency problem.
%
% This diagnostic solves the single-wavenumber QG evolution of volume APV
% and active bottom PV with a direct Legendre--Galerkin inversion. When
% both projected PV gradients are nonzero, it compares the physical
% frequencies with an independent Chebyshev--Lobatto discretization of
% Yassin's eigenvalue-dependent endpoint problem. The resulting
% geostrophic fields are projected into the existing balanced-plus-bottom
% Galerkin coordinates without changing the public coefficient layout.
%
% ```matlab
% audit = problem.auditCoupledPVFrequencyOracle( ...
%     horizontalMode=[1 0],volumePVGradient=[0 2e-11], ...
%     bottomSlope=[0 0.01],polynomialDegree=20);
% ```
%
% - Topic: Audit coupled PV dynamics
% - Declaration: audit = auditCoupledPVFrequencyOracle(options)
% - Parameter horizontalMode: retained nonzero integer mode `[k l]`
% - Parameter volumePVGradient: horizontal background APV gradient
% - Parameter bottomSlope: dimensionless horizontal bottom slope
% - Parameter polynomialDegree: positive vertical polynomial degree
% - Parameter quadratureOrder: optional Gauss--Legendre quadrature order
% - Returns audit: coupled-PV frequencies, invariants, and basis diagnostics
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.horizontalMode (1,2) double {mustBeInteger} = [1 0]
    options.volumePVGradient (1,2) double = [0 2e-11]
    options.bottomSlope (1,2) double = [0 0.01]
    options.polynomialDegree (1,1) double {mustBeInteger,mustBePositive} = 20
    options.quadratureOrder double {mustBeInteger} = []
end

if any(~isfinite(options.horizontalMode)) || all(options.horizontalMode == 0)
    error("WVTerrainEnergyGalerkin:InvalidCoupledPVHorizontalMode", ...
        "horizontalMode must contain a retained nonzero integer mode.")
end
if any(~isfinite(options.volumePVGradient)) || ~isreal(options.volumePVGradient)
    error("WVTerrainEnergyGalerkin:InvalidVolumePVGradient", ...
        "volumePVGradient must contain two real finite values.")
end
if any(~isfinite(options.bottomSlope)) || ~isreal(options.bottomSlope)
    error("WVTerrainEnergyGalerkin:InvalidBottomSlope", ...
        "bottomSlope must contain two real finite values.")
end
if options.polynomialDegree < 3
    error("WVTerrainEnergyGalerkin:InvalidPolynomialDegree", ...
        "polynomialDegree must be at least three.")
end
if ~isempty(options.quadratureOrder) && (~isscalar(options.quadratureOrder) || options.quadratureOrder < options.polynomialDegree+3)
    error("WVTerrainEnergyGalerkin:InsufficientQuadratureOrder", ...
        "quadratureOrder must be empty or at least polynomialDegree+3.")
end

iK = find(self.horizontalLayout.kMode == options.horizontalMode(1) & self.horizontalLayout.lMode == options.horizontalMode(2),1);
if isempty(iK)
    error("WVTerrainEnergyGalerkin:UnknownHorizontalMode", ...
        "The horizontal mode (%d,%d) is not retained.",options.horizontalMode(1),options.horizontalMode(2))
end
if isempty(options.quadratureOrder)
    quadratureOrder = max(2*options.polynomialDegree+7,28);
else
    quadratureOrder = options.quadratureOrder;
end

audit = buildCoupledPVFrequencyOracleAudit(self,iK,options.volumePVGradient, ...
    options.bottomSlope,options.polynomialDegree,quadratureOrder);
end
