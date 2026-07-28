function audit = auditAPVCompatiblePrimitive(self,options)
% Audit APV-compatible local primitive formulations.
%
% The audit compares an energy-weak mixed polynomial descriptor with
% vorticity--divergence and explicit-APV formulations of the frozen
% constant-slope Branch-P equations. It also evaluates the analytic APV
% source left by freezing a cross-slope Fourier component.
%
% ```matlab
% audit = problem.auditAPVCompatiblePrimitive( ...
%     bottomSlope=[0 0.01],horizontalMode=[1 0], ...
%     polynomialDegree=6);
% ```
%
% - Topic: Audit the Branch-P discrete oracle
% - Declaration: audit = auditAPVCompatiblePrimitive(options)
% - Parameter bottomSlope: real two-component local slope `[h_x h_y]`
% - Parameter horizontalMode: retained integer Fourier mode `[kMode lMode]`
% - Parameter polynomialDegree: positive Legendre polynomial degree
% - Parameter quadratureOrder: optional Gauss--Legendre quadrature order
% - Returns audit: candidate diagnostics and mathematical classification
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.bottomSlope (1,2) double = [0 0.01]
    options.horizontalMode (1,2) double {mustBeInteger} = [1 0]
    options.polynomialDegree (1,1) double {mustBeInteger,mustBePositive} = 6
    options.quadratureOrder double {mustBeInteger} = []
end

if ~isreal(options.bottomSlope) || any(~isfinite(options.bottomSlope))
    error("WVTerrainEnergyGalerkin:InvalidBottomSlope", ...
        "bottomSlope must contain two real finite values.")
end
if any(~isfinite(options.horizontalMode))
    error("WVTerrainEnergyGalerkin:InvalidHorizontalMode", ...
        "horizontalMode must contain two finite integer values.")
end
if all(options.horizontalMode == 0)
    error("WVTerrainEnergyGalerkin:ZeroHorizontalMode", ...
        "The APV-compatible local audit requires a nonzero horizontal mode.")
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

iK = find(self.horizontalLayout.kMode == options.horizontalMode(1) ...
    & self.horizontalLayout.lMode == options.horizontalMode(2),1);
if isempty(iK)
    error("WVTerrainEnergyGalerkin:HorizontalModeNotRetained", ...
        "horizontalMode [%d %d] is not retained by the originating transform.", ...
        options.horizontalMode(1),options.horizontalMode(2))
end

audit = buildAPVCompatiblePrimitiveAudit(self,iK,options.bottomSlope, ...
    options.polynomialDegree,quadratureOrder);
end
