function audit = auditBoundaryGreenIdentity(self,options)
% Audit the local constant-slope boundary Green identity.
%
% The audit freezes the terrain scale at the uniform reference depth while
% retaining the supplied bottom slope in the mapped metric and active bottom
% equation. Horizontal Fourier coefficients therefore remain independent.
%
% ```matlab
% audit = problem.auditBoundaryGreenIdentity(bottomSlope=[0.01 0]);
% ```
%
% - Topic: Audit boundary-dynamical evolution
% - Declaration: audit = auditBoundaryGreenIdentity(options)
% - Parameter bottomSlope: real two-component vector `[h_x h_y]`
% - Returns audit: local descriptor, invariant diagnostics, and branch classification
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.bottomSlope (1,2) double
end

if ~isreal(options.bottomSlope) || any(~isfinite(options.bottomSlope))
    error("WVTerrainEnergyGalerkin:InvalidBottomSlope", ...
        "bottomSlope must contain two real finite values.")
end
audit = buildBoundaryGreenIdentityAudit(self,options.bottomSlope);
end
