function audit = auditHybridPrimitivePVOracle(self,options)
% Audit a hybrid primitive, volume-APV, and bottom-coordinate descriptor.
%
% This diagnostic replaces the redundant stationary primitive rows in one
% fixed zonal Fourier block with independent projected volume-APV and
% bottom-evolution rows. It then asks whether that complete coordinate
% system reproduces the independently derived first-order primitive weak
% equations. No closure correction is fitted or applied.
%
% ```matlab
% audit = problem.auditHybridPrimitivePVOracle( ...
%     zonalMode=1,polynomialDegree=4);
% ```
%
% - Topic: Audit a hybrid primitive-PV descriptor
% - Declaration: audit = auditHybridPrimitivePVOracle(options)
% - Parameter zonalMode: nonzero retained integer zonal mode
% - Parameter polynomialDegree: positive Legendre polynomial degree
% - Parameter quadratureOrder: optional Gauss--Legendre quadrature order
% - Parameter tangentStep: positive dimensionless centered-difference step
% - Returns audit: hybrid descriptor forms, diagnostics, and classification
arguments
    self (1,1) WVTerrainEnergyGalerkin
    options.zonalMode (1,1) double {mustBeInteger} = 1
    options.polynomialDegree (1,1) double {mustBeInteger,mustBePositive} = 4
    options.quadratureOrder double {mustBeInteger} = []
    options.tangentStep (1,1) double {mustBePositive} = 1e-3
end

if options.zonalMode == 0
    error("WVTerrainEnergyGalerkin:HybridOracleZeroZonalMode", ...
        "zonalMode must be nonzero.")
end
if options.polynomialDegree < 2
    error("WVTerrainEnergyGalerkin:InvalidPolynomialDegree", ...
        "polynomialDegree must be at least two.")
end
if ~isempty(options.quadratureOrder) ...
        && (~isscalar(options.quadratureOrder) || options.quadratureOrder < options.polynomialDegree+3)
    error("WVTerrainEnergyGalerkin:InsufficientQuadratureOrder", ...
        "quadratureOrder must be empty or at least polynomialDegree+3.")
end
if ~isfinite(options.tangentStep) || options.tangentStep >= 0.1
    error("WVTerrainEnergyGalerkin:InvalidTangentStep", ...
        "tangentStep must be finite and less than 0.1.")
end

xVariation = self.topographicHeight-self.topographicHeight(1,:);
terrainScale = max(abs(self.topographicHeight),[],"all");
if max(abs(xVariation),[],"all") > 1e3*eps*max(terrainScale,1)
    error("WVTerrainEnergyGalerkin:HybridOracleRequiresZonallyInvariantTerrain", ...
        "The fixed-zonal hybrid oracle requires topography independent of x.")
end
if ~any(self.horizontalLayout.kMode == options.zonalMode) ...
        || ~any(self.horizontalLayout.kMode == -options.zonalMode)
    error("WVTerrainEnergyGalerkin:HybridOracleMissingZonalBlock", ...
        "Both requested zonal block and its Fourier conjugate must be retained.")
end

primitive = self.auditGlobalSmallTerrainPrimitive( ...
    polynomialDegree=options.polynomialDegree, ...
    quadratureOrder=options.quadratureOrder, ...
    tangentStep=options.tangentStep);
audit = buildHybridPrimitivePVOracleAudit(self,primitive,options.zonalMode);
end
