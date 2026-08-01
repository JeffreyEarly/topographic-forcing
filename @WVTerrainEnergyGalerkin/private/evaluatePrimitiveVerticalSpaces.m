function spaces = evaluatePrimitiveVerticalSpaces(wvt,degree, ...
    verticalCoordinate,xi)
% Evaluate complete primitive vertical spaces at arbitrary physical depths.
arguments
    wvt (1,1) WVTransformBoussinesq
    degree (1,1) double {mustBeInteger,mustBePositive}
    verticalCoordinate (1,1) string ...
        {mustBeMember(verticalCoordinate,["legendre","wkb"])}
    xi (:,1) double
end
D = wvt.Lz;
if verticalCoordinate == "legendre"
    coordinate = 2*xi/D+1;
    derivative = (2/D)*ones(size(xi));
    [F,Fc] = legendreValuesForPrimitiveAudit(coordinate,degree);
    [Pressure,PressureC] = legendreValuesForPrimitiveAudit( ...
        coordinate,degree+1);
else
    N2 = sampledN2(wvt,xi);
    if any(N2 <= 0)
        error("WVTerrainEnergyGalerkin:WKBPrimitiveStratificationNotPositive", ...
            "The WKB coordinate requires positive N2 at every evaluation depth.")
    end
    N = @(z)sqrt(sampledN2(wvt,z));
    referenceXi = linspace(-D,0,max(129,4*numel(xi)+1)).';
    referenceN2 = sampledN2(wvt,referenceXi);
    scale = max(referenceN2);
    if any(referenceN2 <= 100*eps*scale)
        error("WVTerrainEnergyGalerkin:WKBPrimitiveStratificationNotPositive", ...
            "The WKB coordinate requires N2 resolved away from zero.")
    end
    integralN = integral(N,-D,0,"AbsTol",1e-13*D*sqrt(scale), ...
        "RelTol",1e-13);
    coordinate = zeros(size(xi));
    for iPoint = 1:numel(xi)
        coordinate(iPoint) = 2*integral(N,-D,xi(iPoint), ...
            "AbsTol",1e-13*integralN,"RelTol",1e-13)/integralN-1;
    end
    derivative = 2*sqrt(N2)/integralN;
    [F,Fc] = chebyshevValues(coordinate,degree);
    [Pressure,PressureC] = chebyshevValues(coordinate,degree+1);
end
Fxi = derivative.*Fc;
G = (1-coordinate.^2).*F(:,1:degree);
Gc = -2*coordinate.*F(:,1:degree) ...
    +(1-coordinate.^2).*Fc(:,1:degree);
Gxi = derivative.*Gc;
H = [G (1-coordinate)/2];
Hxi = [Gxi -derivative/2];
spaces = struct("F",F,"Fxi",Fxi,"G",G,"Gxi",Gxi, ...
    "H",H,"Hxi",Hxi,"Pressure",Pressure, ...
    "PressureXi",derivative.*PressureC);
end

function values = sampledN2(wvt,xi)
targetSize = size(xi);
values = wvt.N2Function(xi);
if isscalar(values)
    values = repmat(values,size(xi));
end
if numel(values) ~= numel(xi) || ~isreal(values) ...
        || any(~isfinite(values))
    error("WVTerrainEnergyGalerkin:InvalidWKBPrimitiveStratification", ...
        "N2Function must return one finite real value per requested depth.")
end
values = reshape(values,targetSize);
end

function [T,Tc] = chebyshevValues(s,degree)
T = zeros(numel(s),degree+1);
Tc = zeros(numel(s),degree+1);
T(:,1) = 1;
if degree == 0
    return
end
T(:,2) = s;
Tc(:,2) = 1;
for n = 2:degree
    T(:,n+1) = 2*s.*T(:,n)-T(:,n-1);
    Tc(:,n+1) = 2*T(:,n)+2*s.*Tc(:,n)-Tc(:,n-1);
end
end
