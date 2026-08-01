function [quadrature,spaces,diagnostics] = buildPrimitiveVerticalDiscretization( ...
    wvt,degree,order,verticalCoordinate)
% Build the complete primitive vertical quadrature and mixed spaces.
arguments
    wvt (1,1) WVTransformBoussinesq
    degree (1,1) double {mustBeInteger,mustBePositive}
    order (1,1) double {mustBeInteger,mustBePositive}
    verticalCoordinate (1,1) string ...
        {mustBeMember(verticalCoordinate,["legendre","wkb"])}
end

[r,weightR] = gaussLegendreRule(order);
D = wvt.Lz;
if verticalCoordinate == "legendre"
    xi = D*(r-1)/2;
    coordinateDerivative = (2/D)*ones(size(r));
    weight = D*weightR/2;
    endpointCoordinateDerivative = (2/D)*ones(2,1);
    spaces = legendreSpaces(degree,r,coordinateDerivative, ...
        endpointCoordinateDerivative,D);
    coordinate = struct("r",r,"xi",xi,"weight",weight, ...
        "coordinateDerivative",coordinateDerivative);
    diagnostics = coordinateDiagnostics("legendre",coordinate,D,1, ...
        zeros(size(r)),true,"");
else
    [coordinate,diagnostics] = wkbCoordinate(wvt,r,weightR);
    spaces = chebyshevSpaces(degree,coordinate.r, ...
        coordinate.coordinateDerivative, ...
        coordinate.endpointCoordinateDerivative,D);
end
quadrature = coordinate;
end

function [coordinate,diagnostics] = wkbCoordinate(wvt,s,weightS)
D = wvt.Lz;
qualificationCount = max(129,4*numel(s)+1);
qualificationXi = D*(cos(pi*(0:qualificationCount-1)' ...
    /(qualificationCount-1))-1)/2;
qualificationN2 = sampledN2(wvt,qualificationXi);
scale = max(qualificationN2);
if any(qualificationN2 <= 100*eps*scale)
    error("WVTerrainEnergyGalerkin:WKBPrimitiveStratificationNotPositive", ...
        "The WKB coordinate requires N2 to remain positive and resolved away from zero.")
end

N = @(xi)sqrt(sampledN2(wvt,xi));
integralN = integral(N,-D,0,"AbsTol",1e-13*D*sqrt(scale), ...
    "RelTol",1e-13);
if ~isfinite(integralN) || integralN <= 0
    error("WVTerrainEnergyGalerkin:InvalidWKBPrimitiveIntegral", ...
        "The reference buoyancy-frequency integral must be positive and finite.")
end

targets = [-1;s;1];
odeOptions = odeset("RelTol",1e-13,"AbsTol",1e-13*D);
[~,xiValues] = ode113(@(~,xi)integralN./(2*N(xi)), ...
    targets,-D,odeOptions);
xiValues = xiValues(:);
xiValues(1) = -D;
xiValues(end) = 0;
xi = xiValues(2:end-1);
for iPoint = 1:numel(xi)
    xi(iPoint) = refinedInverse(N,integralN,s(iPoint),xi(iPoint),D);
end

NAtXi = N(xi);
coordinateDerivative = 2*NAtXi/integralN;
jacobian = 1./coordinateDerivative;
weight = weightS.*jacobian;
forward = arrayfun(@(value)2*integral(N,-D,value, ...
    "AbsTol",1e-13*D*sqrt(scale),"RelTol",1e-13)/integralN-1,xi);
forward = forward(:);
forwardDefect = max(abs(forward-s));
jacobianDefect = max(abs(jacobian.*coordinateDerivative-1));
endpointDefect = max(abs([2*integral(N,-D,-D)/integralN; ...
    2*integral(N,-D,0,"AbsTol",1e-13*D*sqrt(scale), ...
    "RelTol",1e-13)/integralN-2]));
quadratureDefect = abs(sum(weight)-D)/D;
conditionNumber = max(NAtXi)/min(NAtXi);
isQualified = forwardDefect <= 1e-12 ...
    && jacobianDefect <= 1e-12 && endpointDefect <= 1e-12 ...
    && all(diff(xi) > 0);
if ~isQualified
    error("WVTerrainEnergyGalerkin:WKBPrimitiveCoordinateFailure", ...
        "The WKB coordinate failed its inverse, Jacobian, endpoint, monotonicity, or quadrature checks.")
end
coordinate = struct("r",s,"xi",xi,"weight",weight, ...
    "coordinateDerivative",coordinateDerivative, ...
    "endpointCoordinateDerivative",2*N([-D;0])/integralN);
diagnostics = coordinateDiagnostics("wkb",coordinate,D,conditionNumber, ...
    forward-s,isQualified,"");
diagnostics.forwardMapDefect = forwardDefect;
diagnostics.inverseMapDefect = forwardDefect;
diagnostics.jacobianDefect = jacobianDefect;
diagnostics.endpointDefect = endpointDefect;
diagnostics.quadratureDefect = quadratureDefect;
diagnostics.integratedBuoyancyFrequency = integralN;
diagnostics.minimumN2 = min(qualificationN2);
diagnostics.maximumN2 = max(qualificationN2);
end

function xi = refinedInverse(N,integralN,target,initial,D)
xi = min(max(initial,-D),0);
lower = -D;
upper = 0;
for iIteration = 1:8
    residual = 2*integral(N,-D,xi,"AbsTol",1e-13*integralN, ...
        "RelTol",1e-13)/integralN-1-target;
    if abs(residual) <= 2e-14
        return
    end
    if residual < 0
        lower = xi;
    else
        upper = xi;
    end
    candidate = xi-residual*integralN/(2*N(xi));
    if candidate <= lower || candidate >= upper || ~isfinite(candidate)
        candidate = (lower+upper)/2;
    end
    xi = candidate;
end
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

function spaces = legendreSpaces(degree,r,coordinateDerivative, ...
    endpointCoordinateDerivative,D)
[F,Fr] = legendreValuesForPrimitiveAudit(r,degree);
Fxi = coordinateDerivative.*Fr;
G = (1-r.^2).*F(:,1:degree);
Gr = -2*r.*F(:,1:degree)+(1-r.^2).*Fr(:,1:degree);
Gxi = coordinateDerivative.*Gr;
H = [G (1-r)/2];
Hxi = [Gxi -coordinateDerivative/2];
[Pressure,PressureR] = legendreValuesForPrimitiveAudit(r,degree+1);
PressureXi = coordinateDerivative.*PressureR;
[Fendpoint,FendpointR] = legendreValuesForPrimitiveAudit([-1;1],degree);
spaces = completeSpaces(F,Fxi,G,Gxi,H,Hxi,Pressure,PressureXi, ...
    Fendpoint,endpointCoordinateDerivative.*FendpointR, ...
    [zeros(2,degree) [1;0]],D);
end

function spaces = chebyshevSpaces(degree,s,coordinateDerivative, ...
    endpointCoordinateDerivative,D)
[F,Fs] = chebyshevValues(s,degree);
Fxi = coordinateDerivative.*Fs;
G = (1-s.^2).*F(:,1:degree);
Gs = -2*s.*F(:,1:degree)+(1-s.^2).*Fs(:,1:degree);
Gxi = coordinateDerivative.*Gs;
H = [G (1-s)/2];
Hxi = [Gxi -coordinateDerivative/2];
[Pressure,PressureS] = chebyshevValues(s,degree+1);
PressureXi = coordinateDerivative.*PressureS;
[Fendpoint,FendpointS] = chebyshevValues([-1;1],degree);
spaces = completeSpaces(F,Fxi,G,Gxi,H,Hxi,Pressure,PressureXi, ...
    Fendpoint,endpointCoordinateDerivative.*FendpointS, ...
    [zeros(2,degree) [1;0]],D);
end

function spaces = completeSpaces(F,Fxi,G,Gxi,H,Hxi,Pressure, ...
    PressureXi,Fendpoint,FendpointXi,Hendpoint,D)
spaces = struct("F",F,"Fxi",Fxi,"G",G,"Gxi",Gxi, ...
    "H",H,"Hxi",Hxi,"Pressure",Pressure, ...
    "PressureXi",PressureXi,"Fendpoint",Fendpoint, ...
    "FendpointXi",FendpointXi,"Hendpoint",Hendpoint, ...
    "depth",D);
end

function diagnostics = coordinateDiagnostics(name,coordinate,D,conditionNumber, ...
    inverseResidual,isQualified,reason)
diagnostics = struct("name",name,"isQualified",isQualified, ...
    "reason",reason,"conditionNumber",conditionNumber, ...
    "forwardMapDefect",max(abs(inverseResidual),[],"all"), ...
    "inverseMapDefect",max(abs(inverseResidual),[],"all"), ...
    "jacobianDefect",0,"endpointDefect",0, ...
    "quadratureDefect",abs(sum(coordinate.weight)-D)/D, ...
    "integratedBuoyancyFrequency",NaN,"minimumN2",NaN, ...
    "maximumN2",NaN);
end

function [r,weight] = gaussLegendreRule(order)
index = (1:order-1)';
offDiagonal = index./sqrt(4*index.^2-1);
[vectors,values] = eig(diag(offDiagonal,1)+diag(offDiagonal,-1),"vector");
[r,permutation] = sort(values);
vectors = vectors(:,permutation);
weight = 2*(vectors(1,:)').^2;
end

function [T,Ts] = chebyshevValues(s,degree)
T = zeros(numel(s),degree+1);
Ts = zeros(numel(s),degree+1);
T(:,1) = 1;
if degree == 0
    return
end
T(:,2) = s;
Ts(:,2) = 1;
for n = 2:degree
    T(:,n+1) = 2*s.*T(:,n)-T(:,n-1);
    Ts(:,n+1) = 2*T(:,n)+2*s.*Ts(:,n)-Ts(:,n-1);
end
end
