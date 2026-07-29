function block = buildCoupledPVInversionBlock(wvt,kappa,degree,quadrature)
% Build one self-adjoint volume--boundary PV inversion block.

W = diag(quadrature.weight);
N2 = wvt.N2Function(quadrature.xi);
N2 = N2(:);
if numel(N2) ~= numel(quadrature.xi) || ~isreal(N2) || any(~isfinite(N2)) || any(N2 <= 0)
    error("WVTerrainEnergyGalerkin:InvalidOracleStratification", ...
        "N2Function must return one positive finite value per quadrature point.")
end
sigma = wvt.f^2./N2;

[Pfull,Prfull] = legendreValues(quadrature.r,degree);
PfullZ = (2/wvt.Lz)*Prfull;
[Pendpoint,PrEndpoint] = legendreValues([-1;1],degree);
topFlux = (wvt.f^2/wvt.N2Function(0))*(2/wvt.Lz)*PrEndpoint(2,:);
bottomFlux = (wvt.f^2/wvt.N2Function(-wvt.Lz))*(2/wvt.Lz)*PrEndpoint(1,:);
Tpsi = deterministicColumns(null(topFlux));
Ppsi = Pfull*Tpsi;
PpsiZ = PfullZ*Tpsi;
psiBottom = Pendpoint(1,:)*Tpsi;
Pq = Pfull(:,1:degree-1);
Mq = Pq'*W*Pq;
C = Ppsi'*W*Pq;
A = PpsiZ'*W*(sigma.*PpsiZ)+kappa^2*(Ppsi'*W*Ppsi);
B = [C psiBottom'];
Gcoefficient = -A\B;
G = Ppsi*Gcoefficient;
Gz = PpsiZ*Gcoefficient;
volumeProjection = Mq\(Pq'*W*Ppsi);
E = Gcoefficient'*A*Gcoefficient;
E = (E+E')/2;

inversionDefect = norm(A*Gcoefficient+B,"fro")/max(norm(B,"fro"),realmin);
greenIdentityDefect = norm(E+B'*Gcoefficient,"fro")/max(norm(E,"fro"),realmin);
topFluxMap = topFlux*Tpsi*Gcoefficient;
bottomFluxMap = bottomFlux*Tpsi*Gcoefficient;
expectedBottomFluxMap = zeros(1,size(Gcoefficient,2));
expectedBottomFluxMap(end) = 1;

block = struct();
block.kappa = kappa;
block.numberOfVolumeCoefficients = size(Pq,2);
block.volumeMassMatrix = Mq;
block.inversionStiffnessMatrix = A;
block.inversionRightHandSide = B;
block.inversionMap = Gcoefficient;
block.energyMatrix = E;
block.stateToPsi = G;
block.stateToPsiXi = Gz;
block.stateToPsiBottom = psiBottom*Gcoefficient;
block.stateToSurfaceFlux = topFluxMap;
block.stateToBottomFlux = bottomFluxMap;
block.volumeProjection = volumeProjection;
block.polynomial = struct("psiTransformation",Tpsi,"psiValues",Ppsi, ...
    "psiDerivative",PpsiZ,"qValues",Pq,"bottomEvaluation",psiBottom);
block.diagnostics = struct("inversionDefect",inversionDefect, ...
    "greenIdentityDefect",greenIdentityDefect, ...
    "surfaceFluxDefect",norm(topFluxMap)/max(norm(Gcoefficient,"fro"),realmin), ...
    "bottomFluxDefect",norm(bottomFluxMap-expectedBottomFluxMap) ...
    /max(norm(expectedBottomFluxMap),realmin));
end

function [P,Pr] = legendreValues(r,degree)
P = zeros(numel(r),degree+1);
Pr = zeros(numel(r),degree+1);
P(:,1) = 1;
if degree == 0
    return
end
P(:,2) = r;
Pr(:,2) = 1;
for n = 2:degree
    P(:,n+1) = ((2*n-1)*r.*P(:,n)-(n-1)*P(:,n-1))/n;
    Pr(:,n+1) = ((2*n-1)*(P(:,n)+r.*Pr(:,n))-(n-1)*Pr(:,n-1))/n;
end
end

function V = deterministicColumns(V)
for column = 1:size(V,2)
    [~,row] = max(abs(V(:,column)));
    if V(row,column) ~= 0
        V(:,column) = V(:,column)/(V(row,column)/abs(V(row,column)));
    end
end
end
