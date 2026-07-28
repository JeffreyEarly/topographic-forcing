function audit = buildCoupledPVFrequencyOracleAudit(problem,iK,volumePVGradient,bottomSlope,degree,quadratureOrder)
% Build the single-wavenumber coupled volume--boundary PV oracle.

wvt = problem.originatingTransform;
k = problem.horizontalLayout.k(iK);
l = problem.horizontalLayout.l(iK);
kappa = hypot(k,l);
betaK = k*volumePVGradient(2)-l*volumePVGradient(1);
sK = wvt.f*(k*bottomSlope(2)-l*bottomSlope(1));
quadrature = legendreQuadrature(quadratureOrder,wvt.Lz);

direct = directOracle(wvt,kappa,betaK,sK,degree,quadrature);
if betaK ~= 0 && sK ~= 0
    endpoint = endpointOracle(wvt,kappa,betaK,sK,degree,quadrature,direct);
else
    endpoint = unavailableEndpoint(degree);
end
basis = balancedBasisProjection(problem,iK,direct);
if endpoint.isAvailable
    resolved = endpoint.resolvedIndices;
else
    [~,frequencyOrder] = sort(abs(direct.frequency),"descend");
    resolved = frequencyOrder(1:min(4,numel(frequencyOrder)));
end
basis.resolvedIndices = resolved;
basis.diagnostics.resolvedReconstructionDefect = max(basis.diagnostics.reconstructionDefect(resolved));

qgPasses = direct.diagnostics.inversionDefect <= 1e-12 ...
    && direct.diagnostics.greenIdentityDefect <= 1e-12 ...
    && direct.diagnostics.energyDefect <= 1e-12 ...
    && direct.diagnostics.frequencyImaginaryDefect <= 1e-11 ...
    && (~direct.pseudoenstrophy.isAvailable || direct.diagnostics.pseudoenstrophyDefect <= 1e-12);
endpointPasses = ~endpoint.isAvailable ...
    || (endpoint.diagnostics.resolvedFrequencyDefect <= 1e-10 ...
    && endpoint.diagnostics.resolvedEigenfunctionDefect <= 1e-9);
basisPasses = basis.diagnostics.rank == basis.diagnostics.expectedRank ...
    && all(isfinite(basis.diagnostics.reconstructionDefect));

if ~qgPasses
    status = "qg-formulation-blocker";
    diagnosis = "The direct coupled volume--boundary inversion does not satisfy its QG energy or Green identities.";
elseif ~endpointPasses
    status = "endpoint-equivalence-blocker";
    diagnosis = "The direct frequency problem and independent Yassin endpoint problem do not agree.";
elseif ~basisPasses
    status = "basis-representation-blocker";
    diagnosis = "The mathematical oracle passes, but its geostrophic fields cannot be represented in the current balanced-plus-bottom coordinates.";
else
    status = "passed";
    diagnosis = "The local coupled-PV frequency problem, Yassin endpoint problem, and current boundary-complete basis are mutually consistent at this truncation.";
end

audit = struct();
audit.scope = "single-horizontal-mode-coupled-pv-frequency-oracle";
audit.status = status;
audit.isCompatible = status == "passed";
audit.diagnosis = diagnosis;
audit.horizontalIndex = iK;
audit.horizontalMode = [problem.horizontalLayout.kMode(iK) problem.horizontalLayout.lMode(iK)];
audit.horizontalWavenumber = [k l];
audit.kappa = kappa;
audit.volumePVGradient = volumePVGradient;
audit.bottomSlope = bottomSlope;
audit.projectedVolumePVGradient = betaK;
audit.projectedBottomPVGradient = sK;
audit.polynomialDegree = degree;
audit.quadratureOrder = quadratureOrder;
audit.direct = direct;
audit.endpoint = endpoint;
audit.basis = basis;
audit.requiredTolerance = struct("structure",1e-12,"frequency",1e-10, ...
    "eigenfunction",1e-9,"frequencyImaginary",1e-11);
end

function direct = directOracle(wvt,kappa,betaK,sK,degree,quadrature)
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
Tendency = [betaK*volumeProjection;sK*psiBottom];
L = -1i*Tendency*Gcoefficient;
E = Gcoefficient'*A*Gcoefficient;
E = (E+E')/2;

[V,lambda] = eig(L,"vector");
omegaComplex = 1i*lambda;
[~,order] = sort(real(omegaComplex));
V = V(:,order);
omegaComplex = omegaComplex(order);
for iMode = 1:size(V,2)
    V(:,iMode) = V(:,iMode)/sqrt(real(V(:,iMode)'*E*V(:,iMode)));
end
omega = real(omegaComplex);
psi = G*V;
psiZ = Gz*V;

inversionDefect = norm(A*Gcoefficient+B,"fro")/max(norm(B,"fro"),realmin);
greenIdentityDefect = norm(E+B'*Gcoefficient,"fro")/max(norm(E,"fro"),realmin);
energyDefect = norm(L'*E+E*L,"fro")/max(norm(E*L,"fro")+norm(L'*E,"fro"),realmin);
frequencyImaginaryDefect = max(abs(imag(omegaComplex)))/max(max(abs(real(omegaComplex))),realmin);
eigenResidual = norm(L*V+1i*V.*omega.',"fro")/max(norm(L*V,"fro")+norm(V.*omega.',"fro"),realmin);
topFluxValue = topFlux*Tpsi*Gcoefficient*V;
topFluxDefect = norm(topFluxValue)/max(norm(V,"fro"),realmin);
bottomFluxValue = bottomFlux*Tpsi*Gcoefficient*V;
bottomFluxDefect = norm(bottomFluxValue-V(end,:))/max(norm(V(end,:)),realmin);

pseudo = pseudoenstrophy(Mq,betaK,sK,L);
qNorm = sqrt(max(real(sum(conj(V(1:end-1,:)).*(Mq*V(1:end-1,:)),1)),0));
stateNorm = sqrt(sum(abs(V).^2,1));
bottomParticipation = abs(V(end,:))./max(stateNorm,realmin);

direct = struct();
direct.volumePVGradient = betaK;
direct.bottomPVGradient = sK;
direct.volumeMassMatrix = Mq;
direct.inversionStiffnessMatrix = A;
direct.inversionRightHandSide = B;
direct.inversionMap = Gcoefficient;
direct.generator = L;
direct.energyMatrix = E;
direct.stateEigenvectors = V;
direct.frequency = omega;
direct.frequencyComplex = omegaComplex;
direct.psiCoefficient = Gcoefficient*V;
direct.psi = psi;
direct.psiZ = psiZ;
direct.qNorm = qNorm;
direct.bottomParticipation = bottomParticipation;
direct.pseudoenstrophy = pseudo;
direct.polynomial = struct("psiTransformation",Tpsi,"psiValues",Ppsi, ...
    "psiDerivative",PpsiZ,"qValues",Pq,"bottomEvaluation",psiBottom);
direct.diagnostics = struct("inversionDefect",inversionDefect, ...
    "greenIdentityDefect",greenIdentityDefect,"energyDefect",energyDefect, ...
    "pseudoenstrophyDefect",pseudo.conservationDefect, ...
    "frequencyImaginaryDefect",frequencyImaginaryDefect, ...
    "eigenResidual",eigenResidual,"surfaceFluxDefect",topFluxDefect, ...
    "bottomFluxDefect",bottomFluxDefect);
end

function pseudo = pseudoenstrophy(Mq,betaK,sK,L)
if betaK == 0 || sK == 0
    pseudo = struct("isAvailable",false,"classification","singular", ...
        "matrix",[],"inertia",[NaN NaN NaN],"conservationDefect",NaN);
    return
end
Z = blkdiag(Mq/betaK,1/sK);
Z = (Z+Z')/2;
eigenvalue = eig(Z);
tolerance = 100*eps*max(max(abs(eigenvalue)),realmin);
inertia = [nnz(eigenvalue > tolerance) nnz(eigenvalue < -tolerance) nnz(abs(eigenvalue) <= tolerance)];
if inertia(1) == 0 || inertia(2) == 0
    classification = "definite";
else
    classification = "pontryagin";
end
defect = norm(L'*Z+Z*L,"fro")/max(norm(Z*L,"fro")+norm(L'*Z,"fro"),realmin);
pseudo = struct("isAvailable",true,"classification",classification, ...
    "matrix",Z,"inertia",inertia,"conservationDefect",defect);
end

function endpoint = endpointOracle(wvt,kappa,betaK,sK,degree,quadrature,direct)
[x,Dx] = chebyshevLobatto(degree);
z = wvt.Lz*(x-1)/2;
Dz = (2/wvt.Lz)*Dx;
N2 = wvt.N2Function(z);
N2 = N2(:);
sigma = wvt.f^2./N2;
operator = Dz*diag(sigma)*Dz;
Gamma = sK/betaK;
n = degree+1;
A = zeros(n);
B = zeros(n);
interior = 2:degree;
A(interior,:) = operator(interior,:);
B(interior,interior) = -eye(degree-1);
A(1,:) = sigma(1)*Dz(1,:);
bottom = zeros(1,n);
bottom(end) = 1;
A(end,:) = sigma(end)*Dz(end,:)+Gamma*kappa^2*bottom;
B(end,:) = -Gamma*bottom;

[V,lambda] = eig(A,B,"vector");
finite = isfinite(lambda) & abs(lambda) < 1/eps;
V = V(:,finite);
lambda = lambda(finite);
omegaComplex = -betaK./(lambda+kappa^2);

[Tnode,~] = chebyshevValues(x,degree);
coefficient = Tnode\V;
[Tq,Tqr] = chebyshevValues(quadrature.r,degree);
psi = Tq*coefficient;
psiZ = (2/wvt.Lz)*Tqr*coefficient;
W = diag(quadrature.weight);
N2q = wvt.N2Function(quadrature.xi);
sigmaq = wvt.f^2./N2q(:);
for iMode = 1:size(psi,2)
    normMode = sqrt(real(psi(:,iMode)'*W*(kappa^2*psi(:,iMode)) ...
        +psiZ(:,iMode)'*W*(sigmaq.*psiZ(:,iMode))));
    psi(:,iMode) = psi(:,iMode)/normMode;
    psiZ(:,iMode) = psiZ(:,iMode)/normMode;
    V(:,iMode) = V(:,iMode)/normMode;
end

[omegaEndpoint,orderEndpoint] = sort(real(omegaComplex));
omegaComplex = omegaComplex(orderEndpoint);
omegaEndpoint = omegaEndpoint(:);
psi = psi(:,orderEndpoint);
psiZ = psiZ(:,orderEndpoint);
V = V(:,orderEndpoint);
lambda = lambda(orderEndpoint);
[omegaDirect,orderDirect] = sort(direct.frequency);
psiDirect = direct.psi(:,orderDirect);
psiZDirect = direct.psiZ(:,orderDirect);
nMatch = min(numel(omegaDirect),numel(omegaEndpoint));
omegaDirect = omegaDirect(1:nMatch);
omegaEndpoint = omegaEndpoint(1:nMatch);
psiDirect = psiDirect(:,1:nMatch);
psiZDirect = psiZDirect(:,1:nMatch);
psi = psi(:,1:nMatch);
psiZ = psiZ(:,1:nMatch);

frequencyScale = max([abs(omegaDirect);abs(omegaEndpoint);realmin]);
frequencyDefect = abs(omegaDirect-omegaEndpoint)./max(abs(omegaDirect)+abs(omegaEndpoint),frequencyScale*eps);
overlap = zeros(nMatch,1);
for iMode = 1:nMatch
    cross = psiDirect(:,iMode)'*W*(kappa^2*psi(:,iMode)) ...
        +psiZDirect(:,iMode)'*W*(sigmaq.*psiZ(:,iMode));
    directNorm = sqrt(real(psiDirect(:,iMode)'*W*(kappa^2*psiDirect(:,iMode)) ...
        +psiZDirect(:,iMode)'*W*(sigmaq.*psiZDirect(:,iMode))));
    endpointNorm = sqrt(real(psi(:,iMode)'*W*(kappa^2*psi(:,iMode)) ...
        +psiZ(:,iMode)'*W*(sigmaq.*psiZ(:,iMode))));
    overlap(iMode) = abs(cross)/(directNorm*endpointNorm);
end
eigenfunctionDefect = max(0,1-overlap);
[~,resolvedOrder] = sort(abs(omegaDirect),"descend");
resolved = resolvedOrder(1:min(4,nMatch));

endpoint = struct();
endpoint.isAvailable = true;
endpoint.lambdaY = lambda;
endpoint.frequency = real(omegaComplex);
endpoint.frequencyComplex = omegaComplex;
endpoint.eigenvectors = V;
endpoint.psi = psi;
endpoint.psiZ = psiZ;
endpoint.directFrequency = omegaDirect;
endpoint.matchedFrequency = omegaEndpoint;
endpoint.frequencyDefect = frequencyDefect;
endpoint.energyOverlap = overlap;
endpoint.eigenfunctionDefect = eigenfunctionDefect;
endpoint.resolvedIndices = resolved;
endpoint.diagnostics = struct( ...
    "resolvedFrequencyDefect",max(frequencyDefect(resolved)), ...
    "resolvedEigenfunctionDefect",max(eigenfunctionDefect(resolved)), ...
    "frequencyImaginaryDefect",max(abs(imag(omegaComplex)))/max(max(abs(real(omegaComplex))),realmin), ...
    "numberOfFiniteModes",numel(lambda),"expectedNumberOfFiniteModes",degree);
end

function endpoint = unavailableEndpoint(degree)
endpoint = struct("isAvailable",false,"lambdaY",[],"frequency",[], ...
    "frequencyComplex",[],"eigenvectors",[],"psi",[],"psiZ",[], ...
    "directFrequency",[],"matchedFrequency",[],"frequencyDefect",[], ...
    "energyOverlap",[],"eigenfunctionDefect",[],"resolvedIndices",[], ...
    "diagnostics",struct("resolvedFrequencyDefect",NaN, ...
    "resolvedEigenfunctionDefect",NaN,"frequencyImaginaryDefect",NaN, ...
    "numberOfFiniteModes",0,"expectedNumberOfFiniteModes",degree));
end

function basis = balancedBasisProjection(problem,iK,direct)
wvt = problem.originatingTransform;
rows = find(problem.stateLayout.horizontalIndex == iK);
components = problem.stateLayout.component(rows);
iBalancedLocal = find(components == "A0");
iBottomLocal = find(components == "etaB",1);
selectedLocal = [iBalancedLocal;iBottomLocal];
block = problem.basisBlocks{iK};
U = block.uHat(:,selectedLocal);
V = block.vHat(:,selectedLocal);
Eta = block.etaHat(:,selectedLocal);
zWeight = wvt.z_int(:);
N2 = wvt.N2(:);
R = [sqrt(zWeight).*U;sqrt(zWeight).*V;sqrt(zWeight.*N2).*Eta];

rNative = 2*wvt.z(:)/wvt.Lz+1;
[Pnative,PrNative] = legendreValues(rNative,size(direct.polynomial.psiTransformation,1)-1);
psiBasis = Pnative*direct.polynomial.psiTransformation;
psiBasisZ = (2/wvt.Lz)*PrNative*direct.polynomial.psiTransformation;
psi = psiBasis*direct.psiCoefficient;
psiZ = psiBasisZ*direct.psiCoefficient;
k = problem.horizontalLayout.k(iK);
l = problem.horizontalLayout.l(iK);
u = -1i*l*psi;
v = 1i*k*psi;
eta = -(wvt.f./N2).*psiZ;
eta(1,:) = -direct.stateEigenvectors(end,:)/wvt.f;
target = [sqrt(zWeight).*u;sqrt(zWeight).*v;sqrt(zWeight.*N2).*eta];

bottomCoefficient = eta(1,:);
balancedCoefficient = R(:,1:end-1)\(target-R(:,end).*bottomCoefficient);
coefficient = [balancedCoefficient;bottomCoefficient];
reconstruction = R*coefficient;
defect = vecnorm(reconstruction-target)./max(vecnorm(target),realmin);
rankR = rank(R,max(size(R))*eps(max(norm(R,2),1)));
bottomValueDefect = abs(bottomCoefficient-eta(1,:))./max(abs(eta(1,:)),1);

basis = struct();
basis.selectedStateRows = rows(selectedLocal);
basis.coefficient = coefficient;
basis.reconstructionDefect = defect;
basis.bottomCoefficient = bottomCoefficient;
basis.targetBottomDisplacement = eta(1,:);
basis.diagnostics = struct("rank",rankR,"expectedRank",size(R,2), ...
    "reconstructionDefect",defect,"maximumReconstructionDefect",max(defect), ...
    "bottomValueDefect",bottomValueDefect, ...
    "maximumBottomValueDefect",max(bottomValueDefect));
end

function quadrature = legendreQuadrature(order,D)
index = (1:order-1)';
offDiagonal = index./sqrt(4*index.^2-1);
[vectors,values] = eig(diag(offDiagonal,1)+diag(offDiagonal,-1),"vector");
[r,permutation] = sort(values);
vectors = vectors(:,permutation);
weightR = 2*(vectors(1,:)').^2;
quadrature = struct("r",r,"xi",D*(r-1)/2,"weight",D*weightR/2);
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

function [x,D] = chebyshevLobatto(degree)
index = (0:degree)';
x = cos(pi*index/degree);
c = [2;ones(degree-1,1);2].*(-1).^index;
dX = x-x';
D = (c*(1./c)')./(dX+eye(degree+1));
D = D-diag(sum(D,2));
end

function [T,Tr] = chebyshevValues(r,degree)
T = zeros(numel(r),degree+1);
Tr = zeros(numel(r),degree+1);
T(:,1) = 1;
if degree == 0
    return
end
T(:,2) = r;
Tr(:,2) = 1;
for n = 2:degree
    T(:,n+1) = 2*r.*T(:,n)-T(:,n-1);
    Tr(:,n+1) = 2*T(:,n)+2*r.*Tr(:,n)-Tr(:,n-1);
end
end
