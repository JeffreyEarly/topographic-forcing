function audit = buildBoundaryDynamicalDescriptorAudit(problem,bottomSlope)
% Build and audit the primitive F--G boundary descriptor.

if nargin < 2
    bottomSlope = [0 0];
end
bottomSlope = reshape(bottomSlope,1,2);
h = problem.topographicHeight;
terrainRange = max(h,[],"all")-min(h,[],"all");
terrainScale = max(1,max(abs(h),[],"all"));
if terrainRange > 20*eps(terrainScale)
    error("WVTerrainEnergyGalerkin:BoundaryDescriptorRequiresUniformDepth", ...
        "Milestone 5 supports only flat or spatially uniform terrain.")
end

gamma = problem.gamma(1);
nHorizontal = height(problem.horizontalLayout);
blocks = cell(nHorizontal,1);
maximum = struct("descriptorResidual",0,"continuityResidual",0, ...
    "bottomResidual",0,"surfaceResidual",0,"pressureGaugeResidual",0, ...
    "publicMapResidual",0,"bottomCoordinateMapResidual",0,"frequencyResidual",0);
for iK = 1:nHorizontal
    blocks{iK} = buildOneBlock(problem,iK,gamma,bottomSlope);
    names = string(fieldnames(maximum));
    for name = names.'
        maximum.(name) = max(maximum.(name),blocks{iK}.diagnostics.(name));
    end
end

residualNames = ["descriptorResidual","continuityResidual","bottomResidual", ...
    "surfaceResidual","pressureGaugeResidual","publicMapResidual", ...
    "bottomCoordinateMapResidual"];
failedConditions = residualNames(arrayfun(@(name)maximum.(name) > 1e-12,residualNames)).';
if maximum.frequencyResidual > 1e-11
    failedConditions(end+1,1) = "frequencyResidual";
end
structuralFailures = strings(6*nHorizontal,1);
nStructuralFailures = 0;
for iK = 1:nHorizontal
    diagnostics = blocks{iK}.diagnostics;
    if diagnostics.numberOfFiniteModes ~= diagnostics.expectedFiniteDimension
        nStructuralFailures = nStructuralFailures+1;
        structuralFailures(nStructuralFailures) = "finiteModeDimension";
    end
    if all(bottomSlope == 0) && diagnostics.stationaryDimension ~= diagnostics.expectedStationaryDimension
        nStructuralFailures = nStructuralFailures+1;
        structuralFailures(nStructuralFailures) = "stationaryDimension";
    end
    if diagnostics.descriptorRank ~= diagnostics.expectedDescriptorRank
        nStructuralFailures = nStructuralFailures+1;
        structuralFailures(nStructuralFailures) = "descriptorRank";
    end
    if diagnostics.constraintRank ~= diagnostics.expectedConstraintRank
        nStructuralFailures = nStructuralFailures+1;
        structuralFailures(nStructuralFailures) = "constraintRank";
    end
    if diagnostics.pressureNullity ~= 0
        nStructuralFailures = nStructuralFailures+1;
        structuralFailures(nStructuralFailures) = "pressureNullspace";
    end
    if diagnostics.publicCoordinateRank ~= diagnostics.expectedFiniteDimension
        nStructuralFailures = nStructuralFailures+1;
        structuralFailures(nStructuralFailures) = "publicCoordinateRank";
    end
end
failedConditions = [failedConditions;structuralFailures(1:nStructuralFailures)];
failedConditions = unique(failedConditions,"stable");

audit = struct;
audit.status = "compatible";
audit.isCompatible = isempty(failedConditions);
if ~audit.isCompatible
    audit.status = "incompatible";
end
audit.failedConditions = failedConditions;
audit.gamma = gamma;
audit.bottomSlope = bottomSlope;
audit.blocks = blocks;
audit.diagnostics = struct("maximum",maximum, ...
    "numberOfHorizontalBlocks",nHorizontal, ...
    "numberOfFiniteModes",sum(cellfun(@(block)numel(block.frequency),blocks)), ...
    "numberOfInfiniteModes",sum(cellfun(@(block)block.diagnostics.numberOfInfiniteModes,blocks)));
end

function block = buildOneBlock(problem,iK,gamma,bottomSlope)
wvt = problem.originatingTransform;
hydro = problem.hydrostaticTransform;
zWeight = wvt.z_int(:);
W = diag(zWeight);
N2 = wvt.N2Function(gamma*wvt.z(:));
if isscalar(N2)
    N2 = repmat(N2,wvt.Nz,1);
end
N2 = N2(:);
xi = wvt.z(:);
D = wvt.Lz;
sx = bottomSlope(1);
sy = bottomSlope(2);

jF = problem.verticalModeIndices;
jG = jF(jF > 0);
iF = arrayfun(@(j)find(wvt.j == j,1),jF);
iG = arrayfun(@(j)find(wvt.j == j,1),jG);
F = hydro.FinvMatrix(:,iF);
G = hydro.GinvMatrix(:,iG);
DzF = -(hydro.N2/hydro.g).*hydro.QG0inv*(squeeze(hydro.Q0./hydro.P0).*hydro.PF0);
DzG = hydro.PF0inv*(squeeze(hydro.P0./(hydro.Q0.*hydro.h_0)).*hydro.QG0);
Fxi = DzF*F;
Gxi = DzG*G;

nF = size(F,2);
nG = size(G,2);
iu = 1:nF;
iv = nF+(1:nF);
iw = 2*nF+(1:nG);
ieta = 2*nF+nG+(1:nG);
ib = 2*nF+2*nG+1;
ip = ib+(1:nF);
nUnknown = ip(end);

k = problem.horizontalLayout.k(iK);
l = problem.horizontalLayout.l(iK);
profile = problem.bottomInversionProfiles{iK};
uBottom = -1i*l*profile.psi;
vBottom = 1i*k*profile.psi;
etaBottom = profile.eta;
pBottom = profile.pressure;
pBottomXi = wvt.rho0*wvt.f*profile.psiXi;

Ru = zeros(wvt.Nz,nUnknown);
Rv = zeros(wvt.Nz,nUnknown);
Rw = zeros(wvt.Nz,nUnknown);
Reta = zeros(wvt.Nz,nUnknown);
Rp = zeros(wvt.Nz,nUnknown);
RpXi = zeros(wvt.Nz,nUnknown);
Ru(:,iu) = F;
Rv(:,iv) = F;
Rw(:,iw) = G;
Reta(:,ieta) = G;
Ru(:,ib) = uBottom;
Rv(:,ib) = vBottom;
Reta(:,ib) = etaBottom;
Rp(:,ib) = pBottom;
Rp(:,ip) = F;
RpXi(:,ib) = pBottomXi;
RpXi(:,ip) = Fxi;

nRows = nUnknown;
M = zeros(nRows);
K = zeros(nRows);
ru = iu;
rv = iv;
rw = iw;
reta = ieta;
rb = ib;
rc = ip;

M(ru,:) = F'*W*Ru;
M(rv,:) = F'*W*Rv;
RuPhysical = Ru/gamma;
RvPhysical = Rv/gamma;
RwPhysical = Rw-(xi/(D*gamma)).*(sx*Ru+sy*Rv);
M(rw,:) = G'*W*RwPhysical;
M(reta,:) = G'*W*diag(N2)*Reta;
M(rb,ib) = 1;

pressureX = 1i*k*Rp+(xi*sx/(D*gamma)).*RpXi;
pressureY = 1i*l*Rp+(xi*sy/(D*gamma)).*RpXi;
K(ru,:) = wvt.f*(F'*W*Rv)-(gamma/wvt.rho0)*(F'*W*pressureX);
K(rv,:) = -wvt.f*(F'*W*Ru)-(gamma/wvt.rho0)*(F'*W*pressureY);
K(rw,:) = -(1/(wvt.rho0*gamma))*(G'*W*RpXi)-G'*W*diag(N2)*Reta;
K(reta,:) = G'*W*diag(N2)*RwPhysical;
K(rb,:) = (sx*Ru(1,:)+sy*Rv(1,:))/gamma;

continuity = F'*W*(1i*k*Ru+1i*l*Rv+Gxi*selector(nUnknown,iw)');
K(rc,:) = continuity;
if hypot(k,l) == 0
    iMean = find(jF == 0,1);
    gauge = zWeight'*Rp;
    K(rc(iMean),:) = gauge;
end

variable = [repmat("uF",nF,1);repmat("vF",nF,1);repmat("wG",nG,1); ...
    repmat("etaG",nG,1);"etaB";repmat("pF",nF,1)];
verticalMode = [jF(:);jF(:);jG(:);jG(:);NaN;jF(:)];
isPressure = variable == "pF";
layout = table((1:nUnknown)',variable,verticalMode,isPressure, ...
    'VariableNames',["index","variable","j","isPressure"]);

continuityRows = rc;
pressureGaugeRow = [];
if hypot(k,l) == 0
    pressureGaugeRow = rc(find(jF == 0,1));
    continuityRows(continuityRows == pressureGaugeRow) = [];
end
publicRows = find(problem.stateLayout.horizontalIndex == iK);
publicBasis = problem.basisBlocks{iK};
sqrtN2 = sqrt(N2);
publicReconstruction = [publicBasis.uHat;publicBasis.vHat;publicBasis.wHat;sqrtN2.*publicBasis.etaHat];
[V,lambda,pencilScaling] = finiteDescriptorModes(K,M,problem.flatModeBlocks{iK}.frequency);
frequency = real(1i*lambda);
[frequency,order] = sort(frequency);
lambda = lambda(order);
V = V(:,order);

descriptorResidual = norm(K*V-M*(V.*lambda.'),"fro")/max(norm(K*V,"fro")+norm(M*(V.*lambda.'),"fro"),realmin);
pressureGaugeResidual = 0;
if ~isempty(pressureGaugeRow)
    pressureGaugeResidual = norm(K(pressureGaugeRow,:)*V)/max(norm(K(pressureGaugeRow,:))*norm(V,"fro"),realmin);
end
continuityResidual = norm(K(continuityRows,:)*V,"fro")/max(norm(K(continuityRows,:),"fro")*norm(V,"fro"),realmin);
bottomResidual = norm(V(ib,:).*lambda.'-K(rb,:)*V,"fro") ...
    /max(norm(lambda)*norm(V(ib,:))+norm(K(rb,:))*norm(V,"fro"),realmin);
surfaceResidual = norm(Rw(end,:)*V,"fro")/max(norm(Rw,"fro")*norm(V,"fro"),realmin);

descriptorReconstruction = [Ru;Rv;Rw;sqrtN2.*Reta];
publicCoordinates = publicReconstruction\(descriptorReconstruction*V);
publicMapResidual = norm(publicReconstruction*publicCoordinates-descriptorReconstruction*V,"fro") ...
    /max(norm(descriptorReconstruction*V,"fro"),realmin);
publicRankTolerance = max(size(publicCoordinates))*eps(norm(publicCoordinates,2));
publicCoordinateRank = rank(publicCoordinates,publicRankTolerance);
stationaryTolerance = 1e3*eps*max(1,max(abs(frequency)));
stationary = abs(frequency) <= stationaryTolerance;
iPublicBottom = find(problem.stateLayout.component(publicRows) == "etaB",1);
bottomCoordinateMapResidual = norm(publicCoordinates(iPublicBottom,:)-V(ib,:))/max(norm(V(ib,:)),realmin);
L = publicCoordinates*diag(lambda)/publicCoordinates;
localForms = constructLocalForms(problem,iK,gamma,bottomSlope,L,Rp(:,ip),RpXi(:,ip));

expectedFrequency = sort(problem.flatModeBlocks{iK}.frequency);
if all(bottomSlope == 0) && gamma == 1 && numel(expectedFrequency) == numel(frequency)
    frequencyResidual = norm(frequency-expectedFrequency)/max(norm(expectedFrequency),1);
elseif all(bottomSlope == 0) && gamma ~= 1
    frequencyResidual = uniformDepthFrequencyResidual(problem,iK,frequency,gamma,jF);
else
    frequencyResidual = 0;
end

rankTolerance = max(size(K))*eps(max([norm(K,2),norm(M,2),1]));
constraintRank = rank(K(rc,:),rankTolerance);
expectedFiniteDimension = 3*nF-1;
if hypot(k,l) == 0
    expectedFiniteDimension = 3*nF;
end

block = struct;
block.horizontalIndex = iK;
block.kMode = problem.horizontalLayout.kMode(iK);
block.lMode = problem.horizontalLayout.lMode(iK);
block.layout = layout;
block.massMatrix = M;
block.operatorMatrix = K;
block.reconstruction = struct("uHat",Ru,"vHat",Rv,"wHat",Rw,"etaHat",Reta,"pressure",Rp);
block.physicalReconstruction = struct("u",RuPhysical,"v",RvPhysical,"w",RwPhysical);
block.continuityRows = rc;
block.bottomRow = rb;
block.pressureColumns = ip;
block.eigenvalues = lambda;
block.frequency = frequency;
block.eigenvectors = V;
block.pencilScaling = pencilScaling;
block.publicRows = publicRows;
block.publicCoordinates = publicCoordinates;
block.publicGenerator = L;
block.localForms = localForms;
block.diagnostics = struct("descriptorResidual",descriptorResidual, ...
    "continuityResidual",continuityResidual,"bottomResidual",bottomResidual, ...
    "surfaceResidual",surfaceResidual,"pressureGaugeResidual",pressureGaugeResidual, ...
    "publicMapResidual",publicMapResidual, ...
    "bottomCoordinateMapResidual",bottomCoordinateMapResidual, ...
    "frequencyResidual",frequencyResidual,"descriptorRank",rank([K M],rankTolerance), ...
    "constraintRank",constraintRank,"pressureNullity",nullity(K(:,ip),rankTolerance), ...
    "stationaryDimension",nnz(stationary), ...
    "expectedStationaryDimension",nF+double(hypot(k,l) > 0), ...
    "publicCoordinateRank",publicCoordinateRank, ...
    "expectedDescriptorRank",nUnknown, ...
    "expectedConstraintRank",nF, ...
    "expectedFiniteDimension",expectedFiniteDimension, ...
    "numberOfFiniteModes",numel(frequency), ...
    "numberOfInfiniteModes",nUnknown-numel(frequency), ...
    "rankTolerance",rankTolerance);
end

function forms = constructLocalForms(problem,iK,gamma,bottomSlope,L,Rp,RpXi)
wvt = problem.originatingTransform;
basis = problem.basisBlocks{iK};
xi = wvt.z(:);
D = wvt.Lz;
W = diag(wvt.z_int(:));
N2 = wvt.N2Function(gamma*xi);
if isscalar(N2)
    N2 = repmat(N2,wvt.Nz,1);
end
N2 = N2(:);
WN = W*diag(N2);
sx = bottomSlope(1);
sy = bottomSlope(2);
k = problem.horizontalLayout.k(iK);
l = problem.horizontalLayout.l(iK);

U = basis.uHat/gamma;
V = basis.vHat/gamma;
UXi = basis.uHatXi/gamma;
VXi = basis.vHatXi/gamma;
WHat = basis.wHat;
Eta = basis.etaHat;
EtaXi = basis.etaHatXi;
WPhysical = WHat-(xi/(D*gamma)).*(sx*basis.uHat+sy*basis.vHat);
WMetric = WPhysical-WHat;

if gamma == 1
    flat = problem.flatModeBlocks{iK};
    E = flat.E+wvt.rho0*(WPhysical'*W*WPhysical-WHat'*W*WHat);
    J = flat.J+wvt.rho0*(Eta'*WN*WMetric-WMetric'*WN*Eta);
    Q = flat.Q+(sx/D)*V-(sy/D)*U+(xi*sx/D).*VXi-(xi*sy/D).*UXi;
else
    E = wvt.rho0*gamma*(U'*W*U+V'*W*V+WPhysical'*W*WPhysical+Eta'*WN*Eta);
    J = wvt.rho0*gamma*(wvt.f*(U'*W*V-V'*W*U) ...
        +Eta'*WN*WPhysical-WPhysical'*WN*Eta);
    Q = 1i*k*V-1i*l*U+(sx/(D*gamma))*V-(sy/(D*gamma))*U ...
        +(xi*sx/(D*gamma)).*VXi-(xi*sy/(D*gamma)).*UXi ...
        -(wvt.f/gamma)*EtaXi;
end
if gamma == 1 && all(bottomSlope == 0)
    L = E\J;
end
Z = Q'*(gamma*W)*Q;

B = zeros(1,size(L,1));
iBottom = find(problem.stateLayout.component(problem.stateLayout.horizontalIndex == iK) == "etaB",1);
B(iBottom) = 1;
R = (sx*basis.uHat(1,:)+sy*basis.vHat(1,:))/gamma;
pressureX = 1i*k*Rp+(xi*sx/(D*gamma)).*RpXi;
pressureY = 1i*l*Rp+(xi*sy/(D*gamma)).*RpXi;
pressureWorkX = gamma*U'*W*pressureX;
pressureWorkY = gamma*V'*W*pressureY;
pressureWorkZ = WPhysical'*W*RpXi;
pressureWork = -(pressureWorkX+pressureWorkY+pressureWorkZ);
pressureWorkScale = norm(pressureWorkX,"fro")+norm(pressureWorkY,"fro")+norm(pressureWorkZ,"fro");
metricPressureX = gamma*U'*W*((xi*sx/(D*gamma)).*RpXi);
metricPressureY = gamma*V'*W*((xi*sy/(D*gamma)).*RpXi);
metricPressureZ = WMetric'*W*RpXi;
metricPressureScale = norm(metricPressureX,"fro")+norm(metricPressureY,"fro")+norm(metricPressureZ,"fro");
metricPressureDefect = norm(metricPressureX+metricPressureY+metricPressureZ,"fro")/max(metricPressureScale,realmin);
continuity = 1i*k*basis.uHat+1i*l*basis.vHat+basis.wHatXi;
continuityDefect = norm(continuity,"fro")/max(norm(1i*k*basis.uHat,"fro") ...
    +norm(1i*l*basis.vHat,"fro")+norm(basis.wHatXi,"fro"),realmin);
boundaryDefect = norm([basis.wHat(1,:);basis.wHat(end,:)],"fro")/max(norm(basis.wHat,"fro"),realmin);

forms = struct;
forms.energyMatrix = E;
forms.exchangeMatrix = J;
forms.apvMatrix = Q;
forms.potentialEnstrophyMatrix = Z;
forms.generator = L;
forms.bottomValue = B;
forms.bottomTendency = R;
forms.pressureWork = pressureWork;
forms.diagnostics = struct( ...
    "energyHermitianDefect",relativeNorm(E-E',E), ...
    "exchangeSkewDefect",relativeNorm(J+J',J), ...
    "weakEvolutionDefect",norm(E*L-J,"fro")/max(norm(E*L,"fro")+norm(J,"fro"),realmin), ...
    "energySkewDefect",norm(L'*E+E*L,"fro")/max(norm(E,"fro")*norm(L,"fro"),realmin), ...
    "apvTendencyDefect",norm(Q*L,"fro")/max(norm(Q,"fro")*norm(L,"fro"),realmin), ...
    "potentialEnstrophyDefect",norm(L'*Z+Z*L,"fro")/max(norm(Z,"fro")*norm(L,"fro"),realmin), ...
    "bottomEvolutionDefect",norm(B*L-R,"fro")/max(norm(B,"fro")*norm(L,"fro")+norm(R,"fro"),realmin), ...
    "pressureGreenDefect",max([metricPressureDefect,continuityDefect,boundaryDefect]), ...
    "directPressureQuadratureDefect",norm(pressureWork,"fro")/max(pressureWorkScale,realmin), ...
    "minimumEnergyEigenvalue",min(real(eig((E+E')/2))), ...
    "energyConditionNumber",cond((E+E')/2), ...
    "maximumRealGrowthRate",max(abs(real(eig(L)))));
end

function value = relativeNorm(numerator,denominator)
value = norm(numerator,"fro")/max(norm(denominator,"fro"),realmin);
end

function S = selector(n,indices)
S = zeros(n,numel(indices));
S(indices,:) = eye(numel(indices));
end

function [V,lambda,scaling] = finiteDescriptorModes(K,M,referenceFrequency)
omegaScale = max(abs(referenceFrequency),[],"all");
omegaScale = max(omegaScale,eps);
Ks = K/omegaScale;
Ms = M;
columnNorm = sqrt(sum(abs(Ks).^2+abs(Ms).^2,1));
columnScale = 1./max(columnNorm,sqrt(realmin));
Ks = Ks.*columnScale;
Ms = Ms.*columnScale;
rowNorm = sqrt(sum(abs(Ks).^2+abs(Ms).^2,2));
rowScale = 1./max(rowNorm,sqrt(realmin));
Ks = rowScale.*Ks;
Ms = rowScale.*Ms;
[scaledVectors,scaledEigenvalues] = eig(Ks,Ms,"vector");
finite = isfinite(scaledEigenvalues) & ~isnan(scaledEigenvalues) & abs(scaledEigenvalues) < 100;
V = columnScale.'.*scaledVectors(:,finite);
lambda = omegaScale*scaledEigenvalues(finite);
V = V./max(vecnorm(V),sqrt(realmin));
scaling = struct("frequency",omegaScale,"row",rowScale,"column",columnScale.');
end

function value = nullity(A,tolerance)
value = size(A,2)-rank(A,tolerance);
end

function residual = uniformDepthFrequencyResidual(problem,iK,frequency,gamma,jF)
wvt = problem.originatingTransform;
kappa = hypot(problem.horizontalLayout.k(iK),problem.horizontalLayout.l(iK));
if kappa == 0
    residual = 0;
    return
end
N2 = wvt.N2Function(0);
if ~isscalar(N2)
    residual = NaN;
    return
end
internalJ = jF(jF > 0);
m = internalJ*pi/(gamma*wvt.Lz);
omega = sqrt((N2*kappa^2+wvt.f^2*m.^2)./(kappa^2+m.^2));
expected = sort([-omega(:);zeros(numel(jF)+1,1);omega(:)]);
if numel(expected) ~= numel(frequency)
    residual = Inf;
else
    residual = norm(frequency-expected)/max(norm(expected),1);
end
end
