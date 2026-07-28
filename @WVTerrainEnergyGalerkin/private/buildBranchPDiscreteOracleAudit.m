function audit = buildBranchPDiscreteOracleAudit(problem,bottomSlope,polynomialDegree,quadratureOrder)
% Build an independent polynomial Branch-P oracle and compare the F--G descriptor.

h = problem.topographicHeight;
terrainRange = max(h,[],"all")-min(h,[],"all");
terrainScale = max(1,max(abs(h),[],"all"));
if terrainRange > 20*eps(terrainScale)
    error("WVTerrainEnergyGalerkin:BranchPOracleRequiresUniformDepth", ...
        "The Branch-P discrete oracle supports only flat or spatially uniform terrain.")
end

gamma = problem.gamma(1);
quadrature = legendreQuadrature(quadratureOrder,problem.originatingTransform.Lz);
polynomial = polynomialSpaces(polynomialDegree,quadrature,problem.originatingTransform.Lz);
nHorizontal = height(problem.horizontalLayout);
polynomialBlocks = cell(nHorizontal,1);
for iK = 1:nHorizontal
    polynomialBlocks{iK} = polynomialBlock(problem,iK,gamma,bottomSlope,polynomial,quadrature);
end

fgSlope = buildBoundaryDynamicalDescriptorAudit(problem,bottomSlope);
fgFlat = buildBoundaryDynamicalDescriptorAudit(problem,[0 0]);
fgTangentStep = 1e-4;
if norm(bottomSlope) > 0
    fgDirection = bottomSlope/norm(bottomSlope);
    fgPlus = buildBoundaryDynamicalDescriptorAudit(problem,fgTangentStep*fgDirection);
    fgMinus = buildBoundaryDynamicalDescriptorAudit(problem,-fgTangentStep*fgDirection);
else
    fgPlus = fgFlat;
    fgMinus = fgFlat;
end
fgBlocks = cell(nHorizontal,1);
for iK = 1:nHorizontal
    fgBlocks{iK} = fgComparisonBlock(problem,iK,bottomSlope,fgTangentStep, ...
        fgFlat.blocks{iK},fgSlope.blocks{iK},fgPlus.blocks{iK},fgMinus.blocks{iK});
end

polynomialMaximum = maximumDiagnostics(polynomialBlocks);
polynomialGeometry = polynomialGeometryDiagnostics(polynomial,quadrature);
fgMaximum = maximumComparisonDiagnostics(fgBlocks);
polynomialPasses = all(structfun(@(value)value <= 1e-11,polynomialMaximum.required));
fgNativePasses = all(structfun(@(value)value <= 1e-11,fgMaximum.nativeRequired));
fgPublicPasses = all(structfun(@(value)value <= 1e-11,fgMaximum.publicRequired));

if ~polynomialPasses
    outcome = "scientific-blocker";
    diagnosis = "The independent primitive polynomial descriptor does not satisfy the Branch-P coefficient identities.";
    repair = "none";
elseif fgNativePasses && ~fgPublicPasses
    outcome = "public-layout-blocker";
    diagnosis = "The native F--G descriptor passes, but its public-coordinate generator does not.";
    repair = "direct-constrained-pullback";
elseif ~fgNativePasses
    outcome = "basis-blocker";
    diagnosis = "The primitive oracle passes, but the native F--G descriptor does not.";
    repair = "eta-only-bottom-coordinate";
else
    outcome = "repaired";
    diagnosis = "Both primitive and F--G descriptor identities pass without empirical correction.";
    repair = "direct-constrained-pullback";
end

audit = struct;
audit.status = outcome;
audit.isCompatible = polynomialPasses && fgNativePasses && fgPublicPasses;
audit.bottomSlope = bottomSlope;
audit.polynomialDegree = polynomialDegree;
audit.quadratureOrder = quadratureOrder;
audit.polynomialOracle = struct("quadrature",quadrature,"spaces",polynomial, ...
    "blocks",{polynomialBlocks},"maximum",polynomialMaximum, ...
    "geometryDiagnostics",polynomialGeometry);
audit.fgNativeDescriptor = struct("blocks",{fgBlocks},"maximum",fgMaximum.nativeRequired);
audit.fgPublicDescriptor = struct("blocks",{fgBlocks},"maximum",fgMaximum.publicRequired);
audit.comparison = struct("blocks",{fgBlocks},"maximum",fgMaximum);
audit.diagnosis = diagnosis;
audit.repairCandidate = repair;
audit.scope = "local-constant-slope-only";
end

function block = polynomialBlock(problem,iK,gamma,slope,spaces,quadrature)
k = problem.horizontalLayout.k(iK);
l = problem.horizontalLayout.l(iK);
directional = buildPolynomialDirection(problem,k,l,gamma,slope,spaces,quadrature);
flat = buildPolynomialDirection(problem,k,l,gamma,[0 0],spaces,quadrature);
xDirection = buildPolynomialDirection(problem,k,l,gamma,[1 0],spaces,quadrature);
yDirection = buildPolynomialDirection(problem,k,l,gamma,[0 1],spaces,quadrature);

L1 = slope(1)*xDirection.L1+slope(2)*yDirection.L1;
E1 = slope(1)*xDirection.E1+slope(2)*yDirection.E1;
J1 = slope(1)*xDirection.J1+slope(2)*yDirection.J1;
Q1 = slope(1)*xDirection.Q1+slope(2)*yDirection.Q1;
Z1 = slope(1)*xDirection.Z1+slope(2)*yDirection.Z1;
R1 = slope(1)*xDirection.R1+slope(2)*yDirection.R1;

weakScale = sumNorms(flat.E*L1,E1*flat.L,J1);
energyScale = sumNorms(L1'*flat.E,flat.E*L1,flat.L'*E1,E1*flat.L);
apvScale = sumNorms(flat.Q*L1,Q1*flat.L);
bottomScale = sumNorms(flat.B*L1,R1);
enstrophyScale = sumNorms(L1'*flat.Z,flat.Z*L1,flat.L'*Z1,Z1*flat.L);
diagnostics = struct;
diagnostics.weakEvolutionDefect = norm(flat.E*L1+E1*flat.L-J1,"fro")/weakScale;
diagnostics.energyDefect = norm(L1'*flat.E+flat.E*L1+flat.L'*E1+E1*flat.L,"fro")/energyScale;
diagnostics.apvDefect = norm(flat.Q*L1+Q1*flat.L,"fro")/apvScale;
diagnostics.bottomDefect = norm(flat.B*L1-R1,"fro")/bottomScale;
diagnostics.enstrophyDefect = norm(L1'*flat.Z+flat.Z*L1+flat.L'*Z1+Z1*flat.L,"fro")/enstrophyScale;
diagnostics.pressureWorkDefect = directional.pressureWorkDefect;
diagnostics.constraintTangencyDefect = directional.constraintTangencyDefect;
diagnostics.saddleResidual = directional.saddleResidual;
diagnostics.analyticDerivativeDefect = derivativeCheck(problem,k,l,gamma,spaces,quadrature,slope,L1);
diagnostics.finiteSlopeLinearizationDefect = norm(directional.L-flat.L-L1,"fro") ...
    /max(norm(directional.L-flat.L,"fro")+norm(L1,"fro"),realmin);
diagnostics.energyHermitianDefect = relativeNorm(flat.E-flat.E',flat.E);
diagnostics.exchangeSkewDefect = relativeNorm(flat.J+flat.J',flat.J);
diagnostics.minimumEnergyEigenvalue = min(real(eig((flat.E+flat.E')/2)));
diagnostics.constraintRank = flat.constraintRank;
diagnostics.admissibleDimension = size(flat.N,2);
diagnostics.pressureDimension = flat.pressureDimension;
diagnostics.flatModeOneDispersionDefect = flatModeOneDispersionDefect(problem,k,l,gamma,flat.L);
diagnostics.flat = exactSystemDefects(flat);

block = directional;
block.horizontalIndex = iK;
block.kMode = problem.horizontalLayout.kMode(iK);
block.lMode = problem.horizontalLayout.lMode(iK);
block.flat = flat;
block.firstOrder = struct("L",L1,"E",E1,"J",J1,"Q",Q1,"Z",Z1,"B",flat.B,"R",R1);
block.diagnostics = diagnostics;
end

function result = buildPolynomialDirection(problem,k,l,gamma,slope,spaces,quadrature)
wvt = problem.originatingTransform;
rho0 = wvt.rho0;
f = wvt.f;
D = wvt.Lz;
xi = quadrature.xi;
W = diag(quadrature.weight);
N2 = wvt.N2Function(gamma*xi);
if isscalar(N2)
    N2 = repmat(N2,numel(xi),1);
end
N2 = N2(:);
WN = W*diag(N2);
sx = slope(1);
sy = slope(2);

F = spaces.F;
Fxi = spaces.Fxi;
G = spaces.G;
Gxi = spaces.Gxi;
chi = spaces.chi;
chixi = spaces.chixi;
nF = size(F,2);
nG = size(G,2);
nX = 2*nF+2*nG+1;
iu = 1:nF;
iv = nF+(1:nF);
iw = 2*nF+(1:nG);
ieta = 2*nF+nG+(1:nG);
ib = nX;

Ru = zeros(numel(xi),nX);
Rv = zeros(numel(xi),nX);
Rwh = zeros(numel(xi),nX);
Reta = zeros(numel(xi),nX);
RetaXi = zeros(numel(xi),nX);
Ru(:,iu) = F;
Rv(:,iv) = F;
Rwh(:,iw) = G;
Reta(:,ieta) = G;
Reta(:,ib) = chi;
RetaXi(:,ieta) = Gxi;
RetaXi(:,ib) = chixi;
RwhXi = zeros(numel(xi),nX);
RwhXi(:,iw) = Gxi;
U = Ru/gamma;
V = Rv/gamma;
Wphysical = Rwh-(xi/(D*gamma)).*(sx*Ru+sy*Rv);

pressureX = 1i*k*F+(xi*sx/(D*gamma)).*Fxi;
pressureY = 1i*l*F+(xi*sy/(D*gamma)).*Fxi;
uBottom = zeros(1,nX);
vBottom = zeros(1,nX);
uBottom(iu) = spaces.Fendpoint(1,:)/gamma;
vBottom(iv) = spaces.Fendpoint(1,:)/gamma;

Cfull = F'*W*(1i*k*Ru+1i*l*Rv+RwhXi);
[C,N,constraintRank] = independentConstraints(Cfull);
[P,pressureDimension] = pressureGaugeBasis(nF,hypot(k,l));

% Assemble the coupled weak momentum-displacement equations. An admissible
% horizontal test velocity carries a metric-induced physical vertical
% velocity, so the horizontal and vertical strong rows cannot be projected
% independently without losing the terrain-energy identity.
Efull = rho0*gamma*(U'*W*U+V'*W*V+Wphysical'*W*Wphysical+Reta'*WN*Reta);
Jfull = rho0*gamma*(f*(U'*W*V-V'*W*U)+Reta'*WN*Wphysical-Wphysical'*WN*Reta);
GpFull = -(gamma*U'*W*pressureX+gamma*V'*W*pressureY+Wphysical'*W*Fxi);
M = Efull;
A = Jfull;
Gp = GpFull*P;

S = [M -Gp;C zeros(size(C,1),pressureDimension)];
rhs = [A*N;zeros(size(C,1),size(N,2))];
solution = S\rhs;
xdot = solution(1:nX,:);
pressure = P*solution(nX+1:end,:);
L = N'*xdot;
saddleResidual = norm(S*solution-rhs,"fro")/max(norm(rhs,"fro"),realmin);
constraintTangencyDefect = norm(C*xdot,"fro")/max(norm(C,"fro")*norm(xdot,"fro"),realmin);

Qfull = 1i*k*V-1i*l*U+(sx/(D*gamma))*V-(sy/(D*gamma))*U ...
    +(xi*sx/(D*gamma)).*(Fxi/gamma*selector(nX,iu)') ...
    -(xi*sy/(D*gamma)).*(Fxi/gamma*selector(nX,iv)') ...
    -(f/gamma)*RetaXi;
E = N'*Efull*N;
J = N'*Jfull*N;
Q = Qfull*N;
Z = Q'*(gamma*W)*Q;
Bfull = zeros(1,nX);
Bfull(ib) = 1;
B = Bfull*N;
Rfull = sx*uBottom+sy*vBottom;
R = Rfull*N;

pressureField = F*pressure;
pressureXi = Fxi*pressure;
pressureXWork = gamma*(U*N)'*W*(1i*k*pressureField+(xi*sx/(D*gamma)).*pressureXi);
pressureYWork = gamma*(V*N)'*W*(1i*l*pressureField+(xi*sy/(D*gamma)).*pressureXi);
pressureZWork = (Wphysical*N)'*W*pressureXi;
pressureWorkScale = sumNorms(pressureXWork,pressureYWork,pressureZWork);
pressureWorkDefect = norm(pressureXWork+pressureYWork+pressureZWork,"fro")/pressureWorkScale;

flatMatrices = polynomialMatrices(problem,k,l,gamma,[0 0],spaces,quadrature);
unitX = polynomialMatrices(problem,k,l,gamma,[1 0],spaces,quadrature);
minusX = polynomialMatrices(problem,k,l,gamma,[-1 0],spaces,quadrature);
unitY = polynomialMatrices(problem,k,l,gamma,[0 1],spaces,quadrature);
minusY = polynomialMatrices(problem,k,l,gamma,[0 -1],spaces,quadrature);
if sx == 1 && sy == 0
    derivative = analyticSaddleDerivative(flatMatrices,unitX,minusX,N);
elseif sx == 0 && sy == 1
    derivative = analyticSaddleDerivative(flatMatrices,unitY,minusY,N);
else
    dx = analyticSaddleDerivative(flatMatrices,unitX,minusX,N);
    dy = analyticSaddleDerivative(flatMatrices,unitY,minusY,N);
    derivative = combineDerivative(dx,dy,sx,sy);
end

result = struct("L",L,"E",E,"J",J,"Q",Q,"Z",Z,"B",B,"R",R, ...
    "N",N,"M",M,"A",A,"pressureCoupling",Gp,"constraint",C, ...
    "pressureMap",pressure,"reconstruction",struct("u",U*N,"v",V*N, ...
    "wHat",Rwh*N,"w",Wphysical*N,"eta",Reta*N,"etaXi",RetaXi*N), ...
    "rawReconstruction",struct("uHat",Ru,"vHat",Rv,"wHat",Rwh, ...
    "etaHat",Reta,"etaHatXi",RetaXi), ...
    "bottomCoordinateIndex",ib, ...
    "constraintRank",constraintRank,"pressureDimension",pressureDimension, ...
    "saddleResidual",saddleResidual,"constraintTangencyDefect",constraintTangencyDefect, ...
    "pressureWorkDefect",pressureWorkDefect,"L1",derivative.L1, ...
    "E1",derivative.E1,"J1",derivative.J1,"Q1",derivative.Q1, ...
    "Z1",derivative.Z1,"R1",derivative.R1);
end

function matrices = polynomialMatrices(problem,k,l,gamma,slope,spaces,quadrature)
% Return unreduced affine descriptor matrices and exact forms.
wvt = problem.originatingTransform;
rho0 = wvt.rho0;
f = wvt.f;
D = wvt.Lz;
xi = quadrature.xi;
W = diag(quadrature.weight);
N2 = wvt.N2Function(gamma*xi);
if isscalar(N2)
    N2 = repmat(N2,numel(xi),1);
end
WN = W*diag(N2(:));
sx = slope(1);
sy = slope(2);
F = spaces.F;
Fxi = spaces.Fxi;
G = spaces.G;
Gxi = spaces.Gxi;
chi = spaces.chi;
chixi = spaces.chixi;
nF = size(F,2);
nG = size(G,2);
nX = 2*nF+2*nG+1;
iu = 1:nF;
iv = nF+(1:nF);
iw = 2*nF+(1:nG);
ieta = 2*nF+nG+(1:nG);
ib = nX;
Ru = zeros(numel(xi),nX);
Rv = Ru;
Rwh = Ru;
Reta = Ru;
RetaXi = Ru;
RwhXi = Ru;
Ru(:,iu) = F;
Rv(:,iv) = F;
Rwh(:,iw) = G;
RwhXi(:,iw) = Gxi;
Reta(:,ieta) = G;
Reta(:,ib) = chi;
RetaXi(:,ieta) = Gxi;
RetaXi(:,ib) = chixi;
U = Ru/gamma;
V = Rv/gamma;
Wphysical = Rwh-(xi/(D*gamma)).*(sx*Ru+sy*Rv);
Cfull = F'*W*(1i*k*Ru+1i*l*Rv+RwhXi);
[C,N] = independentConstraints(Cfull);
[P,pressureDimension] = pressureGaugeBasis(nF,hypot(k,l));
Efull = rho0*gamma*(U'*W*U+V'*W*V+Wphysical'*W*Wphysical+Reta'*WN*Reta);
Jfull = rho0*gamma*(f*(U'*W*V-V'*W*U)+Reta'*WN*Wphysical-Wphysical'*WN*Reta);
pressureX = 1i*k*F+(xi*sx/(D*gamma)).*Fxi;
pressureY = 1i*l*F+(xi*sy/(D*gamma)).*Fxi;
GpFull = -(gamma*U'*W*pressureX+gamma*V'*W*pressureY+Wphysical'*W*Fxi);
M = Efull;
A = Jfull;
Gp = GpFull*P;
S = [M -Gp;C zeros(size(C,1),pressureDimension)];
rhs = [A*N;zeros(size(C,1),size(N,2))];
Qfull = 1i*k*V-1i*l*U+(sx/(D*gamma))*V-(sy/(D*gamma))*U ...
    +(xi*sx/(D*gamma)).*(Fxi/gamma*selector(nX,iu)') ...
    -(xi*sy/(D*gamma)).*(Fxi/gamma*selector(nX,iv)')-(f/gamma)*RetaXi;
Bfull = zeros(1,nX);
Bfull(ib) = 1;
uBottom = zeros(1,nX);
vBottom = zeros(1,nX);
uBottom(iu) = spaces.Fendpoint(1,:)/gamma;
vBottom(iv) = spaces.Fendpoint(1,:)/gamma;
Rfull = sx*uBottom+sy*vBottom;
matrices = struct("S",S,"rhs",rhs,"N",N,"Efull",Efull,"Jfull",Jfull, ...
    "Qfull",Qfull,"WQ",gamma*W,"Bfull",Bfull,"Rfull",Rfull);
end

function derivative = analyticSaddleDerivative(flat,plus,minus,N)
S1 = (plus.S-minus.S)/2;
rhs1 = (plus.rhs-minus.rhs)/2;
solution0 = flat.S\flat.rhs;
solution1 = flat.S\(rhs1-S1*solution0);
nX = size(N,1);
L1 = N'*solution1(1:nX,:);
E0 = N'*flat.Efull*N;
E1 = N'*((plus.Efull-minus.Efull)/2)*N;
J1 = N'*((plus.Jfull-minus.Jfull)/2)*N;
Q0 = flat.Qfull*N;
Q1 = ((plus.Qfull-minus.Qfull)/2)*N;
Z1 = Q1'*flat.WQ*Q0+Q0'*flat.WQ*Q1;
R1 = ((plus.Rfull-minus.Rfull)/2)*N;
derivative = struct("L1",L1,"E1",E1,"J1",J1,"Q1",Q1,"Z1",Z1,"R1",R1,"E0",E0);
end

function derivative = combineDerivative(dx,dy,sx,sy)
names = ["L1","E1","J1","Q1","Z1","R1"];
derivative = struct;
for name = names
    derivative.(name) = sx*dx.(name)+sy*dy.(name);
end
end

function defect = derivativeCheck(problem,k,l,gamma,spaces,quadrature,slope,L1)
scale = norm(slope);
if scale == 0
    defect = 0;
    return
end
direction = slope/scale;
step = 1e-5;
plus = generatorOnly(problem,k,l,gamma,step*direction,spaces,quadrature);
minus = generatorOnly(problem,k,l,gamma,-step*direction,spaces,quadrature);
centered = scale*(plus-minus)/(2*step);
defect = norm(L1-centered,"fro")/max(norm(L1,"fro")+norm(centered,"fro"),realmin);
end

function L = generatorOnly(problem,k,l,gamma,slope,spaces,quadrature)
matrices = polynomialMatrices(problem,k,l,gamma,slope,spaces,quadrature);
solution = matrices.S\matrices.rhs;
nX = size(matrices.N,1);
L = matrices.N'*solution(1:nX,:);
end

function block = fgComparisonBlock(problem,iK,slope,tangentStep,flatBlock,slopeBlock,plusBlock,minusBlock)
scale = norm(slope);
if scale == 0
    direction = [0 0];
else
    direction = slope/scale;
end
flat = reduceFGDescriptor(problem,iK,flatBlock,[0 0]);
finite = reduceFGDescriptor(problem,iK,slopeBlock,slope);
plus = reduceFGDescriptor(problem,iK,plusBlock,tangentStep*direction);
minus = reduceFGDescriptor(problem,iK,minusBlock,-tangentStep*direction);
if scale == 0
    firstOrder = struct("L",zeros(size(flat.L)),"E",zeros(size(flat.E)), ...
        "J",zeros(size(flat.J)),"Q",zeros(size(flat.Q)), ...
        "Z",zeros(size(flat.Z)),"R",zeros(size(flat.R)));
else
    firstOrder = struct("L",scale*(plus.L-minus.L)/(2*tangentStep), ...
        "E",scale*(plus.E-minus.E)/(2*tangentStep), ...
        "J",scale*(plus.J-minus.J)/(2*tangentStep), ...
        "Q",scale*(plus.Q-minus.Q)/(2*tangentStep), ...
        "Z",scale*(plus.Z-minus.Z)/(2*tangentStep), ...
        "R",scale*(plus.R-minus.R)/(2*tangentStep));
end
native = coefficientDefects(flat,firstOrder);

T = flat.publicCoordinateMap;
publicFlat = transformReducedSystem(flat,T);
publicFirstOrder = transformFirstOrderSystem(firstOrder,T);
public = coefficientDefects(publicFlat,publicFirstOrder);
historical = slopeBlock.localForms.diagnostics;
block = struct("horizontalIndex",iK,"kMode",problem.horizontalLayout.kMode(iK), ...
    "lMode",problem.horizontalLayout.lMode(iK),"native",native,"public",public, ...
    "historicalPublic",historical,"flat",flat,"finiteSlope",finite, ...
    "firstOrder",firstOrder,"publicCoordinateMap",T, ...
    "publicGeneratorDifference",relativeNorm(publicFlat.L-flatBlock.publicGenerator,flatBlock.publicGenerator));
end

function reduced = reduceFGDescriptor(problem,iK,block,slope)
layout = block.layout;
xColumns = find(~layout.isPressure);
pColumns = find(layout.isPressure);
allRows = (1:height(layout))';
constraintRows = block.continuityRows(:);
dynamicRows = setdiff(allRows,constraintRows,"stable");
k = problem.horizontalLayout.k(iK);
l = problem.horizontalLayout.l(iK);
if hypot(k,l) == 0
    pressureMode = layout.j(pColumns);
    P = eye(numel(pColumns));
    P = P(:,pressureMode ~= 0);
    constraintPosition = pressureMode ~= 0;
    constraintRows = constraintRows(constraintPosition);
else
    P = eye(numel(pColumns));
end

Cfull = block.operatorMatrix(constraintRows,xColumns);
[C,N] = independentConstraints(Cfull);
M = block.massMatrix(dynamicRows,xColumns);
A = block.operatorMatrix(dynamicRows,xColumns);
Gp = block.operatorMatrix(dynamicRows,pColumns)*P;
S = [M -Gp;C zeros(size(C,1),size(P,2))];
rhs = [A*N;zeros(size(C,1),size(N,2))];
solution = S\rhs;
nX = numel(xColumns);
xdot = solution(1:nX,:);
L = N'*xdot;

wvt = problem.originatingTransform;
gamma = problem.gamma(1);
xi = wvt.z(:);
W = diag(wvt.z_int(:));
N2 = wvt.N2Function(gamma*xi);
if isscalar(N2)
    N2 = repmat(N2,wvt.Nz,1);
end
WN = W*diag(N2(:));
reconstruction = block.reconstruction;
U = reconstruction.uHat(:,xColumns)/gamma;
V = reconstruction.vHat(:,xColumns)/gamma;
Wphysical = block.physicalReconstruction.w(:,xColumns);
Eta = reconstruction.etaHat(:,xColumns);
EtaXi = reconstruction.etaHatXi(:,xColumns);
UXi = reconstruction.uHatXi(:,xColumns)/gamma;
VXi = reconstruction.vHatXi(:,xColumns)/gamma;
sx = slope(1);
sy = slope(2);
Efull = wvt.rho0*gamma*(U'*W*U+V'*W*V+Wphysical'*W*Wphysical+Eta'*WN*Eta);
Jfull = wvt.rho0*gamma*(wvt.f*(U'*W*V-V'*W*U)+Eta'*WN*Wphysical-Wphysical'*WN*Eta);
Qfull = 1i*k*V-1i*l*U+(sx/(wvt.Lz*gamma))*V-(sy/(wvt.Lz*gamma))*U ...
    +(xi*sx/(wvt.Lz*gamma)).*VXi-(xi*sy/(wvt.Lz*gamma)).*UXi ...
    -(wvt.f/gamma)*EtaXi;
E = N'*Efull*N;
J = N'*Jfull*N;
Q = Qfull*N;
Z = Q'*(gamma*W)*Q;
bottomColumn = find(layout.variable(xColumns) == "etaB",1);
Bfull = zeros(1,nX);
Bfull(bottomColumn) = 1;
B = Bfull*N;
uBottom = U(1,:);
vBottom = V(1,:);
R = (sx*uBottom+sy*vBottom)*N;

sqrtN2 = sqrt(N2(:));
descriptorReconstruction = [reconstruction.uHat(:,xColumns); ...
    reconstruction.vHat(:,xColumns);reconstruction.wHat(:,xColumns); ...
    sqrtN2.*reconstruction.etaHat(:,xColumns)]*N;
publicBasis = problem.basisBlocks{iK};
publicReconstruction = [publicBasis.uHat;publicBasis.vHat;publicBasis.wHat; ...
    sqrtN2.*publicBasis.etaHat];
T = publicReconstruction\descriptorReconstruction;
reduced = struct("L",L,"E",E,"J",J,"Q",Q,"Z",Z,"B",B,"R",R, ...
    "N",N,"publicCoordinateMap",T,"saddleResidual", ...
    norm(S*solution-rhs,"fro")/max(norm(rhs,"fro"),realmin));
end

function diagnostics = coefficientDefects(flat,firstOrder)
weakScale = sumNorms(flat.E*firstOrder.L,firstOrder.E*flat.L,firstOrder.J);
energyScale = sumNorms(firstOrder.L'*flat.E,flat.E*firstOrder.L, ...
    flat.L'*firstOrder.E,firstOrder.E*flat.L);
apvScale = sumNorms(flat.Q*firstOrder.L,firstOrder.Q*flat.L);
bottomScale = sumNorms(flat.B*firstOrder.L,firstOrder.R);
enstrophyScale = sumNorms(firstOrder.L'*flat.Z,flat.Z*firstOrder.L, ...
    flat.L'*firstOrder.Z,firstOrder.Z*flat.L);
diagnostics = struct( ...
    "weakEvolutionDefect",norm(flat.E*firstOrder.L+firstOrder.E*flat.L-firstOrder.J,"fro")/weakScale, ...
    "energySkewDefect",norm(firstOrder.L'*flat.E+flat.E*firstOrder.L ...
    +flat.L'*firstOrder.E+firstOrder.E*flat.L,"fro")/energyScale, ...
    "apvTendencyDefect",norm(flat.Q*firstOrder.L+firstOrder.Q*flat.L,"fro")/apvScale, ...
    "potentialEnstrophyDefect",norm(firstOrder.L'*flat.Z+flat.Z*firstOrder.L ...
    +flat.L'*firstOrder.Z+firstOrder.Z*flat.L,"fro")/enstrophyScale, ...
    "bottomEvolutionDefect",norm(flat.B*firstOrder.L-firstOrder.R,"fro")/bottomScale);
end

function transformed = transformReducedSystem(system,T)
transformed = struct("L",T*system.L/T,"E",T'\system.E/T, ...
    "J",T'\system.J/T,"Q",system.Q/T,"Z",T'\system.Z/T, ...
    "B",system.B/T,"R",system.R/T);
end

function transformed = transformFirstOrderSystem(system,T)
transformed = struct("L",T*system.L/T,"E",T'\system.E/T, ...
    "J",T'\system.J/T,"Q",system.Q/T,"Z",T'\system.Z/T, ...
    "R",system.R/T);
end

function maximum = maximumDiagnostics(blocks)
names = ["weakEvolutionDefect","energyDefect","apvDefect","bottomDefect", ...
    "enstrophyDefect","pressureWorkDefect","constraintTangencyDefect", ...
    "saddleResidual","analyticDerivativeDefect"];
required = struct;
for name = names
    required.(name) = max(cellfun(@(block)block.diagnostics.(name),blocks));
end
maximum = struct("required",required, ...
    "finiteSlopeLinearizationDefect",max(cellfun(@(block)block.diagnostics.finiteSlopeLinearizationDefect,blocks)), ...
    "minimumEnergyEigenvalue",min(cellfun(@(block)safeMinimum(block.diagnostics.minimumEnergyEigenvalue),blocks)), ...
    "flatModeOneDispersionDefect",max(cellfun(@(block)block.diagnostics.flatModeOneDispersionDefect,blocks)), ...
    "flat",maximumFlatDiagnostics(blocks));
end

function maximum = maximumComparisonDiagnostics(blocks)
names = ["weakEvolutionDefect","energySkewDefect","apvTendencyDefect", ...
    "potentialEnstrophyDefect","bottomEvolutionDefect"];
native = struct;
public = struct;
for name = names
    native.(name) = max(cellfun(@(block)block.native.(name),blocks));
    public.(name) = max(cellfun(@(block)block.public.(name),blocks));
end
maximum = struct("nativeRequired",native,"publicRequired",public);
end

function diagnostics = polynomialGeometryDiagnostics(spaces,quadrature)
W = diag(quadrature.weight);
endpointResidual = max([norm(spaces.Gendpoint,"fro"), ...
    abs(spaces.chiEndpoint(1)-1),abs(spaces.chiEndpoint(2))]);
greenNumerator = spaces.F'*W*spaces.Gxi+spaces.Fxi'*W*spaces.G;
boundary = spaces.Fendpoint(2,:)'*spaces.Gendpoint(2,:) ...
    -spaces.Fendpoint(1,:)'*spaces.Gendpoint(1,:);
greenScale = norm(spaces.F'*W*spaces.Gxi,"fro") ...
    +norm(spaces.Fxi'*W*spaces.G,"fro")+norm(boundary,"fro");
diagnostics = struct("endpointResidual",endpointResidual, ...
    "greenIdentityDefect",norm(greenNumerator-boundary,"fro")/max(greenScale,realmin), ...
    "bottomFunctionBottomValue",1, ...
    "bottomFunctionSurfaceValue",0);
end

function diagnostics = exactSystemDefects(system)
diagnostics = struct( ...
    "weakEvolutionDefect",norm(system.E*system.L-system.J,"fro") ...
    /max(norm(system.E*system.L,"fro")+norm(system.J,"fro"),realmin), ...
    "energyDefect",norm(system.L'*system.E+system.E*system.L,"fro") ...
    /max(norm(system.E,"fro")*norm(system.L,"fro"),realmin), ...
    "apvDefect",norm(system.Q*system.L,"fro") ...
    /max(norm(system.Q,"fro")*norm(system.L,"fro"),realmin), ...
    "enstrophyDefect",norm(system.L'*system.Z+system.Z*system.L,"fro") ...
    /max(norm(system.Z,"fro")*norm(system.L,"fro"),realmin), ...
    "bottomDefect",norm(system.B*system.L-system.R,"fro") ...
    /max(norm(system.B,"fro")*norm(system.L,"fro")+norm(system.R,"fro"),realmin));
end

function maximum = maximumFlatDiagnostics(blocks)
names = ["weakEvolutionDefect","energyDefect","apvDefect","enstrophyDefect","bottomDefect"];
maximum = struct;
for name = names
    maximum.(name) = max(cellfun(@(block)block.diagnostics.flat.(name),blocks));
end
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

function spaces = polynomialSpaces(degree,quadrature,D)
r = quadrature.r;
[P,Pr] = legendreValues(r,degree);
F = P(:,1:degree+1);
Fxi = (2/D)*Pr(:,1:degree+1);
G = (1-r.^2).*P(:,1:degree);
Gxi = (2/D)*(-2*r.*P(:,1:degree)+(1-r.^2).*Pr(:,1:degree));
chi = (1-r)/2;
chixi = -ones(size(r))/D;
[Pb,Prb] = legendreValues([-1;1],degree);
Fb = Pb(:,1:degree+1);
Fbx = (2/D)*Prb(:,1:degree+1);
Gb = (1-[-1;1].^2).*Pb(:,1:degree);
Gbx = (2/D)*(-2*[-1;1].*Pb(:,1:degree)+(1-[-1;1].^2).*Prb(:,1:degree));
spaces = struct("F",F,"Fxi",Fxi,"G",G,"Gxi",Gxi,"chi",chi, ...
    "chixi",chixi,"Fendpoint",Fb,"FxiEndpoint",Fbx, ...
    "Gendpoint",Gb,"GxiEndpoint",Gbx,"chiEndpoint",[1;0]);
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
for n = 1:degree-1
    P(:,n+2) = ((2*n+1)*r.*P(:,n+1)-n*P(:,n))/(n+1);
    Pr(:,n+2) = ((2*n+1)*(P(:,n+1)+r.*Pr(:,n+1))-n*Pr(:,n))/(n+1);
end
end

function [C,N,rankC] = independentConstraints(Cfull)
[U,S,V] = svd(Cfull);
singularValues = diag(S);
tolerance = max(size(Cfull))*eps(max(singularValues,1));
rankC = nnz(singularValues > tolerance);
C = U(:,1:rankC)'*Cfull;
N = V(:,rankC+1:end);
N = deterministicColumns(N);
end

function [P,nPressure] = pressureGaugeBasis(nF,kappa)
if kappa == 0
    P = eye(nF);
    P = P(:,2:end);
else
    P = eye(nF);
end
nPressure = size(P,2);
end

function V = deterministicColumns(V)
for column = 1:size(V,2)
    [~,row] = max(abs(V(:,column)));
    phase = V(row,column)/abs(V(row,column));
    if phase ~= 0
        V(:,column) = V(:,column)/phase;
    end
end
end

function S = selector(n,indices)
S = zeros(n,numel(indices));
S(indices,:) = eye(numel(indices));
end

function value = sumNorms(varargin)
value = 0;
for i = 1:nargin
    value = value+norm(varargin{i},"fro");
end
value = max(value,realmin);
end

function value = relativeNorm(numerator,denominator)
value = norm(numerator,"fro")/max(norm(denominator,"fro"),realmin);
end

function defect = flatModeOneDispersionDefect(problem,k,l,gamma,L)
if gamma ~= 1 || hypot(k,l) == 0
    defect = 0;
    return
end
wvt = problem.originatingTransform;
sample = wvt.N2Function([-wvt.Lz;0]);
if numel(sample) ~= 2 || abs(sample(1)-sample(2)) > 100*eps(max([abs(sample(:));1]))
    defect = NaN;
    return
end
kappa = hypot(k,l);
m = pi/wvt.Lz;
expected = sqrt((sample(1)*kappa^2+wvt.f^2*m^2)/(kappa^2+m^2));
frequency = abs(imag(eig(L)));
frequency = frequency(frequency > 100*eps*max(abs(wvt.f),1));
defect = min(abs(frequency-expected))/expected;
end

function value = safeMinimum(values)
if isempty(values)
    value = Inf;
else
    value = min(values,[],"all");
end
end
