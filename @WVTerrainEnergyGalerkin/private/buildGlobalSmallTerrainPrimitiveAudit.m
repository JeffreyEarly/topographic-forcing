function audit = buildGlobalSmallTerrainPrimitiveAudit(problem,polynomialDegree,quadratureOrder,tangentStep)
% Build and classify the global first-order primitive terrain oracle.

context = buildContext(problem,polynomialDegree,quadratureOrder);
flat = descriptorAtScale(context,0);
analytic = analyticTangent(context,flat);
steps = tangentStep./(2.^(0:2));
for iStep = 1:numel(steps)
    plus = descriptorAtScale(context,steps(iStep));
    minus = descriptorAtScale(context,-steps(iStep));
    if iStep == 1
        centered = repmat(centeredTangent(plus,minus,steps(iStep)),numel(steps),1);
    else
        centered(iStep) = centeredTangent(plus,minus,steps(iStep));
    end
end

tangentAgreement = tangentAgreementDiagnostics(analytic,centered);
compatibility = compatibilityDiagnostics(context,flat,analytic);
strongAPV = strongAPVDiagnostics(context,flat,analytic);
fourier = fourierDiagnostics(context,analytic,compatibility.apvCancellation);
flatDiagnostics = flatReferenceDiagnostics(context,flat);
structure = structureDiagnostics(context,flat,analytic);

requiredTolerance = struct("tangent",1e-9,"weak",1e-11,"energy",1e-11, ...
    "apv",1e-10,"enstrophy",1e-10,"bottom",1e-11,"conjugacy",1e-11);
tangentPasses = tangentAgreement.maximumRelativeDefect <= requiredTolerance.tangent;
identityPasses = compatibility.weakEvolutionDefect <= requiredTolerance.weak ...
    && compatibility.energyDefect <= requiredTolerance.energy ...
    && compatibility.apvDefect <= requiredTolerance.apv ...
    && compatibility.enstrophyDefect <= requiredTolerance.enstrophy ...
    && compatibility.bottomDefect <= requiredTolerance.bottom ...
    && compatibility.conjugacyDefect <= requiredTolerance.conjugacy;
structurePasses = structure.maximumHermitianDefect <= 1e-12 ...
    && structure.maximumSkewHermitianDefect <= 1e-12 ...
    && flat.diagnostics.pressureNullity == 1 ...
    && flat.diagnostics.gaugedSaddleResidual <= 1e-11;

if ~tangentPasses
    status = "implementation-unresolved";
    diagnosis = "The analytic and independently centered global terrain tangents do not agree.";
elseif structurePasses && identityPasses
    status = "compatible-global-first-order-oracle";
    diagnosis = "Global terrain coupling restores simultaneous first-order energy, APV, enstrophy, and bottom compatibility.";
elseif structurePasses && compatibility.weakEvolutionDefect <= requiredTolerance.weak ...
        && compatibility.energyDefect <= requiredTolerance.energy ...
        && compatibility.bottomDefect <= requiredTolerance.bottom
    status = "representation-blocker";
    diagnosis = "The global strong equations preserve first-order APV, but the finite primitive Galerkin projection does not.";
else
    status = "scientific-blocker";
    diagnosis = "The independently verified global first-order primitive system fails at least one mandatory physical identity.";
end

audit = struct;
audit.status = status;
audit.isCompatible = status == "compatible-global-first-order-oracle";
audit.scope = "global-small-terrain-first-order-only";
audit.polynomialDegree = polynomialDegree;
audit.quadratureOrder = quadratureOrder;
audit.tangentSteps = steps;
audit.layout = context.layout;
audit.flatReference = flat;
audit.analyticTangent = analytic;
audit.centeredTangents = centered;
audit.tangentAgreement = tangentAgreement;
audit.compatibility = compatibility;
audit.apvCancellation = compatibility.apvCancellation;
audit.strongAPV = strongAPV;
audit.fourier = fourier;
audit.flatDiagnostics = flatDiagnostics;
audit.structure = structure;
audit.requiredTolerance = requiredTolerance;
audit.diagnosis = diagnosis;
audit.nextScope = "finite-amplitude-periodic-terrain-only-if-compatible";
end

function context = buildContext(problem,degree,quadratureOrder)
wvt = problem.originatingTransform;
quadrature = legendreQuadrature(quadratureOrder,wvt.Lz);
spaces = polynomialSpaces(degree,quadrature,wvt.Lz);
nK = height(problem.horizontalLayout);
nF = size(spaces.F,2);
nG = size(spaces.G,2);
nH = size(spaces.H,2);
nP = size(spaces.Pressure,2);
nXBlock = 2*nF+nG+nH;
nX = nK*nXBlock;
nPressure = nK*nP;

oversampling = problem.horizontalOversamplingFactor;
Nx = oversampling*wvt.Nx;
Ny = oversampling*wvt.Ny;
nXY = Nx*Ny;
nZ = numel(quadrature.xi);
[x,y] = ndgrid((0:Nx-1)'*wvt.Lx/Nx,(0:Ny-1)'*wvt.Ly/Ny);
phase = exp(1i*(x(:)*problem.horizontalLayout.k.'+y(:)*problem.horizontalLayout.l.'));
phaseGramDefect = norm(phase'*phase/nXY-eye(nK),"fro")/sqrt(nK);
if phaseGramDefect > 1e-12
    error("WVTerrainEnergyGalerkin:GlobalPrimitivePhaseFailure", ...
        "The oversampled signed Fourier phases are not orthonormal; defect %.3g.",phaseGramDefect)
end

h = interpft(interpft(problem.topographicHeight,Nx,1),Ny,2);
if max(abs(imag(h)),[],"all") > 1e-12*max(1,max(abs(h),[],"all"))
    error("WVTerrainEnergyGalerkin:TerrainInterpolationFailure", ...
        "Periodic terrain interpolation produced a material imaginary component.")
end
h = real(h);
[hX,hY] = horizontalDerivatives(h,wvt.Lx,wvt.Ly);
H = h/wvt.Lz;
HX = hX/wvt.Lz;
HY = hY/wvt.Lz;

nGrid = nXY*nZ;
Ru = zeros(nGrid,nX);
Rv = zeros(nGrid,nX);
Rwh = zeros(nGrid,nX);
RwhXi = zeros(nGrid,nX);
RuXi = zeros(nGrid,nX);
RvXi = zeros(nGrid,nX);
Reta = zeros(nGrid,nX);
RetaXi = zeros(nGrid,nX);
TF = zeros(nGrid,nK*nF);
TH = zeros(nGrid,nK*nH);
RP = zeros(nGrid,nPressure);
RPXi = zeros(nGrid,nPressure);
RPX = zeros(nGrid,nPressure);
RPY = zeros(nGrid,nPressure);
uBottom = zeros(nXY,nX);
vBottom = zeros(nXY,nX);
bottomValueFull = zeros(nK,nX);

for iK = 1:nK
    xCols = (iK-1)*nXBlock+(1:nXBlock);
    iu = xCols(1:nF);
    iv = xCols(nF+(1:nF));
    iw = xCols(2*nF+(1:nG));
    ieta = xCols(2*nF+nG+(1:nH));
    fCols = (iK-1)*nF+(1:nF);
    hCols = (iK-1)*nH+(1:nH);
    pCols = (iK-1)*nP+(1:nP);
    phaseK = phase(:,iK);
    Ru(:,iu) = kron(spaces.F,phaseK);
    Rv(:,iv) = kron(spaces.F,phaseK);
    RuXi(:,iu) = kron(spaces.Fxi,phaseK);
    RvXi(:,iv) = kron(spaces.Fxi,phaseK);
    Rwh(:,iw) = kron(spaces.G,phaseK);
    RwhXi(:,iw) = kron(spaces.Gxi,phaseK);
    Reta(:,ieta) = kron(spaces.H,phaseK);
    RetaXi(:,ieta) = kron(spaces.Hxi,phaseK);
    TF(:,fCols) = kron(spaces.F,phaseK);
    TH(:,hCols) = kron(spaces.H,phaseK);
    RP(:,pCols) = kron(spaces.Pressure,phaseK);
    RPXi(:,pCols) = kron(spaces.PressureXi,phaseK);
    RPX(:,pCols) = 1i*problem.horizontalLayout.k(iK)*RP(:,pCols);
    RPY(:,pCols) = 1i*problem.horizontalLayout.l(iK)*RP(:,pCols);
    uBottom(:,iu) = phaseK*spaces.Fendpoint(1,:);
    vBottom(:,iv) = phaseK*spaces.Fendpoint(1,:);
    bottomValueFull(iK,ieta) = spaces.Hendpoint(1,:);
end

kByState = repelem(problem.horizontalLayout.k,nXBlock);
lByState = repelem(problem.horizontalLayout.l,nXBlock);
RuX = Ru.*(1i*kByState.');
RuY = Ru.*(1i*lByState.');
RvX = Rv.*(1i*kByState.');
RvY = Rv.*(1i*lByState.');
divergence = RuX+RvY+RwhXi;
volumeWeight = kron(quadrature.weight,ones(nXY,1)/nXY);
continuity = TF'*(volumeWeight.*divergence);
[N,admissibleDimensions,admissibleRanges,continuityDiagnostics] = admissibleBasis(problem,continuity,nXBlock,nF);

zeroK = find(problem.horizontalLayout.kMode == 0 & problem.horizontalLayout.lMode == 0,1);
if isempty(zeroK)
    error("WVTerrainEnergyGalerkin:MissingMeanHorizontalMode", ...
        "The retained signed layout must include the zero horizontal mode.")
end
pressureGauge = zeros(nX+nPressure,1);
pressureGauge(nX+(zeroK-1)*nP+1) = 1;

horizontalConjugateIndex = conjugateHorizontalIndex(problem);
coordinateConjugateIndex = zeros(size(N,2),1);
for iK = 1:nK
    partner = horizontalConjugateIndex(iK);
    rows = admissibleRanges{iK};
    partnerRows = admissibleRanges{partner};
    coordinateConjugateIndex(rows) = partnerRows;
end

N20 = wvt.N2Function(quadrature.xi);
if isscalar(N20)
    N20 = repmat(N20,nZ,1);
end
N20 = N20(:);
if numel(N20) ~= nZ || ~isreal(N20) || any(~isfinite(N20)) || any(N20 <= 0)
    error("WVTerrainEnergyGalerkin:InvalidPrimitiveStratification", ...
        "N2Function must return finite positive values on the primitive quadrature.")
end
[Pq,Pqr] = legendreValues(quadrature.r,nZ-1);
N2Xi = (2/wvt.Lz)*(Pqr/Pq)*N20;

layout = struct("numberOfHorizontalModes",nK,"numberOfStateCoefficients",nX, ...
    "numberOfAdmissibleCoefficients",size(N,2),"numberOfPressureCoefficients",nPressure, ...
    "stateBlockSize",nXBlock,"admissibleBlockSizes",admissibleDimensions, ...
    "admissibleRanges",{admissibleRanges}, ...
    "nF",nF,"nG",nG,"nH",nH,"nP",nP,"oversampledSize",[Nx Ny nZ]);
context = struct("problem",problem,"wvt",wvt,"quadrature",quadrature,"spaces",spaces, ...
    "layout",layout,"nK",nK,"nF",nF,"nG",nG,"nH",nH,"nP",nP, ...
    "nXBlock",nXBlock,"nX",nX,"nPressure",nPressure,"nXY",nXY,"nZ",nZ, ...
    "phase",phase,"H",H(:),"HX",HX(:),"HY",HY(:),"hX",hX(:),"hY",hY(:), ...
    "xiGrid",kron(quadrature.xi,ones(nXY,1)),"volumeWeight",volumeWeight, ...
    "Ru",Ru,"Rv",Rv,"Rwh",Rwh,"RwhXi",RwhXi,"Reta",Reta,"RetaXi",RetaXi, ...
    "RuXi",RuXi,"RvXi",RvXi,"RuX",RuX,"RuY",RuY,"RvX",RvX,"RvY",RvY, ...
    "TF",TF,"TH",TH,"RP",RP,"RPXi",RPXi,"RPX",RPX,"RPY",RPY, ...
    "uBottom",uBottom,"vBottom",vBottom,"bottomValueFull",bottomValueFull, ...
    "continuity",continuity,"N",N,"pressureGauge",pressureGauge, ...
    "N20",N20,"N2Xi",N2Xi,"coordinateConjugateIndex",coordinateConjugateIndex, ...
    "continuityDiagnostics",continuityDiagnostics,"phaseGramDefect",phaseGramDefect);
end

function direction = descriptorAtScale(c,delta)
coefficients = exactCoefficients(c,delta);
[S,F,matrices] = descriptorMatrices(c,coefficients);
[Y,solveDiagnostics] = solveGaugedSystem(S,F,c.pressureGauge);
xDot = Y(1:c.nX,:);
pressureMap = Y(c.nX+(1:c.nPressure),:);
L = c.N'*xDot;
R = coefficients.bottomTendencyFull*c.N;
Q = coefficients.Q*c.N;
E = c.N'*coefficients.Efull*c.N;
J = c.N'*coefficients.Jfull*c.N;
Z = c.N'*(coefficients.Q'*(coefficients.apvWeight.*coefficients.Q))*c.N;
B = c.bottomValueFull*c.N;
direction = struct("scale",delta,"L",L,"E",E,"J",J,"Q",Q,"Z",Z,"B",B,"R",R, ...
    "pressureMap",pressureMap,"primitiveTendency",xDot,"matrices",matrices);
direction.diagnostics = solveDiagnostics;
direction.diagnostics.continuityTangencyDefect = norm(c.continuity*xDot,"fro") ...
    /max(norm(c.continuity,"fro")*norm(xDot,"fro"),realmin);
direction.diagnostics.weakEvolutionDefect = productDefect(E*L-J,{E*L,J});
direction.diagnostics.energyDefect = productDefect(L'*E+E*L,{L'*E,E*L});
direction.diagnostics.apvDefect = norm(Q*L,"fro") ...
    /max(norm(Q,"fro")*norm(L,"fro"),realmin);
direction.diagnostics.bottomDefect = productDefect(B*L-R,{B*L,R});
direction.diagnostics.enstrophyDefect = norm(L'*Z+Z*L,"fro") ...
    /max(norm(Z,"fro")*norm(L,"fro"),realmin);
end

function analytic = analyticTangent(c,flat)
first = firstOrderCoefficients(c);
[S1,F1,matrices1] = descriptorMatrices(c,first);
S0 = flat.matrices.S;
F0 = flat.matrices.F;
Y0 = [flat.primitiveTendency;flat.pressureMap];
[~,solve0] = solveGaugedSystem(S0,F0,c.pressureGauge);
leftNull = solve0.leftNullVector;
bordered = [S0 leftNull;c.pressureGauge' 0];
rhs1 = F1-S1*Y0;
Y1Augmented = bordered\[rhs1;zeros(1,size(rhs1,2))];
Y1 = Y1Augmented(1:end-1,:);
xDot1 = Y1(1:c.nX,:);
pressureMap1 = Y1(c.nX+(1:c.nPressure),:);
L1 = c.N'*xDot1;
E1 = c.N'*first.Efull*c.N;
J1 = c.N'*first.Jfull*c.N;
Q1 = first.Q*c.N;
R1 = first.bottomTendencyFull*c.N;
W0 = c.volumeWeight;
W1 = -repmat(c.H,c.nZ,1).*W0;
Q0Full = flat.Q;
Z1 = Q1'*(W0.*Q0Full)+Q0Full'*(W0.*Q1)+Q0Full'*(W1.*Q0Full);

analytic = struct("L",L1,"E",E1,"J",J1,"Q",Q1,"Z",Z1,"R",R1, ...
    "pressureMap",pressureMap1,"primitiveTendency",xDot1,"matrices",matrices1);
analytic.diagnostics = struct;
analytic.diagnostics.tangentSaddleResidual = norm(S0*Y1-rhs1,"fro")/max(norm(rhs1,"fro"),realmin);
analytic.diagnostics.gaugeDefect = norm(c.pressureGauge'*Y1,"fro")/max(norm(Y1,"fro"),realmin);
end

function coefficients = exactCoefficients(c,delta)
H3 = repmat(c.H,c.nZ,1);
HX3 = repmat(c.HX,c.nZ,1);
HY3 = repmat(c.HY,c.nZ,1);
gamma = 1-delta*H3;
if any(gamma <= 0)
    error("WVTerrainEnergyGalerkin:NonpositivePrimitiveGamma", ...
        "The centered primitive terrain scale produced gamma<=0.")
end
gradLnGammaX = -delta*HX3./gamma;
gradLnGammaY = -delta*HY3./gamma;
U = c.Ru./gamma;
V = c.Rv./gamma;
Wphysical = c.Rwh+c.xiGrid.*(gradLnGammaX.*c.Ru+gradLnGammaY.*c.Rv);
mappedDepth = gamma.*c.xiGrid;
N2 = c.wvt.N2Function(mappedDepth);
if isscalar(N2)
    N2 = repmat(N2,size(mappedDepth));
end
N2 = N2(:);
gammaN2 = gamma.*N2;
Dxp = c.RPX-c.xiGrid.*gradLnGammaX.*c.RPXi;
Dyp = c.RPY-c.xiGrid.*gradLnGammaY.*c.RPXi;
UXi = c.RuXi./gamma;
VXi = c.RvXi./gamma;
DxV = horizontalDerivative(V,c,"x")-c.xiGrid.*gradLnGammaX.*VXi;
DyU = horizontalDerivative(U,c,"y")-c.xiGrid.*gradLnGammaY.*UXi;
Q = DxV-DyU-(c.wvt.f./gamma).*c.RetaXi;
apvWeight = c.volumeWeight.*gamma;
gammaBottom = 1-delta*c.H;
bottomTendencyGrid = delta*(c.hX.*(c.uBottom./gammaBottom)+c.hY.*(c.vBottom./gammaBottom));
bottomTendencyFull = c.phase'*bottomTendencyGrid/c.nXY;
[Efull,Jfull] = physicalForms(c,gamma,gammaN2,U,V,Wphysical);
coefficients = struct("gamma",gamma,"gammaN2",gammaN2,"U",U,"V",V, ...
    "Wphysical",Wphysical,"gammaWphysical",gamma.*Wphysical, ...
    "gammaN2Wphysical",gammaN2.*Wphysical,"Dxp",Dxp,"Dyp",Dyp, ...
    "Q",Q,"apvWeight",apvWeight, ...
    "bottomTendencyFull",bottomTendencyFull,"Efull",Efull,"Jfull",Jfull);
end

function coefficients = firstOrderCoefficients(c)
H3 = repmat(c.H,c.nZ,1);
HX3 = repmat(c.HX,c.nZ,1);
HY3 = repmat(c.HY,c.nZ,1);
N20 = kron(c.N20,ones(c.nXY,1));
N2Xi = kron(c.N2Xi,ones(c.nXY,1));
U1 = H3.*c.Ru;
V1 = H3.*c.Rv;
W1 = -c.xiGrid.*(HX3.*c.Ru+HY3.*c.Rv);
gammaN21 = -H3.*(N20+c.xiGrid.*N2Xi);
gammaW1 = W1-H3.*c.Rwh;
gammaN2W1 = gammaN21.*c.Rwh+N20.*W1;
Dxp1 = c.xiGrid.*HX3.*c.RPXi;
Dyp1 = c.xiGrid.*HY3.*c.RPXi;
Q1 = H3.*(c.RvX-c.RuY-c.wvt.f*c.RetaXi) ...
    +HX3.*(c.Rv+c.xiGrid.*c.RvXi)-HY3.*(c.Ru+c.xiGrid.*c.RuXi);
bottomTendencyFull = c.phase'*(c.hX.*c.uBottom+c.hY.*c.vBottom)/c.nXY;
[Efull,Jfull] = firstOrderPhysicalForms(c,H3,N20,gammaN21,U1,V1,W1);
coefficients = struct("gamma", -H3,"gammaN2",gammaN21,"U",U1,"V",V1, ...
    "Wphysical",W1,"gammaWphysical",gammaW1, ...
    "gammaN2Wphysical",gammaN2W1,"Dxp",Dxp1,"Dyp",Dyp1,"Q",Q1, ...
    "apvWeight",-H3.*c.volumeWeight,"bottomTendencyFull",bottomTendencyFull, ...
    "Efull",Efull,"Jfull",Jfull,"isFirstOrder",true);
end

function [S,F,matrices] = descriptorMatrices(c,a)
w = c.volumeWeight;
Mu = c.TF'*(w.*a.U);
Mv = c.TF'*(w.*a.V);
Meta = c.TH'*(w.*(a.gammaN2.*c.Reta));
Au = c.wvt.f*c.TF'*(w.*a.V);
Av = -c.wvt.f*c.TF'*(w.*a.U);
Aw = -c.TH'*(w.*(a.gammaN2.*c.Reta));
Gu = -(c.TF'*(w.*a.Dxp))/c.wvt.rho0;
Gv = -(c.TF'*(w.*a.Dyp))/c.wvt.rho0;
if isfield(a,"isFirstOrder")
    Gw = zeros(c.nK*c.nH,c.nPressure);
else
    Gw = -(c.TH'*(w.*c.RPXi))/c.wvt.rho0;
end
Geta = zeros(c.nK*c.nH,c.nPressure);
G = [Gu;Gv;Gw;Geta];
Mw = c.TH'*(w.*a.gammaWphysical);
Aeta = c.TH'*(w.*a.gammaN2Wphysical);
M = [Mu;Mv;Mw;Meta];
A = [Au;Av;Aw;Aeta];
if isfield(a,"isFirstOrder")
    continuity = zeros(size(c.continuity));
else
    continuity = c.continuity;
end
S = [M -G;continuity zeros(size(continuity,1),c.nPressure)];
F = [A*c.N;zeros(size(continuity,1),size(c.N,2))];
matrices = struct("M",M,"A",A,"G",G,"S",S,"F",F);
end

function [E,J] = physicalForms(c,gamma,gammaN2,U,V,Wphysical)
w = c.volumeWeight;
E = c.wvt.rho0*(c.Ru'*(w.*U)+c.Rv'*(w.*V) ...
    +Wphysical'*(w.*(gamma.*Wphysical))+c.Reta'*(w.*(gammaN2.*c.Reta)));
exchange = c.wvt.f*(c.Ru'*(w.*V)-c.Rv'*(w.*U)) ...
    +c.Reta'*(w.*(gammaN2.*Wphysical))-Wphysical'*(w.*(gammaN2.*c.Reta));
J = c.wvt.rho0*exchange;
end

function [E1,J1] = firstOrderPhysicalForms(c,H3,N20,gammaN21,U1,V1,W1)
w = c.volumeWeight;
verticalEnergy = W1'*(w.*c.Rwh)+c.Rwh'*(w.*W1)-c.Rwh'*(w.*(H3.*c.Rwh));
E1 = c.wvt.rho0*(c.Ru'*(w.*U1)+c.Rv'*(w.*V1)+verticalEnergy ...
    +c.Reta'*(w.*(gammaN21.*c.Reta)));
horizontalExchange = c.wvt.f*(c.Ru'*(w.*V1)-c.Rv'*(w.*U1));
verticalExchange = c.Reta'*(w.*(gammaN21.*c.Rwh+N20.*W1)) ...
    -(gammaN21.*c.Rwh+N20.*W1)'*(w.*c.Reta);
J1 = c.wvt.rho0*(horizontalExchange+verticalExchange);
end

function [Y,diagnostics] = solveGaugedSystem(S,F,gauge)
[U,Sigma,~] = svd(S,"econ");
singularValues = diag(Sigma);
tolerance = max(size(S))*eps(max(singularValues))*100;
nullity = nnz(singularValues <= tolerance);
if nullity ~= 1
    error("WVTerrainEnergyGalerkin:UnexpectedPrimitivePressureNullity", ...
        "The ungauged global primitive saddle has nullity %d instead of one.",nullity)
end
leftNull = U(:,end);
bordered = [S leftNull;gauge' 0];
solution = bordered\[F;zeros(1,size(F,2))];
Y = solution(1:end-1,:);
diagnostics = struct("pressureNullity",nullity,"singularValues",singularValues, ...
    "nullTolerance",tolerance,"leftNullVector",leftNull, ...
    "gaugedSaddleResidual",norm(S*Y-F,"fro")/max(norm(F,"fro"),realmin), ...
    "pressureGaugeDefect",norm(gauge'*Y,"fro")/max(norm(Y,"fro"),realmin), ...
    "borderedReciprocalConditionNumber",rcond(bordered));
end

function tangent = centeredTangent(plus,minus,step)
names = ["L","E","J","Q","Z","R","pressureMap","primitiveTendency"];
tangent = struct("step",step);
for name = names
    tangent.(name) = (plus.(name)-minus.(name))/(2*step);
end
end

function diagnostics = tangentAgreementDiagnostics(analytic,centered)
names = ["L","E","J","Q","Z","R","pressureMap","primitiveTendency"];
relative = struct;
maximum = 0;
for name = names
    values = zeros(numel(centered),1);
    for iStep = 1:numel(centered)
        values(iStep) = norm(centered(iStep).(name)-analytic.(name),"fro") ...
            /max(norm(analytic.(name),"fro"),realmin);
    end
    relative.(name) = values;
    maximum = max(maximum,values(end));
end
orders = log2(max(relative.L(1:end-1),realmin)./max(relative.L(2:end),realmin));
diagnostics = struct("relativeDefects",relative,"observedOrders",orders, ...
    "maximumRelativeDefect",maximum);
end

function diagnostics = compatibilityDiagnostics(c,flat,first)
E0 = flat.E;
L0 = flat.L;
J0 = flat.J;
Q0 = flat.Q;
Z0 = flat.Z;
B = flat.B;
weakTerms = {E0*first.L,first.E*L0,-first.J};
energyTerms = {first.L'*E0,E0*first.L,L0'*first.E,first.E*L0};
apvTerms = {Q0*first.L,first.Q*L0};
bottomTerms = {B*first.L,-first.R};
enstrophyTerms = {first.L'*Z0,Z0*first.L,L0'*first.Z,first.Z*L0};
C = sparse((1:size(L0,1))',c.coordinateConjugateIndex,1,size(L0,1),size(L0,1));
conjugacyTerms = {first.L*C,-C*conj(first.L)};
diagnostics = struct( ...
    "weakEvolutionDefect",identityDefect(weakTerms,norm(E0,"fro")*norm(first.L,"fro") ...
        +norm(first.E,"fro")*norm(L0,"fro")+norm(first.J,"fro")), ...
    "energyDefect",identityDefect(energyTerms,2*norm(E0,"fro")*norm(first.L,"fro") ...
        +2*norm(first.E,"fro")*norm(L0,"fro")), ...
    "apvDefect",identityDefect(apvTerms,norm(Q0,"fro")*norm(first.L,"fro") ...
        +norm(first.Q,"fro")*norm(L0,"fro")), ...
    "bottomDefect",identityDefect(bottomTerms,norm(B,"fro")*norm(first.L,"fro") ...
        +norm(first.R,"fro")), ...
    "enstrophyDefect",identityDefect(enstrophyTerms,2*norm(Z0,"fro")*norm(first.L,"fro") ...
        +2*norm(first.Z,"fro")*norm(L0,"fro")), ...
    "conjugacyDefect",identityDefect(conjugacyTerms,2*norm(first.L,"fro")));
diagnostics.apvCancellation = struct("tendencyCorrection",apvTerms{1}, ...
    "mapCorrection",apvTerms{2},"total",apvTerms{1}+apvTerms{2}, ...
    "relativeDefect",diagnostics.apvDefect);
diagnostics.flatWeakDefect = productDefect(E0*L0-J0,{E0*L0,J0});
end

function diagnostics = strongAPVDiagnostics(c,flat,first)
H3 = repmat(c.H,c.nZ,1);
HX3 = repmat(c.HX,c.nZ,1);
HY3 = repmat(c.HY,c.nZ,1);
N2Xi = kron(c.N2Xi,ones(c.nXY,1));
state = c.N;
u = c.Ru*state;
v = c.Rv*state;
eta = c.Reta*state;
pX = c.RPX*flat.pressureMap;
pY = c.RPY*flat.pressureMap;
pXi = c.RPXi*flat.pressureMap;
uDot0 = c.Ru*flat.primitiveTendency;
vDot0 = c.Rv*flat.primitiveTendency;
uDot1 = (H3.*pX-c.xiGrid.*HX3.*pXi)/c.wvt.rho0;
vDot1 = (H3.*pY-c.xiGrid.*HY3.*pXi)/c.wvt.rho0;
wHatDot1 = c.xiGrid.*(HX3.*uDot0+HY3.*vDot0) ...
    -H3.*pXi/c.wvt.rho0+H3.*c.xiGrid.*N2Xi.*eta;
etaDot1 = -c.xiGrid.*(HX3.*u+HY3.*v);
etaDot1Xi = -(HX3.*(u+c.xiGrid.*(c.RuXi*state)) ...
    +HY3.*(v+c.xiGrid.*(c.RvXi*state)));
q0Tendency = horizontalDerivative(vDot1,c,"x") ...
    -horizontalDerivative(uDot1,c,"y")-c.wvt.f*etaDot1Xi;
q1Tendency = first.Q*flat.L;
total = q0Tendency+q1Tendency;
scale = norm(q0Tendency,"fro")+norm(q1Tendency,"fro");
diagnostics = struct("q0Tendency",q0Tendency,"q1Tendency",q1Tendency, ...
    "total",total,"identityDefect",norm(total,"fro")/max(scale,realmin), ...
    "wHatTendency",wHatDot1,"etaTendency",etaDot1);
end

function diagnostics = fourierDiagnostics(c,first,apvCancellation)
terrainModes = terrainFourierModes(c);
allowed = false(c.nK);
for iOut = 1:c.nK
    for iIn = 1:c.nK
        difference = [c.problem.horizontalLayout.kMode(iOut)-c.problem.horizontalLayout.kMode(iIn), ...
            c.problem.horizontalLayout.lMode(iOut)-c.problem.horizontalLayout.lMode(iIn)];
        allowed(iOut,iIn) = any(all(terrainModes == difference,2));
    end
end
allowedCoordinates = false(size(first.L));
for iOut = 1:c.nK
    for iIn = 1:c.nK
        if allowed(iOut,iIn)
            allowedCoordinates(c.layout.admissibleRanges{iOut},c.layout.admissibleRanges{iIn}) = true;
        end
    end
end
disallowed = first.L;
disallowed(allowedCoordinates) = 0;
leakage = norm(disallowed,"fro")/max(norm(first.L,"fro"),realmin);
interiorModes = false(c.nK,1);
for iIn = 1:c.nK
    destinations = [c.problem.horizontalLayout.kMode(iIn)+terrainModes(:,1), ...
        c.problem.horizontalLayout.lMode(iIn)+terrainModes(:,2)];
    retained = false(size(destinations,1),1);
    for iDestination = 1:size(destinations,1)
        retained(iDestination) = any(c.problem.horizontalLayout.kMode == destinations(iDestination,1) ...
            & c.problem.horizontalLayout.lMode == destinations(iDestination,2));
    end
    interiorModes(iIn) = all(retained);
end
interiorColumns = horzcat(c.layout.admissibleRanges{interiorModes});
edgeColumns = setdiff(1:size(first.L,2),interiorColumns);
diagnostics = struct("terrainModes",terrainModes,"allowedHorizontalBlocks",allowed, ...
    "couplingLeakage",leakage,"interiorHorizontalModes",interiorModes, ...
    "numberOfInteriorHorizontalModes",nnz(interiorModes), ...
    "interiorAPVDefect",cancellationDefect(apvCancellation,interiorColumns), ...
    "edgeAPVDefect",cancellationDefect(apvCancellation,edgeColumns));
end

function value = cancellationDefect(cancellation,columns)
if isempty(columns)
    value = NaN;
    return
end
first = cancellation.tendencyCorrection(:,columns);
second = cancellation.mapCorrection(:,columns);
value = norm(first+second,"fro")/max(norm(first,"fro")+norm(second,"fro"),realmin);
end

function diagnostics = flatReferenceDiagnostics(c,flat)
N2Variation = max(c.N20)-min(c.N20);
modeOneDispersionDefect = NaN;
if N2Variation <= 1e-12*max(c.N20)
    defects = zeros(c.nK,1);
    for iK = 1:c.nK
        kappa = hypot(c.problem.horizontalLayout.k(iK),c.problem.horizontalLayout.l(iK));
        if kappa == 0
            defects(iK) = 0;
            continue
        end
        m = pi/c.wvt.Lz;
        omega = sqrt((c.N20(1)*kappa^2+c.wvt.f^2*m^2)/(kappa^2+m^2));
        rows = c.layout.admissibleRanges{iK};
        eigenvalues = eig(flat.L(rows,rows));
        defects(iK) = min(abs(abs(imag(eigenvalues))-omega))/omega;
    end
    modeOneDispersionDefect = max(defects);
end
diagnostics = struct("modeOneDispersionDefect",modeOneDispersionDefect, ...
    "continuity",c.continuityDiagnostics,"phaseGramDefect",c.phaseGramDefect);
end

function diagnostics = structureDiagnostics(~,flat,first)
hermitian = [matrixDefect(flat.E-flat.E',flat.E),matrixDefect(first.E-first.E',first.E)];
skew = [matrixDefect(flat.J+flat.J',flat.J),matrixDefect(first.J+first.J',first.J)];
diagnostics = struct("energyHermitianDefects",hermitian, ...
    "exchangeSkewHermitianDefects",skew, ...
    "maximumHermitianDefect",max(hermitian), ...
    "maximumSkewHermitianDefect",max(skew));
end

function values = horizontalDerivative(values,c,direction)
Nx = c.layout.oversampledSize(1);
Ny = c.layout.oversampledSize(2);
nColumns = size(values,2);
array = reshape(values,[Nx Ny c.nZ nColumns]);
if direction == "x"
    k = 2*pi*[0:floor((Nx-1)/2) -floor(Nx/2):-1]'/c.wvt.Lx;
    multiplier = reshape(1i*k,[Nx 1 1 1]);
    dimension = 1;
else
    l = 2*pi*[0:floor((Ny-1)/2) -floor(Ny/2):-1]/c.wvt.Ly;
    multiplier = reshape(1i*l,[1 Ny 1 1]);
    dimension = 2;
end
values = reshape(ifft(fft(array,[],dimension).*multiplier,[],dimension),size(values));
end

function [xDerivative,yDerivative] = horizontalDerivatives(values,Lx,Ly)
[Nx,Ny] = size(values);
k = 2*pi*[0:floor((Nx-1)/2) -floor(Nx/2):-1]'/Lx;
l = 2*pi*[0:floor((Ny-1)/2) -floor(Ny/2):-1]/Ly;
spectrum = fft2(values);
xDerivative = real(ifft2((1i*k).*spectrum));
yDerivative = real(ifft2(spectrum.*(1i*l)));
end

function [N,dimensions,ranges,diagnostics] = admissibleBasis(problem,continuity,nXBlock,nF)
nK = height(problem.horizontalLayout);
horizontalConjugateIndex = conjugateHorizontalIndex(problem);
local = cell(nK,1);
rankValues = zeros(nK,1);
for iK = 1:nK
    rows = (iK-1)*nF+(1:nF);
    cols = (iK-1)*nXBlock+(1:nXBlock);
    if isempty(local{iK})
        local{iK} = null(continuity(rows,cols));
    end
    partner = horizontalConjugateIndex(iK);
    if partner ~= iK
        local{partner} = conj(local{iK});
    else
        local{iK} = real(local{iK});
    end
    rankValues(iK) = rank(continuity(rows,cols));
end
dimensions = cellfun(@(value)size(value,2),local);
ranges = cell(nK,1);
first = 1;
for iK = 1:nK
    ranges{iK} = first:first+dimensions(iK)-1;
    first = first+dimensions(iK);
end
N = zeros(nK*nXBlock,sum(dimensions));
for iK = 1:nK
    rows = (iK-1)*nXBlock+(1:nXBlock);
    cols = ranges{iK};
    N(rows,cols) = local{iK};
end
diagnostics = struct("rankByHorizontalMode",rankValues, ...
    "expectedRank",nF,"nullspaceDefect",norm(continuity*N,"fro") ...
    /max(norm(continuity,"fro"),realmin),"orthogonalityDefect", ...
    norm(N'*N-eye(size(N,2)),"fro")/sqrt(size(N,2)));
end

function conjugateIndex = conjugateHorizontalIndex(problem)
nK = height(problem.horizontalLayout);
conjugateIndex = zeros(nK,1);
for iK = 1:nK
    partner = find(problem.horizontalLayout.kMode == -problem.horizontalLayout.kMode(iK) ...
        & problem.horizontalLayout.lMode == -problem.horizontalLayout.lMode(iK),1);
    if isempty(partner)
        error("WVTerrainEnergyGalerkin:IncompletePrimitiveHorizontalConjugacy", ...
            "The retained horizontal layout is missing a Fourier conjugate.")
    end
    conjugateIndex(iK) = partner;
end
end

function modes = terrainFourierModes(c)
Nx = c.wvt.Nx;
Ny = c.wvt.Ny;
spectrum = fft2(c.problem.topographicHeight)/(Nx*Ny);
[kMode,lMode] = ndgrid(c.wvt.kMode_dft,c.wvt.lMode_dft);
mask = abs(spectrum) > 1e-12*max(abs(spectrum),[],"all");
modes = [kMode(mask) lMode(mask)];
end

function value = identityDefect(terms,naturalScale)
residual = zeros(size(terms{1}));
scale = 0;
for iTerm = 1:numel(terms)
    residual = residual+terms{iTerm};
    scale = scale+norm(terms{iTerm},"fro");
end
scale = max(scale,100*eps*naturalScale);
value = norm(residual,"fro")/max(scale,realmin);
end

function value = productDefect(residual,terms)
scale = 0;
if iscell(terms)
    for iTerm = 1:numel(terms)
        scale = scale+norm(terms{iTerm},"fro");
    end
else
    scale = norm(terms,"fro");
end
value = norm(residual,"fro")/max(scale,realmin);
end

function value = matrixDefect(residual,reference)
value = norm(residual,"fro")/max(norm(reference,"fro"),realmin);
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
[F,Fr] = legendreValues(r,degree);
Fxi = (2/D)*Fr;
G = (1-r.^2).*F(:,1:degree);
Gxi = (2/D)*(-2*r.*F(:,1:degree)+(1-r.^2).*Fr(:,1:degree));
chi = (1-r)/2;
H = [G chi];
Hxi = [Gxi -ones(size(r))/D];
[Pressure,PressureR] = legendreValues(r,degree+1);
PressureXi = (2/D)*PressureR;
[Fendpoint,~] = legendreValues([-1;1],degree);
Gendpoint = (1-[-1;1].^2).*Fendpoint(:,1:degree);
Hendpoint = [Gendpoint [1;0]];
spaces = struct("F",F,"Fxi",Fxi,"G",G,"Gxi",Gxi,"H",H,"Hxi",Hxi, ...
    "Pressure",Pressure,"PressureXi",PressureXi,"Fendpoint",Fendpoint, ...
    "Hendpoint",Hendpoint);
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
