function audit = buildAPVCompatiblePrimitiveAudit(problem,iK,bottomSlope,polynomialDegree,quadratureOrder)
% Compare local Branch-P primitive formulations and classify APV compatibility.

h = problem.topographicHeight;
if max(abs(h),[],"all") > 20*eps(max(1,max(abs(h),[],"all")))
    error("WVTerrainEnergyGalerkin:APVPrimitiveAuditRequiresFlatReference", ...
        "The local APV-compatible primitive audit requires topographicHeight=0.")
end

wvt = problem.originatingTransform;
k = problem.horizontalLayout.k(iK);
l = problem.horizontalLayout.l(iK);
quadrature = legendreQuadrature(quadratureOrder,wvt.Lz);
spaces = polynomialSpaces(polynomialDegree,quadrature,wvt.Lz);
geometry = buildGeometry(wvt,k,l,bottomSlope,spaces,quadrature);

energyWeak = buildCandidate(geometry,"energy-weak");
vorticityDivergence = buildCandidate(geometry,"vorticity-divergence");
explicitAPV = explicitAPVCoordinates(vorticityDivergence,geometry);
flat = buildCandidate(buildGeometry(wvt,k,l,[0 0],spaces,quadrature),"energy-weak");
flat.diagnostics.modeOneDispersionDefect = flatModeOneDispersionDefect(wvt,k,l,flat.L);

crossSlopeWavenumber = k*bottomSlope(2)-l*bottomSlope(1);
sourceCoefficient = 1i*crossSlopeWavenumber/(wvt.rho0*wvt.Lz);
analyticSource = sourceCoefficient*energyWeak.pressureField;
algebraicIdentityDefect = analyticFrozenAPVIdentity(geometry,sourceCoefficient);
sourceScale = norm(energyWeak.apvTendency,"fro")+norm(analyticSource,"fro");
if sourceScale <= 100*eps
    sourceComparisonDefect = 0;
else
    sourceComparisonDefect = norm(energyWeak.apvTendency-analyticSource,"fro")/sourceScale;
end

obstructionTolerance = 100*eps*max([abs(k)*norm(bottomSlope),abs(l)*norm(bottomSlope),realmin]);
hasCrossSlopeObstruction = abs(crossSlopeWavenumber) > obstructionTolerance;
requiredTolerance = 1e-11;
energyWeakPasses = energyWeak.diagnostics.energyDefect <= requiredTolerance ...
    && energyWeak.diagnostics.bottomDefect <= requiredTolerance ...
    && energyWeak.diagnostics.weakEvolutionDefect <= requiredTolerance ...
    && energyWeak.diagnostics.apvDefect <= requiredTolerance;
if hasCrossSlopeObstruction
    status = "mathematical-blocker";
    isCompatible = false;
    diagnosis = "The frozen Branch-P equations generate APV for a cross-slope Fourier component before discretization.";
else
    status = "compatible-aligned-special-case";
    isCompatible = energyWeakPasses;
    diagnosis = "The analytic APV source vanishes because the horizontal wavenumber and local slope are parallel.";
end

audit = struct;
audit.status = status;
audit.isCompatible = isCompatible;
audit.scope = "local-frozen-constant-slope-only";
audit.horizontalIndex = iK;
audit.horizontalMode = [problem.horizontalLayout.kMode(iK) problem.horizontalLayout.lMode(iK)];
audit.horizontalWavenumber = [k l];
audit.bottomSlope = bottomSlope;
audit.polynomialDegree = polynomialDegree;
audit.quadratureOrder = quadratureOrder;
audit.spaces = spaces;
audit.quadrature = quadrature;
audit.commutingDiagnostics = commutingDiagnostics(geometry);
audit.energyWeak = energyWeak;
audit.vorticityDivergence = vorticityDivergence;
audit.explicitAPV = explicitAPV;
audit.flatReference = flat;
audit.analyticObstruction = struct( ...
    "crossSlopeWavenumber",crossSlopeWavenumber, ...
    "sourceCoefficient",sourceCoefficient, ...
    "sourceFormula","q_t=i(k s_y-l s_x)p/(rho_0 D)", ...
    "sourceMatrix",analyticSource, ...
    "algebraicIdentityDefect",algebraicIdentityDefect, ...
    "sourceComparisonDefect",sourceComparisonDefect, ...
    "hasCrossSlopeObstruction",hasCrossSlopeObstruction);
audit.candidateVerdicts = struct( ...
    "compatibleMixedSpace","preserves energy and bottom evolution but reproduces the nonzero analytic APV source", ...
    "vorticityDivergence","enforces APV and bottom evolution but changes the Branch-P weak equations and loses physical-energy conservation", ...
    "explicitAPV","makes q_t=0 explicit but is a coordinate form of the same non-energy-conserving APV-enforced system");
audit.diagnosis = diagnosis;
audit.nextScope = "global-small-terrain-or-covariant-local-formulation";
end

function geometry = buildGeometry(wvt,k,l,slope,spaces,quadrature)
rho0 = wvt.rho0;
f = wvt.f;
D = wvt.Lz;
xi = quadrature.xi;
W = diag(quadrature.weight);
N2 = wvt.N2Function(xi);
if isscalar(N2)
    N2 = repmat(N2,numel(xi),1);
end
N2 = N2(:);
WN = W*diag(N2);
F = spaces.F;
Fxi = spaces.Fxi;
Fxxi = spaces.Fxxi;
G = spaces.G;
Gxi = spaces.Gxi;
H = spaces.H;
Hxi = spaces.Hxi;
nF = size(F,2);
nG = size(G,2);
nX = 2*nF+nG+nF;
iu = 1:nF;
iv = nF+(1:nF);
iw = 2*nF+(1:nG);
ieta = 2*nF+nG+(1:nF);

Ru = zeros(numel(xi),nX);
Rv = Ru;
Rwh = Ru;
Reta = Ru;
RwhXi = Ru;
Ru(:,iu) = F;
Rv(:,iv) = F;
Rwh(:,iw) = G;
RwhXi(:,iw) = Gxi;
Reta(:,ieta) = H;
sx = slope(1);
sy = slope(2);
X = diag(xi);
Wphysical = Rwh-X*(sx*Ru+sy*Rv)/D;
Qfull = 1i*k*Rv-1i*l*Ru+(sx/D)*(Rv+X*(Fxi*selector(nX,iv)')) ...
    -(sy/D)*(Ru+X*(Fxi*selector(nX,iu)'))-f*(Hxi*selector(nX,ieta)');
continuityFull = F'*W*(1i*k*Ru+1i*l*Rv+RwhXi);
[continuity,N,constraintRank] = independentConstraints(continuityFull);

Efull = rho0*(Ru'*W*Ru+Rv'*W*Rv+Wphysical'*W*Wphysical+Reta'*WN*Reta);
Jfull = rho0*(f*(Ru'*W*Rv-Rv'*W*Ru)+Reta'*WN*Wphysical-Wphysical'*WN*Reta);
bottomValueFull = zeros(1,nX);
bottomValueFull(ieta(end)) = 1;
uBottomFull = zeros(1,nX);
vBottomFull = zeros(1,nX);
uBottomFull(iu) = spaces.Fendpoint(1,:);
vBottomFull(iv) = spaces.Fendpoint(1,:);
bottomRhsFull = sx*uBottomFull+sy*vBottomFull;

geometry = struct("rho0",rho0,"f",f,"D",D,"k",k,"l",l,"slope",slope, ...
    "xi",xi,"W",W,"WN",WN,"N2",N2,"F",F,"Fxi",Fxi,"Fxxi",Fxxi,"G",G, ...
    "Gxi",Gxi,"H",H,"Hxi",Hxi,"X",X,"nF",nF,"nG",nG,"nX",nX, ...
    "iu",iu,"iv",iv,"iw",iw,"ieta",ieta,"Ru",Ru,"Rv",Rv,"Rwh",Rwh, ...
    "Reta",Reta,"RwhXi",RwhXi,"Wphysical",Wphysical,"Qfull",Qfull, ...
    "continuity",continuity,"N",N,"constraintRank",constraintRank, ...
    "Efull",Efull,"Jfull",Jfull,"bottomValueFull",bottomValueFull, ...
    "bottomRhsFull",bottomRhsFull);
end

function candidate = buildCandidate(g,kind)
if kind == "energy-weak"
    [P,Pr] = legendreValues(g.xi*2/g.D+1,g.nF);
    pressure = P;
    pressureXi = (2/g.D)*Pr;
    nPressure = g.nF+1;
    pressureX = 1i*g.k*pressure+g.X*(g.slope(1)*pressureXi)/g.D;
    pressureY = 1i*g.l*pressure+g.X*(g.slope(2)*pressureXi)/g.D;
    M = [g.F'*g.W*g.Ru;g.F'*g.W*g.Rv;g.H'*g.W*g.Wphysical;g.H'*g.WN*g.Reta];
    A = [g.f*g.F'*g.W*g.Rv;-g.f*g.F'*g.W*g.Ru; ...
        -g.H'*g.WN*g.Reta;g.H'*g.WN*g.Wphysical];
    Gp = [-(g.F'*g.W*pressureX)/g.rho0;-(g.F'*g.W*pressureY)/g.rho0; ...
        -(g.H'*g.W*pressureXi)/g.rho0;zeros(g.nF,nPressure)];
else
    pressure = g.F;
    pressureXi = g.Fxi;
    nPressure = g.nF;
    pressureX = 1i*g.k*pressure+g.X*(g.slope(1)*pressureXi)/g.D;
    pressureY = 1i*g.l*pressure+g.X*(g.slope(2)*pressureXi)/g.D;
    Mu = g.F'*g.W*g.Ru;
    Mv = g.F'*g.W*g.Rv;
    Au = g.f*g.F'*g.W*g.Rv;
    Av = -g.f*g.F'*g.W*g.Ru;
    Gpu = -(g.F'*g.W*pressureX)/g.rho0;
    Gpv = -(g.F'*g.W*pressureY)/g.rho0;
    M = [1i*g.k*Mu+1i*g.l*Mv;g.F'*g.W*g.Qfull; ...
        g.G'*g.W*g.Wphysical;g.H'*g.WN*g.Reta];
    A = [1i*g.k*Au+1i*g.l*Av;zeros(g.nF,g.nX); ...
        -g.G'*g.WN*g.Reta;g.H'*g.WN*g.Wphysical];
    Gp = [1i*g.k*Gpu+1i*g.l*Gpv;zeros(g.nF,nPressure); ...
        -(g.G'*g.W*pressureXi)/g.rho0;zeros(g.nF,nPressure)];
end

S = [M -Gp;g.continuity zeros(size(g.continuity,1),nPressure)];
rhs = [A*g.N;zeros(size(g.continuity,1),size(g.N,2))];
solution = S\rhs;
xdot = solution(1:g.nX,:);
pressureMap = solution(g.nX+1:end,:);
L = g.N'*xdot;
E = g.N'*g.Efull*g.N;
J = g.N'*g.Jfull*g.N;
Q = g.Qfull*g.N;
Z = Q'*g.W*Q;
B = g.bottomValueFull*g.N;
R = g.bottomRhsFull*g.N;
apvTendency = Q*L;

candidate = struct;
candidate.name = kind;
candidate.L = L;
candidate.E = E;
candidate.J = J;
candidate.Q = Q;
candidate.Z = Z;
candidate.B = B;
candidate.R = R;
candidate.pressureMap = pressureMap;
candidate.pressureField = pressure*pressureMap;
candidate.apvTendency = apvTendency;
candidate.reconstruction = struct("u",g.Ru*g.N,"v",g.Rv*g.N, ...
    "wHat",g.Rwh*g.N,"w",g.Wphysical*g.N,"eta",g.Reta*g.N);
candidate.diagnostics = candidateDiagnostics(candidate,S,solution,rhs,g);
end

function diagnostics = candidateDiagnostics(candidate,S,solution,rhs,g)
energyScale = max(norm(candidate.E,"fro")*norm(candidate.L,"fro"),realmin);
apvScale = max(norm(candidate.Q,"fro")*norm(candidate.L,"fro"),realmin);
enstrophyScale = max(norm(candidate.Z,"fro")*norm(candidate.L,"fro"),realmin);
bottomScale = max(norm(candidate.B,"fro")*norm(candidate.L,"fro")+norm(candidate.R,"fro"),realmin);
diagnostics = struct( ...
    "weakEvolutionDefect",norm(candidate.E*candidate.L-candidate.J,"fro")/energyScale, ...
    "energyDefect",norm(candidate.L'*candidate.E+candidate.E*candidate.L,"fro")/energyScale, ...
    "apvDefect",norm(candidate.apvTendency,"fro")/apvScale, ...
    "enstrophyDefect",norm(candidate.L'*candidate.Z+candidate.Z*candidate.L,"fro")/enstrophyScale, ...
    "bottomDefect",norm(candidate.B*candidate.L-candidate.R,"fro")/bottomScale, ...
    "saddleResidual",norm(S*solution-rhs,"fro")/max(norm(rhs,"fro"),realmin), ...
    "continuityTangencyDefect",norm(g.continuity*(g.N*candidate.L),"fro") ...
    /max(norm(g.continuity,"fro")*norm(candidate.L,"fro"),realmin), ...
    "minimumEnergyEigenvalue",min(real(eig((candidate.E+candidate.E')/2))), ...
    "maximumRealEigenvalue",max(abs(real(eig(candidate.L)))));
end

function coordinates = explicitAPVCoordinates(candidate,g)
MF = g.F'*g.W*g.F;
qCoefficientFull = MF\(g.F'*g.W*g.Qfull);
zetaCoefficientFull = zeros(g.nF,g.nX);
zetaCoefficientFull(:,g.iu) = -1i*g.l*eye(g.nF);
zetaCoefficientFull(:,g.iv) = 1i*g.k*eye(g.nF);
verticalVelocityFull = selector(g.nX,g.iw)';
displacementFull = selector(g.nX,g.ieta)';
apvMap = [qCoefficientFull;verticalVelocityFull;displacementFull]*g.N;
vorticityDivergenceMap = [zetaCoefficientFull;verticalVelocityFull;displacementFull]*g.N;
apvRank = rankWithTolerance(apvMap);
vorticityDivergenceRank = rankWithTolerance(vorticityDivergenceMap);
if apvRank == size(apvMap,1)
    Lapv = apvMap*candidate.L/apvMap;
    Eapv = apvMap'\candidate.E/apvMap;
    qRowDefect = norm(Lapv(1:g.nF,:),"fro")/max(norm(Lapv,"fro"),realmin);
    energyDefect = norm(Lapv'*Eapv+Eapv*Lapv,"fro") ...
        /max(norm(Eapv,"fro")*norm(Lapv,"fro"),realmin);
else
    Lapv = [];
    Eapv = [];
    qRowDefect = Inf;
    energyDefect = Inf;
end
coordinates = struct( ...
    "apvMap",apvMap, ...
    "vorticityDivergenceMap",vorticityDivergenceMap, ...
    "apvCoordinateRank",apvRank, ...
    "vorticityDivergenceRank",vorticityDivergenceRank, ...
    "coordinateDimension",size(apvMap,1), ...
    "apvMapConditionNumber",cond(apvMap), ...
    "vorticityDivergenceConditionNumber",cond(vorticityDivergenceMap), ...
    "generator",Lapv, ...
    "energyMatrix",Eapv, ...
    "qRowTendencyDefect",qRowDefect, ...
    "energyDefect",energyDefect);
end

function diagnostics = commutingDiagnostics(g)
MF = g.F'*g.W*g.F;
MH = g.H'*g.W*g.H;
DHtoF = MF\(g.F'*g.W*g.Hxi);
XFtoH = MH\(g.H'*g.W*g.X*g.F);
XDFtoF = MF\(g.F'*g.W*g.X*g.Fxi);
identityDefect = norm(DHtoF*XFtoH-XDFtoF-eye(g.nF),"fro") ...
    /max(norm(eye(g.nF),"fro"),realmin);
diagnostics = struct( ...
    "productRuleDefect",identityDefect, ...
    "continuityRank",g.constraintRank, ...
    "expectedContinuityRank",g.nF, ...
    "admissibleDimension",size(g.N,2), ...
    "expectedAdmissibleDimension",3*g.nF-1);
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
[P,Pr,Prr] = legendreValues(r,degree);
F = P;
Fxi = (2/D)*Pr;
Fxxi = (2/D)^2*Prr;
G = (1-r.^2).*P(:,1:degree);
Gxi = (2/D)*(-2*r.*P(:,1:degree)+(1-r.^2).*Pr(:,1:degree));
chi = (1-r)/2;
chixi = -ones(size(r))/D;
H = [G chi];
Hxi = [Gxi chixi];
[Pendpoint,~] = legendreValues([-1;1],degree);
spaces = struct("F",F,"Fxi",Fxi,"Fxxi",Fxxi,"G",G,"Gxi",Gxi,"H",H,"Hxi",Hxi, ...
    "chi",chi,"chixi",chixi,"Fendpoint",Pendpoint, ...
    "bottomValue",H(1,:),"surfaceValue",H(end,:));
end

function [P,Pr,Prr] = legendreValues(r,degree)
P = zeros(numel(r),degree+1);
Pr = zeros(numel(r),degree+1);
Prr = zeros(numel(r),degree+1);
P(:,1) = 1;
if degree == 0
    return
end
P(:,2) = r;
Pr(:,2) = 1;
for n = 1:degree-1
    P(:,n+2) = ((2*n+1)*r.*P(:,n+1)-n*P(:,n))/(n+1);
    Pr(:,n+2) = ((2*n+1)*(P(:,n+1)+r.*Pr(:,n+1))-n*Pr(:,n))/(n+1);
    Prr(:,n+2) = ((2*n+1)*(2*Pr(:,n+1)+r.*Prr(:,n+1))-n*Prr(:,n))/(n+1);
end
end

function [C,N,rankC] = independentConstraints(Cfull)
[U,S,V] = svd(Cfull);
singularValues = diag(S);
tolerance = max(size(Cfull))*eps(max(singularValues,1));
rankC = nnz(singularValues > tolerance);
C = U(:,1:rankC)'*Cfull;
N = deterministicColumns(V(:,rankC+1:end));
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

function rankA = rankWithTolerance(A)
rankA = rank(A,max(size(A))*eps(max(norm(A,2),1)));
end

function S = selector(n,indices)
S = zeros(n,numel(indices));
S(indices,:) = eye(numel(indices));
end

function defect = flatModeOneDispersionDefect(wvt,k,l,L)
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

function defect = analyticFrozenAPVIdentity(g,sourceCoefficient)
pressureX = 1i*g.k*g.F+g.X*(g.slope(1)*g.Fxi)/g.D;
pressureY = 1i*g.l*g.F+g.X*(g.slope(2)*g.Fxi)/g.D;
pressureXxi = 1i*g.k*g.Fxi+(g.slope(1)/g.D)*(g.Fxi+g.X*g.Fxxi);
pressureYxi = 1i*g.l*g.Fxi+(g.slope(2)/g.D)*(g.Fxi+g.X*g.Fxxi);
uTendency = -pressureX/g.rho0;
vTendency = -pressureY/g.rho0;
uTendencyXi = -pressureXxi/g.rho0;
vTendencyXi = -pressureYxi/g.rho0;
apvTendency = 1i*g.k*vTendency-1i*g.l*uTendency ...
    +(g.slope(1)/g.D)*(vTendency+g.X*vTendencyXi) ...
    -(g.slope(2)/g.D)*(uTendency+g.X*uTendencyXi);
expected = sourceCoefficient*g.F;
defect = norm(apvTendency-expected,"fro") ...
    /max(norm(apvTendency,"fro")+norm(expected,"fro"),realmin);
end
