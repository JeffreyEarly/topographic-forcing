function space = constructCompleteStationarySpace(c,direction,trustedBounds, ...
    fullDegree,trustedDegree,terrainScale)
% Construct the complete tangent stationary space and bottom complement.

fullRaw = rawStationaryScalarSpace( ...
    c,fullDegree,true(c.nK,1),terrainScale);
[fullScalar,fullNonTangentScalar] = splitByBottomTangency(fullRaw);
[GfullRaw,fullRepresentation] = projectGeostrophicState(c,fullScalar);
[Gfull,fullScalar,fullOrthogonalization] = energyOrthonormalize( ...
    GfullRaw,fullScalar,direction.E);

if isempty(fullNonTangentScalar.coefficientBasis)
    bottomSeed = zeros(size(direction.E,1),0);
    bottomSeedOrthogonalization = struct( ...
        "energyOrthogonalityDefect",0, ...
        "stationaryOrthogonalityDefect",0, ...
        "gramReciprocalConditionNumber",NaN);
    nonTangentRepresentation = emptyRepresentation;
else
    [bottomSeedRaw,nonTangentRepresentation] = ...
        projectGeostrophicState(c,fullNonTangentScalar);
    bottomSeedRaw = bottomSeedRaw-Gfull*(Gfull'*direction.E*bottomSeedRaw);
    [bottomSeed,bottomSeedOrthogonalization] = ...
        orthonormalizeComplement(bottomSeedRaw,Gfull,direction.E);
end

trustedHorizontal = abs(c.horizontalLayout.kMode) <= trustedBounds(1) ...
    & abs(c.horizontalLayout.lMode) <= trustedBounds(2);
trustedRaw = rawStationaryScalarSpace( ...
    c,trustedDegree,trustedHorizontal,terrainScale);
[trustedScalar,trustedNonTangentScalar] = splitByBottomTangency(trustedRaw);
[GtrustedRaw,trustedRepresentation] = ...
    projectGeostrophicState(c,trustedScalar);
[Gtrusted,trustedScalar,trustedOrthogonalization] = ...
    energyOrthonormalize(GtrustedRaw,trustedScalar,direction.E);
if isempty(trustedNonTangentScalar.coefficientBasis)
    trustedBottomSeed = zeros(size(direction.E,1),0);
    trustedBottomSeedOrthogonalization = struct( ...
        "energyOrthogonalityDefect",0, ...
        "stationaryOrthogonalityDefect",0, ...
        "gramReciprocalConditionNumber",NaN);
else
    [trustedBottomSeedRaw,~] = ...
        projectGeostrophicState(c,trustedNonTangentScalar);
    trustedBottomSeedRaw = trustedBottomSeedRaw ...
        -Gtrusted*(Gtrusted'*direction.E*trustedBottomSeedRaw);
    [trustedBottomSeed,trustedBottomSeedOrthogonalization] = ...
        orthonormalizeComplement( ...
        trustedBottomSeedRaw,Gtrusted,direction.E);
end

space = struct;
space.fullBasis = Gfull;
space.trustedBasis = Gtrusted;
space.bottomSeedBasis = bottomSeed;
space.trustedBottomSeedBasis = trustedBottomSeed;
space.fullScalar = fullScalar;
space.trustedScalar = trustedScalar;
space.fullRawScalar = fullRaw;
space.trustedRawScalar = trustedRaw;
space.fullNonTangentScalar = fullNonTangentScalar;
space.trustedNonTangentScalar = trustedNonTangentScalar;
space.fullRepresentation = fullRepresentation;
space.trustedRepresentation = trustedRepresentation;
space.nonTangentRepresentation = nonTangentRepresentation;
space.fullOrthogonalization = fullOrthogonalization;
space.trustedOrthogonalization = trustedOrthogonalization;
space.bottomSeedOrthogonalization = bottomSeedOrthogonalization;
space.trustedBottomSeedOrthogonalization = ...
    trustedBottomSeedOrthogonalization;
space.diagnostics = struct( ...
    "numberOfFullStationaryStates",size(Gfull,2), ...
    "numberOfTrustedStationaryStates",size(Gtrusted,2), ...
    "numberOfBottomSeedStates",size(bottomSeed,2), ...
    "fullStationaryRowDefect",normalizedProduct( ...
        direction.J*Gfull,direction.J,Gfull), ...
    "trustedStationaryRowDefect",normalizedProduct( ...
        direction.J*Gtrusted,direction.J,Gtrusted), ...
    "fullBottomTangencyDefect",normalizedProduct( ...
        direction.R*Gfull,direction.R,Gfull), ...
    "trustedBottomTangencyDefect",normalizedProduct( ...
        direction.R*Gtrusted,direction.R,Gtrusted), ...
    "bottomSeedStationaryOrthogonalityDefect", ...
        bottomSeedOrthogonalization.stationaryOrthogonalityDefect);
end

function raw = rawStationaryScalarSpace(c,degree,horizontalMask,terrainScale)
n = (0:degree)';
topDerivative = n.*(n+1)/c.wvt.Lz;
[~,singularMatrix,V] = svd(topDerivative.');
singularValues = diag(singularMatrix);
tolerance = max(size(topDerivative.'))*eps(max(singularValues));
rankTop = nnz(singularValues > tolerance);
topBasis = canonicalColumnPhases(V(:,rankTop+1:end));
verticalCoefficients = zeros(c.nF,size(topBasis,2));
verticalCoefficients(1:degree+1,:) = topBasis;
S = c.spaces.F*verticalCoefficients;
Sxi = c.spaces.Fxi*verticalCoefficients;
bottomValue = (-1).^(0:degree);
Sb = bottomValue*topBasis;
nVertical = size(S,2);
selectedHorizontal = find(horizontalMask);
nSelected = numel(selectedHorizontal);
nRaw = nSelected*nVertical;
nGrid = c.nXY*c.nZ;
values = zeros(nGrid,nRaw);
valuesXi = zeros(nGrid,nRaw);
valuesX = zeros(nGrid,nRaw);
valuesY = zeros(nGrid,nRaw);
bottomValues = zeros(c.nXY,nRaw);
bottomX = zeros(c.nXY,nRaw);
bottomY = zeros(c.nXY,nRaw);
for iSelected = 1:nSelected
    iK = selectedHorizontal(iSelected);
    columns = (iSelected-1)*nVertical+(1:nVertical);
    values(:,columns) = kron(S,c.phase(:,iK));
    valuesXi(:,columns) = kron(Sxi,c.phase(:,iK));
    valuesX(:,columns) = 1i*c.horizontalLayout.k(iK)*values(:,columns);
    valuesY(:,columns) = 1i*c.horizontalLayout.l(iK)*values(:,columns);
    bottomValues(:,columns) = c.phase(:,iK)*Sb;
    bottomX(:,columns) = 1i*c.horizontalLayout.k(iK)*bottomValues(:,columns);
    bottomY(:,columns) = 1i*c.horizontalLayout.l(iK)*bottomValues(:,columns);
end

tangency = terrainScale*(bottomX.*c.hY-bottomY.*c.hX);
projectedTangency = c.Pxy*tangency;
zeroK = find(c.horizontalLayout.kMode == 0 ...
    & c.horizontalLayout.lMode == 0,1);
zeroSelected = find(selectedHorizontal == zeroK,1);
if isempty(zeroSelected)
    error("WVTerrainEnergyGalerkin:StationaryScalarMissingGauge", ...
        "The stationary scalar support must include the zero horizontal mode.")
end
constantVertical = topBasis\eye(degree+1,1);
gauge = zeros(nRaw,1);
gauge((zeroSelected-1)*nVertical+(1:nVertical)) = constantVertical;
gauge = gauge/norm(gauge);

selectedLayout = c.horizontalLayout(selectedHorizontal,:);
horizontalConjugate = conjugateHorizontalIndex(selectedLayout);
conjugacy = kron(sparse((1:nSelected)',horizontalConjugate,1, ...
    nSelected,nSelected),eye(nVertical));
raw = struct("values",values,"valuesXi",valuesXi, ...
    "valuesX",valuesX,"valuesY",valuesY, ...
    "bottomValues",bottomValues,"bottomX",bottomX,"bottomY",bottomY, ...
    "tangency",tangency,"projectedTangency",projectedTangency, ...
    "gauge",gauge,"conjugacy",conjugacy, ...
    "Pxy",c.Pxy,"hX",terrainScale*c.hX, ...
    "hY",terrainScale*c.hY,"terrainScale",terrainScale, ...
    "selectedHorizontal",selectedHorizontal, ...
    "numberOfVerticalFunctions",nVertical, ...
    "topDerivative",topDerivative,"topBasis",topBasis, ...
    "scalarDegree",degree,"rankTop",rankTop);
end

function [stationary,nonTangent] = splitByBottomTangency(raw)
if norm(raw.projectedTangency,"fro") > 0
    constraint = [raw.projectedTangency ...
        /norm(raw.projectedTangency,2);raw.gauge'];
else
    constraint = raw.gauge';
end
[~,singularMatrix,V] = svd(constraint);
singularValues = diag(singularMatrix);
nominalTolerance = max(size(constraint))*eps(max(singularValues))*100;
rankValues = arrayfun(@(factor) ...
    nnz(singularValues > factor*nominalTolerance),[0.1 1 10]);
rankConstraint = rankValues(2);
stationaryBasis = canonicalColumnPhases(V(:,rankConstraint+1:end));

rowBasis = V(:,1:rankConstraint);
gaugeCoordinate = rowBasis'*raw.gauge;
if rankConstraint <= 1
    nonTangentBasis = zeros(size(rowBasis,1),0);
else
    [gaugeQ,~] = qr(gaugeCoordinate);
    nonTangentBasis = canonicalColumnPhases(rowBasis*gaugeQ(:,2:end));
end

stationary = applyScalarBasis(raw,stationaryBasis);
nonTangent = applyScalarBasis(raw,nonTangentBasis);
stationary.diagnostics = scalarDiagnostics(raw,stationaryBasis, ...
    rankConstraint,rankValues,singularValues,nominalTolerance);
nonTangent.diagnostics = struct( ...
    "numberOfNonTangentScalarCoordinates",size(nonTangentBasis,2), ...
    "projectedTangencyRank",max(rankConstraint-1,0));
end

function scalar = applyScalarBasis(raw,basis)
scalar = struct( ...
    "values",raw.values*basis, ...
    "valuesXi",raw.valuesXi*basis, ...
    "valuesX",raw.valuesX*basis, ...
    "valuesY",raw.valuesY*basis, ...
    "bottomValues",raw.bottomValues*basis, ...
    "bottomValuesByHorizontalMode",raw.Pxy*raw.bottomValues*basis, ...
    "coefficientBasis",basis,"terrainScale",raw.terrainScale);
end

function diagnostics = scalarDiagnostics(raw,basis,rankConstraint, ...
    rankValues,singularValues,nominalTolerance)
tangencyResidual = raw.tangency*basis;
projectedTangencyResidual = raw.projectedTangency*basis;
bottomScale = norm(raw.bottomX,"fro")*norm(raw.hY) ...
    +norm(raw.bottomY,"fro")*norm(raw.hX);
if bottomScale == 0
    bottomScale = norm(raw.bottomX,"fro")+norm(raw.bottomY,"fro");
end
P = basis*basis';
conjugacyDefect = norm(P*raw.conjugacy ...
    -raw.conjugacy*conj(P),"fro")/max(2*norm(P,"fro"),realmin);
diagnostics = struct( ...
    "scalarDegree",raw.scalarDegree, ...
    "numberOfHorizontalModes",numel(raw.selectedHorizontal), ...
    "numberOfRawScalarCoordinates",size(raw.values,2), ...
    "numberOfStationaryScalarCoordinates",size(basis,2), ...
    "topConstraintRank",raw.rankTop, ...
    "bottomTangencyRank",rankConstraint-1, ...
    "constraintSingularValues",singularValues, ...
    "constraintTolerance",nominalTolerance, ...
    "constraintRanks",rankValues, ...
    "rankStable",all(rankValues == rankConstraint), ...
    "surfaceDerivativeDefect", ...
    norm(raw.topDerivative.'*raw.topBasis,"fro") ...
    /max(norm(raw.topDerivative)*norm(raw.topBasis,"fro"),realmin), ...
    "bottomTangencyDefect",norm(tangencyResidual,"fro") ...
    /max(bottomScale,realmin), ...
    "projectedBottomTangencyDefect",norm(projectedTangencyResidual,"fro") ...
    /max(norm(raw.projectedTangency,"fro")*norm(basis,"fro"),realmin), ...
    "scalarConjugacyDefect",conjugacyDefect);
end

function [G,diagnostics] = projectGeostrophicState(c,scalar)
H3 = scalar.terrainScale*repmat(c.H,c.nZ,1);
HX3 = scalar.terrainScale*repmat(c.HX,c.nZ,1);
HY3 = scalar.terrainScale*repmat(c.HY,c.nZ,1);
gamma = 1-H3;
gammaX = -HX3;
gammaY = -HY3;
N2 = c.wvt.N2Function(gamma.*c.xiGrid);
if isscalar(N2)
    N2 = repmat(N2,size(c.xiGrid));
end
N2 = N2(:);
u = -gamma.*scalar.valuesY+c.xiGrid.*gammaY.*scalar.valuesXi;
v = gamma.*scalar.valuesX-c.xiGrid.*gammaX.*scalar.valuesXi;
w = c.xiGrid.*(gammaX.*scalar.valuesY-gammaY.*scalar.valuesX);
eta = -c.wvt.f*scalar.valuesXi./(gamma.*N2);

[uColumns,vColumns,wColumns,etaColumns] = componentColumns(c);
raw = zeros(c.nX,size(u,2));
raw(uColumns,:) = weightedProject(c.Ru(:,uColumns),c.volumeWeight,u);
raw(vColumns,:) = weightedProject(c.Rv(:,vColumns),c.volumeWeight,v);
raw(wColumns,:) = weightedProject(c.Rwh(:,wColumns),c.volumeWeight,w);
raw(etaColumns,:) = weightedProject(c.Reta(:,etaColumns),c.volumeWeight,eta);
G = c.N'*raw;
representedRaw = c.N*G;
represented = [c.Ru*representedRaw,c.Rv*representedRaw, ...
    c.Rwh*representedRaw,c.Reta*representedRaw];
targets = [u,v,w,eta];
diagnostics = struct( ...
    "rawCoordinateDefect",blockDefect(representedRaw,raw), ...
    "maximumDefect",blockDefect(represented,targets), ...
    "continuityDefect",norm(c.continuity*representedRaw,"fro") ...
    /max(norm(c.continuity,"fro")*norm(representedRaw,"fro"),realmin));
end

function coefficients = weightedProject(reconstruction,weight,values)
mass = reconstruction'*(weight.*reconstruction);
coefficients = mass\(reconstruction'*(weight.*values));
end

function [G,scalar,diagnostics] = energyOrthonormalize(G,scalar,E)
gram = G'*E*G;
gram = (gram+gram')/2;
[R,flag] = chol(gram);
if flag ~= 0
    error("WVTerrainEnergyGalerkin:StationaryEnergyRankFailure", ...
        "The represented stationary geostrophic states are not independent in physical energy.")
end
G = G/R;
scalar.values = scalar.values/R;
scalar.valuesXi = scalar.valuesXi/R;
scalar.valuesX = scalar.valuesX/R;
scalar.valuesY = scalar.valuesY/R;
scalar.bottomValues = scalar.bottomValues/R;
diagnostics = struct("energyOrthogonalityDefect", ...
    norm(G'*E*G-eye(size(G,2)),"fro")/sqrt(size(G,2)), ...
    "gramReciprocalConditionNumber",rcond(gram));
end

function [D,diagnostics] = orthonormalizeComplement(D,G,E)
gram = D'*E*D;
gram = (gram+gram')/2;
[R,flag] = chol(gram);
if flag ~= 0
    error("WVTerrainEnergyGalerkin:BottomSeedEnergyRankFailure", ...
        "The non-tangent bottom seed states are not independent in physical energy.")
end
D = D/R;
diagnostics = struct( ...
    "energyOrthogonalityDefect", ...
    norm(D'*E*D-eye(size(D,2)),"fro")/sqrt(size(D,2)), ...
    "stationaryOrthogonalityDefect", ...
    norm(G'*E*D,"fro")/max(norm(G,"fro")*norm(E,"fro") ...
    *norm(D,"fro"),realmin), ...
    "gramReciprocalConditionNumber",rcond(gram));
end

function [uColumns,vColumns,wColumns,etaColumns] = componentColumns(c)
uColumns = zeros(c.nK*c.nF,1);
vColumns = zeros(c.nK*c.nF,1);
wColumns = zeros(c.nK*c.nG,1);
etaColumns = zeros(c.nK*c.nH,1);
for iK = 1:c.nK
    block = (iK-1)*c.nXBlock;
    uColumns((iK-1)*c.nF+(1:c.nF)) = block+(1:c.nF);
    vColumns((iK-1)*c.nF+(1:c.nF)) = block+c.nF+(1:c.nF);
    wColumns((iK-1)*c.nG+(1:c.nG)) = block+2*c.nF+(1:c.nG);
    etaColumns((iK-1)*c.nH+(1:c.nH)) = ...
        block+2*c.nF+c.nG+(1:c.nH);
end
end

function conjugateIndex = conjugateHorizontalIndex(horizontalLayout)
nK = height(horizontalLayout);
conjugateIndex = zeros(nK,1);
for iK = 1:nK
    partner = find(horizontalLayout.kMode ...
        == -horizontalLayout.kMode(iK) ...
        & horizontalLayout.lMode == -horizontalLayout.lMode(iK),1);
    if isempty(partner)
        error("WVTerrainEnergyGalerkin:IncompletePrimitiveHorizontalConjugacy", ...
            "The retained horizontal layout is missing a Fourier conjugate.")
    end
    conjugateIndex(iK) = partner;
end
end

function values = canonicalColumnPhases(values)
for iColumn = 1:size(values,2)
    [~,pivot] = max(abs(values(:,iColumn)));
    if values(pivot,iColumn) ~= 0
        values(:,iColumn) = values(:,iColumn) ...
            *conj(values(pivot,iColumn))/abs(values(pivot,iColumn));
    end
end
end

function value = normalizedProduct(residual,left,right)
value = norm(residual,"fro") ...
    /max(norm(left,"fro")*norm(right,"fro"),realmin);
end

function value = blockDefect(first,second)
value = norm(first-second,"fro") ...
    /max(norm(first,"fro")+norm(second,"fro"),realmin);
end

function diagnostics = emptyRepresentation
diagnostics = struct("rawCoordinateDefect",0, ...
    "maximumDefect",0,"continuityDefect",0);
end
