function audit = buildCompleteStationaryBalancedSpaceAudit(problem,trustedBounds,supportBounds,degrees,paddingFactors,quadratureOrder)
% Build the complete finite-terrain stationary balanced-space oracle.

nRefinement = numel(degrees);
nPadding = numel(paddingFactors);
details = cell(nRefinement,nPadding);
refinement = repmat(emptySummary,0,1);
for iRefinement = 1:nRefinement
    layout = retainedLayout(problem,supportBounds(iRefinement,:));
    terrain = terrainSupport(problem);
    validateGuard(layout,trustedBounds,terrain);
    degree = degrees(iRefinement);
    if isempty(quadratureOrder)
        order = max(2*degree+7,18);
    else
        order = quadratureOrder;
    end
    for iPadding = 1:nPadding
        [primitive,context] = buildGlobalSmallTerrainPrimitiveAudit(problem, ...
            degree,order,1e-3,horizontalLayout=layout, ...
            paddingFactor=paddingFactors(iPadding),trustedModeBounds=trustedBounds, ...
            rejectTerrainNyquist=true,evaluationScales=1);
        result = stationarySpaceAudit(context,primitive.evaluatedDirections(1), ...
            trustedBounds,degrees(1));
        result.projection = primitive.projection;
        details{iRefinement,iPadding} = result;
        refinement(end+1,1) = summarize(result,degree,order, ...
            supportBounds(iRefinement,:),paddingFactors(iPadding)); %#ok<AGROW>
    end
end

finest = details{end,1};
padding = paddingDiagnostics(details(end,:));
convergence = convergenceDiagnostics(details(:,1));
tolerance = struct("projection",1e-12,"structure",1e-12, ...
    "representation",1e-10,"green",1e-10,"stationary",1e-10, ...
    "bottomTangency",1e-10,"projector",1e-9,"conjugacy",1e-11, ...
    "padding",1e-10,"rankGap",100);

allResults = [details{:}];
structurePasses = max([allResults.projectionAdjointDefect]) <= tolerance.projection ...
    && max([allResults.exactConvolutionDefect]) <= tolerance.projection ...
    && max([allResults.energyHermitianDefect]) <= tolerance.structure ...
    && max([allResults.exchangeSkewHermitianDefect]) <= tolerance.structure ...
    && min([allResults.energyRcond]) > 1e-12;
stationaryPasses = finest.stateRepresentationDefect <= tolerance.representation ...
    && finest.greenIdentityDefect <= tolerance.green ...
    && finest.stationaryRowDefect <= tolerance.stationary ...
    && finest.trustedStationaryRowDefect <= tolerance.stationary ...
    && finest.projectedBottomTangencyDefect <= tolerance.bottomTangency ...
    && finest.strongBottomTangencyDefect <= tolerance.bottomTangency ...
    && finest.stationaryNullspaceGap >= tolerance.rankGap ...
    && finest.scalarRankStable ...
    && finest.stationaryNullityStable ...
    && finest.conjugacyDefect <= tolerance.conjugacy ...
    && padding.maximumDefect <= tolerance.padding;
slopeScale = max(hypot(problem.originatingTransform.diffX(problem.topographicHeight), ...
    problem.originatingTransform.diffY(problem.topographicHeight)),[],"all");
bottomComplementPasses = slopeScale == 0 || finest.nonstationaryBottomDimension > 0;

if ~structurePasses
    status = "implementation-unresolved";
    diagnosis = "The common projection or primitive energy/exchange structure failed before stationary-space classification.";
elseif ~bottomComplementPasses
    status = "stationary-state-representation-blocker";
    diagnosis = "The stationary construction absorbed every bottom direction over sloping terrain instead of preserving a dynamical bottom complement.";
elseif stationaryPasses
    status = "complete-stationary-space-oracle";
    diagnosis = "The tangent geostrophic construction, complete APV Green identity, and zero-frequency primitive weak subspace agree at finite terrain amplitude.";
elseif finest.stateRepresentationDefect > tolerance.representation ...
        || finest.greenIdentityDefect > tolerance.green ...
        || finest.strongBottomTangencyDefect > tolerance.bottomTangency
    status = "stationary-state-representation-blocker";
    diagnosis = "The primitive trial space does not represent the complete tangent geostrophic family or its boundary Green identity to the required tolerance.";
else
    status = "stationary-nullspace-blocker";
    diagnosis = "The represented tangent geostrophic family does not equal the zero-frequency nullspace of the unmodified primitive exchange form.";
end

audit = struct;
audit.scope = "milestone-8-complete-finite-terrain-stationary-balanced-space";
audit.status = status;
audit.isCompatible = status == "complete-stationary-space-oracle";
audit.diagnosis = diagnosis;
audit.trustedModeBounds = trustedBounds;
audit.supportModeBounds = supportBounds;
audit.primitivePolynomialDegrees = degrees;
audit.paddingFactors = paddingFactors;
audit.quadratureOrder = quadratureOrder;
audit.refinement = refinement;
audit.convergence = convergence;
audit.padding = padding;
audit.finest = finest;
audit.details = details;
audit.requiredTolerance = tolerance;
audit.nextScope = "milestone-9-only-if-compatible-and-separately-authorized";
end

function result = stationarySpaceAudit(c,direction,trustedBounds,scalarDegree)
fullScalar = stationaryScalarSpace(c,scalarDegree,true(c.nK,1));
[GfullRaw,fullRepresentation] = projectGeostrophicState(c,fullScalar);
[~,fullScalar,fullOrthogonalization] = energyOrthonormalize( ...
    GfullRaw,fullScalar,direction.E);
trustedHorizontal = abs(c.horizontalLayout.kMode) <= trustedBounds(1) ...
    & abs(c.horizontalLayout.lMode) <= trustedBounds(2);
scalar = stationaryScalarSpace(c,scalarDegree,trustedHorizontal);
[GtestRaw,representation] = projectGeostrophicState(c,scalar);
[Gtest,scalar,orthogonalization] = energyOrthonormalize( ...
    GtestRaw,scalar,direction.E);

H3 = repmat(c.H,c.nZ,1);
gamma = 1-H3;
weakRow = Gtest'*direction.E;
volumeRow = -c.wvt.rho0*(scalar.values'*(c.volumeWeight.*gamma.*direction.Q));
gammaBottom = 1-c.H;
uBottom = (c.uBottom*c.N)./gammaBottom;
vBottom = (c.vBottom*c.N)./gammaBottom;
etaBottom = c.phase*direction.B;
cBottom = c.wvt.f*etaBottom+vBottom.*c.hX-uBottom.*c.hY;
boundaryRow = c.wvt.rho0*(scalar.bottomValues'*cBottom)/c.nXY;
strongRow = volumeRow+boundaryRow;

energyFactor = chol((direction.E+direction.E')/2);
energyExchange = energyFactor'\(direction.J/energyFactor);
[~,singularMatrix,V] = svd(energyExchange);
singularValues = diag(singularMatrix);
[gapRatio,gapIndex] = max(singularValues(1:end-1) ...
    ./max(singularValues(2:end),realmin));
nWeakNull = numel(singularValues)-gapIndex;
Ynull = V(:,gapIndex+1:end);
Gnull = energyFactor\Ynull;
complementSingularValue = singularValues(gapIndex);
stationarySingularValues = singularValues(end-nWeakNull+1:end);
stationarySingularValue = max(stationarySingularValues,[],"all");
gapScale = max(stationarySingularValue,eps(max(singularValues)));
nullspaceGap = complementSingularValue/gapScale;
projectorNull = Ynull*Ynull';
Ytest = energyFactor*Gtest;
trustedProjectorDefect = norm((eye(size(projectorNull))-projectorNull)*Ytest,"fro") ...
    /max(norm(Ytest,"fro"),realmin);
principalCosines = svd(Ytest'*Ynull);
rankCutoff = sqrt(max(stationarySingularValue,realmin) ...
    *max(complementSingularValue,realmin));
rankCutoffs = rankCutoff*[0.1 1 10];
nullities = arrayfun(@(cutoff)nnz(singularValues <= cutoff),rankCutoffs);

n = size(direction.J,1);
C = sparse((1:n)',c.coordinateConjugateIndex,1,n,n);
geostrophicTestProjector = Gtest*(Gtest'*direction.E);
conjugacyResidual = geostrophicTestProjector*C ...
    -C*conj(geostrophicTestProjector);
conjugacyDefect = norm(conjugacyResidual,"fro") ...
    /max(2*norm(geostrophicTestProjector,"fro"),realmin);

stationaryResidual = direction.J*Gtest;
projectedBottomResidual = direction.R*Gtest;
fullSupportBottomResidual = direction.R*Gnull;
trustedStationaryResidual = direction.J*Gtest;
QG = direction.QProjected*Gtest;
bottomValues = direction.B*Gtest;
fullBottomRank = numericalRank(direction.B);
stationaryBottomRank = numericalRank(bottomValues);
bottomStreamfunction = scalar.bottomValuesByHorizontalMode;
stationaryBottomStreamfunctionRank = numericalRank(bottomStreamfunction);
expectedZeroAPVDimension = max(stationaryBottomStreamfunctionRank-1,0);
[~,apvSingularMatrix,apvVectors] = svd(QG,"econ");
apvSingularValues = diag(apvSingularMatrix);
if expectedZeroAPVDimension == 0
    zeroAPV = zeros(size(Gtest,2),0);
else
    zeroAPV = apvVectors(:,end-expectedZeroAPVDimension+1:end);
end
qRank = size(Gtest,2)-expectedZeroAPVDimension;
zeroAPVBottomRank = numericalRank(bottomValues*zeroAPV);
zeroK = find(c.horizontalLayout.kMode == 0 & c.horizontalLayout.lMode == 0,1);
nonmeanRows = horzcat(c.layout.admissibleRanges{setdiff(1:c.nK,zeroK)});
expectedMeanDimension = max(scalar.diagnostics.scalarDegree-1,0);
if isempty(nonmeanRows) || expectedMeanDimension == 0
    meanSubspaceDefect = 0;
else
    [~,~,meanVectors] = svd(Gtest(nonmeanRows,:),"econ");
    meanBasis = meanVectors(:,end-expectedMeanDimension+1:end);
    meanSubspaceDefect = norm(Gtest(nonmeanRows,:)*meanBasis,"fro") ...
        /max(norm(Gtest(nonmeanRows,:),"fro"),realmin);
end

result = struct;
result.stationaryBasis = Gtest;
result.geostrophicTestBasis = Gtest;
result.stationaryEnergyProjector = geostrophicTestProjector;
result.geostrophicTestEnergyProjector = geostrophicTestProjector;
result.weakNullspaceBasisEnergyCoordinates = Ynull;
result.scalar = scalar.diagnostics;
result.fullScalar = fullScalar.diagnostics;
result.stateRepresentation = representation;
result.stateRepresentationDefect = representation.rawCoordinateDefect;
result.unprojectedStateRepresentationDefect = representation.maximumDefect;
result.fullSupportStateRepresentation = fullRepresentation;
result.orthogonalization = fullOrthogonalization;
result.trustedStateRepresentation = representation;
result.trustedOrthogonalization = orthogonalization;
result.numberOfStationaryStates = size(Gtest,2);
result.numberOfTrustedGeostrophicTests = size(Gtest,2);
result.numberOfProjectedGeostrophicStates = ...
    fullScalar.diagnostics.numberOfStationaryScalarCoordinates;
result.numberOfUnclassifiedWeakNullStates = nWeakNull ...
    -fullScalar.diagnostics.numberOfStationaryScalarCoordinates;
result.numberOfWeakNullStates = nWeakNull;
result.stationaryVolumeAPVRank = qRank;
result.stationaryZeroVolumeAPVDimension = expectedZeroAPVDimension;
result.stationaryZeroVolumeAPVCandidateDimension = expectedZeroAPVDimension;
result.stationaryAPVSingularValues = apvSingularValues;
result.zeroVolumeAPVBottomRank = zeroAPVBottomRank;
result.zeroVolumeAPVBottomCandidateRank = zeroAPVBottomRank;
result.meanDensityAnomalyDimension = expectedMeanDimension;
result.meanDensityAnomalySubspaceDefect = meanSubspaceDefect;
result.baselineBottomInversionMaximumAPVResidual = ...
    c.problem.constructionDiagnostics.bottomInversion.maximumAPVResidual;
result.baselineBottomInversionMaximumResidual = ...
    c.problem.constructionDiagnostics.bottomInversion.maximumInversionResidual;
result.fullBottomDimension = fullBottomRank;
result.stationaryBottomDimension = stationaryBottomRank;
result.stationaryBottomStreamfunctionDimension = stationaryBottomStreamfunctionRank;
result.nonstationaryBottomDimension = scalar.diagnostics.numberOfHorizontalModes ...
    -stationaryBottomStreamfunctionRank;
result.greenIdentityDefect = pairDefect(weakRow,strongRow);
result.greenVolumeRowNorm = norm(volumeRow,"fro");
result.greenBoundaryRowNorm = norm(boundaryRow,"fro");
result.stationaryRowDefect = norm(stationaryResidual,"fro") ...
    /max(norm(direction.J,"fro")*norm(Gtest,"fro"),realmin);
result.trustedStationaryRowDefect = norm(trustedStationaryResidual,"fro") ...
    /max(norm(direction.J,"fro")*norm(Gtest,"fro"),realmin);
result.projectedBottomTangencyDefect = norm(projectedBottomResidual,"fro") ...
    /max(norm(direction.R,"fro")*norm(Gtest,"fro"),realmin);
result.fullSupportBottomTangencyDefect = norm(fullSupportBottomResidual,"fro") ...
    /max(norm(direction.R,"fro")*norm(Gnull,"fro"),realmin);
result.strongBottomTangencyDefect = scalar.diagnostics.bottomTangencyDefect;
result.stationaryProjectorDefect = trustedProjectorDefect;
result.trustedStationaryContainmentDefect = trustedProjectorDefect;
result.minimumPrincipalCosine = min(principalCosines);
result.stationaryNullspaceGap = nullspaceGap;
result.stationarySpectralGapRatio = gapRatio;
result.stationaryProjectorRoundoffScale = ...
    sqrt(numel(singularValues))*eps(max(singularValues)) ...
    /max(complementSingularValue,realmin);
result.stationarySingularValue = stationarySingularValue;
result.complementSingularValue = complementSingularValue;
result.weakSingularValues = singularValues;
result.stationaryNullityCutoffs = rankCutoffs;
result.stationaryNullities = nullities;
result.stationaryNullityStable = all(nullities == nWeakNull);
result.scalarRankStable = scalar.diagnostics.rankStable ...
    && fullScalar.diagnostics.rankStable;
result.conjugacyDefect = conjugacyDefect;
result.energyHermitianDefect = norm(direction.E-direction.E',"fro") ...
    /max(norm(direction.E,"fro"),realmin);
result.exchangeSkewHermitianDefect = norm(direction.J+direction.J',"fro") ...
    /max(norm(direction.J,"fro"),realmin);
result.energyRcond = rcond((direction.E+direction.E')/2);
result.projectionAdjointDefect = c.projectionDiagnostics.adjointDefect;
result.exactConvolutionDefect = c.projectionDiagnostics.maximumExactConvolutionDefect;
result.bottomStreamfunctionByHorizontalMode = bottomStreamfunction;
result.bottomStreamfunctionProjector = subspaceProjector(bottomStreamfunction);
result.horizontalLayout = c.horizontalLayout;
end

function scalar = stationaryScalarSpace(c,degree,horizontalMask)
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

tangency = bottomX.*c.hY-bottomY.*c.hX;
projectedTangency = c.Pxy*tangency;
zeroK = find(c.horizontalLayout.kMode == 0 & c.horizontalLayout.lMode == 0,1);
zeroSelected = find(selectedHorizontal == zeroK,1);
if isempty(zeroSelected)
    error("WVTerrainEnergyGalerkin:StationaryScalarMissingGauge", ...
        "The stationary scalar support must include the zero horizontal mode.")
end
constantVertical = topBasis\eye(degree+1,1);
gauge = zeros(nRaw,1);
gauge((zeroSelected-1)*nVertical+(1:nVertical)) = constantVertical;
gauge = gauge/norm(gauge);
if norm(projectedTangency,"fro") > 0
    constraint = [projectedTangency/norm(projectedTangency,2);gauge'];
else
    constraint = gauge';
end
[~,singularMatrix,V] = svd(constraint);
singularValues = diag(singularMatrix);
nominalTolerance = max(size(constraint))*eps(max(singularValues))*100;
rankValues = arrayfun(@(factor)nnz(singularValues > factor*nominalTolerance), ...
    [0.1 1 10]);
rankConstraint = rankValues(2);
stationaryBasis = canonicalColumnPhases(V(:,rankConstraint+1:end));

values = values*stationaryBasis;
valuesXi = valuesXi*stationaryBasis;
valuesX = valuesX*stationaryBasis;
valuesY = valuesY*stationaryBasis;
bottomValues = bottomValues*stationaryBasis;
tangencyResidual = tangency*stationaryBasis;
projectedTangencyResidual = projectedTangency*stationaryBasis;
bottomScale = norm(bottomX,"fro")*norm(c.hY) ...
    +norm(bottomY,"fro")*norm(c.hX);

selectedLayout = c.horizontalLayout(selectedHorizontal,:);
horizontalConjugate = conjugateHorizontalIndex(selectedLayout);
C = kron(sparse((1:nSelected)',horizontalConjugate,1,nSelected,nSelected),eye(nVertical));
P = stationaryBasis*stationaryBasis';
conjugacyDefect = norm(P*C-C*conj(P),"fro")/max(2*norm(P,"fro"),realmin);
diagnostics = struct("scalarDegree",degree, ...
    "numberOfHorizontalModes",nSelected, ...
    "numberOfRawScalarCoordinates",nRaw, ...
    "numberOfStationaryScalarCoordinates",size(stationaryBasis,2), ...
    "topConstraintRank",rankTop,"bottomTangencyRank",rankConstraint-1, ...
    "constraintSingularValues",singularValues, ...
    "constraintTolerance",nominalTolerance,"constraintRanks",rankValues, ...
    "rankStable",all(rankValues == rankConstraint), ...
    "surfaceDerivativeDefect",norm(topDerivative.'*topBasis,"fro") ...
        /max(norm(topDerivative)*norm(topBasis,"fro"),realmin), ...
    "bottomTangencyDefect",norm(tangencyResidual,"fro") ...
        /max(bottomScale,realmin), ...
    "projectedBottomTangencyDefect",norm(projectedTangencyResidual,"fro") ...
        /max(norm(projectedTangency,"fro")*norm(stationaryBasis,"fro"),realmin), ...
    "scalarConjugacyDefect",conjugacyDefect);
scalar = struct("values",values,"valuesXi",valuesXi,"valuesX",valuesX, ...
    "valuesY",valuesY,"bottomValues",bottomValues, ...
    "bottomValuesByHorizontalMode",c.Pxy*bottomValues, ...
    "coefficientBasis",stationaryBasis,"diagnostics",diagnostics);
end

function [G,diagnostics] = projectGeostrophicState(c,scalar)
H3 = repmat(c.H,c.nZ,1);
HX3 = repmat(c.HX,c.nZ,1);
HY3 = repmat(c.HY,c.nZ,1);
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
diagnostics = struct("rawCoordinateDefect",blockDefect(representedRaw,raw), ...
    "maximumDefect",blockDefect(represented,targets), ...
    "continuityDefect",norm(c.continuity*representedRaw,"fro") ...
        /max(norm(c.continuity,"fro")*norm(representedRaw,"fro"),realmin));
end

function coefficients = weightedProject(reconstruction,weight,values)
mass = reconstruction'*(weight.*reconstruction);
coefficients = mass\(reconstruction'*(weight.*values));
end

function [G,scalar,diagnostics] = energyOrthonormalize(G,scalar,E)
gram = (G'*E*G);
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
    etaColumns((iK-1)*c.nH+(1:c.nH)) = block+2*c.nF+c.nG+(1:c.nH);
end
end

function diagnostics = convergenceDiagnostics(results)
fields = ["stateRepresentationDefect","greenIdentityDefect", ...
    "stationaryRowDefect","projectedBottomTangencyDefect", ...
    "strongBottomTangencyDefect","stationaryProjectorDefect"];
diagnostics = struct;
for field = fields
    values = reshape(cellfun(@(result)result.(field),results),1,[]);
    diagnostics.(field) = values;
    diagnostics.(field+"Reduction") = values(1:end-1)./max(values(2:end),realmin);
end
diagnostics.numberOfStationaryStates = cellfun( ...
    @(result)result.numberOfStationaryStates,results);
diagnostics.numberOfWeakNullStates = cellfun( ...
    @(result)result.numberOfWeakNullStates,results);
diagnostics.nonstationaryBottomDimension = cellfun( ...
    @(result)result.nonstationaryBottomDimension,results);
end

function diagnostics = paddingDiagnostics(results)
reference = results{1};
fields = ["geostrophicTestEnergyProjector","bottomStreamfunctionProjector"];
diagnostics = struct;
maximum = 0;
for field = fields
    values = zeros(numel(results)-1,1);
    for iResult = 2:numel(results)
        values(iResult-1) = blockDefect(results{iResult}.(field),reference.(field));
    end
    diagnostics.(field+"Defect") = values;
    maximum = max(maximum,max(values));
end
diagnostics.maximumDefect = maximum;
end

function summary = summarize(result,degree,order,support,padding)
summary = emptySummary;
summary.primitivePolynomialDegree = degree;
summary.quadratureOrder = order;
summary.supportModeBounds = support;
summary.paddingFactor = padding;
summary.numberOfStationaryStates = result.numberOfStationaryStates;
summary.numberOfWeakNullStates = result.numberOfWeakNullStates;
summary.nonstationaryBottomDimension = result.nonstationaryBottomDimension;
summary.stateRepresentationDefect = result.stateRepresentationDefect;
summary.greenIdentityDefect = result.greenIdentityDefect;
summary.stationaryRowDefect = result.stationaryRowDefect;
summary.projectedBottomTangencyDefect = result.projectedBottomTangencyDefect;
summary.strongBottomTangencyDefect = result.strongBottomTangencyDefect;
summary.stationaryProjectorDefect = result.stationaryProjectorDefect;
summary.conjugacyDefect = result.conjugacyDefect;
end

function summary = emptySummary
summary = struct("primitivePolynomialDegree",0,"quadratureOrder",0, ...
    "supportModeBounds",[0 0],"paddingFactor",0, ...
    "numberOfStationaryStates",0,"numberOfWeakNullStates",0, ...
    "nonstationaryBottomDimension",0,"stateRepresentationDefect",NaN, ...
    "greenIdentityDefect",NaN,"stationaryRowDefect",NaN, ...
    "projectedBottomTangencyDefect",NaN,"strongBottomTangencyDefect",NaN, ...
    "stationaryProjectorDefect",NaN,"conjugacyDefect",NaN);
end

function value = pairDefect(first,second)
value = norm(first-second,"fro") ...
    /max(norm(first,"fro")+norm(second,"fro"),realmin);
end

function value = blockDefect(first,second)
value = norm(first-second,"fro") ...
    /max(norm(first,"fro")+norm(second,"fro"),realmin);
end

function rankValue = numericalRank(values)
singularValues = svd(values);
if isempty(singularValues)
    rankValue = 0;
    return
end
tolerance = max(size(values))*eps(max(singularValues))*100;
rankValue = nnz(singularValues > tolerance);
end

function projector = subspaceProjector(values)
[U,S,~] = svd(values,"econ");
singularValues = diag(S);
if isempty(singularValues)
    projector = zeros(size(values,1));
    return
end
tolerance = max(size(values))*eps(max(singularValues))*100;
rankValue = nnz(singularValues > tolerance);
projector = U(:,1:rankValue)*U(:,1:rankValue)';
end

function values = canonicalColumnPhases(values)
for iColumn = 1:size(values,2)
    [~,pivot] = max(abs(values(:,iColumn)));
    if values(pivot,iColumn) ~= 0
        values(:,iColumn) = values(:,iColumn)*conj(values(pivot,iColumn)) ...
            /abs(values(pivot,iColumn));
    end
end
end

function conjugateIndex = conjugateHorizontalIndex(horizontalLayout)
nK = height(horizontalLayout);
conjugateIndex = zeros(nK,1);
for iK = 1:nK
    partner = find(horizontalLayout.kMode == -horizontalLayout.kMode(iK) ...
        & horizontalLayout.lMode == -horizontalLayout.lMode(iK),1);
    if isempty(partner)
        error("WVTerrainEnergyGalerkin:IncompleteStationaryHorizontalConjugacy", ...
            "The retained stationary scalar layout is missing a Fourier conjugate.")
    end
    conjugateIndex(iK) = partner;
end
end

function layout = retainedLayout(problem,bounds)
source = problem.horizontalLayout;
mask = abs(source.kMode) <= bounds(1) & abs(source.lMode) <= bounds(2);
layout = source(mask,:);
expected = (2*bounds(1)+1)*(2*bounds(2)+1);
if height(layout) ~= expected
    error("WVTerrainEnergyGalerkin:UnavailableStationarySupport", ...
        "The originating transform does not retain the complete signed rectangle [%d %d].", ...
        bounds(1),bounds(2))
end
end

function terrain = terrainSupport(problem)
wvt = problem.originatingTransform;
spectrum = fft2(problem.topographicHeight)/(wvt.Nx*wvt.Ny);
nyquist = logical(WVGeometryDoublyPeriodic.maskForNyquistModes(wvt.Nx,wvt.Ny));
scale = max(abs(spectrum),[],"all");
tolerance = 100*eps*max(scale,1);
if max(abs(spectrum(nyquist)),[],"all") > tolerance
    error("WVTerrainEnergyGalerkin:StationaryTerrainNyquist", ...
        "The stationary-space oracle requires terrain with no material Nyquist coefficient.")
end
[kMode,lMode] = ndgrid(wvt.kMode_dft,wvt.lMode_dft);
active = abs(spectrum) > tolerance & ~nyquist;
terrain = struct("kMode",kMode(active),"lMode",lMode(active), ...
    "coefficient",spectrum(active),"spectralTolerance",tolerance);
end

function validateGuard(layout,trustedBounds,terrain)
for iK = find(abs(layout.kMode) <= trustedBounds(1) ...
        & abs(layout.lMode) <= trustedBounds(2)).'
    destinations = [layout.kMode(iK)+terrain.kMode layout.lMode(iK)+terrain.lMode];
    for iDestination = 1:size(destinations,1)
        if ~any(layout.kMode == destinations(iDestination,1) ...
                & layout.lMode == destinations(iDestination,2))
            error("WVTerrainEnergyGalerkin:InsufficientStationaryGuard", ...
                "supportModeBounds must retain every first terrain sideband of the trusted band.")
        end
    end
end
end
