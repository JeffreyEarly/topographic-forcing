function factor = factorGlobalTangentStationarySpace(c,direction,exact,scope)
% Factor a polynomial stationary space into zero-trace and trace lifts.
arguments
    c (1,1) struct
    direction (1,1) struct
    exact (1,1) struct
    scope (1,1) string {mustBeMember(scope,["full","trusted"])} = "full"
end

if scope == "full"
    raw = exact.fullRawScalar;
    exactBasis = exact.fullBasis;
else
    raw = exact.trustedRawScalar;
    exactBasis = exact.trustedBasis;
end

trace = traceSpaces(c,raw);
% A zero bottom trace already excludes the constant streamfunction gauge.
% Adding the separately assembled gauge row here would duplicate that
% constraint at finite precision and can spuriously remove one direction.
zeroBasis = rightNullspace(trace.rawTrace);
zeroScalar = applyScalarBasis(raw,zeroBasis);
[zeroStates,zeroRepresentation] = projectGeostrophicState(c,zeroScalar);

[tangentSeedBasis,tangentRightInverse] = traceSeedBasis( ...
    trace.rawTrace,trace.tangentBasis);
tangentSeedScalar = applyScalarBasis(raw,tangentSeedBasis);
[tangentSeedStates,tangentSeedRepresentation] = ...
    projectGeostrophicState(c,tangentSeedScalar);
[tangentStates,tangentScalar,tangentSolve] = energyMinimalLifts( ...
    zeroStates,zeroScalar,tangentSeedStates,tangentSeedScalar,direction.E);

[complementSeedBasis,complementRightInverse] = traceSeedBasis( ...
    trace.rawTrace,trace.complementBasis);
complementSeedScalar = applyScalarBasis(raw,complementSeedBasis);
[complementSeedStates,complementSeedRepresentation] = ...
    projectGeostrophicState(c,complementSeedScalar);
[complementStates,complementScalar,complementSolve] = energyMinimalLifts( ...
    zeroStates,zeroScalar,complementSeedStates, ...
    complementSeedScalar,direction.E);

stationaryRaw = [zeroStates tangentStates];
stationaryScalar = concatenateScalar(zeroScalar,tangentScalar);
[stationaryBasis,stationaryScalar,stationaryOrthonormalization] = ...
    energyOrthonormalize(stationaryRaw,stationaryScalar,direction.E);

if isempty(complementStates)
    complementBasis = zeros(size(stationaryBasis,1),0);
    complementOrthogonality = 0;
else
    complementStates = complementStates-stationaryBasis ...
        *(stationaryBasis'*direction.E*complementStates);
    [complementBasis,~,~] = ...
        energyOrthonormalize(complementStates,complementScalar,direction.E);
    complementOrthogonality = norm( ...
        stationaryBasis'*direction.E*complementBasis,"fro") ...
        /max(norm(stationaryBasis,"fro")*norm(direction.E,"fro") ...
        *norm(complementBasis,"fro"),realmin);
end

projectorDefect = principalProjectorDefect( ...
    stationaryBasis,exactBasis,direction.E);
traceScale = max(norm(trace.rawTrace,"fro"),realmin);
factor = struct;
factor.basis = stationaryBasis;
factor.scalar = stationaryScalar;
factor.zeroTraceBasis = zeroStates;
factor.tangentLiftBasis = tangentStates;
factor.nonTangentLiftBasis = complementBasis;
factor.trace = trace;
factor.solvers = struct("tangent",tangentSolve, ...
    "nonTangent",complementSolve, ...
    "tangentRightInverse",tangentRightInverse, ...
    "nonTangentRightInverse",complementRightInverse, ...
    "orthonormalization",stationaryOrthonormalization);
factor.representation = struct("zeroTrace",zeroRepresentation, ...
    "tangentSeed",tangentSeedRepresentation, ...
    "nonTangentSeed",complementSeedRepresentation);
factor.diagnostics = struct( ...
    "numberOfRawCoordinates",size(raw.values,2), ...
    "numberOfZeroTraceCoordinates",size(zeroStates,2), ...
    "numberOfTangentTraceCoordinates",size(tangentStates,2), ...
    "numberOfNonTangentTraceCoordinates",size(complementBasis,2), ...
    "numberOfStationaryCoordinates",size(stationaryBasis,2), ...
    "numberOfExactStationaryCoordinates",size(exactBasis,2), ...
    "dimensionDefect",abs(size(stationaryBasis,2)-size(exactBasis,2)), ...
    "projectorDefect",projectorDefect, ...
    "zeroTraceDefect",norm(raw.bottomValues*zeroBasis,"fro")/traceScale, ...
    "tangentTraceDefect",trace.tangentDefect, ...
    "traceConjugacyDefect",trace.conjugacyDefect, ...
    "rankStable",trace.rankStable, ...
    "tangentLiftTraceDefect",tangentSolve.traceDefect, ...
    "tangentLiftEnergyDefect",tangentSolve.energyDefect, ...
    "tangentLiftAPVDefect",weakAPVMomentDefect( ...
        zeroStates,tangentStates,direction.E), ...
    "stationaryRowDefect",normalizedProduct( ...
        direction.J,stationaryBasis), ...
    "stationaryBottomTangencyDefect",normalizedProduct( ...
        direction.R,stationaryBasis), ...
    "stationaryConjugacyDefect",subspaceConjugacyDefect( ...
        stationaryBasis,direction.E,c.coordinateConjugateIndex), ...
    "complementStationaryOrthogonalityDefect", ...
        complementOrthogonality);
end

function value = weakAPVMomentDefect(zeroStates,liftStates,E)
value = norm(zeroStates'*E*liftStates,"fro") ...
    /max(norm(zeroStates,"fro")*norm(E,"fro") ...
    *norm(liftStates,"fro"),realmin);
end

function trace = traceSpaces(c,raw)
selected = raw.selectedHorizontal;
k = reshape(c.horizontalLayout.k(selected),1,[]);
l = reshape(c.horizontalLayout.l(selected),1,[]);
phase = c.phase(:,selected);
traceX = phase.*(1i*k);
traceY = phase.*(1i*l);
tangency = traceX.*raw.hY-traceY.*raw.hX;
operator = raw.Pxy*tangency;
scale = norm(operator,2);
if scale == 0
    scaledOperator = operator;
else
    scaledOperator = operator/scale;
end
zeroK = find(c.horizontalLayout.kMode(selected) == 0 ...
    & c.horizontalLayout.lMode(selected) == 0,1);
if isempty(zeroK)
    error("WVTerrainEnergyGalerkin:GlobalTangentTraceMissingGauge", ...
        "The bottom-trace support must contain the zero horizontal mode.")
end
gauge = zeros(numel(selected),1);
gauge(zeroK) = 1;
[~,singularMatrix,V] = svd([scaledOperator;gauge']);
singularValues = diag(singularMatrix);
nominalTolerance = max(size(scaledOperator,1)+1,numel(selected)) ...
    *eps(max([singularValues;1]))*100;
rankValues = arrayfun(@(value)nnz( ...
    singularValues > value*nominalTolerance),[0.1 1 10]);
rankConstraint = rankValues(2);
tangentBasis = canonicalColumnPhases(V(:,rankConstraint+1:end));

rowBasis = V(:,1:rankConstraint);
gaugeCoordinate = rowBasis'*gauge;
if rankConstraint <= 1
    complementBasis = zeros(numel(selected),0);
else
    [gaugeQ,~] = qr(gaugeCoordinate);
    complementBasis = canonicalColumnPhases(rowBasis*gaugeQ(:,2:end));
end

rawTrace = raw.Pxy*raw.bottomValues;
rawTrace = rawTrace(selected,:);
horizontalConjugate = conjugateHorizontalIndex( ...
    c.horizontalLayout(selected,:));
conjugacy = sparse((1:numel(selected))',horizontalConjugate,1, ...
    numel(selected),numel(selected));
projector = tangentBasis*tangentBasis';
conjugacyDefect = norm(projector*conjugacy ...
    -conjugacy*conj(projector),"fro") ...
    /max(2*norm(projector,"fro"),realmin);
trace = struct("selectedHorizontal",selected, ...
    "operator",operator,"rawTrace",rawTrace,"gauge",gauge, ...
    "tangentBasis",tangentBasis,"complementBasis",complementBasis, ...
    "singularValues",singularValues,"tolerance",nominalTolerance, ...
    "rank",rankConstraint-1,"rankValues",rankValues-1, ...
    "rankStable",all(rankValues == rankConstraint), ...
    "tangentDefect",norm(operator*tangentBasis,"fro") ...
        /max(norm(operator,"fro")*norm(tangentBasis,"fro"),realmin), ...
    "conjugacyDefect",conjugacyDefect);
end

function [basis,diagnostics] = traceSeedBasis(traceMap,target)
if isempty(target)
    basis = zeros(size(traceMap,2),0);
    diagnostics = struct("rank",size(traceMap,1), ...
        "minimumSingularValue",NaN,"conditionNumber",NaN, ...
        "residual",0);
    return
end
[U,S,V] = svd(traceMap,"econ");
singularValues = diag(S);
tolerance = max(size(traceMap))*eps(max(singularValues))*100;
rankValue = nnz(singularValues > tolerance);
if rankValue ~= size(traceMap,1)
    error("WVTerrainEnergyGalerkin:GlobalTangentTraceRankFailure", ...
        "The polynomial scalar space does not span every retained bottom trace.")
end
basis = V*((U'*target)./singularValues);
diagnostics = struct("rank",rankValue, ...
    "minimumSingularValue",singularValues(end), ...
    "conditionNumber",singularValues(1)/singularValues(end), ...
    "residual",norm(traceMap*basis-target,"fro") ...
        /max(norm(target,"fro"),realmin));
end

function basis = rightNullspace(constraint)
[~,S,V] = svd(constraint);
singularValues = diag(S);
tolerance = max(size(constraint))*eps(max([singularValues;1]))*100;
rankValue = nnz(singularValues > tolerance);
basis = canonicalColumnPhases(V(:,rankValue+1:end));
end

function scalar = applyScalarBasis(raw,basis)
scalar = struct("values",raw.values*basis, ...
    "valuesXi",raw.valuesXi*basis, ...
    "valuesX",raw.valuesX*basis, ...
    "valuesY",raw.valuesY*basis, ...
    "bottomValues",raw.bottomValues*basis, ...
    "coefficientBasis",basis,"terrainScale",raw.terrainScale);
end

function scalar = concatenateScalar(first,second)
fields = ["values","valuesXi","valuesX","valuesY", ...
    "bottomValues","coefficientBasis"];
scalar = first;
for field = fields
    scalar.(field) = [first.(field) second.(field)];
end
end

function [states,scalar,diagnostics] = energyMinimalLifts( ...
    zeroStates,zeroScalar,seedStates,seedScalar,E)
if isempty(seedStates)
    states = seedStates;
    scalar = seedScalar;
    diagnostics = struct("rank",size(zeroStates,2), ...
        "minimumSingularValue",NaN,"conditionNumber",NaN, ...
        "traceDefect",0,"energyDefect",0);
    return
end
R = chol((E+E')/2);
whitenedZero = R*zeroStates;
if isempty(zeroStates)
    correction = zeros(0,size(seedStates,2));
    singularValues = zeros(0,1);
    rankValue = 0;
else
    [U,S,V] = svd(whitenedZero,"econ");
    singularValues = diag(S);
    tolerance = max(size(whitenedZero)) ...
        *eps(max(singularValues))*100;
    rankValue = nnz(singularValues > tolerance);
    if rankValue ~= size(zeroStates,2)
        error("WVTerrainEnergyGalerkin:GlobalTangentLiftRankFailure", ...
            "The zero-bottom-trace states are rank deficient in physical energy.")
    end
    correction = -V*((U'*(R*seedStates))./singularValues);
end
states = seedStates+zeroStates*correction;
scalar = seedScalar;
fields = ["values","valuesXi","valuesX","valuesY", ...
    "bottomValues","coefficientBasis"];
for field = fields
    scalar.(field) = seedScalar.(field)+zeroScalar.(field)*correction;
end
traceScale = norm(seedScalar.bottomValues,"fro");
energyScale = norm(zeroStates,"fro")*norm(E,"fro") ...
    *norm(states,"fro");
diagnostics = struct("rank",rankValue, ...
    "minimumSingularValue",minimumValue(singularValues), ...
    "conditionNumber",conditionNumber(singularValues), ...
    "traceDefect",norm(scalar.bottomValues ...
        -seedScalar.bottomValues,"fro")/max(traceScale,realmin), ...
    "energyDefect",norm(zeroStates'*E*states,"fro") ...
        /max(energyScale,realmin));
end

function [states,scalar,diagnostics] = energyOrthonormalize(states,scalar,E)
if isempty(states)
    diagnostics = struct("rank",0,"minimumSingularValue",NaN, ...
        "conditionNumber",NaN,"defect",0);
    return
end
R = chol((E+E')/2);
[~,S,V] = svd(R*states,"econ");
singularValues = diag(S);
tolerance = max(size(states))*eps(max(singularValues))*100;
rankValue = nnz(singularValues > tolerance);
if rankValue ~= size(states,2)
    error("WVTerrainEnergyGalerkin:GlobalTangentStationaryRankFailure", ...
        "The factored stationary states are rank deficient in physical energy.")
end
transform = V./singularValues.';
states = states*transform;
fields = ["values","valuesXi","valuesX","valuesY", ...
    "bottomValues","coefficientBasis"];
for field = fields
    scalar.(field) = scalar.(field)*transform;
end
diagnostics = struct("rank",rankValue, ...
    "minimumSingularValue",singularValues(end), ...
    "conditionNumber",singularValues(1)/singularValues(end), ...
    "defect",norm(states'*E*states-eye(size(states,2)),"fro") ...
        /sqrt(size(states,2)));
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
raw(uColumns,:) = weightedProject( ...
    c.Ru(:,uColumns),c.volumeWeight,u);
raw(vColumns,:) = weightedProject( ...
    c.Rv(:,vColumns),c.volumeWeight,v);
raw(wColumns,:) = weightedProject( ...
    c.Rwh(:,wColumns),c.volumeWeight,w);
raw(etaColumns,:) = weightedProject( ...
    c.Reta(:,etaColumns),c.volumeWeight,eta);
G = c.N'*raw;
representedRaw = c.N*G;
represented = [c.Ru*representedRaw,c.Rv*representedRaw, ...
    c.Rwh*representedRaw,c.Reta*representedRaw];
targets = [u,v,w,eta];
diagnostics = struct("rawCoordinateDefect", ...
    blockDefect(representedRaw,raw), ...
    "maximumDefect",blockDefect(represented,targets), ...
    "continuityDefect",norm(c.continuity*representedRaw,"fro") ...
        /max(norm(c.continuity,"fro") ...
        *norm(representedRaw,"fro"),realmin));
end

function coefficients = weightedProject(reconstruction,weight,values)
mass = reconstruction'*(weight.*reconstruction);
coefficients = mass\(reconstruction'*(weight.*values));
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

function value = principalProjectorDefect(first,second,E)
if size(first,2) ~= size(second,2)
    value = Inf;
    return
end
R = chol((E+E')/2);
firstWhitened = R*first;
secondWhitened = R*second;
firstProjector = firstWhitened*firstWhitened';
secondProjector = secondWhitened*secondWhitened';
value = norm(firstProjector-secondProjector,"fro") ...
    /max(norm(firstProjector,"fro") ...
    +norm(secondProjector,"fro"),realmin);
end

function value = subspaceConjugacyDefect(X,H,conjugateIndex)
n = size(H,1);
C = sparse((1:n)',conjugateIndex,1,n,n);
P = X*(X'*H);
value = norm(P*C-C*conj(P),"fro") ...
    /max(2*norm(P,"fro"),realmin);
end

function conjugateIndex = conjugateHorizontalIndex(horizontalLayout)
nK = height(horizontalLayout);
conjugateIndex = zeros(nK,1);
for iK = 1:nK
    partner = find(horizontalLayout.kMode ...
        == -horizontalLayout.kMode(iK) ...
        & horizontalLayout.lMode == -horizontalLayout.lMode(iK),1);
    if isempty(partner)
        error("WVTerrainEnergyGalerkin:IncompleteTraceConjugacy", ...
            "The retained bottom-trace layout is missing a conjugate.")
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

function value = normalizedProduct(operator,state)
value = norm(operator*state,"fro") ...
    /max(norm(operator,"fro")*norm(state,"fro"),realmin);
end

function value = minimumValue(values)
if isempty(values)
    value = NaN;
else
    value = values(end);
end
end

function value = conditionNumber(values)
if isempty(values)
    value = NaN;
else
    value = values(1)/max(values(end),realmin);
end
end

function value = blockDefect(first,second)
value = norm(first-second,"fro") ...
    /max(norm(first,"fro")+norm(second,"fro"),realmin);
end
