function audit = buildHybridPrimitivePVOracleAudit(problem,primitive,zonalMode)
% Build and classify the fixed-zonal hybrid primitive-PV oracle.

positiveHorizontal = find(problem.horizontalLayout.kMode == zonalMode);
negativeHorizontal = conjugateHorizontalModes(problem,positiveHorizontal);
positive = hybridBlock(primitive,positiveHorizontal);
negative = hybridBlock(primitive,negativeHorizontal);

centeredCell = cell(1,numel(primitive.centeredTangents));
for iStep = 1:numel(centeredCell)
    centeredCell{iStep} = hybridTangent(positive,primitive.centeredTangents(iStep));
end
centered = [centeredCell{:}];
analyticCenteredDefect = zeros(numel(centered),1);
for iStep = 1:numel(centered)
    analyticCenteredDefect(iStep) = matrixDefect( ...
        centered(iStep).L1-positive.L1,positive.L1);
end

negativeColumnsInPositiveOrder = primitive.layout.admissibleRanges(negativeHorizontal);
negativeColumnsInPositiveOrder = horzcat(negativeColumnsInPositiveOrder{:});
conjugacyDefect = matrixDefect( ...
    negative.L1-conj(positive.L1),positive.L1);
if ~isequal(negative.columns,negativeColumnsInPositiveOrder)
    error("WVTerrainEnergyGalerkin:HybridOracleConjugateOrdering", ...
        "The conjugate fixed-zonal block does not retain the declared coordinate ordering.")
end

publicBridge = publicFlatWaveBridge(problem,positive,positiveHorizontal);
qgBridge = projectedQGBridge(problem,positive,positiveHorizontal,primitive.polynomialDegree);

requiredTolerance = struct("coordinateRcond",1e-10,"flatRecovery",1e-10, ...
    "tangent",1e-9,"selectedRows",1e-12,"primitive",1e-10, ...
    "qgClosure",1e-12,"publicFrequency",1e-8,"conjugacy",1e-10);
coordinatePasses = positive.numberOfRows == positive.numberOfCoordinates ...
    && positive.descriptorRank == positive.numberOfCoordinates ...
    && positive.descriptorReciprocalConditionNumber > requiredTolerance.coordinateRcond;
flatPasses = positive.diagnostics.flatRecoveryDefect <= requiredTolerance.flatRecovery;
tangentPasses = analyticCenteredDefect(end) <= requiredTolerance.tangent;
selectedRowsPass = positive.diagnostics.waveWeakDefect <= requiredTolerance.selectedRows ...
    && positive.diagnostics.projectedAPVDefect <= requiredTolerance.selectedRows ...
    && positive.diagnostics.bottomDefect <= requiredTolerance.selectedRows;
qgPasses = qgBridge.periodicOracleStatus == "passed" ...
    || qgBridge.periodicOracleStatus == "basis-bridge-blocker";
qgPasses = qgPasses ...
    && qgBridge.periodicCoreMaximumDefect <= requiredTolerance.qgClosure ...
    && qgBridge.projectedClosureDefect <= requiredTolerance.qgClosure;
publicPasses = publicBridge.frequencyDefect <= requiredTolerance.publicFrequency;
primitivePasses = positive.diagnostics.fullWeakDefect <= requiredTolerance.primitive ...
    && positive.diagnostics.energyDefect <= requiredTolerance.primitive ...
    && positive.diagnostics.fullAPVDefect <= requiredTolerance.primitive ...
    && positive.diagnostics.enstrophyDefect <= requiredTolerance.primitive;

if ~coordinatePasses
    status = "coordinate-rank-blocker";
    diagnosis = "The wave, projected-volume-APV, and bottom rows do not form a complete well-conditioned coordinate system.";
elseif ~flatPasses
    status = "flat-recovery-blocker";
    diagnosis = "The hybrid descriptor does not recover the independently derived flat primitive generator.";
elseif ~tangentPasses
    status = "tangent-implementation-blocker";
    diagnosis = "Analytic and independently centered hybrid terrain tangents do not agree.";
elseif ~selectedRowsPass
    status = "hybrid-row-blocker";
    diagnosis = "The hybrid generator does not satisfy one of its defining wave, projected-APV, or bottom rows.";
elseif ~primitivePasses
    status = "primitive-equivalence-blocker";
    diagnosis = "The complete hybrid coordinates close, but their generator is not equivalent to the unmodified primitive weak equations and therefore fails physical energy and full sampled APV.";
elseif ~qgPasses
    status = "qg-bridge-blocker";
    diagnosis = "The hybrid projected-PV rows do not retain the verified periodic QG closure.";
elseif ~publicPasses
    status = "public-basis-bridge-blocker";
    diagnosis = "The primitive flat wave rows are not resolved by the existing public wave basis.";
elseif conjugacyDefect > requiredTolerance.conjugacy
    status = "conjugacy-blocker";
    diagnosis = "The positive and negative fixed-zonal hybrid generators violate Fourier conjugacy.";
else
    status = "passed";
    diagnosis = "The hybrid primitive-PV descriptor reproduces the primitive weak equations and all required first-order identities.";
end

audit = struct;
audit.scope = "fixed-zonal-first-order-hybrid-primitive-pv-oracle";
audit.status = status;
audit.isCompatible = status == "passed";
audit.diagnosis = diagnosis;
audit.zonalMode = zonalMode;
audit.polynomialDegree = primitive.polynomialDegree;
audit.quadratureOrder = primitive.quadratureOrder;
audit.primitive = primitive;
audit.positiveBlock = positive;
audit.negativeBlock = negative;
audit.centeredTangents = centered;
audit.tangentAgreement = struct("relativeDefect",analyticCenteredDefect, ...
    "maximumFinalDefect",analyticCenteredDefect(end));
audit.conjugacyDefect = conjugacyDefect;
audit.publicBasisBridge = publicBridge;
audit.qgBridge = qgBridge;
audit.requiredTolerance = requiredTolerance;
audit.nextScope = "do-not-proceed-to-periodic-primitive-terrain";
end

function block = hybridBlock(primitive,horizontalIndices)
columns = primitive.layout.admissibleRanges(horizontalIndices);
columns = horzcat(columns{:});
E0 = primitive.flatReference.E(columns,columns);
J0 = primitive.flatReference.J(columns,columns);
L0Primitive = primitive.flatReference.L(columns,columns);
Q0 = primitive.flatReference.Q(:,columns);
B = primitive.flatReference.B(horizontalIndices,columns);
first = restrictTangent(primitive.analyticTangent,horizontalIndices,columns);

[Uq,Sq,~] = svd(Q0,"econ");
singularValues = diag(Sq);
rankTolerance = max(size(Q0))*eps(max(singularValues))*100;
qRank = nnz(singularValues > rankTolerance);
Uq = Uq(:,1:qRank);
projectedQ0 = Uq'*Q0;
projectedQ1 = Uq'*first.Q;

E0h = (E0+E0')/2;
J0h = (J0-J0')/2;
[waveVectors,frequency] = eig(1i*J0h,E0h,"vector");
frequency = real(frequency);
frequencyTolerance = 1e3*eps*max([abs(frequency);1]);
waveIndex = find(abs(frequency) > frequencyTolerance);
waveVectors = waveVectors(:,waveIndex);
waveFrequency = frequency(waveIndex);
[waveFrequency,order] = sort(waveFrequency);
waveVectors = waveVectors(:,order);
for iWave = 1:size(waveVectors,2)
    waveVectors(:,iWave) = waveVectors(:,iWave) ...
        /sqrt(real(waveVectors(:,iWave)'*E0h*waveVectors(:,iWave)));
end

n = numel(columns);
nBottom = size(B,1);
M0 = [waveVectors'*E0;projectedQ0;B];
K0 = [waveVectors'*J0;zeros(qRank,n);zeros(nBottom,n)];
M1 = [waveVectors'*first.E;projectedQ1;zeros(size(B))];
K1 = [waveVectors'*first.J;zeros(qRank,n);first.R];
descriptorRank = rank(M0,max(size(M0))*eps(max(norm(M0,2),1))*100);
if size(M0,1) == n && descriptorRank == n
    L0 = M0\K0;
    L1 = M0\(K1-M1*L0);
    reciprocalConditionNumber = rcond(M0);
else
    L0 = NaN(n);
    L1 = NaN(n);
    reciprocalConditionNumber = 0;
end

weakResidual = E0*L1+first.E*L0-first.J;
energyResidual = L1'*E0+E0*L1+L0'*first.E+first.E*L0;
projectedAPVResidual = projectedQ0*L1+projectedQ1*L0;
fullAPVResidual = Q0*L1+first.Q*L0;
bottomResidual = B*L1-first.R;
enstrophyResidual = L1'*primitive.flatReference.Z(columns,columns) ...
    +primitive.flatReference.Z(columns,columns)*L1 ...
    +L0'*first.Z+first.Z*L0;
waveWeakResidual = waveVectors'*weakResidual;

diagnostics = struct( ...
    "flatRecoveryDefect",productDefect(L0-L0Primitive,{L0,L0Primitive}), ...
    "waveWeakDefect",productDefect(waveWeakResidual, ...
        {waveVectors'*E0*L1,waveVectors'*first.E*L0,waveVectors'*first.J}), ...
    "fullWeakDefect",productDefect(weakResidual,{E0*L1,first.E*L0,first.J}), ...
    "energyDefect",productDefect(energyResidual, ...
        {L1'*E0,E0*L1,L0'*first.E,first.E*L0}), ...
    "projectedAPVDefect",productDefect(projectedAPVResidual, ...
        {projectedQ0*L1,projectedQ1*L0}), ...
    "fullAPVDefect",productDefect(fullAPVResidual,{Q0*L1,first.Q*L0}), ...
    "bottomDefect",productDefect(bottomResidual,{B*L1,first.R}), ...
    "enstrophyDefect",productDefect(enstrophyResidual, ...
        {L1'*primitive.flatReference.Z(columns,columns), ...
        primitive.flatReference.Z(columns,columns)*L1,L0'*first.Z,first.Z*L0}));
localRanges = cell(size(horizontalIndices));
iFirst = 1;
for i = 1:numel(horizontalIndices)
    count = numel(primitive.layout.admissibleRanges{horizontalIndices(i)});
    localRanges{i} = iFirst:(iFirst+count-1);
    iFirst = iFirst+count;
end
apvDefectByHorizontalMode = zeros(numel(horizontalIndices),1);
for i = 1:numel(horizontalIndices)
    current = localRanges{i};
    apvDefectByHorizontalMode(i) = productDefect(fullAPVResidual(:,current), ...
        {Q0*L1(:,current),first.Q*L0(:,current)});
end

block = struct("horizontalIndices",horizontalIndices,"columns",columns, ...
    "numberOfCoordinates",n,"numberOfRows",size(M0,1), ...
    "numberOfWaveRows",size(waveVectors,2),"numberOfProjectedAPVRows",qRank, ...
    "numberOfBottomRows",nBottom,"descriptorRank",descriptorRank, ...
    "descriptorReciprocalConditionNumber",reciprocalConditionNumber, ...
    "volumeAPVRangeBasis",Uq,"volumeAPVSingularValues",singularValues, ...
    "volumeAPVRankTolerance",rankTolerance,"projectedQ0",projectedQ0, ...
    "projectedQ1",projectedQ1,"bottomSelector",B, ...
    "waveVectors",waveVectors,"waveFrequency",waveFrequency, ...
    "M0",M0,"K0",K0,"M1",M1,"K1",K1,"L0",L0,"L1",L1, ...
    "primitiveFirstOrder",first, ...
    "residuals",struct("weak",weakResidual,"energy",energyResidual, ...
    "projectedAPV",projectedAPVResidual,"fullAPV",fullAPVResidual, ...
    "bottom",bottomResidual,"enstrophy",enstrophyResidual), ...
    "diagnostics",diagnostics);
block.diagnostics.fullAPVDefectByHorizontalMode = apvDefectByHorizontalMode;
end

function tangent = hybridTangent(block,centered)
first = restrictTangent(centered,block.horizontalIndices,block.columns);
M1 = [block.waveVectors'*first.E;block.volumeAPVRangeBasis'*first.Q; ...
    zeros(size(block.bottomSelector))];
K1 = [block.waveVectors'*first.J; ...
    zeros(block.numberOfProjectedAPVRows,block.numberOfCoordinates);first.R];
tangent = struct("step",centered.step,"M1",M1,"K1",K1, ...
    "L1",block.M0\(K1-M1*block.L0));
end

function first = restrictTangent(source,horizontalIndices,columns)
first = struct("E",source.E(columns,columns),"J",source.J(columns,columns), ...
    "Q",source.Q(:,columns),"Z",source.Z(columns,columns), ...
    "R",source.R(horizontalIndices,columns));
end

function partner = conjugateHorizontalModes(problem,indices)
partner = zeros(size(indices));
for i = 1:numel(indices)
    partner(i) = find(problem.horizontalLayout.kMode == -problem.horizontalLayout.kMode(indices(i)) ...
        & problem.horizontalLayout.lMode == -problem.horizontalLayout.lMode(indices(i)),1);
end
end

function bridge = publicFlatWaveBridge(problem,block,horizontalIndices)
publicFrequency = [];
for i = 1:numel(horizontalIndices)
    flat = problem.flatModeBlocks{horizontalIndices(i)};
    publicFrequency = [publicFrequency;flat.frequency(flat.waveIndices)]; %#ok<AGROW>
end
publicFrequency = sort(publicFrequency);
primitiveFrequency = sort(block.waveFrequency);
if numel(publicFrequency) == numel(primitiveFrequency)
    frequencyDefect = norm(publicFrequency-primitiveFrequency) ...
        /max(norm(publicFrequency),realmin);
else
    matched = zeros(size(primitiveFrequency));
    available = true(size(publicFrequency));
    for i = 1:numel(primitiveFrequency)
        candidates = find(available & sign(publicFrequency) == sign(primitiveFrequency(i)));
        if isempty(candidates)
            frequencyDefect = Inf;
            matched = [];
            break
        end
        [~,nearest] = min(abs(publicFrequency(candidates)-primitiveFrequency(i)));
        matched(i) = publicFrequency(candidates(nearest));
        available(candidates(nearest)) = false;
    end
    if ~isempty(matched)
        frequencyDefect = norm(matched-primitiveFrequency)/max(norm(matched),realmin);
    end
end
bridge = struct("publicFrequency",publicFrequency, ...
    "primitiveFrequency",primitiveFrequency,"frequencyDefect",frequencyDefect, ...
    "publicCoefficientLayoutUnchanged",true);
end

function bridge = projectedQGBridge(problem,block,horizontalIndices,primitiveDegree)
qg = problem.auditPeriodicCoupledPVOracle(polynomialDegree=primitiveDegree+2);
qgHorizontal = zeros(size(horizontalIndices));
for i = 1:numel(horizontalIndices)
    qgHorizontal(i) = find(qg.layout.kMode == problem.horizontalLayout.kMode(horizontalIndices(i)) ...
        & qg.layout.lMode == problem.horizontalLayout.lMode(horizontalIndices(i)),1);
end
qgColumns = qg.layout.ranges(qgHorizontal);
qgColumns = horzcat(qgColumns{:});
qgGenerator = qg.dynamics.generator(qgColumns,qgColumns);

[~,Sigma,V] = svd(block.L0);
singularValues = diag(Sigma);
rankTolerance = max(size(block.L0))*eps(max(singularValues))*100;
dynamicRank = nnz(singularValues > rankTolerance);
stationaryBasis = V(:,dynamicRank+1:end);
pvMap = [block.projectedQ0;block.bottomSelector];
pvCoordinateMap = pvMap*stationaryBasis;
if size(pvCoordinateMap,1) == size(pvCoordinateMap,2) ...
        && rank(pvCoordinateMap,max(size(pvCoordinateMap))*eps(max(norm(pvCoordinateMap,2),1))*100) == size(pvCoordinateMap,1)
    projectedGenerator = pvMap*block.L1*stationaryBasis/pvCoordinateMap;
    projectedClosureDefect = productDefect( ...
        block.projectedQ0*block.L1+block.projectedQ1*block.L0, ...
        {block.projectedQ0*block.L1,block.projectedQ1*block.L0});
else
    projectedGenerator = NaN(size(pvCoordinateMap,1));
    projectedClosureDefect = Inf;
end

hybridFrequency = sort(real(1i*eig(projectedGenerator)));
qgFrequency = sort(real(1i*eig(qgGenerator)));
if numel(hybridFrequency) == numel(qgFrequency)
    frequencyDefect = norm(hybridFrequency-qgFrequency)/max(norm(qgFrequency),realmin);
else
    frequencyDefect = Inf;
end
bridge = struct("periodicOracleStatus",qg.status, ...
    "projectedGenerator",projectedGenerator,"qgGenerator",qgGenerator, ...
    "hybridFrequency",hybridFrequency,"qgFrequency",qgFrequency, ...
    "frequencyDefect",frequencyDefect, ...
    "projectedClosureDefect",projectedClosureDefect, ...
    "stationaryDimension",size(stationaryBasis,2), ...
    "pvCoordinateDimension",size(pvCoordinateMap,1), ...
    "pvCoordinateRank",rank(pvCoordinateMap));
bridge.periodicCoreMaximumDefect = max([ ...
    qg.forms.diagnostics.maximumInversionDefect, ...
    qg.forms.diagnostics.maximumGreenIdentityDefect, ...
    qg.dynamics.diagnostics.energyDefect, ...
    qg.dynamics.diagnostics.apvDefect, ...
    qg.dynamics.diagnostics.enstrophyDefect, ...
    qg.dynamics.diagnostics.bottomDefect, ...
    qg.dynamics.diagnostics.conjugacyDefect]);
end

function value = productDefect(residual,terms)
scale = 0;
for i = 1:numel(terms)
    scale = scale+norm(terms{i},"fro");
end
value = norm(residual,"fro")/max(scale,realmin);
end

function value = matrixDefect(residual,reference)
value = norm(residual,"fro")/max(norm(reference,"fro"),realmin);
end
