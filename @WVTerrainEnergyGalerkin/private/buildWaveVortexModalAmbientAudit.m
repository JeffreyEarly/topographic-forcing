function audit = buildWaveVortexModalAmbientAudit(problem,trustedBounds, ...
    targetWaveIndices,waveCounts,geostrophicCounts,orders, ...
    referenceDegrees,paddingFactors,terrainScales,providerOrders, ...
    shouldRunTwoDimensionalControl,cacheDirectory,quadratureOrder,experiment)
% Build the Milestone-10.2.3 flat wave-vortex modal ambient oracle.

if ~isempty(experiment)
    audit = buildTerrainDressedRobinAudit(problem,trustedBounds, ...
        targetWaveIndices,waveCounts,geostrophicCounts,orders, ...
        referenceDegrees,paddingFactors,terrainScales,providerOrders, ...
        shouldRunTwoDimensionalControl,cacheDirectory,quadratureOrder, ...
        experiment);
    return
end

cache = emptyCacheDiagnostics(cacheDirectory);
terrain = terrainSupport(problem);
shells = scatteringShells(problem,terrain,orders,trustedBounds);
maximumLayout = shells.layouts{end};

[primary,record] = cachedAuditValue(problem,cacheDirectory, ...
    "zonal-primary",caseConfiguration("zonal-primary", ...
    referenceDegrees(end),maximumLayout,paddingFactors(1), ...
    terrainScales,quadratureOrder,max(waveCounts), ...
    max(geostrophicCounts),providerOrders), ...
    @() modalCase(problem,maximumLayout,trustedBounds, ...
    targetWaveIndices,waveCounts,geostrophicCounts, ...
    referenceDegrees(end),paddingFactors(1),terrainScales, ...
    providerOrders,quadratureOrder));
cache = appendCacheRecord(cache,record);

degreeCases = repmat(emptyComparisonCase,numel(referenceDegrees),1);
for iDegree = 1:numel(referenceDegrees)
    degree = referenceDegrees(iDegree);
    if degree == referenceDegrees(end)
        degreeCases(iDegree) = primary;
        continue
    end
    [degreeCases(iDegree),record] = cachedAuditValue(problem, ...
        cacheDirectory,"zonal-degree",caseConfiguration( ...
        "zonal-degree",degree,maximumLayout,paddingFactors(1), ...
        terrainScales,quadratureOrder,max(waveCounts), ...
        max(geostrophicCounts),providerOrders), ...
        @() modalCase(problem,maximumLayout,trustedBounds, ...
        targetWaveIndices,max(waveCounts),max(geostrophicCounts), ...
        degree,paddingFactors(1),terrainScales,providerOrders, ...
        quadratureOrder));
    cache = appendCacheRecord(cache,record);
end

paddingCases = repmat(emptyComparisonCase,numel(paddingFactors),1);
paddingCases(1) = primary;
for iPadding = 2:numel(paddingFactors)
    [paddingCases(iPadding),record] = cachedAuditValue(problem, ...
        cacheDirectory,"zonal-padding",caseConfiguration( ...
        "zonal-padding",referenceDegrees(end),maximumLayout, ...
        paddingFactors(iPadding),[0;1],quadratureOrder, ...
        max(waveCounts),max(geostrophicCounts),providerOrders), ...
        @() modalCase(problem,maximumLayout,trustedBounds, ...
        targetWaveIndices,max(waveCounts),max(geostrophicCounts), ...
        referenceDegrees(end),paddingFactors(iPadding),[0;1], ...
        providerOrders,quadratureOrder));
    cache = appendCacheRecord(cache,record);
end

supportCases = repmat(emptyComparisonCase,numel(shells.layouts),1);
for iSupport = 1:numel(shells.layouts)
    if iSupport == numel(shells.layouts)
        supportCases(iSupport) = primary;
        continue
    end
    layout = shells.layouts{iSupport};
    [supportCases(iSupport),record] = cachedAuditValue(problem, ...
        cacheDirectory,"zonal-support",caseConfiguration( ...
        "zonal-support",referenceDegrees(end),layout, ...
        paddingFactors(1),terrainScales,quadratureOrder, ...
        max(waveCounts),max(geostrophicCounts),providerOrders), ...
        @() modalCase(problem,layout,trustedBounds, ...
        targetWaveIndices,max(waveCounts),max(geostrophicCounts), ...
        referenceDegrees(end),paddingFactors(1),terrainScales, ...
        providerOrders,quadratureOrder));
    cache = appendCacheRecord(cache,record);
end

convergence = summarizeConvergence(primary,degreeCases,paddingCases, ...
    supportCases,waveCounts,geostrophicCounts);
tolerance = struct("provider",2e-11,"embedding",1e-11, ...
    "structure",1e-12,"stationary",1e-10,"ritz",1e-10, ...
    "projector",1e-8,"frequency",1e-8,"apv",1e-8, ...
    "bottom",1e-10,"strong",1e-5,"padding",1e-8, ...
    "support",1e-8,"reference",1e-8,"compression",2);
[classification,diagnosis,gates] = classify(primary,convergence,tolerance);

twoDimensional = struct("wasRequested", ...
    shouldRunTwoDimensionalControl,"wasAttempted",false, ...
    "status","not-authorized-by-zonal-result");
if shouldRunTwoDimensionalControl && classification == "wave-vortex-modal-acceleration"
    controlTrusted = [1 1];
    controlShells = scatteringShells(problem,terrain,orders,controlTrusted);
    controlLayout = controlShells.layouts{end};
    [control,record] = cachedAuditValue(problem,cacheDirectory, ...
        "two-dimensional-control",caseConfiguration( ...
        "two-dimensional-control",referenceDegrees(end), ...
        controlLayout,paddingFactors(1),terrainScales, ...
        quadratureOrder,max(waveCounts),max(geostrophicCounts), ...
        providerOrders),@() modalCase(problem,controlLayout, ...
        controlTrusted,targetWaveIndices(1),max(waveCounts), ...
        max(geostrophicCounts),referenceDegrees(end), ...
        paddingFactors(1),terrainScales,providerOrders,quadratureOrder));
    cache = appendCacheRecord(cache,record);
    twoDimensional = control;
    twoDimensional.wasRequested = true;
    twoDimensional.wasAttempted = true;
    if control.primary.maximumProjectorDefect > tolerance.projector ...
            || control.primary.maximumPhysicalDefect > tolerance.strong
        classification = "wave-modal-blocker";
        diagnosis = "The zonal modal ambient passes, but the required two-dimensional control does not reproduce the primitive internal-wave block.";
    end
end

audit = struct;
audit.scope = "milestone-10.2.3-flat-wave-vortex-modal-ambient";
audit.status = classification;
audit.classification = classification;
audit.isCompatible = ismember(classification, ...
    ["wave-vortex-modal-acceleration","wave-vortex-modal-equivalent"]);
audit.diagnosis = diagnosis;
audit.gates = gates;
audit.primary = primary;
audit.referenceDegreeConvergence = degreeCases;
audit.paddingConvergence = paddingCases;
audit.supportConvergence = supportCases;
audit.convergence = convergence;
audit.twoDimensionalControl = twoDimensional;
audit.shells = shells;
audit.cache = finalizeCacheDiagnostics(cache);
audit.trustedModeBounds = trustedBounds;
audit.targetWaveModeIndices = targetWaveIndices;
audit.waveGuardModeCounts = waveCounts;
audit.geostrophicModeCounts = geostrophicCounts;
audit.scatteringOrders = orders;
audit.primitiveReferenceDegrees = referenceDegrees;
audit.paddingFactors = paddingFactors;
audit.terrainScales = terrainScales;
audit.internalModesEVPOrders = providerOrders;
audit.internalModesEVPCommit = ...
    "df86687e91faa31bf65941299062d125a96904b1";
audit.requiredTolerance = tolerance;
audit.nextScope = "milestone-10.2.3-stop-for-review";
end

function audit = buildTerrainDressedRobinAudit(problem,trustedBounds, ...
    targetWaveIndices,waveCounts,geostrophicCounts,~,referenceDegrees, ...
    paddingFactors,terrainScales,providerOrders, ...
    shouldRunTwoDimensionalControl,cacheDirectory,quadratureOrder,experiment)
% Build the Milestone-10.2.4 terrain-dressed Robin stationary oracle.

cache = emptyCacheDiagnostics(cacheDirectory);
layouts = supportLayouts(problem,experiment.supportModeBounds);
trainingDegree = referenceDegrees(ceil(numel(referenceDegrees)/2));
trainingConfiguration = dressedCaseConfiguration("training", ...
    trainingDegree,layouts{end},2,terrainScales,quadratureOrder, ...
    waveCounts,geostrophicCounts,providerOrders, ...
    experiment.robinLengthRatios,experiment.trainingGeostrophicModeCount, ...
    experiment.tangentStep);
[training,record] = cachedAuditValue(problem,cacheDirectory, ...
    "dressed-robin-training",trainingConfiguration,@() trainingSweep( ...
    problem,layouts{end},trustedBounds,trainingDegree,2,terrainScales, ...
    geostrophicCounts,providerOrders,quadratureOrder,experiment));
cache = appendCacheRecord(cache,record);

selectedRatio = training.selectedRobinLengthRatio;
if ~training.hasEligibleCandidate
    audit = stoppedDressedRobinAudit(training,cache,trustedBounds, ...
        targetWaveIndices,waveCounts,geostrophicCounts,referenceDegrees, ...
        paddingFactors,terrainScales,providerOrders,experiment, ...
        "dressed-stationary-blocker", ...
        "No Robin candidate passes the provider, endpoint, rank, conjugacy, tangent-inclusion, and Green-identity training gates.");
    return
end

primaryConfiguration = dressedCaseConfiguration("zonal-primary", ...
    referenceDegrees(end),layouts{end},paddingFactors(1),terrainScales, ...
    quadratureOrder,waveCounts,geostrophicCounts,providerOrders, ...
    selectedRatio,experiment.trainingGeostrophicModeCount, ...
    experiment.tangentStep);
[primary,record] = cachedAuditValue(problem,cacheDirectory, ...
    "dressed-robin-primary",primaryConfiguration,@() dressedModalCase( ...
    problem,layouts{end},trustedBounds,targetWaveIndices,waveCounts, ...
    geostrophicCounts,referenceDegrees(end),paddingFactors(1), ...
    terrainScales,providerOrders,quadratureOrder,selectedRatio, ...
    experiment.tangentStep));
cache = appendCacheRecord(cache,record);

degreeCases = repmat(emptyDressedComparisonCase,numel(referenceDegrees),1);
for iDegree = 1:numel(referenceDegrees)
    degree = referenceDegrees(iDegree);
    if degree == referenceDegrees(end)
        degreeCases(iDegree) = primary;
        continue
    end
    configuration = dressedCaseConfiguration("zonal-degree",degree, ...
        layouts{end},paddingFactors(1),terrainScales,quadratureOrder, ...
        waveCounts,geostrophicCounts,providerOrders,selectedRatio, ...
        experiment.trainingGeostrophicModeCount,experiment.tangentStep);
    [degreeCases(iDegree),record] = cachedAuditValue(problem, ...
        cacheDirectory,"dressed-robin-degree",configuration, ...
        @() dressedModalCase(problem,layouts{end},trustedBounds, ...
        targetWaveIndices,waveCounts,geostrophicCounts,degree, ...
        paddingFactors(1),terrainScales,providerOrders,quadratureOrder, ...
        selectedRatio,experiment.tangentStep));
    cache = appendCacheRecord(cache,record);
end

paddingCases = repmat(emptyDressedComparisonCase,numel(paddingFactors),1);
paddingCases(1) = primary;
for iPadding = 2:numel(paddingFactors)
    configuration = dressedCaseConfiguration("zonal-padding", ...
        referenceDegrees(end),layouts{end},paddingFactors(iPadding), ...
        [0;1],quadratureOrder,waveCounts,geostrophicCounts,providerOrders, ...
        selectedRatio,experiment.trainingGeostrophicModeCount, ...
        experiment.tangentStep);
    [paddingCases(iPadding),record] = cachedAuditValue(problem, ...
        cacheDirectory,"dressed-robin-padding",configuration, ...
        @() dressedModalCase(problem,layouts{end},trustedBounds, ...
        targetWaveIndices,waveCounts,geostrophicCounts, ...
        referenceDegrees(end),paddingFactors(iPadding),[0;1], ...
        providerOrders,quadratureOrder,selectedRatio, ...
        experiment.tangentStep));
    cache = appendCacheRecord(cache,record);
end

supportCases = repmat(emptyDressedComparisonCase,numel(layouts),1);
for iSupport = 1:numel(layouts)
    if iSupport == numel(layouts)
        supportCases(iSupport) = primary;
        continue
    end
    configuration = dressedCaseConfiguration("zonal-support", ...
        referenceDegrees(end),layouts{iSupport},paddingFactors(1), ...
        terrainScales,quadratureOrder,waveCounts,geostrophicCounts, ...
        providerOrders,selectedRatio,experiment.trainingGeostrophicModeCount, ...
        experiment.tangentStep);
    [supportCases(iSupport),record] = cachedAuditValue(problem, ...
        cacheDirectory,"dressed-robin-support",configuration, ...
        @() dressedModalCase(problem,layouts{iSupport},trustedBounds, ...
        targetWaveIndices,waveCounts,geostrophicCounts, ...
        referenceDegrees(end),paddingFactors(1),terrainScales, ...
        providerOrders,quadratureOrder,selectedRatio, ...
        experiment.tangentStep));
    cache = appendCacheRecord(cache,record);
end

convergence = summarizeDressedConvergence(primary,degreeCases, ...
    paddingCases,supportCases,waveCounts,geostrophicCounts);
tolerance = struct("provider",1e-10,"centered",1e-9, ...
    "structure",1e-11,"green",1e-10,"stationary",1e-10, ...
    "stationaryRepresentation",1e-8,"projector",1e-8, ...
    "frequency",1e-8,"ritz",1e-10,"apv",1e-8, ...
    "bottom",1e-10,"strong",1e-5,"convergence",1e-8, ...
    "compression",2);
[classification,diagnosis,gates] = classifyDressedRobin( ...
    training,primary,convergence,tolerance);

twoDimensional = struct("wasRequested", ...
    shouldRunTwoDimensionalControl,"wasAttempted",false, ...
    "status","not-authorized-by-zonal-result");
if shouldRunTwoDimensionalControl && ismember(classification, ...
        ["dressed-robin-acceleration","dressed-flat-acceleration"])
    controlTrusted = [1 1];
    configuration = dressedCaseConfiguration("two-dimensional-control", ...
        referenceDegrees(end),layouts{end},paddingFactors(1),terrainScales, ...
        quadratureOrder,waveCounts(1),geostrophicCounts,providerOrders, ...
        selectedRatio,experiment.trainingGeostrophicModeCount, ...
        experiment.tangentStep);
    [control,record] = cachedAuditValue(problem,cacheDirectory, ...
        "dressed-robin-two-dimensional",configuration, ...
        @() dressedModalCase(problem,layouts{end},controlTrusted, ...
        targetWaveIndices(1),waveCounts(1),geostrophicCounts, ...
        referenceDegrees(end),paddingFactors(1),terrainScales, ...
        providerOrders,quadratureOrder,selectedRatio, ...
        experiment.tangentStep));
    cache = appendCacheRecord(cache,record);
    twoDimensional = control;
    twoDimensional.wasRequested = true;
    twoDimensional.wasAttempted = true;
    if ~dressedResultPasses(control.primary,tolerance)
        classification = "dressed-stationary-blocker";
        diagnosis = "The zonal dressed stationary ambient passes, but the required two-dimensional control does not reproduce the primitive physical block.";
    end
end

audit = struct;
audit.scope = "milestone-10.2.4-terrain-dressed-robin-stationary-ambient";
audit.status = classification;
audit.classification = classification;
audit.isCompatible = ismember(classification, ...
    ["dressed-robin-acceleration","dressed-flat-acceleration", ...
    "dressed-stationary-equivalent"]);
audit.diagnosis = diagnosis;
audit.gates = gates;
audit.training = training;
audit.selectedRobinLengthRatio = selectedRatio;
audit.primary = primary;
audit.referenceDegreeConvergence = degreeCases;
audit.paddingConvergence = paddingCases;
audit.supportConvergence = supportCases;
audit.convergence = convergence;
audit.twoDimensionalControl = twoDimensional;
audit.cache = finalizeCacheDiagnostics(cache);
audit.trustedModeBounds = trustedBounds;
audit.supportModeBounds = experiment.supportModeBounds;
audit.targetWaveModeIndices = targetWaveIndices;
audit.waveGuardModeCounts = waveCounts;
audit.geostrophicModeCounts = geostrophicCounts;
audit.robinLengthRatios = experiment.robinLengthRatios;
audit.trainingGeostrophicModeCount = ...
    experiment.trainingGeostrophicModeCount;
audit.primitiveReferenceDegrees = referenceDegrees;
audit.paddingFactors = paddingFactors;
audit.terrainScales = terrainScales;
audit.internalModesEVPOrders = providerOrders;
audit.internalModesEVPCommit = ...
    "df86687e91faa31bf65941299062d125a96904b1";
audit.tangentStep = experiment.tangentStep;
audit.requiredTolerance = tolerance;
audit.nextScope = "milestone-10.3-only-after-dressed-acceleration";
end

function audit = stoppedDressedRobinAudit(training,cache,trustedBounds, ...
    targetWaveIndices,waveCounts,geostrophicCounts,referenceDegrees, ...
    paddingFactors,terrainScales,providerOrders,experiment,status,diagnosis)
audit = struct("scope", ...
    "milestone-10.2.4-terrain-dressed-robin-stationary-ambient", ...
    "status",status,"classification",status,"isCompatible",false, ...
    "diagnosis",diagnosis,"gates",struct,"training",training, ...
    "selectedRobinLengthRatio",training.selectedRobinLengthRatio, ...
    "primary",emptyDressedComparisonCase, ...
    "referenceDegreeConvergence",repmat(emptyDressedComparisonCase,0,1), ...
    "paddingConvergence",repmat(emptyDressedComparisonCase,0,1), ...
    "supportConvergence",repmat(emptyDressedComparisonCase,0,1), ...
    "convergence",struct,"twoDimensionalControl",struct( ...
    "wasRequested",false,"wasAttempted",false,"status","not-attempted"), ...
    "cache",finalizeCacheDiagnostics(cache), ...
    "trustedModeBounds",trustedBounds, ...
    "supportModeBounds",experiment.supportModeBounds, ...
    "targetWaveModeIndices",targetWaveIndices, ...
    "waveGuardModeCounts",waveCounts, ...
    "geostrophicModeCounts",geostrophicCounts, ...
    "robinLengthRatios",experiment.robinLengthRatios, ...
    "trainingGeostrophicModeCount", ...
    experiment.trainingGeostrophicModeCount, ...
    "primitiveReferenceDegrees",referenceDegrees, ...
    "paddingFactors",paddingFactors,"terrainScales",terrainScales, ...
    "internalModesEVPOrders",providerOrders, ...
    "internalModesEVPCommit", ...
    "df86687e91faa31bf65941299062d125a96904b1", ...
    "tangentStep",experiment.tangentStep,"requiredTolerance",struct, ...
    "nextScope","milestone-10.2.4-stop-for-review");
end

function training = trainingSweep(problem,layout,trusted,degree,padding, ...
    scales,~,providerOrders,quadratureOrder,experiment)
[primitive,c] = buildGlobalSmallTerrainPrimitiveAudit(problem,degree, ...
    quadratureOrder,experiment.tangentStep,horizontalLayout=layout, ...
    paddingFactor=padding,trustedModeBounds=trusted, ...
    rejectTerrainNyquist=true,evaluationScales=scales, ...
    shouldAuditTangent=true);
direction = primitive.evaluatedDirections(end);
count = experiment.trainingGeostrophicModeCount;
ratios = experiment.robinLengthRatios;
results = repmat(emptyStationaryTrainingResult,numel(ratios),3);
exact = constructCompleteStationarySpace(c,direction,trusted,degree, ...
    min(4,degree),1);
ordinaryCatalog = verticalCatalog(problem,c,1,count,providerOrders,Inf);
ordinaryCandidate = stationaryCandidateData( ...
    c,ordinaryCatalog,count,experiment.tangentStep);
ordinaryResult = stationaryCandidateResult(c,primitive,direction, ...
    ordinaryCandidate,exact,1,"ordinary-flat",experiment.tangentStep);
for iRatio = 1:numel(ratios)
    robinLength = ratios(iRatio)*c.wvt.Lz;
    catalog = verticalCatalog(problem,c,1,count,providerOrders,robinLength);
    candidate = stationaryCandidateData( ...
        c,catalog,count,experiment.tangentStep);
    results(iRatio,1) = ordinaryResult;
    results(iRatio,2) = stationaryCandidateResult(c,primitive,direction, ...
        candidate,exact,1,"robin-undressed",experiment.tangentStep);
    results(iRatio,3) = stationaryCandidateResult(c,primitive,direction, ...
        candidate,exact,1,"robin-dressed",experiment.tangentStep);
end

eligible = false(numel(ratios),1);
score = Inf(numel(ratios),3);
for iRatio = 1:numel(ratios)
    value = results(iRatio,3);
    eligible(iRatio) = value.providerPasses && value.fullRank ...
        && value.maximumCenteredDerivativeDefect <= 1e-9 ...
        && value.maximumGreenIdentityDefect <= 1e-10 ...
        && value.dressedTangentDefect <= 1e-10 ...
        && value.bottomTangencyDefect <= 1e-10 ...
        && value.conjugacyDefect <= 1e-11;
    score(iRatio,:) = [value.stationaryRepresentationDefect, ...
        value.gramConditionNumber,iRatio];
end
indices = find(eligible);
if isempty(indices)
    [~,order] = sortrows(score,[1 2 3]);
    selected = order(1);
    selectedRatio = NaN;
    hasEligible = false;
else
    [~,order] = sortrows(score(indices,:),[1 2 3]);
    selected = indices(order(1));
    selectedRatio = ratios(selected);
    hasEligible = true;
end
training = struct("ratios",ratios,"results",results, ...
    "eligible",eligible,"score",score,"selectedIndex",selected, ...
    "selectedRobinLengthRatio",selectedRatio, ...
    "bestUnqualifiedIndex",selected, ...
    "bestUnqualifiedRobinLengthRatio",ratios(selected), ...
    "hasEligibleCandidate",hasEligible, ...
    "geostrophicModeCount",count,"primitiveReferenceDegree",degree, ...
    "paddingFactor",padding,"terrainScale",1, ...
    "horizontalModes",[layout.kMode layout.lMode]);
end

function candidate = stationaryCandidateData(c,catalog,count,tangentStep)
[G0,G1,scalar,dressedDiagnostics] = dressedGeostrophicColumns( ...
    c,catalog,count,[Inf Inf],tangentStep);
candidate = struct("catalog",catalog,"G0",G0,"G1",G1, ...
    "scalar",scalar,"dressedDiagnostics",dressedDiagnostics, ...
    "bottom",bottomInversionColumns(c,catalog,[Inf Inf]), ...
    "mean",stationaryMeanColumns(c,catalog,count));
end

function result = stationaryCandidateResult(c,primitive,direction,candidate, ...
    exact,terrainScale,ablation,tangentStep)
G0 = candidate.G0;
G1 = candidate.G1;
scalar = candidate.scalar;
dressedDiagnostics = candidate.dressedDiagnostics;
catalogUsed = candidate.catalog;
[dressed,tangentDiagnostics] = dressedStationaryCoordinates( ...
    G0,G1,primitive.analyticTangent,terrainScale);
switch ablation
    case "ordinary-flat"
        geostrophic = G0;
    case "robin-undressed"
        geostrophic = G0;
    case "robin-dressed"
        geostrophic = dressed;
    otherwise
        error("WVTerrainEnergyGalerkin:UnknownDressedRobinAblation", ...
            "Unknown stationary-coordinate ablation '%s'.",ablation)
end
raw = [geostrophic candidate.bottom candidate.mean];
[X,rankDiagnostics] = energyIndependentColumns(raw,direction.E);
H = X'*direction.E*X;
J = X'*direction.J*X;
[G,stationaryDiagnostics] = stationarySubspace( ...
    exact.trustedBasis,direction,X,H,J);
conjugacy = subspaceConjugacyDefect(X,direction.E, ...
    c.coordinateConjugateIndex);
green = tangentRobinGreenIdentity(c,primitive,G0,G1,scalar, ...
    tangentDiagnostics.tangentBasis,tangentStep);
singularValues = rankDiagnostics.singularValues;
if isempty(singularValues)
    conditionNumber = Inf;
else
    conditionNumber = singularValues(1) ...
        /max(singularValues(end),realmin);
end
provider = catalogUsed.diagnostics;
result = emptyStationaryTrainingResult;
result.ablation = ablation;
result.robinLengthRatio = catalogUsed.diagnostics.robinLength/c.wvt.Lz;
result.numberOfRawCoordinates = size(raw,2);
result.numberOfCoordinates = size(X,2);
result.numberOfStationaryCoordinates = size(G,2);
result.fullRank = rankDiagnostics.rank == rankDiagnostics.nominalDimension;
result.gramConditionNumber = conditionNumber^2;
result.stationaryRepresentationDefect = ...
    stationaryDiagnostics.representationDefect;
result.stationaryRowDefect = stationaryDiagnostics.rowDefect;
result.bottomTangencyDefect = stationaryDiagnostics.bottomTangencyDefect;
result.dressedTangentDefect = tangentDiagnostics.bottomTangencyDefect;
result.numberOfTangentCoordinates = ...
    tangentDiagnostics.numberOfTangentCoordinates;
result.conjugacyDefect = conjugacy;
result.maximumCenteredDerivativeDefect = ...
    dressedDiagnostics.maximumCenteredDerivativeDefect;
result.maximumEmbeddingDefect = dressedDiagnostics.maximumEmbeddingDefect;
result.maximumGreenIdentityDefect = green.maximumDefect;
result.greenIdentity = green;
result.geostrophicDiagnostics = dressedDiagnostics;
result.geostrophicDiagnostics.tangent = tangentDiagnostics;
result.provider = provider;
result.providerPasses = provider.providerDifference <= 1e-10 ...
    && provider.maximumEndpointDefect <= 1e-10 ...
    && provider.negativeRobinModeCountDefect == 0;
result.rankDiagnostics = rankDiagnostics;
end

function result = dressedModalCase(problem,layout,trusted,targetWaveIndices, ...
    waveCounts,geostrophicCounts,degree,padding,scales,providerOrders, ...
    quadratureOrder,robinLengthRatio,tangentStep)
[primitive,c] = buildGlobalSmallTerrainPrimitiveAudit(problem,degree, ...
    quadratureOrder,tangentStep,horizontalLayout=layout, ...
    paddingFactor=padding,trustedModeBounds=trusted, ...
    rejectTerrainNyquist=true,evaluationScales=scales, ...
    shouldAuditTangent=true);
catalog = verticalCatalog(problem,c,max(waveCounts), ...
    max(geostrophicCounts),providerOrders,robinLengthRatio*c.wvt.Lz);
seed = waveColumns(c,catalog,targetWaveIndices,trusted);
flatSeed = subspaceRitz(seed,primitive.evaluatedDirections(1));
reference = continueReference(c,primitive,seed,scales,trusted,degree);
[G0All,G1All,~,geostrophicDiagnostics] = dressedGeostrophicColumns( ...
    c,catalog,max(geostrophicCounts),[Inf Inf],tangentStep);
exactEnd = constructCompleteStationarySpace(c, ...
    primitive.evaluatedDirections(end),trusted,degree,min(4,degree), ...
    scales(end));

grid = repmat(emptyModalResult,numel(waveCounts), ...
    numel(geostrophicCounts));
for iWave = 1:numel(waveCounts)
    for iGeostrophic = 1:numel(geostrophicCounts)
        grid(iWave,iGeostrophic) = dressedReducedModalResult(c, ...
            primitive.evaluatedDirections(end),reference,catalog, ...
            waveCounts(iWave),geostrophicCounts(iGeostrophic),trusted, ...
            scales(end),primitive.analyticTangent,G0All,G1All, ...
            geostrophicDiagnostics,exactEnd);
    end
end
primary = selectDressedPrimaryResult(grid);
amplitude = repmat(emptyModalResult,numel(scales),1);
for iScale = 1:numel(scales)
    exactScale = constructCompleteStationarySpace(c, ...
        primitive.evaluatedDirections(iScale),trusted,degree, ...
        min(4,degree),scales(iScale));
    amplitude(iScale) = dressedReducedModalResult(c, ...
        primitive.evaluatedDirections(iScale),reference.history(iScale), ...
        catalog,primary.waveModeCount,primary.geostrophicModeCount, ...
        trusted,scales(iScale),primitive.analyticTangent, ...
        G0All,G1All,geostrophicDiagnostics,exactScale);
end
exactFlat = constructCompleteStationarySpace(c, ...
    primitive.evaluatedDirections(1),trusted,degree,min(4,degree),0);
flatControl = dressedReducedModalResult(c, ...
    primitive.evaluatedDirections(1),reference.history(1),catalog, ...
    max(targetWaveIndices),max(1,min(geostrophicCounts)),trusted,0, ...
    primitive.analyticTangent,G0All,G1All, ...
    geostrophicDiagnostics,exactFlat);

result = emptyDressedComparisonCase;
result.degree = degree;
result.paddingFactor = padding;
result.horizontalModes = [layout.kMode layout.lMode];
result.numberOfPrimitiveCoordinates = size(primitive.evaluatedDirections(end).E,1);
result.provider = catalog.diagnostics;
result.reference = reference;
result.flatSeed = flatSeed;
result.grid = grid;
result.primary = primary;
result.amplitude = amplitude;
result.flatControl = flatControl;
result.contextSummary = contextSummary(c);
result.robinLengthRatio = robinLengthRatio;
end

function result = dressedReducedModalResult(c,direction,reference,catalog, ...
    waveCount,geostrophicCount,trusted,terrainScale,first,G0All, ...
    G1All,geostrophicDiagnostics,exactStationary)
waves = waveColumns(c,catalog,1:waveCount,[Inf Inf]);
keepGeostrophic = geostrophicDiagnostics.modeIndex <= geostrophicCount;
G0 = G0All(:,keepGeostrophic);
G1 = G1All(:,keepGeostrophic);
[geostrophic,tangentDiagnostics] = dressedStationaryCoordinates( ...
    G0,G1,first,terrainScale);
bottom = bottomInversionColumns(c,catalog,[Inf Inf]);
mean = meanSectorColumns(c,catalog,waveCount,geostrophicCount);
raw = [waves geostrophic bottom mean];
[X,rankDiagnostics] = energyIndependentColumns(raw,direction.E);
H = X'*direction.E*X;
J = X'*direction.J*X;
[G,stationaryDiagnostics] = stationarySubspace( ...
    exactStationary.trustedBasis,direction,X,H,J);
dense = reducedPhysicalEigensystem(X,H,J,G,direction.E);
selected = selectByOverlap(dense.vectors,reference.basis, ...
    direction.E,size(reference.basis,2));
selected = closePhysicalBlock(selected,dense,direction.E, ...
    c.coordinateConjugateIndex);
basis = energyNormalizeColumns(dense.vectors(:,selected),direction.E);
frequency = columnRayleighQuotients(basis,direction.E,1i*direction.J);
physical = physicalDiagnostics(c,direction,basis,frequency);

result = emptyModalResult;
result.waveModeCount = waveCount;
result.geostrophicModeCount = geostrophicCount;
result.numberOfModalCoordinates = size(X,2);
result.numberOfPrimitiveCoordinates = size(direction.E,1);
result.compressionFactor = size(direction.E,1)/size(X,2);
result.numberOfStationaryCoordinates = size(G,2);
result.numberOfDynamicalCoordinates = size(X,2)-size(G,2);
result.numberOfSelectedCoordinates = size(basis,2);
result.basis = basis;
result.frequency = frequency;
result.maximumProjectorDefect = maximumPrincipalSine( ...
    basis,reference.basis,direction.E);
result.maximumFrequencyDefect = sortedFrequencyDefect( ...
    frequency,reference.frequency);
result.maximumRitzResidual = physical.maximumRitzResidual;
result.maximumAPVDefect = physical.maximumAPVDefect;
result.maximumBottomDefect = physical.maximumBottomDefect;
result.maximumStrongResidual = physical.maximumStrongResidual;
result.maximumPhysicalDefect = max([physical.maximumRitzResidual, ...
    physical.maximumAPVDefect,physical.maximumBottomDefect, ...
    physical.maximumStrongResidual]);
result.energyHermitianDefect = norm(H-H',"fro") ...
    /max(norm(H,"fro"),realmin);
result.exchangeSkewHermitianDefect = norm(J+J',"fro") ...
    /max(norm(J,"fro"),realmin);
result.stationaryRowDefect = stationaryDiagnostics.rowDefect;
result.stationaryBottomTangencyDefect = ...
    stationaryDiagnostics.bottomTangencyDefect;
result.stationaryRepresentationDefect = ...
    stationaryDiagnostics.representationDefect;
result.maximumEmbeddingDefect = max( ...
    [geostrophicDiagnostics.maximumEmbeddingDefect, ...
    rankDiagnostics.maximumRejectedRelativeSingularValue]);
result.rankDiagnostics = rankDiagnostics;
result.geostrophicDiagnostics = geostrophicDiagnostics;
result.geostrophicDiagnostics.tangent = tangentDiagnostics;
result.trustedGeostrophicDiagnostics = ...
    exactStationary.trustedRepresentation;
result.trustedModeBounds = trusted;
result.terrainScale = terrainScale;
end

function [G0,G1,scalar,diagnostics] = dressedGeostrophicColumns( ...
    c,catalog,count,bounds,tangentStep)
selected = find(abs(c.horizontalLayout.kMode) <= bounds(1) ...
    & abs(c.horizontalLayout.lMode) <= bounds(2) ...
    & ~(c.horizontalLayout.kMode == 0 & c.horizontalLayout.lMode == 0));
nColumns = numel(selected)*count;
nGrid = c.nXY*c.nZ;
values = zeros(nGrid,nColumns);
bottomValues = zeros(c.nXY,nColumns);
modeIndex = zeros(nColumns,1);
G0 = zeros(size(c.N,2),nColumns);
G1 = zeros(size(c.N,2),nColumns);
centeredNumerator = 0;
centeredDenominator = 0;
embeddingDefect = 0;
column = 0;
N20 = c.N20(:);
H3 = repmat(c.H,c.nZ,1);
HX3 = repmat(c.HX,c.nZ,1);
HY3 = repmat(c.HY,c.nZ,1);
N2 = kron(c.N20,ones(c.nXY,1));
N2Xi = kron(c.N2Xi,ones(c.nXY,1));
embedding = globalEmbeddingOperators(c);
for iHorizontal = reshape(selected,1,[])
    k = c.horizontalLayout.k(iHorizontal);
    l = c.horizontalLayout.l(iHorizontal);
    kappa = hypot(k,l);
    family = catalog.family{find(abs(catalog.kappa-kappa) ...
        <=100*eps*max(kappa,1),1)};
    for iMode = 1:count
        column = column+1;
        modeIndex(column) = iMode;
        F = family.robin.F(:,iMode);
        Fxi = -N20.*family.robin.eta(:,iMode)/c.wvt.f;
        value = kron(F,c.phase(:,iHorizontal));
        valueXi = kron(Fxi,c.phase(:,iHorizontal));
        valueX = 1i*k*value;
        valueY = 1i*l*value;
        values(:,column) = value;
        [~,iBottom] = min(c.quadrature.xi);
        bottomValues(:,column) = c.phase(:,iHorizontal)*F(iBottom);

        eta0 = -c.wvt.f*valueXi./N2;
        fields0 = struct("u",-valueY,"v",valueX, ...
            "w",zeros(size(value)),"eta",eta0);
        fields1 = struct( ...
            "u",H3.*valueY-c.xiGrid.*HY3.*valueXi, ...
            "v",-H3.*valueX+c.xiGrid.*HX3.*valueXi, ...
            "w",c.xiGrid.*(-HX3.*valueY+HY3.*valueX), ...
            "eta",H3.*(1+c.xiGrid.*N2Xi./N2).*eta0);
        [G0(:,column),flatEmbedding] = embedGlobalState( ...
            c,fields0,embedding);
        [G1(:,column),terrainEmbedding] = embedGlobalState( ...
            c,fields1,embedding);
        scalarColumn = struct("values",value,"valuesXi",valueXi, ...
            "valuesX",valueX,"valuesY",valueY);
        [plus,plusEmbedding] = embedGlobalState(c, ...
            exactMappedGeostrophicFields(c,scalarColumn,tangentStep), ...
            embedding);
        [minus,minusEmbedding] = embedGlobalState(c, ...
            exactMappedGeostrophicFields(c,scalarColumn,-tangentStep), ...
            embedding);
        centered = (plus-minus)/(2*tangentStep);
        centeredNumerator = centeredNumerator ...
            +norm(centered-G1(:,column))^2;
        centeredDenominator = centeredDenominator ...
            +(norm(centered)+norm(G1(:,column)))^2;
        embeddingDefect = max([embeddingDefect, ...
            flatEmbedding.maximumDefect,terrainEmbedding.maximumDefect, ...
            plusEmbedding.maximumDefect,minusEmbedding.maximumDefect]);
    end
end
scalar = struct("values",values,"bottomValues",bottomValues);
centeredDefect = sqrt(centeredNumerator) ...
    /max(sqrt(centeredDenominator),realmin);
diagnostics = struct("maximumCenteredDerivativeDefect",centeredDefect, ...
    "maximumEmbeddingDefect",embeddingDefect,"modeIndex",modeIndex);
end

function fields = exactMappedGeostrophicFields(c,scalar,delta)
H3 = repmat(c.H,c.nZ,1);
HX3 = repmat(c.HX,c.nZ,1);
HY3 = repmat(c.HY,c.nZ,1);
gamma = 1-delta*H3;
gammaX = -delta*HX3;
gammaY = -delta*HY3;
z = gamma.*c.xiGrid;
N2 = c.wvt.N2Function(z);
if isscalar(N2)
    N2 = repmat(N2,size(z));
end
fields = struct( ...
    "u",-gamma.*scalar.valuesY+c.xiGrid.*gammaY.*scalar.valuesXi, ...
    "v",gamma.*scalar.valuesX-c.xiGrid.*gammaX.*scalar.valuesXi, ...
    "w",c.xiGrid.*(gammaX.*scalar.valuesY-gammaY.*scalar.valuesX), ...
    "eta",-c.wvt.f*scalar.valuesXi./(gamma.*N2));
end

function [coordinates,diagnostics] = embedGlobalState(c,fields,embedding)
raw = zeros(c.nX,size(fields.u,2));
raw(embedding.uColumns,:) = embedding.projectU*fields.u;
raw(embedding.vColumns,:) = embedding.projectV*fields.v;
raw(embedding.wColumns,:) = embedding.projectW*fields.w;
raw(embedding.etaColumns,:) = embedding.projectEta*fields.eta;
coordinates = c.N'*raw;
representedRaw = c.N*coordinates;
represented = {c.Ru*representedRaw,c.Rv*representedRaw, ...
    c.Rwh*representedRaw,c.Reta*representedRaw};
targets = {fields.u,fields.v,fields.w,fields.eta};
defects = zeros(4,1);
for iField = 1:4
    numerator = norm(represented{iField}-targets{iField},"fro");
    denominator = norm(targets{iField},"fro") ...
        +norm(represented{iField},"fro");
    if denominator == 0
        defects(iField) = 0;
    else
        defects(iField) = numerator/denominator;
    end
end

diagnostics = struct("fieldDefects",defects, ...
    "maximumDefect",norm(vertcat(represented{:})-vertcat(targets{:}),"fro") ...
        /max(norm(vertcat(represented{:}),"fro") ...
        +norm(vertcat(targets{:}),"fro"),realmin), ...
    "continuityDefect",norm(c.continuity*representedRaw,"fro") ...
    /max(norm(c.continuity,"fro")*norm(representedRaw,"fro"),realmin));
end

function embedding = globalEmbeddingOperators(c)
[uColumns,vColumns,wColumns,etaColumns] = globalComponentColumns(c);
embedding = struct("uColumns",uColumns,"vColumns",vColumns, ...
    "wColumns",wColumns,"etaColumns",etaColumns, ...
    "projectU",weightedProjectionMatrix( ...
    c.Ru(:,uColumns),c.volumeWeight), ...
    "projectV",weightedProjectionMatrix( ...
    c.Rv(:,vColumns),c.volumeWeight), ...
    "projectW",weightedProjectionMatrix( ...
    c.Rwh(:,wColumns),c.volumeWeight), ...
    "projectEta",weightedProjectionMatrix( ...
    c.Reta(:,etaColumns),c.volumeWeight));
end

function matrix = weightedProjectionMatrix(reconstruction,weight)
matrix = (reconstruction'*(weight.*reconstruction)) ...
    \(reconstruction'.*weight.');
end

function [uColumns,vColumns,wColumns,etaColumns] = globalComponentColumns(c)
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

function [columns,diagnostics] = dressedStationaryCoordinates( ...
    G0,G1,first,terrainScale)
bottomMap = first.R*G0;
scale = norm(bottomMap,2);
if scale == 0
    tangent = eye(size(G0,2));
    complement = zeros(size(G0,2),0);
    singularValues = zeros(0,1);
    rankValue = 0;
else
    tolerance = max(size(bottomMap))*eps(scale)*100;
    tangent = null(bottomMap,tolerance);
    complement = orth(bottomMap');
    singularValues = svd(bottomMap);
    rankValue = nnz(singularValues > tolerance);
end
columns = [(G0+terrainScale*G1)*tangent G0*complement];
tangentResidual = first.R*(G0*tangent);
diagnostics = struct("tangentBasis",tangent, ...
    "complementBasis",complement, ...
    "numberOfTangentCoordinates",size(tangent,2), ...
    "numberOfComplementCoordinates",size(complement,2), ...
    "bottomMapRank",rankValue,"bottomMapSingularValues",singularValues, ...
    "bottomTangencyDefect",norm(tangentResidual,"fro") ...
        /max(norm(first.R,"fro")*norm(G0*tangent,"fro"),realmin));
end

function diagnostics = tangentRobinGreenIdentity(c,primitive,G0,G1, ...
    scalar,tangent,tangentStep)
flat = primitive.flatReference;
first = primitive.analyticTangent;
weak0 = tangent'*(G0'*flat.E);
weak1 = tangent'*(G1'*flat.E+G0'*first.E);
strong0 = robinGreenStrongRow(c,flat,scalar,0);
plus = tangentDirection(flat,first,tangentStep);
minus = tangentDirection(flat,first,-tangentStep);
strongPlus = robinGreenStrongRow(c,plus,scalar,tangentStep);
strongMinus = robinGreenStrongRow(c,minus,scalar,-tangentStep);
strong1 = (strongPlus-strongMinus)/(2*tangentStep);
strong0 = tangent'*strong0;
strong1 = tangent'*strong1;
flatDefect = matrixRelativeDifference(weak0,strong0);
tangentDefect = matrixRelativeDifference(weak1,strong1);
diagnostics = struct("flatDefect",flatDefect, ...
    "tangentDefect",tangentDefect, ...
    "maximumDefect",max(flatDefect,tangentDefect));
end

function direction = tangentDirection(flat,first,scale)
direction = flat;
fields = ["E","J","Q","QProjected","R","L"];
for field = fields
    if isfield(flat,field) && isfield(first,field)
        direction.(field) = flat.(field)+scale*first.(field);
    end
end
end

function row = robinGreenStrongRow(c,direction,scalar,terrainScale)
H3 = repmat(c.H,c.nZ,1);
gamma = 1-terrainScale*H3;
volume = -c.wvt.rho0*(scalar.values' ...
    *(c.volumeWeight.*gamma.*direction.Q));
gammaBottom = 1-terrainScale*c.H;
uBottom = (c.uBottom*c.N)./gammaBottom;
vBottom = (c.vBottom*c.N)./gammaBottom;
etaBottom = c.phase*direction.B;
cBottom = c.wvt.f*etaBottom ...
    +terrainScale*(vBottom.*c.hX-uBottom.*c.hY);
boundary = c.wvt.rho0*(scalar.bottomValues'*cBottom)/c.nXY;
row = volume+boundary;
end

function value = matrixRelativeDifference(first,second)
value = norm(first-second,"fro") ...
    /max(norm(first,"fro")+norm(second,"fro"),realmin);
end

function columns = bottomInversionColumns(c,catalog,bounds)
columns = zeros(size(c.N,2),0);
for iHorizontal = 1:c.nK
    if abs(c.horizontalLayout.kMode(iHorizontal)) > bounds(1) ...
            || abs(c.horizontalLayout.lMode(iHorizontal)) > bounds(2)
        continue
    end
    k = c.horizontalLayout.k(iHorizontal);
    l = c.horizontalLayout.l(iHorizontal);
    kappa = hypot(k,l);
    if kappa == 0
        continue
    end
    family = catalog.family{find(abs(catalog.kappa-kappa) ...
        <=100*eps*max(kappa,1),1)};
    state = struct("u",-1i*l*family.bottom.psi, ...
        "v",1i*k*family.bottom.psi, ...
        "w",zeros(size(family.bottom.psi)), ...
        "eta",family.bottom.eta);
    columns(:,end+1) = embedState(c,iHorizontal,state); %#ok<AGROW>
end
end

function columns = stationaryMeanColumns(c,catalog,geostrophicCount)
iHorizontal = find(c.horizontalLayout.kMode == 0 ...
    & c.horizontalLayout.lMode == 0,1);
if isempty(iHorizontal)
    columns = zeros(size(c.N,2),0);
    return
end
family = catalog.family{1};
constant = ones(size(c.quadrature.xi));
etaProfiles = family.robin.eta(:,1:min(geostrophicCount, ...
    size(family.robin.eta,2)));
columns = zeros(size(c.N,2),0);
for iProfile = 1:size(etaProfiles,2)
    if norm(etaProfiles(:,iProfile)) <=100*eps
        continue
    end
    state = struct("u",zeros(size(constant)), ...
        "v",zeros(size(constant)),"w",zeros(size(constant)), ...
        "eta",etaProfiles(:,iProfile));
    columns(:,end+1) = embedState(c,iHorizontal,state); %#ok<AGROW>
end
chi = -c.quadrature.xi/c.wvt.Lz;
state = struct("u",zeros(size(chi)),"v",zeros(size(chi)), ...
    "w",zeros(size(chi)),"eta",chi);
columns(:,end+1) = embedState(c,iHorizontal,state);
end

function value = subspaceConjugacyDefect(X,H,conjugateIndex)
n = size(H,1);
C = sparse((1:n)',conjugateIndex,1,n,n);
P = X*(X'*H);
value = norm(P*C-C*conj(P),"fro")/max(2*norm(P,"fro"),realmin);
end

function result = selectDressedPrimaryResult(grid)
result = grid(end,end);
eligible = false(size(grid));
for iValue = 1:numel(grid)
    eligible(iValue) = grid(iValue).maximumProjectorDefect <=1e-8 ...
        && grid(iValue).maximumFrequencyDefect <=1e-8 ...
        && grid(iValue).maximumRitzResidual <=1e-10 ...
        && grid(iValue).maximumAPVDefect <=1e-8 ...
        && grid(iValue).maximumBottomDefect <=1e-10 ...
        && grid(iValue).maximumStrongResidual <=1e-5 ...
        && grid(iValue).stationaryRepresentationDefect <=1e-8 ...
        && grid(iValue).stationaryRowDefect <=1e-10 ...
        && grid(iValue).stationaryBottomTangencyDefect <=1e-10;
end
indices = find(eligible);
if isempty(indices)
    return
end
compression = arrayfun(@(value)value.compressionFactor,grid(indices));
projector = arrayfun(@(value)value.maximumProjectorDefect,grid(indices));
[~,order] = sortrows([-compression(:) projector(:) indices(:)],[1 2 3]);
result = grid(indices(order(1)));
end

function summary = summarizeDressedConvergence(primary,degrees,padding, ...
    support,waveCounts,geostrophicCounts)
grid = primary.grid;
summary = struct;
summary.waveCounts = waveCounts;
summary.geostrophicCounts = geostrophicCounts;
summary.projectorGrid = arrayfun(@(value)value.maximumProjectorDefect,grid);
summary.stationaryGrid = arrayfun( ...
    @(value)value.stationaryRepresentationDefect,grid);
summary.compressionGrid = arrayfun(@(value)value.compressionFactor,grid);
summary.waveProjectorDefect = min(summary.projectorGrid,[],2);
summary.geostrophicProjectorDefect = min(summary.projectorGrid,[],1).';
summary.geostrophicStationaryDefect = min(summary.stationaryGrid,[],1).';
summary.referenceDegreeProjectorDrift = caseDriftDressed(degrees);
summary.paddingProjectorDrift = caseDriftDressed(padding);
summary.supportProjectorDrift = caseDriftDressed(support);
end

function values = caseDriftDressed(cases)
values = arrayfun(@(value)value.primary.maximumProjectorDefect,cases).';
end

function [classification,diagnosis,gates] = classifyDressedRobin( ...
    training,primary,convergence,tolerance)
result = primary.primary;
selected = training.results(training.selectedIndex,3);
ordinary = training.results(find(isinf(training.ratios),1),1);
undressed = training.results(training.selectedIndex,2);
gates = struct;
gates.training = training.hasEligibleCandidate;
gates.provider = selected.providerPasses;
gates.centered = selected.maximumCenteredDerivativeDefect ...
    <=tolerance.centered;
gates.green = selected.maximumGreenIdentityDefect <= tolerance.green;
gates.structure = result.energyHermitianDefect <= tolerance.structure ...
    && result.exchangeSkewHermitianDefect <= tolerance.structure;
gates.stationary = result.stationaryRepresentationDefect ...
    <=tolerance.stationaryRepresentation ...
    && result.stationaryRowDefect <= tolerance.stationary ...
    && result.stationaryBottomTangencyDefect <= tolerance.stationary;
gates.physical = dressedResultPasses(result,tolerance);
gates.waveConvergence = finalTwoPass( ...
    convergence.waveProjectorDefect,tolerance.convergence);
gates.geostrophicConvergence = finalTwoPass( ...
    convergence.geostrophicProjectorDefect,tolerance.convergence) ...
    && finalTwoPass(convergence.geostrophicStationaryDefect, ...
    tolerance.stationaryRepresentation);
gates.reference = finalTwoPass( ...
    convergence.referenceDegreeProjectorDrift,tolerance.convergence);
gates.padding = finalTwoPass( ...
    convergence.paddingProjectorDrift,tolerance.convergence);
gates.support = finalTwoPass( ...
    convergence.supportProjectorDrift,tolerance.convergence);
gates.economy = result.compressionFactor >= tolerance.compression;
bestBaseline = min(ordinary.stationaryRepresentationDefect, ...
    undressed.stationaryRepresentationDefect);
gates.materialImprovement = bestBaseline ...
    /max(selected.stationaryRepresentationDefect,realmin) >= 4;
physicalPass = gates.training && gates.provider && gates.centered ...
    && gates.green && gates.structure && gates.stationary ...
    && gates.physical && gates.waveConvergence ...
    && gates.geostrophicConvergence && gates.reference ...
    && gates.padding && gates.support;
if physicalPass && gates.economy
    if isinf(training.selectedRobinLengthRatio)
        classification = "dressed-flat-acceleration";
        diagnosis = "Global terrain dressing makes the ordinary-boundary geostrophic family an economical stationary ambient space; finite Robin localization is unnecessary for this oracle.";
    else
        classification = "dressed-robin-acceleration";
        diagnosis = "A trained-and-frozen finite Robin length combined with global terrain dressing reproduces the stationary and internal-wave projectors economically.";
    end
elseif physicalPass
    classification = "dressed-stationary-equivalent";
    diagnosis = "The trained-and-frozen dressed stationary coordinates reproduce the physical projectors but do not retain factor-two compression.";
elseif gates.materialImprovement && gates.training && gates.provider ...
        && gates.centered && gates.green
    classification = "dressed-stationary-seed";
    diagnosis = "Global dressing materially improves the stationary representation but does not yet pass every finite-amplitude physical and economy gate.";
else
    classification = "dressed-stationary-blocker";
    diagnosis = "Neither Robin localization nor the global first terrain correction gives a convergent economical stationary representation in the declared modal ambient space.";
end
end

function passes = dressedResultPasses(result,tolerance)
passes = result.maximumProjectorDefect <= tolerance.projector ...
    && result.maximumFrequencyDefect <= tolerance.frequency ...
    && result.maximumRitzResidual <= tolerance.ritz ...
    && result.maximumAPVDefect <= tolerance.apv ...
    && result.maximumBottomDefect <= tolerance.bottom ...
    && result.maximumStrongResidual <= tolerance.strong ...
    && result.stationaryRepresentationDefect ...
    <=tolerance.stationaryRepresentation ...
    && result.stationaryRowDefect <= tolerance.stationary ...
    && result.stationaryBottomTangencyDefect <= tolerance.stationary;
end

function layouts = supportLayouts(problem,bounds)
source = problem.horizontalLayout;
layouts = cell(size(bounds,1),1);
for iBound = 1:size(bounds,1)
    keep = abs(source.kMode) <= bounds(iBound,1) ...
        & abs(source.lMode) <= bounds(iBound,2);
    layouts{iBound} = source(keep,:);
end
end

function value = dressedCaseConfiguration(role,degree,layout,padding, ...
    scales,quadrature,waveCounts,geostrophicCounts,providerOrders, ...
    robinRatios,trainingCount,tangentStep)
value = struct("role",role,"degree",degree,"kMode",layout.kMode, ...
    "lMode",layout.lMode,"padding",padding,"scales",scales, ...
    "quadrature",quadrature,"waveCounts",waveCounts, ...
    "geostrophicCounts",geostrophicCounts, ...
    "providerOrders",providerOrders,"robinRatios",robinRatios, ...
    "trainingCount",trainingCount,"tangentStep",tangentStep);
end

function value = emptyStationaryTrainingResult
value = struct("ablation","","robinLengthRatio",NaN, ...
    "numberOfRawCoordinates",0,"numberOfCoordinates",0, ...
    "numberOfStationaryCoordinates",0,"fullRank",false, ...
    "gramConditionNumber",Inf,"stationaryRepresentationDefect",Inf, ...
    "stationaryRowDefect",Inf,"bottomTangencyDefect",Inf, ...
    "dressedTangentDefect",Inf,"numberOfTangentCoordinates",0, ...
    "conjugacyDefect",Inf,"maximumCenteredDerivativeDefect",Inf, ...
    "maximumEmbeddingDefect",Inf,"maximumGreenIdentityDefect",Inf, ...
    "greenIdentity",struct,"geostrophicDiagnostics",struct, ...
    "provider",struct,"providerPasses",false, ...
    "rankDiagnostics",struct);
end

function value = emptyDressedComparisonCase
value = struct("degree",0,"paddingFactor",0, ...
    "horizontalModes",zeros(0,2),"numberOfPrimitiveCoordinates",0, ...
    "provider",struct,"reference",struct,"flatSeed",struct, ...
    "grid",repmat(emptyModalResult,0,0),"primary",emptyModalResult, ...
    "amplitude",repmat(emptyModalResult,0,1), ...
    "flatControl",emptyModalResult,"contextSummary",struct, ...
    "robinLengthRatio",NaN);
end

function result = modalCase(problem,layout,trusted,targetWaveIndices, ...
    waveCounts,geostrophicCounts,degree,padding,scales,providerOrders, ...
    quadratureOrder)
[primitive,c] = buildGlobalSmallTerrainPrimitiveAudit(problem,degree, ...
    quadratureOrder,1e-3,horizontalLayout=layout, ...
    paddingFactor=padding,trustedModeBounds=trusted, ...
    rejectTerrainNyquist=true,evaluationScales=scales, ...
    shouldAuditTangent=false);
catalog = verticalCatalog(problem,c,max(waveCounts), ...
    max(geostrophicCounts),providerOrders);
seed = waveColumns(c,catalog,targetWaveIndices,trusted);
flatSeed = subspaceRitz(seed,primitive.evaluatedDirections(1));
reference = continueReference(c,primitive,seed,scales,trusted,degree);

pairedCount = min(numel(waveCounts),numel(geostrophicCounts));
paired = repmat(emptyModalResult,pairedCount,1);
for iCount = 1:pairedCount
    paired(iCount) = reducedModalResult(c, ...
        primitive.evaluatedDirections(end),reference,catalog, ...
        waveCounts(iCount),geostrophicCounts(iCount),trusted,1);
end
waveSweep = repmat(emptyModalResult,numel(waveCounts),1);
for iWave = 1:numel(waveCounts)
    waveSweep(iWave) = reducedModalResult(c, ...
        primitive.evaluatedDirections(end),reference,catalog, ...
        waveCounts(iWave),max(geostrophicCounts),trusted,1);
end
geostrophicSweep = repmat(emptyModalResult,numel(geostrophicCounts),1);
for iGeostrophic = 1:numel(geostrophicCounts)
    geostrophicSweep(iGeostrophic) = reducedModalResult(c, ...
        primitive.evaluatedDirections(end),reference,catalog, ...
        max(waveCounts),geostrophicCounts(iGeostrophic),trusted,1);
end
flatReference = struct("basis",reference.history(1).basis, ...
    "frequency",reference.history(1).frequency);
flatControl = reducedModalResult(c,primitive.evaluatedDirections(1), ...
    flatReference,catalog,max(waveCounts),max(geostrophicCounts), ...
    trusted,0);

result = emptyComparisonCase;
result.degree = degree;
result.paddingFactor = padding;
result.horizontalModes = [layout.kMode layout.lMode];
result.numberOfPrimitiveCoordinates = size(primitive.evaluatedDirections(end).E,1);
result.provider = catalog.diagnostics;
result.reference = reference;
result.flatSeed = flatSeed;
result.paired = paired;
result.waveSweep = waveSweep;
result.geostrophicSweep = geostrophicSweep;
result.primary = selectPrimaryModalResult(paired);
result.flatControl = flatControl;
result.contextSummary = contextSummary(c);
end

function result = selectPrimaryModalResult(values)
result = values(end);
for iValue = 1:numel(values)
    value = values(iValue);
    if value.maximumProjectorDefect <=1e-8 ...
            && value.maximumFrequencyDefect <=1e-8 ...
            && value.maximumRitzResidual <=1e-10 ...
            && value.maximumAPVDefect <=1e-8 ...
            && value.maximumBottomDefect <=1e-10 ...
            && value.maximumStrongResidual <=1e-5 ...
            && value.stationaryRowDefect <=1e-10 ...
            && value.stationaryBottomTangencyDefect <=1e-10
        result = value;
        return
    end
end
end

function catalog = verticalCatalog(problem,c,waveCount,geostrophicCount,orders,robinLength)
if nargin < 6
    robinLength = Inf;
end
kappa = hypot(c.horizontalLayout.k,c.horizontalLayout.l);
uniqueKappa = unique(kappa(kappa > 0),"sorted");
count = max(waveCount,geostrophicCount+1);
high = cell(numel(uniqueKappa),1);
low = cell(numel(uniqueKappa),1);
difference = zeros(numel(uniqueKappa),1);
endpoint = zeros(numel(uniqueKappa),1);
for iKappa = 1:numel(uniqueKappa)
    high{iKappa} = buildBoundaryCompleteVerticalReferenceFamily( ...
        problem,uniqueKappa(iKappa),count,robinLength,orders(end), ...
        c.quadrature.xi);
    low{iKappa} = buildBoundaryCompleteVerticalReferenceFamily( ...
        problem,uniqueKappa(iKappa),count,robinLength,orders(end-1), ...
        c.quadrature.xi);
    difference(iKappa) = familyDifference(low{iKappa},high{iKappa}, ...
        c.quadrature.weight);
    diagnostics = high{iKappa}.diagnostics;
    endpoint(iKappa) = max([diagnostics.waveEndpointDefect, ...
        diagnostics.surfaceRobinDefect, ...
        diagnostics.bottomRobinDefect, ...
        diagnostics.bottomAPVDefect, ...
        diagnostics.bottomValueDefect, ...
        diagnostics.bottomSurfaceValueDefect]);
end
catalog = struct("kappa",uniqueKappa,"family",{high}, ...
    "waveCount",waveCount,"geostrophicCount",geostrophicCount, ...
    "diagnostics",struct("providerDifference",max(difference), ...
    "maximumEndpointDefect",max(endpoint), ...
    "negativeRobinModeCountDefect",max(cellfun(@(value) ...
    value.diagnostics.negativeRobinModeCountDefect,high)), ...
    "orders",orders,"robinLength",robinLength));
end

function result = reducedModalResult(c,direction,reference,catalog, ...
    waveCount,geostrophicCount,trusted,terrainScale)
waves = waveColumns(c,catalog,1:waveCount,[Inf Inf]);
[geostrophic,geostrophicDiagnostics] = geostrophicColumns( ...
    c,catalog,geostrophicCount,[Inf Inf]);
mean = meanSectorColumns(c,catalog,waveCount,geostrophicCount);
raw = [waves geostrophic mean];
[X,rankDiagnostics] = energyIndependentColumns(raw,direction.E);
H = X'*direction.E*X;
J = X'*direction.J*X;

exactStationary = constructCompleteStationarySpace(c,direction, ...
    trusted,c.nG,min(4,c.nG),terrainScale);
[G,stationaryDiagnostics] = stationarySubspace( ...
    exactStationary.trustedBasis,direction,X,H,J);
dense = reducedPhysicalEigensystem(X,H,J,G,direction.E);
selected = selectByOverlap(dense.vectors,reference.basis, ...
    direction.E,size(reference.basis,2));
selected = closePhysicalBlock(selected,dense,direction.E, ...
    c.coordinateConjugateIndex);
basis = energyNormalizeColumns(dense.vectors(:,selected),direction.E);
frequency = columnRayleighQuotients( ...
    basis,direction.E,1i*direction.J);
physical = physicalDiagnostics(c,direction,basis,frequency);

result = emptyModalResult;
result.waveModeCount = waveCount;
result.geostrophicModeCount = geostrophicCount;
result.numberOfModalCoordinates = size(X,2);
result.numberOfPrimitiveCoordinates = size(direction.E,1);
result.compressionFactor = size(direction.E,1)/size(X,2);
result.numberOfStationaryCoordinates = size(G,2);
result.numberOfDynamicalCoordinates = size(X,2)-size(G,2);
result.numberOfSelectedCoordinates = size(basis,2);
result.basis = basis;
result.frequency = frequency;
result.maximumProjectorDefect = maximumPrincipalSine( ...
    basis,reference.basis,direction.E);
result.maximumFrequencyDefect = sortedFrequencyDefect( ...
    frequency,reference.frequency);
result.maximumRitzResidual = physical.maximumRitzResidual;
result.maximumAPVDefect = physical.maximumAPVDefect;
result.maximumBottomDefect = physical.maximumBottomDefect;
result.maximumStrongResidual = physical.maximumStrongResidual;
result.maximumPhysicalDefect = max([physical.maximumRitzResidual, ...
    physical.maximumAPVDefect,physical.maximumBottomDefect, ...
    physical.maximumStrongResidual]);
result.energyHermitianDefect = norm(H-H',"fro") ...
    /max(norm(H,"fro"),realmin);
result.exchangeSkewHermitianDefect = norm(J+J',"fro") ...
    /max(norm(J,"fro"),realmin);
result.stationaryRowDefect = stationaryDiagnostics.rowDefect;
result.stationaryBottomTangencyDefect = ...
    stationaryDiagnostics.bottomTangencyDefect;
result.stationaryRepresentationDefect = ...
    stationaryDiagnostics.representationDefect;
result.maximumEmbeddingDefect = max( ...
    [geostrophicDiagnostics.maximumEmbeddingDefect, ...
    rankDiagnostics.maximumRejectedRelativeSingularValue]);
result.rankDiagnostics = rankDiagnostics;
result.geostrophicDiagnostics = geostrophicDiagnostics;
result.trustedGeostrophicDiagnostics = ...
    exactStationary.trustedRepresentation;
result.trustedModeBounds = trusted;
result.terrainScale = terrainScale;
end

function [columns,diagnostics] = geostrophicColumns(c,catalog,count,bounds)
columns = zeros(size(c.N,2),0);
embeddingDefect = zeros(0,1);
for iHorizontal = 1:c.nK
    if abs(c.horizontalLayout.kMode(iHorizontal)) > bounds(1) ...
            || abs(c.horizontalLayout.lMode(iHorizontal)) > bounds(2)
        continue
    end
    k = c.horizontalLayout.k(iHorizontal);
    l = c.horizontalLayout.l(iHorizontal);
    kappa = hypot(k,l);
    if kappa == 0
        continue
    end
    family = catalog.family{find(abs(catalog.kappa-kappa) ...
        <=100*eps*max(kappa,1),1)};
    for iMode = 1:count
        F = family.robin.F(:,iMode);
        state = struct("u",-1i*l*F,"v",1i*k*F, ...
            "w",zeros(size(F)),"eta",family.robin.eta(:,iMode));
        [coordinate,defect] = embedState(c,iHorizontal,state);
        columns(:,end+1) = coordinate; %#ok<AGROW>
        embeddingDefect(end+1,1) = defect; %#ok<AGROW>
    end
    state = struct("u",-1i*l*family.bottom.psi, ...
        "v",1i*k*family.bottom.psi, ...
        "w",zeros(size(family.bottom.psi)), ...
        "eta",family.bottom.eta);
    [coordinate,defect] = embedState(c,iHorizontal,state);
    columns(:,end+1) = coordinate; %#ok<AGROW>
    embeddingDefect(end+1,1) = defect; %#ok<AGROW>
end
diagnostics = struct("maximumEmbeddingDefect", ...
    max(embeddingDefect,[],'all'), ...
    "numberOfRawColumns",size(columns,2));
end

function columns = waveColumns(c,catalog,modeIndices,bounds)
columns = zeros(size(c.N,2),0);
for iHorizontal = 1:c.nK
    if abs(c.horizontalLayout.kMode(iHorizontal)) > bounds(1) ...
            || abs(c.horizontalLayout.lMode(iHorizontal)) > bounds(2)
        continue
    end
    k = c.horizontalLayout.k(iHorizontal);
    l = c.horizontalLayout.l(iHorizontal);
    kappa = hypot(k,l);
    if kappa == 0
        continue
    end
    family = catalog.family{find(abs(catalog.kappa-kappa) ...
        <=100*eps*max(kappa,1),1)};
    for iMode = reshape(modeIndices,1,[])
        for sigma = [-1 1]
            omega = family.wave.omega(iMode);
            h = family.wave.h(iMode);
            F = family.wave.F(:,iMode);
            G = family.wave.G(:,iMode);
            state = struct( ...
                "u",(k*omega-sigma*1i*c.wvt.f*l) ...
                /(omega*kappa)*F, ...
                "v",(l*omega+sigma*1i*c.wvt.f*k) ...
                /(omega*kappa)*F, ...
                "w",-1i*kappa*h*G, ...
                "eta",-sigma*kappa*h/omega*G);
            columns(:,end+1) = embedState(c,iHorizontal,state); %#ok<AGROW>
        end
    end
end
end

function columns = meanSectorColumns(c,catalog,waveCount,geostrophicCount)
iHorizontal = find(c.horizontalLayout.kMode == 0 ...
    & c.horizontalLayout.lMode == 0,1);
if isempty(iHorizontal)
    columns = zeros(size(c.N,2),0);
    return
end
family = catalog.family{1};
constant = ones(size(c.quadrature.xi));
profiles = [constant family.robin.F(:,1:waveCount)];
columns = zeros(size(c.N,2),0);
for iProfile = 1:size(profiles,2)
    F = profiles(:,iProfile);
    for sigma = [-1 1]
        state = struct("u",F,"v",sigma*1i*F, ...
            "w",zeros(size(F)),"eta",zeros(size(F)));
        columns(:,end+1) = embedState(c,iHorizontal,state); %#ok<AGROW>
    end
end
etaProfiles = family.robin.eta(:,1:min(geostrophicCount, ...
    size(family.robin.eta,2)));
for iProfile = 1:size(etaProfiles,2)
    if norm(etaProfiles(:,iProfile)) <=100*eps
        continue
    end
    state = struct("u",zeros(size(constant)), ...
        "v",zeros(size(constant)),"w",zeros(size(constant)), ...
        "eta",etaProfiles(:,iProfile));
    columns(:,end+1) = embedState(c,iHorizontal,state); %#ok<AGROW>
end
chi = -c.quadrature.xi/c.wvt.Lz;
state = struct("u",zeros(size(chi)),"v",zeros(size(chi)), ...
    "w",zeros(size(chi)),"eta",chi);
columns(:,end+1) = embedState(c,iHorizontal,state);
end

function reference = continueReference(c,primitive,seed,scales,trusted,degree)
previous = energyOrthonormalize(seed,primitive.evaluatedDirections(1).E);
history = repmat(struct("basis",zeros(0),"frequency",zeros(0,1)), ...
    numel(scales),1);
for iScale = 1:numel(scales)
    direction = primitive.evaluatedDirections(iScale);
    if iScale == 1 && scales(iScale) == 0
        flat = subspaceRitz(previous,direction);
        previous = flat.basis;
        history(iScale).basis = flat.basis;
        history(iScale).frequency = flat.frequency;
        continue
    end
    stationary = constructCompleteStationarySpace(c,direction, ...
        trusted,degree,min(4,degree),scales(iScale));
    previous = projectOut(previous,stationary.trustedBasis,direction.E);
    previous = energyOrthonormalize(previous,direction.E);
    dense = physicalEigensystem(direction,stationary.trustedBasis);
    selected = selectByOverlap(dense.vectors,previous,direction.E, ...
        size(previous,2));
    selected = closePhysicalBlock(selected,dense,direction.E, ...
        c.coordinateConjugateIndex);
    previous = energyNormalizeColumns( ...
        dense.vectors(:,selected),direction.E);
    history(iScale).basis = previous;
    history(iScale).frequency = columnRayleighQuotients( ...
        previous,direction.E,1i*direction.J);
end
reference = history(end);
reference.history = history;
reference.numberOfCoordinates = size(reference.basis,2);
end

function result = subspaceRitz(basis,direction)
H = direction.E;
basis = energyOrthonormalize(basis,H);
K = basis'*(1i*direction.J)*basis;
[vectors,form] = schur((K+K')/2,"complex");
frequency = real(diag(form));
[frequency,order] = sort(frequency);
basis = basis*vectors(:,order);
R = chol((H+H')/2);
residual = 1i*direction.J*basis-(H*basis).*frequency.';
result = struct("basis",basis,"frequency",frequency, ...
    "maximumResidual",max(vecnorm(R'\residual,2,1),[],'all'));
end

function [G,diagnostics] = stationarySubspace(exact,direction,X,H,J)
if isempty(exact)
    G = zeros(size(X,1),0);
    diagnostics = struct("rowDefect",0,"bottomTangencyDefect",0);
    return
end
coordinate = H\(X'*direction.E*exact);
G = X*coordinate;
representationResidual = exact-G;
representationDefect = norm(chol(direction.E)*representationResidual,"fro") ...
    /max(norm(chol(direction.E)*exact,"fro"),realmin);
G = energyOrthonormalize(G,direction.E);
coordinate = H\(X'*direction.E*G);
rowDefect = norm(J*coordinate,"fro") ...
    /max(norm(J,"fro")*norm(coordinate,"fro"),realmin);
bottom = direction.R*G;
bottomDefect = norm(bottom,"fro") ...
    /max(norm(direction.R,"fro")*norm(G,"fro"),realmin);
diagnostics = struct("rowDefect",rowDefect, ...
    "bottomTangencyDefect",bottomDefect, ...
    "representationDefect",representationDefect, ...
    "constraintRank",size(G,2), ...
    "constraintSingularValues",zeros(0,1));
end

function dense = reducedPhysicalEigensystem(X,H,J,G,E)
Gcoordinate = H\(X'*E*G);
R = chol((H+H')/2);
Yg = R*Gcoordinate;
if isempty(Yg)
    complement = eye(size(H));
else
    [complete,~] = qr(Yg);
    complement = complete(:,size(Yg,2)+1:end);
end
K = R'\(1i*J/R);
reduced = complement'*K*complement;
[vectors,form] = schur((reduced+reduced')/2,"complex");
frequency = real(diag(form));
[frequency,order] = sort(frequency);
coordinates = R\(complement*vectors(:,order));
residual = 1i*J*coordinates-(H*coordinates).*frequency.';
uncertainty = vecnorm(R'\residual,2,1).';
vectors = energyNormalizeColumns(X*coordinates,E);
uncertainty = max(uncertainty,10*eps(max(norm(K,2),realmin)));
dense = struct("frequency",frequency,"vectors",vectors, ...
    "absoluteFrequencyUncertainty",uncertainty);
end

function [coordinate,defect] = embedState(c,iHorizontal,state)
weight = c.quadrature.weight;
uCoefficient = weightedFit(c.spaces.F,state.u,weight);
vCoefficient = weightedFit(c.spaces.F,state.v,weight);
wCoefficient = weightedFit(c.spaces.G,state.w,weight);
etaCoefficient = weightedFit(c.spaces.H,state.eta,weight);
rawBlock = [uCoefficient;vCoefficient;wCoefficient;etaCoefficient];
rawRows = (iHorizontal-1)*c.nXBlock+(1:c.nXBlock);
admissible = c.layout.admissibleRanges{iHorizontal};
localBasis = c.N(rawRows,admissible);
coordinate = zeros(size(c.N,2),1);
coordinate(admissible) = localBasis\rawBlock;
represented = localBasis*coordinate(admissible);
defect = norm(represented-rawBlock) ...
    /max(norm(rawBlock),realmin);
end

function coefficient = weightedFit(reconstruction,values,weight)
coefficient = (reconstruction'*(weight.*reconstruction)) ...
    \(reconstruction'*(weight.*values));
end

function [X,diagnostics] = energyIndependentColumns(X,H)
R = chol((H+H')/2);
[~,singularMatrix,right] = svd(R*X,"econ");
singularValues = diag(singularMatrix);
tolerance = max(size(X))*eps(max(singularValues,[],'all'))*100;
rankValue = nnz(singularValues > tolerance);
X = X*right(:,1:rankValue) ...
    /singularMatrix(1:rankValue,1:rankValue);
rejected = singularValues(rankValue+1:end);
diagnostics = struct("rank",rankValue, ...
    "nominalDimension",numel(singularValues), ...
    "singularValues",singularValues,"tolerance",tolerance, ...
    "maximumRejectedRelativeSingularValue", ...
    max([rejected;0])/max(singularValues(1),realmin));
end

function diagnostics = physicalDiagnostics(c,direction,basis,frequency)
H = direction.E;
A = 1i*direction.J;
R = chol((H+H')/2);
residual = A*basis-(H*basis).*frequency.';
ritz = vecnorm(R'\residual,2,1).' ...
    ./max(norm(R'\(A/R),2)+abs(frequency),realmin);
trustedRows = find(repmat(c.layout.trustedHorizontalModes,c.nZ,1));
weightedQ = sqrt(c.apvVerticalWeight(trustedRows)) ...
    .*direction.QProjected(trustedRows,:);
apv = vecnorm(weightedQ*basis,2,1).' ...
    /max(norm(weightedQ/R,2),realmin);
bottomResidual = -1i*(direction.B*basis).*frequency.' ...
    -direction.R*basis;
bottomScale = abs(frequency)*norm(direction.B/R,2) ...
    +norm(direction.R/R,2);
bottom = vecnorm(bottomResidual,2,1).'./max(bottomScale,realmin);
strong = strongModeDiagnostics(c,direction,basis,frequency);
diagnostics = struct("maximumRitzResidual",max(ritz,[],'all'), ...
    "maximumAPVDefect",max(apv,[],'all'), ...
    "maximumBottomDefect",max(bottom,[],'all'), ...
    "maximumStrongResidual",strong.maximumResidual, ...
    "ritzResidual",ritz,"apvDefect",apv, ...
    "bottomDefect",bottom,"strong",strong);
end

function diagnostics = strongModeDiagnostics(c,direction,C,frequency)
rawState = c.N*C;
pressure = direction.pressureMap*C;
stateTendency = -1i*rawState.*frequency.';
residual = direction.matrices.S*[stateTendency;pressure] ...
    -direction.matrices.F*C;
terms = {direction.matrices.S*[stateTendency;pressure], ...
    -direction.matrices.F*C};
defect = columnRelativeResidual(residual,terms);
diagnostics = struct("maximumByMode",defect, ...
    "maximumResidual",max(defect,[],'all'));
end

function value = columnRelativeResidual(residual,terms)
scale = zeros(size(residual,2),1);
for iTerm = 1:numel(terms)
    scale = scale+vecnorm(terms{iTerm},2,1).';
end
value = vecnorm(residual,2,1).'./max(scale,realmin);
end

function summary = summarizeConvergence(primary,degrees,padding,support, ...
    waveCounts,geostrophicCounts)
summary = struct;
summary.waveCounts = waveCounts;
summary.waveProjectorDefect = arrayfun( ...
    @(value)value.maximumProjectorDefect,primary.waveSweep);
summary.geostrophicCounts = geostrophicCounts;
summary.geostrophicProjectorDefect = arrayfun( ...
    @(value)value.maximumProjectorDefect,primary.geostrophicSweep);
summary.pairedProjectorDefect = arrayfun( ...
    @(value)value.maximumProjectorDefect,primary.paired);
summary.pairedCompressionFactor = arrayfun( ...
    @(value)value.compressionFactor,primary.paired);
summary.referenceDegreeProjectorDrift = caseDrift(degrees);
summary.paddingProjectorDrift = caseDrift(padding);
summary.supportProjectorDrift = caseDrift(support);
end

function values = caseDrift(cases)
% The case bases generally have different primitive coordinate layouts.
% Their physically relevant convergence is reported by each case's direct
% modal-versus-primitive projector defect rather than by subtracting raw
% coordinate projectors.
values = arrayfun(@(value)value.primary.maximumProjectorDefect,cases).';
end

function [classification,diagnosis,gates] = classify(primary,convergence,tolerance)
result = primary.primary;
provider = primary.provider;
gates = struct;
gates.provider = provider.providerDifference <= tolerance.provider ...
    && provider.maximumEndpointDefect <= tolerance.provider;
gates.flat = primary.flatControl.maximumProjectorDefect ...
    <=tolerance.projector ...
    && primary.flatControl.maximumRitzResidual <= tolerance.ritz ...
    && primary.flatControl.maximumAPVDefect <= tolerance.apv ...
    && primary.flatControl.maximumBottomDefect <= tolerance.bottom ...
    && primary.flatControl.maximumStrongResidual <= tolerance.strong;
gates.structure = result.energyHermitianDefect <= tolerance.structure ...
    && result.exchangeSkewHermitianDefect <= tolerance.structure;
gates.stationary = result.stationaryRowDefect <= tolerance.stationary ...
    && result.stationaryBottomTangencyDefect <= tolerance.stationary ...
    && result.stationaryRepresentationDefect <= tolerance.stationary;
gates.physical = result.maximumRitzResidual <= tolerance.ritz ...
    && result.maximumProjectorDefect <= tolerance.projector ...
    && result.maximumFrequencyDefect <= tolerance.frequency ...
    && result.maximumAPVDefect <= tolerance.apv ...
    && result.maximumBottomDefect <= tolerance.bottom ...
    && result.maximumStrongResidual <= tolerance.strong;
gates.waveConvergence = finalTwoPass( ...
    convergence.waveProjectorDefect,tolerance.projector);
gates.geostrophicConvergence = finalTwoPass( ...
    convergence.geostrophicProjectorDefect,tolerance.projector);
gates.reference = finalTwoPass( ...
    convergence.referenceDegreeProjectorDrift,tolerance.reference);
gates.padding = finalTwoPass( ...
    convergence.paddingProjectorDrift,tolerance.padding);
gates.support = finalTwoPass( ...
    convergence.supportProjectorDrift,tolerance.support);
gates.economy = result.compressionFactor >= tolerance.compression;
if ~gates.provider
    classification = "provider-blocker";
    diagnosis = "The pinned InternalModesEVP wave or ordinary geostrophic family fails its independent order or endpoint qualification.";
elseif ~gates.flat
    classification = "wave-modal-blocker";
    diagnosis = "The fixed-wavenumber wave coordinates do not recover the independent flat primitive wave block.";
elseif ~gates.geostrophicConvergence || ~gates.stationary
    classification = "geostrophic-modal-blocker";
    diagnosis = "The ordinary flat geostrophic modes plus explicit bottom inversion do not converge to a stationary-compatible finite-terrain ambient space.";
elseif ~gates.waveConvergence
    classification = "wave-modal-blocker";
    diagnosis = "The fixed-wavenumber internal-wave guard sequence does not converge to the primitive internal-wave projector.";
elseif ~gates.reference
    classification = "primitive-reference-blocker";
    diagnosis = "The independent polynomial validation projector has not converged under the declared reference-degree sequence.";
elseif ~gates.structure || ~gates.physical || ~gates.padding || ~gates.support
    classification = "wave-modal-blocker";
    diagnosis = "The modal ambient retains physical energy but does not pass every internal-wave projector, APV, bottom, strong, padding, or support gate.";
elseif gates.economy
    classification = "wave-vortex-modal-acceleration";
    diagnosis = "The flat wave-vortex modal ambient reproduces the primitive internal-wave projector with at least factor-two compression.";
else
    classification = "wave-vortex-modal-equivalent";
    diagnosis = "The flat wave-vortex modal ambient reproduces the primitive internal-wave projector but does not reduce the retained dimension by a factor of two.";
end
end

function value = familyDifference(first,second,weight)
fields = {first.wave.F,second.wave.F;first.wave.G,second.wave.G; ...
    first.robin.F,second.robin.F; ...
    first.robin.eta,second.robin.eta; ...
    first.bottom.psi,second.bottom.psi; ...
    first.bottom.eta,second.bottom.eta};
value = 0;
for iField = 1:size(fields,1)
    value = max(value,columnSubspaceDifference( ...
        fields{iField,1},fields{iField,2},weight));
end
end

function value = columnSubspaceDifference(first,second,weight)
if isvector(first)
    first = first(:);
    second = second(:);
end
first = weightedOrthonormalize(first,weight);
second = weightedOrthonormalize(second,weight);
value = norm(first-second*(second'*(weight.*first)),2);
end

function values = weightedOrthonormalize(values,weight)
[left,singularMatrix,~] = svd(sqrt(weight).*values,"econ");
singularValues = diag(singularMatrix);
tolerance = max(size(values))*eps(max(singularValues,[],'all'))*100;
rankValue = nnz(singularValues > tolerance);
values = left(:,1:rankValue)./sqrt(weight);
end

function dense = physicalEigensystem(direction,G)
H = direction.E;
R = chol((H+H')/2);
Yg = R*G;
if isempty(Yg)
    complement = eye(size(H));
else
    [complete,~] = qr(Yg);
    complement = complete(:,size(Yg,2)+1:end);
end
K = R'\(1i*direction.J/R);
reduced = complement'*K*complement;
[vectors,form] = schur((reduced+reduced')/2,"complex");
frequency = real(diag(form));
[frequency,order] = sort(frequency);
vectors = R\(complement*vectors(:,order));
vectors = energyNormalizeColumns(vectors,H);
residual = 1i*direction.J*vectors-(H*vectors).*frequency.';
uncertainty = vecnorm(R'\residual,2,1).';
uncertainty = max(uncertainty,10*eps(max(norm(K,2),realmin)));
dense = struct("frequency",frequency,"vectors",vectors, ...
    "absoluteFrequencyUncertainty",uncertainty);
end

function indices = selectByOverlap(vectors,target,H,count)
projector = energyProjector(target,H);
participation = real(sum(conj(vectors).*(H*projector*vectors),1)).';
[~,indices] = maxk(participation,min(count,numel(participation)));
indices = sort(indices);
end

function indices = closePhysicalBlock(indices,dense,H,conjugateIndex)
% The declared wave projector has fixed dimension. Backward-error-nearby
% directions remain in the unresolved complement and are reported through
% the isolation audit; they are not silently appended to the physical block.
indices = closeConjugacy(indices,dense.vectors,H,conjugateIndex);
end

function indices = closeConjugacy(indices,vectors,H,conjugateIndex)
n = size(H,1);
conjugacy = sparse((1:n)',conjugateIndex,1,n,n);
for index = reshape(indices,1,[])
    partner = conjugacy*conj(vectors(:,index));
    overlap = abs(vectors'*(H*partner));
    [~,iPartner] = max(overlap);
    indices(end+1,1) = iPartner; %#ok<AGROW>
end
indices = sort(unique(indices));
end

function X = projectOut(X,G,H)
if ~isempty(G)
    X = X-G*(G'*H*X);
end
end

function X = energyOrthonormalize(X,H)
if isempty(X)
    return
end
R = chol((H+H')/2);
[~,singularMatrix,right] = svd(R*X,"econ");
singularValues = diag(singularMatrix);
tolerance = max(size(X))*eps(max(singularValues,[],'all'))*100;
rankValue = nnz(singularValues > tolerance);
X = X*right(:,1:rankValue) ...
    /singularMatrix(1:rankValue,1:rankValue);
end

function X = energyNormalizeColumns(X,H)
energy = real(sum(conj(X).*(H*X),1));
X = X./sqrt(energy);
end

function value = energyProjector(basis,H)
value = basis*((basis'*H*basis)\(basis'*H));
end

function value = maximumPrincipalSine(first,second,H)
if size(first,2) ~= size(second,2)
    value = Inf;
    return
end
R = chol((H+H')/2);
first = orth(R*first);
second = orth(R*second);
singularValues = svd(first'*second);
value = sqrt(max(0,1-min(singularValues,[],'all')^2));
end

function value = sortedFrequencyDefect(first,second)
first = sort(first(:));
second = sort(second(:));
if numel(first) ~= numel(second)
    value = Inf;
else
    value = max(abs(first-second)) ...
        /max([abs(first);abs(second);realmin]);
end
end

function values = columnRayleighQuotients(C,H,A)
values = real(sum(conj(C).*(A*C),1).' ...
    ./sum(conj(C).*(H*C),1).');
end

function passes = finalTwoPass(values,tolerance)
values = values(:);
passes = numel(values) >= 2 && all(values(end-1:end) <= tolerance);
end

function summary = contextSummary(c)
summary = struct("numberOfHorizontalModes",c.nK, ...
    "numberOfPrimitiveCoordinates",size(c.N,2), ...
    "numberOfQuadraturePoints",c.nZ, ...
    "horizontalModes",[c.horizontalLayout.kMode ...
    c.horizontalLayout.lMode]);
end

function terrain = terrainSupport(problem)
wvt = problem.originatingTransform;
spectrum = fft2(problem.topographicHeight)/(wvt.Nx*wvt.Ny);
nyquist = logical(WVGeometryDoublyPeriodic.maskForNyquistModes( ...
    wvt.Nx,wvt.Ny));
scale = max(abs(spectrum),[],'all');
tolerance = 100*eps*max(scale,1);
[kMode,lMode] = ndgrid(wvt.kMode_dft,wvt.lMode_dft);
active = abs(spectrum) > tolerance & ~nyquist ...
    & ~(kMode == 0 & lMode == 0);
terrain = struct("kMode",kMode(active),"lMode",lMode(active), ...
    "coefficient",spectrum(active));
end

function shells = scatteringShells(problem,terrain,orders,trusted)
source = problem.horizontalLayout;
trustedMask = abs(source.kMode) <= trusted(1) ...
    & abs(source.lMode) <= trusted(2);
seed = [source.kMode(trustedMask) source.lMode(trustedMask)];
seed = seed(any(seed ~= 0,2),:);
terrainModes = unique([terrain.kMode terrain.lMode],"rows");
reachable = cell(max(orders)+1,1);
reachable{1} = seed;
for iOrder = 1:max(orders)
    sums = zeros(0,2);
    for iTerrain = 1:size(terrainModes,1)
        sums = [sums;reachable{iOrder}+terrainModes(iTerrain,:)]; %#ok<AGROW>
    end
    reachable{iOrder+1} = unique([reachable{iOrder};sums],"rows");
end
layouts = cell(numel(reachable),1);
for iOrder = 1:numel(reachable)
    requested = unique([[0 0];seed;reachable{iOrder}],"rows");
    [present,index] = ismember(requested, ...
        [source.kMode source.lMode],"rows");
    if ~all(present)
        error("WVTerrainEnergyGalerkin:UnavailableModalAmbientSupport", ...
            "The originating transform does not retain every requested modal scattering sideband.")
    end
    keep = false(height(source),1);
    keep(index) = true;
    layouts{iOrder} = source(keep,:);
end
shells = struct("orders",(0:max(orders)).', ...
    "terrainModes",terrainModes,"reachableModes",{reachable}, ...
    "layouts",{layouts});
end

function value = emptyModalResult
value = struct("waveModeCount",0,"geostrophicModeCount",0, ...
    "numberOfModalCoordinates",0,"numberOfPrimitiveCoordinates",0, ...
    "compressionFactor",0,"numberOfStationaryCoordinates",0, ...
    "numberOfDynamicalCoordinates",0, ...
    "numberOfSelectedCoordinates",0,"basis",zeros(0), ...
    "frequency",zeros(0,1),"maximumProjectorDefect",Inf, ...
    "maximumFrequencyDefect",Inf,"maximumRitzResidual",Inf, ...
    "maximumAPVDefect",Inf,"maximumBottomDefect",Inf, ...
    "maximumStrongResidual",Inf,"maximumPhysicalDefect",Inf, ...
    "energyHermitianDefect",Inf, ...
    "exchangeSkewHermitianDefect",Inf,"stationaryRowDefect",Inf, ...
    "stationaryBottomTangencyDefect",Inf, ...
    "stationaryRepresentationDefect",Inf, ...
    "maximumEmbeddingDefect",Inf,"rankDiagnostics",struct, ...
    "geostrophicDiagnostics",struct, ...
    "trustedGeostrophicDiagnostics",struct, ...
    "trustedModeBounds",[0 0],"terrainScale",0);
end

function value = emptyComparisonCase
value = struct("degree",0,"paddingFactor",0, ...
    "horizontalModes",zeros(0,2),"numberOfPrimitiveCoordinates",0, ...
    "provider",struct,"reference",struct,"flatSeed",struct, ...
    "paired", ...
    repmat(emptyModalResult,0,1),"waveSweep", ...
    repmat(emptyModalResult,0,1),"geostrophicSweep", ...
    repmat(emptyModalResult,0,1),"primary",emptyModalResult, ...
    "flatControl",emptyModalResult,"contextSummary",struct);
end

function value = caseConfiguration(role,degree,layout,padding,scales, ...
    quadrature,waveCount,geostrophicCount,providerOrders)
value = struct("role",role,"degree",degree, ...
    "kMode",layout.kMode,"lMode",layout.lMode, ...
    "padding",padding,"scales",scales,"quadrature",quadrature, ...
    "waveCount",waveCount,"geostrophicCount",geostrophicCount, ...
    "providerOrders",providerOrders);
end

function cache = emptyCacheDiagnostics(directory)
cache = struct("enabled",directory ~= "","directory",directory, ...
    "records",repmat(emptyCacheRecord,0,1));
end

function record = emptyCacheRecord
record = struct("kind","","key","","path","","status","", ...
    "bytes",0,"buildSeconds",0,"secondsSaved",0);
end

function cache = appendCacheRecord(cache,record)
cache.records(end+1,1) = record;
end

function cache = finalizeCacheDiagnostics(cache)
status = string({cache.records.status});
cache.numberOfHits = nnz(status == "hit");
cache.numberOfMisses = nnz(status == "miss" | status == "disabled");
cache.numberOfInvalidations = nnz(status == "invalidated");
cache.bytes = sum([cache.records.bytes]);
cache.secondsSaved = sum([cache.records.secondsSaved]);
end

function [value,record] = cachedAuditValue(problem,directory,kind,configuration,builder)
record = emptyCacheRecord;
record.kind = kind;
start = tic;
if directory == ""
    value = builder();
    record.status = "disabled";
    record.buildSeconds = toc(start);
    return
end
if ~isfolder(directory)
    mkdir(directory);
end
schema = "wave-vortex-modal-ambient-cache-v2";
key = configurationKey(problem,kind,configuration,schema);
path = fullfile(directory,kind+"-"+key+".mat");
record.key = key;
record.path = string(path);
if isfile(path)
    try
        loaded = load(path,"entry");
        entry = loaded.entry;
        if entry.schema == schema && entry.key == key
            value = restoreCachedValue(entry.value,problem);
            information = dir(path);
            record.status = "hit";
            record.bytes = information.bytes;
            record.buildSeconds = entry.buildSeconds;
            record.secondsSaved = entry.buildSeconds;
            return
        end
        record.status = "invalidated";
    catch
        record.status = "invalidated";
    end
end
buildStart = tic;
value = builder();
buildSeconds = toc(buildStart);
entry = struct("schema",schema,"key",key, ...
    "configuration",configuration,"buildSeconds",buildSeconds, ...
    "created",string(datetime("now","TimeZone","UTC")), ...
    "value",stripCachedValue(value));
temporary = string(tempname(directory))+".mat";
cleanup = onCleanup(@() deleteIfPresent(temporary));
save(temporary,"entry","-v7.3");
movefile(temporary,path,"f");
clear cleanup
information = dir(path);
if record.status ~= "invalidated"
    record.status = "miss";
end
record.bytes = information.bytes;
record.buildSeconds = buildSeconds;
end

function value = stripCachedValue(value)
if isfield(value,"reference") && isfield(value.reference,"history")
    for iHistory = 1:numel(value.reference.history)
        if isfield(value.reference.history(iHistory),"direction")
            value.reference.history(iHistory) = rmfield( ...
                value.reference.history(iHistory),"direction");
        end
    end
end
end

function value = restoreCachedValue(value,~)
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end

function key = configurationKey(problem,kind,configuration,schema)
engine = java.security.MessageDigest.getInstance("SHA-256");
updateDigest(engine,uint8(char(schema)));
updateDigest(engine,uint8(char(kind)));
updateDigest(engine,uint8(char(jsonencode(configuration))));
privateRoot = string(fileparts(mfilename("fullpath")));
classRoot = string(fileparts(privateRoot));
files = [string(mfilename("fullpath"))+".m"; ...
    fullfile(classRoot,"auditWaveVortexModalAmbient.m"); ...
    fullfile(classRoot,"auditTerrainDressedRobinStationaryAmbient.m"); ...
    fullfile(privateRoot, ...
    "buildBoundaryCompleteVerticalReferenceFamily.m")];
for iFile = 1:numel(files)
    updateDigest(engine,uint8(fileread(files(iFile))));
end
wvt = problem.originatingTransform;
constants = [wvt.Lx wvt.Ly wvt.Lz wvt.Nx wvt.Ny wvt.Nz ...
    wvt.f wvt.rho0 wvt.g double(wvt.shouldAntialias)];
updateDigest(engine,typecast(double(constants(:)),"uint8"));
updateDigest(engine,typecast(double(problem.topographicHeight(:)),"uint8"));
z = linspace(-wvt.Lz,0,257).';
N2 = wvt.N2Function(z);
if isscalar(N2)
    N2 = repmat(N2,size(z));
end
updateDigest(engine,typecast(double(N2(:)),"uint8"));
digest = typecast(int8(engine.digest()),"uint8");
key = string(lower(reshape(dec2hex(digest,2).',1,[])));
end

function updateDigest(engine,bytes)
engine.update(typecast(uint8(bytes(:)),"int8"));
end
