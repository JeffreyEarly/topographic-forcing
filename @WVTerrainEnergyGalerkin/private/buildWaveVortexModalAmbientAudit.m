function audit = buildWaveVortexModalAmbientAudit(problem,trustedBounds, ...
    targetWaveIndices,waveCounts,geostrophicCounts,orders, ...
    referenceDegrees,paddingFactors,terrainScales,providerOrders, ...
    shouldRunTwoDimensionalControl,cacheDirectory,quadratureOrder)
% Build the Milestone-10.2.3 flat wave-vortex modal ambient oracle.

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

function catalog = verticalCatalog(problem,c,waveCount,geostrophicCount,orders)
kappa = hypot(c.horizontalLayout.k,c.horizontalLayout.l);
uniqueKappa = unique(kappa(kappa > 0),"sorted");
count = max(waveCount,geostrophicCount+1);
high = cell(numel(uniqueKappa),1);
low = cell(numel(uniqueKappa),1);
difference = zeros(numel(uniqueKappa),1);
endpoint = zeros(numel(uniqueKappa),1);
for iKappa = 1:numel(uniqueKappa)
    high{iKappa} = buildBoundaryCompleteVerticalReferenceFamily( ...
        problem,uniqueKappa(iKappa),count,Inf,orders(end), ...
        c.quadrature.xi);
    low{iKappa} = buildBoundaryCompleteVerticalReferenceFamily( ...
        problem,uniqueKappa(iKappa),count,Inf,orders(end-1), ...
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
    "orders",orders,"robinLength",Inf));
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
schema = "wave-vortex-modal-ambient-cache-v1";
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
