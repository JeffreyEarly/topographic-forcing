function audit = buildConvergedPhysicalSubspaceClassificationAudit( ...
    dense,flat,guard)
% Build the Milestone-9.1 converged physical-subspace oracle.

nRefinement = numel(dense.primitivePolynomialDegrees);
nScale = numel(dense.terrainScales);
iFullScale = nScale;

fullScale = cell(nRefinement,1);
for iRefinement = 1:nRefinement
    fullScale{iRefinement} = dense.details{iRefinement,1,iFullScale};
end

transfers = cell(nRefinement-1,1);
for iRefinement = 1:nRefinement-1
    transfers{iRefinement} = nestedTransfer( ...
        fullScale{iRefinement},fullScale{iRefinement+1});
end

guardResults = squeeze(guard.details(:,1,end));
guardTransfer = nestedTransfer(guardResults{1},guardResults{2});
stationary = stationaryClassification(fullScale,transfers,dense.details, ...
    guardResults,guardTransfer);
internal = internalClassification(dense,flat,transfers);
boundary = boundaryClassification(dense,transfers);
internal = applyGuardAndPaddingDiagnostics( ...
    internal,guardResults,guardTransfer,dense.details);
boundary = applyBoundaryGuardAndPaddingDiagnostics( ...
    boundary,guardResults,guardTransfer,dense.details);
classification = finestClassification(fullScale{end},stationary, ...
    internal,boundary);

tolerance = struct( ...
    "transfer",1e-12, ...
    "stationaryPrincipalSine",1e-8, ...
    "stationaryPadding",1e-8, ...
    "projector",1e-10, ...
    "verticalPrincipalSine",1e-3, ...
    "principalSineReduction",4, ...
    "frequency",1e-8, ...
    "frequencyPath",1e-3, ...
    "homotopy",1e-1, ...
    "apv",1e-8, ...
    "bottom",1e-8, ...
    "strong",1e-5, ...
    "strongReduction",4, ...
    "padding",1e-8, ...
    "minimumFrequencyResolutionRatio",1e3, ...
    "bottomParticipationChange",1e-2, ...
    "terrainFrequencyFit",1e-2);

transferPasses = max(cellfun( ...
    @(value)value.maximumDefect,transfers),[],"all") ...
    <= tolerance.transfer;
stationaryPasses = stationary.isConverged ...
    && stationary.principalSine(end) ...
    <= tolerance.stationaryPrincipalSine ...
    && stationary.guardDefect <= tolerance.stationaryPrincipalSine ...
    && stationary.paddingDefect <= tolerance.stationaryPadding;
dynamicPasses = any([internal.tracks.isValidated]) ...
    || boundary.isValidated;
classificationPasses = classification.maximumProjectorDefect ...
    <= tolerance.projector ...
    && classification.conjugacyDefect <= tolerance.projector ...
    && classification.dimensionDefect == 0;

if ~transferPasses || ~classificationPasses
    status = "implementation-unresolved";
    diagnosis = "The nested transfer or physical-energy projector decomposition failed before scientific subspace classification.";
elseif ~stationaryPasses
    status = "stationary-subspace-convergence-blocker";
    diagnosis = "The fixed trusted stationary balanced space did not converge as one physical-energy subspace under nested support and vertical refinement.";
elseif ~dynamicPasses
    status = "dynamical-subspace-convergence-blocker";
    diagnosis = "The stationary space converged, but no dynamical spectral subspace passed the independent frequency, APV, bottom, strong-equation, and projector gates.";
else
    status = "physical-subspace-classification-oracle";
    if boundary.isValidated
        diagnosis = "The stationary, internal-wave, and topographic-boundary-wave projectors converge while the unresolved algebraic completion remains explicitly retained.";
    else
        diagnosis = "The stationary and internal-wave projectors converge; the topographic-boundary candidate remains in the explicitly retained unresolved algebraic completion.";
    end
end

audit = struct;
audit.scope = "milestone-9.1-converged-physical-subspaces";
audit.status = status;
audit.isCompatible = status == "physical-subspace-classification-oracle";
audit.diagnosis = diagnosis;
audit.trustedModeBounds = dense.trustedModeBounds;
audit.supportModeBounds = dense.supportModeBounds;
audit.stationaryPolynomialDegree = dense.stationaryPolynomialDegree;
audit.primitivePolynomialDegrees = dense.primitivePolynomialDegrees;
audit.paddingFactors = dense.paddingFactors;
audit.terrainScales = dense.terrainScales;
audit.quadratureOrder = dense.quadratureOrder;
audit.dense = dense;
audit.flat = flat;
audit.guard = guard;
audit.transfers = transfers;
audit.stationary = stationary;
audit.internal = internal;
audit.boundary = boundary;
audit.classification = classification;
audit.finest = fullScale{end};
audit.requiredTolerance = tolerance;
audit.nextScope = "milestone-10-only-if-compatible-and-separately-authorized";
end

function diagnostics = nestedTransfer(coarse,fine)
coarseRepresentation = coarse.nestedRepresentation;
fineRepresentation = fine.nestedRepresentation;
rawMap = rawNestedMap(coarseRepresentation,fineRepresentation);
embeddedRaw = rawMap*coarseRepresentation.admissibleBasis;
transfer = fineRepresentation.admissibleBasis\embeddedRaw;
representedRaw = fineRepresentation.admissibleBasis*transfer;
rawDefect = norm(representedRaw-embeddedRaw,"fro") ...
    /max(norm(embeddedRaw,"fro"),realmin);

horizontalMap = horizontalNestedMap( ...
    coarseRepresentation.horizontalLayout, ...
    fineRepresentation.horizontalLayout);
bottomResidual = fineRepresentation.bottomMap*transfer ...
    -horizontalMap*coarseRepresentation.bottomMap;
bottomDefect = norm(bottomResidual,"fro") ...
    /max(norm(horizontalMap*coarseRepresentation.bottomMap,"fro"),realmin);

coarseConjugacy = sparse( ...
    (1:numel(coarse.coordinateConjugateIndex))', ...
    coarse.coordinateConjugateIndex,1, ...
    numel(coarse.coordinateConjugateIndex), ...
    numel(coarse.coordinateConjugateIndex));
fineConjugacy = sparse( ...
    (1:numel(fine.coordinateConjugateIndex))', ...
    fine.coordinateConjugateIndex,1, ...
    numel(fine.coordinateConjugateIndex), ...
    numel(fine.coordinateConjugateIndex));
conjugacyResidual = transfer*coarseConjugacy ...
    -fineConjugacy*conj(transfer);
conjugacyDefect = norm(conjugacyResidual,"fro") ...
    /max(2*norm(transfer,"fro"),realmin);

diagnostics = struct( ...
    "matrix",transfer, ...
    "rawMap",rawMap, ...
    "horizontalMap",horizontalMap, ...
    "rawReconstructionDefect",rawDefect, ...
    "bottomValueDefect",bottomDefect, ...
    "conjugacyDefect",conjugacyDefect, ...
    "maximumDefect",max([rawDefect bottomDefect conjugacyDefect]));
end

function rawMap = rawNestedMap(coarse,fine)
nCoarse = coarse.numberOfRawStateCoefficients;
nFine = fine.numberOfRawStateCoefficients;
rawMap = sparse(nFine,nCoarse);
for iCoarse = 1:height(coarse.horizontalLayout)
    mode = [coarse.horizontalLayout.kMode(iCoarse) ...
        coarse.horizontalLayout.lMode(iCoarse)];
    iFine = find(fine.horizontalLayout.kMode == mode(1) ...
        & fine.horizontalLayout.lMode == mode(2),1);
    if isempty(iFine)
        error("WVTerrainEnergyGalerkin:NonNestedPhysicalSubspaceSupport", ...
            "Every coarse signed Fourier mode must be retained by the finer physical-subspace audit.")
    end
    coarseBlock = (iCoarse-1)*coarse.rawStateBlockSize;
    fineBlock = (iFine-1)*fine.rawStateBlockSize;
    rawMap = insertIdentity(rawMap, ...
        fineBlock+(1:coarse.numberOfFCoordinates), ...
        coarseBlock+(1:coarse.numberOfFCoordinates));
    rawMap = insertIdentity(rawMap, ...
        fineBlock+fine.numberOfFCoordinates ...
        +(1:coarse.numberOfFCoordinates), ...
        coarseBlock+coarse.numberOfFCoordinates ...
        +(1:coarse.numberOfFCoordinates));
    rawMap = insertIdentity(rawMap, ...
        fineBlock+2*fine.numberOfFCoordinates ...
        +(1:coarse.numberOfGCoordinates), ...
        coarseBlock+2*coarse.numberOfFCoordinates ...
        +(1:coarse.numberOfGCoordinates));
    fineEta = fineBlock+2*fine.numberOfFCoordinates ...
        +fine.numberOfGCoordinates;
    coarseEta = coarseBlock+2*coarse.numberOfFCoordinates ...
        +coarse.numberOfGCoordinates;
    rawMap = insertIdentity(rawMap, ...
        fineEta+(1:coarse.numberOfGCoordinates), ...
        coarseEta+(1:coarse.numberOfGCoordinates));
    rawMap(fineEta+fine.numberOfHCoordinates, ...
        coarseEta+coarse.numberOfHCoordinates) = 1;
end
end

function values = insertIdentity(values,rows,columns)
values(rows,columns) = speye(numel(rows),numel(columns));
end

function values = horizontalNestedMap(coarse,fine)
row = zeros(height(coarse),1);
for iCoarse = 1:height(coarse)
    row(iCoarse) = find(fine.kMode == coarse.kMode(iCoarse) ...
        & fine.lMode == coarse.lMode(iCoarse),1);
end
values = sparse(row,(1:height(coarse))',1,height(fine),height(coarse));
end

function result = stationaryClassification( ...
    results,transfers,allDetails,guardResults,guardTransfer)
nRefinement = numel(results);
principalSine = zeros(nRefinement-1,1);
dimension = zeros(nRefinement,1);
for iRefinement = 1:nRefinement
    dimension(iRefinement) = ...
        size(results{iRefinement}.stationaryBasis,2);
end
for iRefinement = 1:nRefinement-1
    principalSine(iRefinement) = maximumPrincipalSine( ...
        transfers{iRefinement}.matrix ...
        *results{iRefinement}.stationaryBasis, ...
        results{iRefinement+1}.stationaryBasis, ...
        results{iRefinement+1}.energyMatrix);
end

finest = results{end};
[apvBearing,zeroAPV,split] = splitStationaryByAPV(finest);
paddingDefect = 0;
reference = allDetails{end,1,end};
for iPadding = 2:size(allDetails,2)
    current = allDetails{end,iPadding,end};
    paddingDefect = max(paddingDefect,maximumPrincipalSine( ...
        current.stationaryBasis,reference.stationaryBasis, ...
        reference.energyMatrix));
end
guardDefect = maximumPrincipalSine( ...
    guardTransfer.matrix*guardResults{1}.stationaryBasis, ...
    guardResults{2}.stationaryBasis, ...
    guardResults{2}.energyMatrix);

result = struct;
result.basis = finest.stationaryBasis;
result.apvBearingBasis = apvBearing;
result.zeroVolumeAPVBasis = zeroAPV;
result.apvSplit = split;
result.dimension = dimension;
result.principalSine = principalSine;
result.maximumPrincipalSine = max(principalSine,[],"all");
result.paddingDefect = paddingDefect;
result.guardDefect = guardDefect;
result.isConverged = all(dimension == dimension(1)) ...
    && split.rankStable ...
    && all(diff(principalSine) <= 0);
end

function [apvBearing,zeroAPV,diagnostics] = splitStationaryByAPV(result)
representation = result.nestedRepresentation;
Q = representation.projectedAPVMap;
weight = representation.projectedAPVWeight;
G = result.stationaryBasis;
weightedAPV = sqrt(weight).*Q*G;
[~,singularMatrix,V] = svd(weightedAPV,"econ");
singularValues = diag(singularMatrix);
if isempty(singularValues)
    nominalTolerance = 0;
else
    nominalTolerance = max(size(weightedAPV)) ...
        *eps(max(singularValues))*100;
end
rankValues = arrayfun(@(factor) ...
    nnz(singularValues > factor*nominalTolerance),[0.1 1 10]);
rankValue = rankValues(2);
apvBearing = G*V(:,1:rankValue);
zeroAPV = G*V(:,rankValue+1:end);
diagnostics = struct( ...
    "singularValues",singularValues, ...
    "balancedEnstrophyEigenvalues",singularValues.^2, ...
    "rankTolerance",nominalTolerance, ...
    "rankValues",rankValues, ...
    "rankStable",all(rankValues == rankValue), ...
    "apvBearingDimension",size(apvBearing,2), ...
    "zeroVolumeAPVDimension",size(zeroAPV,2), ...
    "zeroVolumeAPVDefect",norm(weightedAPV*V(:,rankValue+1:end),"fro") ...
    /max(norm(weightedAPV,"fro"),realmin));
end

function result = internalClassification(dense,flat,transfers)
nRefinement = numel(dense.primitivePolynomialDegrees);
nScale = numel(dense.terrainScales);
candidates = cell(nRefinement,1);
for iRefinement = 1:nRefinement
    seedBlocks = spectralBlocks(flat.details{iRefinement,1,1});
    homotopyDefect = zeros(numel(seedBlocks),nScale);
    active = seedBlocks;
    for iScale = 1:nScale
        terrainResult = dense.details{iRefinement,1,iScale};
        terrainBlocks = spectralBlocks(terrainResult);
        [active,defect] = matchBlocks(active,terrainBlocks, ...
            terrainResult.energyMatrix,[]);
        homotopyDefect(:,iScale) = defect;
    end
    for iBlock = 1:numel(active)
        active(iBlock).homotopyDefect = ...
            max(homotopyDefect(iBlock,:),[],"all");
    end
    candidates{iRefinement} = active;
end

tracks = refinementTracks(candidates,transfers, ...
    cellfun(@(value)value.energyMatrix, ...
    squeeze(dense.details(:,1,end)),"UniformOutput",false));
for iTrack = 1:numel(tracks)
    tracks(iTrack) = evaluateInternalTrack(tracks(iTrack), ...
        squeeze(dense.details(:,1,end)),dense.inertialFrequency);
end
result = struct("candidates",{candidates},"tracks",tracks, ...
    "numberOfValidatedTracks",nnz([tracks.isValidated]));
end

function result = boundaryClassification(dense,transfers)
nRefinement = numel(dense.primitivePolynomialDegrees);
nScale = numel(dense.terrainScales);
candidates = cell(nRefinement,1);
for iRefinement = 1:nRefinement
    first = dense.details{iRefinement,1,1};
    blocks = spectralBlocks(first);
    seed = first.stationary.trustedBottomSeedBasis;
    capture = arrayfun(@(block)subspaceCapture(seed,block.basis, ...
        first.energyMatrix),blocks);
    [maximumCapture,index] = max(capture);
    candidate = blocks(index);
    homotopyDefect = 0;
    for iScale = 2:nScale
        current = dense.details{iRefinement,1,iScale};
        currentBlocks = spectralBlocks(current);
        [matched,defect] = matchBlocks(candidate,currentBlocks, ...
            current.energyMatrix,[]);
        candidate = matched(1);
        homotopyDefect = max(homotopyDefect,defect(1));
    end
    candidate.seedCapture = maximumCapture;
    candidate.homotopyDefect = homotopyDefect;
    candidates{iRefinement} = candidate;
end

energy = cellfun(@(value)value.energyMatrix, ...
    squeeze(dense.details(:,1,end)),"UniformOutput",false);
tracks = refinementTracks(candidates,transfers,energy);
track = tracks(1);
track = evaluateBoundaryTrack(track,squeeze(dense.details(:,1,end)), ...
    dense.terrainScales,squeeze(dense.details(end,1,:)));
result = struct("candidates",{candidates},"track",track, ...
    "isValidated",track.isValidated);
end

function result = applyGuardAndPaddingDiagnostics( ...
    result,guardResults,guardTransfer,allDetails)
[guardMatched,guardDefect] = matchBlocks( ...
    spectralBlocks(guardResults{1}), ...
    spectralBlocks(guardResults{2}), ...
    guardResults{2}.energyMatrix,guardTransfer.matrix);
fineGuardBlocks = spectralBlocks(guardResults{2});
for iTrack = 1:numel(result.tracks)
    [fineMatch,fineMatchDefect] = matchBlocks( ...
        result.tracks(iTrack).blocks{end},fineGuardBlocks, ...
        guardResults{2}.energyMatrix,[]);
    result.tracks(iTrack).guardDefect = Inf;
    if isfinite(fineMatchDefect) && ~isempty(fineMatch.indices)
        for iGuard = 1:numel(guardMatched)
            if sameIndices(guardMatched(iGuard),fineMatch)
                result.tracks(iTrack).guardDefect = guardDefect(iGuard);
                break
            end
        end
    end
    result.tracks(iTrack).paddingDefect = paddingBlockDefect( ...
        result.tracks(iTrack).blocks{end},allDetails);
    result.tracks(iTrack).isValidated = ...
        result.tracks(iTrack).isValidated ...
        &&result.tracks(iTrack).guardDefect <= 1e-8 ...
        &&result.tracks(iTrack).paddingDefect <= 1e-8;
end
result.numberOfValidatedTracks = nnz([result.tracks.isValidated]);
end

function result = applyBoundaryGuardAndPaddingDiagnostics( ...
    result,guardResults,guardTransfer,allDetails)
container = struct("tracks",result.track);
container = applyGuardAndPaddingDiagnostics( ...
    container,guardResults,guardTransfer,allDetails);
result.track = container.tracks;
result.isValidated = result.track.isValidated;
end

function value = paddingBlockDefect(referenceBlock,allDetails)
value = 0;
for iPadding = 2:size(allDetails,2)
    current = allDetails{end,iPadding,end};
    [~,defect] = matchBlocks(referenceBlock, ...
        spectralBlocks(current),current.energyMatrix,[]);
    value = max(value,defect);
end
end

function tf = sameIndices(first,second)
tf = isequal(sort(first.indices),sort(second.indices));
end

function blocks = spectralBlocks(result)
frequency = result.frequency;
uncertainty = result.absoluteFrequencyUncertainty;
positive = find(frequency > uncertainty);
used = false(size(positive));
blocks = repmat(emptyBlock,0,1);
for iPositive = 1:numel(positive)
    if used(iPositive)
        continue
    end
    group = iPositive;
    changed = true;
    while changed
        changed = false;
        for jPositive = 1:numel(positive)
            if used(jPositive) || ismember(jPositive,group)
                continue
            end
            first = positive(group);
            overlaps = abs(frequency(positive(jPositive))-frequency(first)) ...
                <= 10*(uncertainty(positive(jPositive))+uncertainty(first));
            if any(overlaps)
                group(end+1) = jPositive; %#ok<AGROW>
                changed = true;
            end
        end
    end
    used(group) = true;
    positiveIndices = positive(group);
    negative = find(frequency < -uncertainty);
    partnerMask = false(size(negative));
    for iMember = 1:numel(positiveIndices)
        partnerMask = partnerMask ...
            |abs(frequency(negative)+frequency(positiveIndices(iMember))) ...
            <=10*(uncertainty(negative) ...
            +uncertainty(positiveIndices(iMember)));
    end
    partnerIndices = negative(partnerMask);
    indices = sort(unique([positiveIndices;partnerIndices]));
    block = emptyBlock;
    block.indices = indices;
    block.basis = result.modeVectors(:,indices);
    block.dimension = numel(indices);
    block.positiveFrequency = sort(frequency(positiveIndices));
    block.maximumEigenResidual = max(result.eigenResidual(indices));
    blocks(end+1,1) = block; %#ok<AGROW>
end
end

function block = emptyBlock
block = struct("indices",zeros(0,1),"basis",zeros(0), ...
    "dimension",0,"positiveFrequency",zeros(0,1), ...
    "maximumEigenResidual",NaN,"homotopyDefect",NaN, ...
    "seedCapture",NaN);
end

function [matched,defect] = matchBlocks(source,target,H,transfer)
nSource = numel(source);
matched = repmat(emptyBlock,nSource,1);
defect = Inf(nSource,1);
cost = Inf(nSource,numel(target));
[energyFactor,flag] = chol((H+H')/2);
if flag ~= 0
    error("WVTerrainEnergyGalerkin:PhysicalSubspaceEnergyNotPositive", ...
        "The physical-subspace comparison requires a positive energy matrix.")
end
sourceOrthonormal = cell(nSource,1);
for iSource = 1:nSource
    sourceBasis = source(iSource).basis;
    if isempty(sourceBasis)
        continue
    end
    if ~isempty(transfer)
        if size(transfer,2) ~= size(sourceBasis,1)
            continue
        end
        sourceBasis = transfer*sourceBasis;
    end
    sourceOrthonormal{iSource} = orthonormalizeWithFactor( ...
        sourceBasis,energyFactor);
end
targetOrthonormal = cell(numel(target),1);
for iTarget = 1:numel(target)
    targetBasis = target(iTarget).basis;
    if ~isempty(targetBasis) && size(targetBasis,1) == size(H,1)
        targetOrthonormal{iTarget} = orthonormalizeWithFactor( ...
            targetBasis,energyFactor);
    end
end
for iSource = 1:nSource
    sourceBasis = sourceOrthonormal{iSource};
    if isempty(sourceBasis)
        continue
    end
    for iTarget = 1:numel(target)
        if source(iSource).dimension == target(iTarget).dimension
            targetBasis = targetOrthonormal{iTarget};
            if size(sourceBasis,2) == size(targetBasis,2)
                residual = sourceBasis ...
                    -targetBasis*(targetBasis'*H*sourceBasis);
                cost(iSource,iTarget) = ...
                    norm(energyFactor*residual,2);
            end
        end
    end
end
finiteCost = cost(isfinite(cost));
if isempty(finiteCost)
    return
end
unmatchedCost = max(2*max(finiteCost),1);
assignments = matchpairs(cost,unmatchedCost,"min");
for iAssignment = 1:size(assignments,1)
    iSource = assignments(iAssignment,1);
    iTarget = assignments(iAssignment,2);
    matched(iSource) = target(iTarget);
    defect(iSource) = cost(iSource,iTarget);
end
end

function tracks = refinementTracks(candidates,transfers,energy)
tracks = repmat(emptyTrack,0,1);
coarse = candidates{1};
for iCandidate = 1:numel(coarse)
    track = emptyTrack;
    track.blocks = cell(numel(candidates),1);
    track.blocks{1} = coarse(iCandidate);
    track.principalSine = NaN(numel(candidates)-1,1);
    tracks(end+1,1) = track; %#ok<AGROW>
end
for iRefinement = 1:numel(candidates)-1
    source = repmat(emptyBlock,numel(tracks),1);
    for iTrack = 1:numel(tracks)
        source(iTrack) = tracks(iTrack).blocks{iRefinement};
    end
    [matched,defect] = matchBlocks(source,candidates{iRefinement+1}, ...
        energy{iRefinement+1},transfers{iRefinement}.matrix);
    for iTrack = 1:numel(tracks)
        tracks(iTrack).blocks{iRefinement+1} = matched(iTrack);
        tracks(iTrack).principalSine(iRefinement) = defect(iTrack);
    end
end
end

function track = emptyTrack
track = struct("blocks",{{}},"principalSine",zeros(0,1), ...
    "frequencyDefect",zeros(0,1),"apvDefect",zeros(0,1), ...
    "bottomDefect",zeros(0,1),"strongResidual",zeros(0,1), ...
    "homotopyDefect",zeros(0,1), ...
    "frequencyResolutionRatio",zeros(0,1), ...
    "bottomParticipation",zeros(0,1), ...
    "guardDefect",Inf,"paddingDefect",Inf, ...
    "isValidated",false);
end

function track = evaluateInternalTrack(track,results,inertialFrequency)
n = numel(track.blocks);
track = populatePhysicalDiagnostics(track,results);
track.frequencyDefect = consecutiveFrequencyDefect(track.blocks);
track.homotopyDefect = cellfun( ...
    @(block)block.homotopyDefect,track.blocks);
principalReduction = track.principalSine(1:end-1) ...
    ./max(track.principalSine(2:end),realmin);
strongReduction = track.strongResidual(1:end-1) ...
    ./max(track.strongResidual(2:end),realmin);
track.isValidated = all(isfinite(track.principalSine)) ...
    && track.principalSine(end) <= 1e-3 ...
    && all(principalReduction >= 4) ...
    && track.frequencyDefect(end) <= 1e-8 ...
    && max(track.frequencyDefect) <= 1e-3 ...
    && max(track.homotopyDefect) <= 1e-1 ...
    && min(track.blocks{end}.positiveFrequency) ...
    >=(1-1e-8)*inertialFrequency ...
    && max(track.bottomParticipation) <= 1e-8 ...
    && track.apvDefect(end) <= 1e-8 ...
    && track.bottomDefect(end) <= 1e-8 ...
    && all(strongReduction >= 4) ...
    && track.strongResidual(end) <= 1e-5 ...
    && all(cellfun(@(value)~isempty(value.indices),track.blocks));
if n < 3
    track.isValidated = false;
end
end

function track = evaluateBoundaryTrack(track,results,scales,amplitudeResults)
track = populatePhysicalDiagnostics(track,results);
track.frequencyDefect = consecutiveFrequencyDefect(track.blocks);
track.homotopyDefect = cellfun( ...
    @(block)block.homotopyDefect,track.blocks);
apvDecreases = all(diff(track.apvDefect) < 0);
bottomDecreases = all(diff(track.bottomDefect) < 0);
strongDecreases = all(diff(track.strongResidual) < 0);
participationChange = abs(track.bottomParticipation(end) ...
    -track.bottomParticipation(end-1)) ...
    /max(track.bottomParticipation(end),realmin);

amplitudeFrequency = zeros(numel(scales),1);
for iScale = 1:numel(scales)
    result = amplitudeResults{iScale};
    blocks = spectralBlocks(result);
    seed = result.stationary.trustedBottomSeedBasis;
    capture = arrayfun(@(block)subspaceCapture( ...
        seed,block.basis,result.energyMatrix),blocks);
    [~,index] = max(capture);
    amplitudeFrequency(iScale) = mean(blocks(index).positiveFrequency);
end
design = [scales scales.^2];
coefficient = design\amplitudeFrequency;
fit = design*coefficient;
fitDefect = norm(fit-amplitudeFrequency) ...
    /max(norm(amplitudeFrequency),realmin);
track.terrainFrequency = amplitudeFrequency;
track.terrainFrequencyFit = fit;
track.terrainFrequencyFitDefect = fitDefect;
track.bottomParticipationChange = participationChange;
track.isValidated = all(isfinite(track.principalSine)) ...
    && track.principalSine(end) <= 1e-8 ...
    && max(track.homotopyDefect) <= 1e-1 ...
    && min(track.frequencyResolutionRatio) >= 1e3 ...
    && apvDecreases && bottomDecreases && strongDecreases ...
    && track.bottomParticipation(end) > 100*eps ...
    && participationChange <= 1e-2 ...
    && fitDefect <= 1e-2;
end

function track = populatePhysicalDiagnostics(track,results)
n = numel(track.blocks);
track.apvDefect = Inf(n,1);
track.bottomDefect = Inf(n,1);
track.strongResidual = Inf(n,1);
track.frequencyResolutionRatio = zeros(n,1);
track.bottomParticipation = zeros(n,1);
for iResult = 1:n
    indices = track.blocks{iResult}.indices;
    if isempty(indices)
        continue
    end
    result = results{iResult};
    track.apvDefect(iResult) = ...
        max(result.trustedAPV.directDefect(indices));
    track.bottomDefect(iResult) = max(result.bottomDefect(indices));
    track.strongResidual(iResult) = ...
        max(result.strong.maximumByMode(indices));
    track.frequencyResolutionRatio(iResult) = ...
        min(result.frequencyResolutionRatio(indices));
    track.bottomParticipation(iResult) = ...
        mean(result.bottomSeedParticipation(indices));
end
end

function defect = consecutiveFrequencyDefect(blocks)
defect = NaN(numel(blocks)-1,1);
for iBlock = 1:numel(blocks)-1
    first = blocks{iBlock}.positiveFrequency;
    second = blocks{iBlock+1}.positiveFrequency;
    if numel(first) ~= numel(second) || isempty(first)
        continue
    end
    scale = max([abs(first);abs(second);realmin]);
    defect(iBlock) = max(abs(sort(first)-sort(second)))/scale;
end
end

function tracks = finestClassification(result,stationary,internal,boundary)
H = result.energyMatrix;
G = stationary.basis;
internalIndices = zeros(0,1);
for iTrack = find([internal.tracks.isValidated])
    internalIndices = [internalIndices; ...
        internal.tracks(iTrack).blocks{end}.indices]; %#ok<AGROW>
end
boundaryIndices = zeros(0,1);
if boundary.isValidated
    boundaryIndices = boundary.track.blocks{end}.indices;
end
waveIndices = sort(unique([internalIndices;boundaryIndices]));
allDynamic = (1:size(result.modeVectors,2))';
unresolvedZero = intersect(result.unresolvedZeroIndices,allDynamic);
unresolvedNonzero = setdiff(allDynamic, ...
    [waveIndices;unresolvedZero],"stable");
dynamicLabels = repmat("unresolved-nonzero",numel(allDynamic),1);
dynamicLabels(unresolvedZero) = "unresolved-zero";
dynamicLabels(internalIndices) = "internal-wave";
dynamicLabels(boundaryIndices) = "topographic-boundary-wave";

stationaryProjector = energyProjector(G,H);
waveProjector = energyProjector(result.modeVectors(:,waveIndices),H);
remainderProjector = eye(size(H))-stationaryProjector-waveProjector;
identity = eye(size(H));
projectorDefects = [ ...
    projectorDefect(stationaryProjector,H), ...
    projectorDefect(waveProjector,H), ...
    projectorDefect(remainderProjector,H), ...
    norm(stationaryProjector*waveProjector,"fro"), ...
    norm(stationaryProjector*remainderProjector,"fro"), ...
    norm(waveProjector*remainderProjector,"fro"), ...
    norm(stationaryProjector+waveProjector+remainderProjector ...
    -identity,"fro")/sqrt(size(identity,1))];

n = size(H,1);
conjugacy = sparse((1:n)',result.coordinateConjugateIndex,1,n,n);
projectors = {stationaryProjector,waveProjector,remainderProjector};
conjugacyDefect = 0;
for iProjector = 1:numel(projectors)
    residual = projectors{iProjector}*conjugacy ...
        -conjugacy*conj(projectors{iProjector});
    conjugacyDefect = max(conjugacyDefect,norm(residual,"fro") ...
        /max(2*norm(projectors{iProjector},"fro"),realmin));
end

tracks = struct;
tracks.stationaryProjector = stationaryProjector;
tracks.internalWaveProjector = energyProjector( ...
    result.modeVectors(:,internalIndices),H);
tracks.topographicBoundaryProjector = energyProjector( ...
    result.modeVectors(:,boundaryIndices),H);
tracks.waveProjector = waveProjector;
tracks.remainderProjector = remainderProjector;
tracks.energyForms = struct( ...
    "stationary",H*stationaryProjector, ...
    "internalWave",H*tracks.internalWaveProjector, ...
    "topographicBoundary",H*tracks.topographicBoundaryProjector, ...
    "wave",H*waveProjector, ...
    "remainder",H*remainderProjector);
tracks.stationaryAPVBearingDimension = ...
    size(stationary.apvBearingBasis,2);
tracks.stationaryZeroVolumeAPVDimension = ...
    size(stationary.zeroVolumeAPVBasis,2);
tracks.internalWaveIndices = internalIndices;
tracks.topographicBoundaryIndices = boundaryIndices;
tracks.unresolvedZeroIndices = unresolvedZero;
tracks.unresolvedNonzeroIndices = unresolvedNonzero;
tracks.demonstratedNumericalIndices = zeros(0,1);
tracks.dynamicLabels = dynamicLabels;
tracks.dimension = struct( ...
    "state",size(H,1), ...
    "stationary",size(G,2), ...
    "wave",numel(waveIndices), ...
    "remainder",size(H,1)-size(G,2)-numel(waveIndices));
tracks.dimensionDefect = tracks.dimension.state ...
    -tracks.dimension.stationary-tracks.dimension.wave ...
    -tracks.dimension.remainder;
tracks.projectorDefects = projectorDefects;
tracks.maximumProjectorDefect = max(projectorDefects);
tracks.conjugacyDefect = conjugacyDefect;
end

function value = energyProjector(C,H)
if isempty(C)
    value = zeros(size(H));
else
    gram = (C'*H*C+C'*H'*C)/2;
    value = C*(gram\(C'*H));
end
end

function value = projectorDefect(P,H)
value = max( ...
    norm(P^2-P,"fro")/max(norm(P,"fro"),realmin), ...
    norm(P'*H-H*P,"fro")/max(norm(H*P,"fro"),realmin));
end

function value = maximumPrincipalSine(first,second,H)
if size(first,2) ~= size(second,2)
    value = Inf;
    return
end
if isempty(first)
    value = 0;
    return
end
first = orthonormalize(first,H);
second = orthonormalize(second,H);
if size(first,2) ~= size(second,2)
    value = Inf;
    return
end
[energyFactor,flag] = chol((H+H')/2);
if flag ~= 0
    error("WVTerrainEnergyGalerkin:PhysicalSubspaceEnergyNotPositive", ...
        "The physical-subspace comparison requires a positive energy matrix.")
end
residual = first-second*(second'*H*first);
value = norm(energyFactor*residual,2);
end

function values = orthonormalize(values,H)
[energyFactor,flag] = chol((H+H')/2);
if flag ~= 0
    error("WVTerrainEnergyGalerkin:PhysicalSubspaceEnergyNotPositive", ...
        "The physical-subspace comparison requires a positive energy matrix.")
end
values = orthonormalizeWithFactor(values,energyFactor);
end

function values = orthonormalizeWithFactor(values,energyFactor)
[~,singularMatrix,right] = svd(energyFactor*values,"econ");
singularValues = diag(singularMatrix);
if isempty(singularValues)
    values = zeros(size(values,1),0);
    return
end
tolerance = max(size(values))*eps(max(singularValues))*100;
rankValue = nnz(singularValues > tolerance);
values = values*right(:,1:rankValue) ...
    /singularMatrix(1:rankValue,1:rankValue);
end

function value = subspaceCapture(seed,candidate,H)
if isempty(seed) || isempty(candidate)
    value = 0;
    return
end
seed = orthonormalize(seed,H);
candidate = orthonormalize(candidate,H);
singularValues = svd(seed'*H*candidate);
value = sum(singularValues.^2)/size(seed,2);
end
