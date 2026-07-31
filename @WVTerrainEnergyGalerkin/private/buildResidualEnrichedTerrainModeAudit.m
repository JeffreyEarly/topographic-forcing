function audit = buildResidualEnrichedTerrainModeAudit(problem,trustedBounds,supportBounds,stationaryDegree,degrees,comparisonDegree,paddingFactors,terrainScales,tangentStep,maximumIterations,quadratureOrder)
% Build the Milestone-10 exact-residual terrain-mode enrichment oracle.

[seed,setup] = buildGlobalFirstOrderTerrainDressingAudit(problem, ...
    trustedBounds,supportBounds,stationaryDegree,degrees,comparisonDegree, ...
    paddingFactors,terrainScales,tangentStep,quadratureOrder);
if seed.classification == "global-dressing-incompatible"
    error("WVTerrainEnergyGalerkin:ResidualEnrichmentIncompatibleSeed", ...
        "Milestone 10 requires a compatible Milestone-9.4 global dressing seed.")
end

nDegree = numel(degrees);
nPadding = numel(paddingFactors);
details = cell(nDegree,nPadding);
for iDegree = 1:nDegree
    for iPadding = 1:nPadding
        details{iDegree,iPadding} = enrichmentCase(setup.details{iDegree,iPadding}, ...
            trustedBounds,stationaryDegree,degrees(iDegree),1,maximumIterations);
    end
end
comparison = cell(nPadding,1);
for iPadding = 1:nPadding
    comparison{iPadding} = enrichmentCase(setup.comparison{iPadding}, ...
        trustedBounds,stationaryDegree,comparisonDegree,1,maximumIterations);
end

primary = comparison{1};
crossPaddingDefect = maximumPrincipalSine(primary.internalModeBasis, ...
    comparison{2}.internalModeBasis,primary.energyMatrix);
denseCrossPaddingDefect = maximumPrincipalSine(primary.denseTargetBasis, ...
    comparison{2}.denseTargetBasis,primary.energyMatrix);
paddingDefect = max([comparison{1}.internalProjectorDefect ...
    comparison{2}.internalProjectorDefect ...
    abs(crossPaddingDefect-denseCrossPaddingDefect)]);
nestedDefect = zeros(nDegree,1);
denseNestedDefect = zeros(nDegree,1);
for iDegree = 1:nDegree-1
    transfer = seed.reference.transfers{iDegree}.matrix;
    nestedDefect(iDegree) = maximumPrincipalSine( ...
        transfer*details{iDegree,1}.internalModeBasis, ...
        details{iDegree+1,1}.internalModeBasis, ...
        details{iDegree+1,1}.energyMatrix);
    denseNestedDefect(iDegree) = maximumPrincipalSine( ...
        transfer*details{iDegree,1}.denseTargetBasis, ...
        details{iDegree+1,1}.denseTargetBasis, ...
        details{iDegree+1,1}.energyMatrix);
end
nestedDefect(end) = maximumPrincipalSine( ...
    seed.comparison{1}.referenceTransfer*details{end,1}.internalModeBasis, ...
    primary.internalModeBasis,primary.energyMatrix);
denseNestedDefect(end) = maximumPrincipalSine( ...
    seed.comparison{1}.referenceTransfer*details{end,1}.denseTargetBasis, ...
    primary.denseTargetBasis,primary.energyMatrix);
nestedExcessDefect = max([cellfun(@(value)value.internalProjectorDefect,details(:,1)); ...
    primary.internalProjectorDefect;abs(nestedDefect-denseNestedDefect)]);

tolerance = struct("structure",1e-12,"stationary",1e-10, ...
    "residual",1e-10,"projector",1e-9,"frequency",1e-9, ...
    "apv",1e-8,"bottom",1e-8,"strong",1e-5, ...
    "padding",1e-8,"nested",1e-8,"classification",1e-10);
allCases = [details(:);comparison(:)];
structurePasses = maximumCellField(allCases,"energyHermitianDefect") <= tolerance.structure ...
    && maximumCellField(allCases,"exchangeSkewHermitianDefect") <= tolerance.structure;
stationaryPasses = maximumCellField(allCases,"stationaryRowDefect") <= tolerance.stationary ...
    && maximumCellField(allCases,"stationaryOrthogonalityDefect") <= tolerance.stationary;
physicalPasses = primary.maximumInternalResidual <= tolerance.residual ...
    && primary.internalProjectorDefect <= tolerance.projector ...
    && primary.internalFrequencyDefect <= tolerance.frequency ...
    && primary.maximumInternalAPVDefect <= tolerance.apv ...
    && primary.maximumInternalBottomDefect <= tolerance.bottom ...
    && primary.maximumInternalStrongResidual <= tolerance.strong ...
    && paddingDefect <= tolerance.padding ...
    && nestedExcessDefect <= tolerance.nested;
classificationPasses = primary.decomposition.maximumOrthogonalityDefect <= tolerance.classification ...
    && primary.decomposition.completenessDefect <= tolerance.classification ...
    && primary.decomposition.conjugacyDefect <= tolerance.classification;

if structurePasses && stationaryPasses && physicalPasses && classificationPasses
    if primary.dimension.compressionFactor >= 2
        classification = "residual-enrichment-acceleration";
        diagnosis = "Exact-residual enrichment reproduces the validated dense internal-wave subspace with at least factor-two compression.";
    else
        classification = "residual-enrichment-equivalent";
        diagnosis = "Exact-residual enrichment reproduces the validated dense internal-wave subspace, but with less than factor-two compression.";
    end
else
    classification = "residual-enrichment-blocker";
    diagnosis = blockerDiagnosis(structurePasses,stationaryPasses,physicalPasses,classificationPasses,primary,paddingDefect,nestedExcessDefect,tolerance);
end

audit = struct;
audit.scope = "milestone-10-residual-enriched-terrain-modes";
audit.status = classification;
audit.classification = classification;
audit.isCompatible = classification ~= "residual-enrichment-blocker";
audit.diagnosis = diagnosis;
audit.seed = seed;
audit.details = details;
audit.comparison = comparison;
audit.primary = primary;
audit.paddingDefect = paddingDefect;
audit.crossPaddingProjectorDefect = crossPaddingDefect;
audit.denseCrossPaddingProjectorDefect = denseCrossPaddingDefect;
audit.nestedProjectorDefect = nestedDefect;
audit.denseNestedProjectorDefect = denseNestedDefect;
audit.nestedExcessDefect = nestedExcessDefect;
audit.maximumIterations = maximumIterations;
audit.requiredTolerance = tolerance;
audit.nextScope = "milestone-11-only-if-separately-authorized";
end

function result = enrichmentCase(setup,trustedBounds,stationaryDegree,degree,terrainScale,maximumIterations)
c = setup.context;
primitive = setup.primitive;
direction = primitive.evaluatedDirections(find(primitive.evaluationScales == terrainScale,1));
if isempty(direction)
    error("WVTerrainEnergyGalerkin:ResidualEnrichmentMissingTerrainScale", ...
        "The exact finite-amplitude terrain direction was not constructed.")
end
H = direction.E;
stationary = constructCompleteStationarySpace(c,direction,trustedBounds,degree,stationaryDegree,terrainScale);
G = stationary.trustedBasis;
internalSeed = internalColumns(setup.internal,terrainScale);
zeroSeed = zeroColumns(setup.zeroFrequency,terrainScale);
X = [internalSeed zeroSeed];
X = projectOut(X,G,H);
X = energyOrthonormalize(X,H);
targetParts = cellfun(@(block)block.basis,setup.targetBlocks,"UniformOutput",false);
targetSeed = energyOrthonormalize(horzcat(targetParts{:}),H);
dense = exactEigensystem(direction);
targetIndices = selectMostCaptured(dense.vectors,targetSeed,H,size(targetSeed,2),zeros(0,1));
target = dense.vectors(:,targetIndices);
flat = flatEigensystem(primitive.flatReference);

history = repmat(emptyIteration,maximumIterations+1,1);
for iIteration = 0:maximumIterations
    iteration = reducedIteration(c,direction,X,G,target,zeroSeed);
    history(iIteration+1) = iteration;
    if iteration.passesMandatoryGates || iIteration == maximumIterations
        break
    end
    selected = unique([iteration.internalModeIndices;iteration.zeroCandidateModeIndices]);
    selected = expandInseparable(selected,iteration.frequency,iteration.absoluteFrequencyUncertainty);
    correction = residualCorrection(iteration.modeVectors(:,selected), ...
        iteration.frequency(selected),iteration.residual(:,selected), ...
        iteration.absoluteFrequencyUncertainty(selected), ...
        iteration.zeroCandidateModeIndices,selected,flat,setup.zeroFrequency, ...
        primitive.flatReference,terrainScale);
    correction = [correction conjugateColumns(c,correction)]; %#ok<AGROW>
    correction = projectOut(correction,G,H);
    correction = projectOut(correction,X,H);
    correction = independentColumns(correction,H);
    if isempty(correction)
        break
    end
    X = energyOrthonormalize([X correction],H);
end
history = history(1:iIteration+1);
final = history(end);

stationaryProjector = energyProjector(G,H);
physicalProjector = energyProjector(final.internalModeBasis,H);
remainderProjector = eye(size(H))-stationaryProjector-physicalProjector;
decomposition = decompositionDiagnostics(stationaryProjector,physicalProjector, ...
    remainderProjector,H,c.coordinateConjugateIndex);
dimension = struct("ambient",size(H,1),"stationary",size(G,2), ...
    "dynamicalTrial",size(X,2),"validatedInternal",size(final.internalModeBasis,2), ...
    "unresolved",size(H,1)-size(G,2)-size(final.internalModeBasis,2), ...
    "compressionFactor",size(H,1)/max(size(G,2)+size(X,2),1));

result = final;
result.history = history;
result.numberOfIterations = numel(history)-1;
result.trialBasis = X;
result.stationaryBasis = G;
result.stationaryProjector = stationaryProjector;
result.physicalProjector = physicalProjector;
result.unresolvedProjector = remainderProjector;
result.decomposition = decomposition;
result.dimension = dimension;
result.energyMatrix = H;
result.denseTargetBasis = target;
result.denseTargetFrequency = dense.frequency(targetIndices);
result.energyHermitianDefect = norm(H-H',"fro")/max(norm(H,"fro"),realmin);
result.exchangeSkewHermitianDefect = norm(direction.J+direction.J',"fro")/max(norm(direction.J,"fro"),realmin);
result.stationaryRowDefect = norm(direction.J*G,"fro")/max(norm(direction.J,"fro")*norm(G,"fro"),realmin);
result.stationaryOrthogonalityDefect = norm(G'*H*final.internalModeBasis,"fro")/max(norm(G,"fro")*norm(H*final.internalModeBasis,"fro"),realmin);
result.polynomialDegree = degree;
result.paddingFactor = c.layout.paddingFactor;
result.supportModeBounds = [max(abs(c.horizontalLayout.kMode)) max(abs(c.horizontalLayout.lMode))];
end

function iteration = reducedIteration(c,direction,X,G,target,zeroSeed)
H = direction.E;
A = 1i*direction.J;
Hr = X'*H*X;
Ar = X'*A*X;
[R,flag] = chol((Hr+Hr')/2);
if flag ~= 0
    error("WVTerrainEnergyGalerkin:ResidualEnrichmentReducedEnergyNotPositive", ...
        "The enriched reduced energy matrix must remain positive definite.")
end
K = R'\(Ar/R);
[U,T] = schur(K,"complex");
frequency = real(diag(T));
[frequency,order] = sort(frequency);
vectors = X*(R\U(:,order));
for iMode = 1:size(vectors,2)
    vectors(:,iMode) = vectors(:,iMode)/sqrt(real(vectors(:,iMode)'*H*vectors(:,iMode)));
end
residual = A*vectors-(H*vectors).*frequency.';
energyFactor = chol((H+H')/2);
energyResidual = energyFactor'\residual;
absoluteUncertainty = vecnorm(energyResidual,2,1).'./max(vecnorm(energyFactor*vectors,2,1).',realmin);
operatorNorm = norm(energyFactor'\(A/energyFactor),2);
relativeResidual = absoluteUncertainty./max(operatorNorm+abs(frequency),realmin);

internalCount = size(target,2);
internalIndices = selectMostCaptured(vectors,target,H,internalCount,zeros(0,1));
zeroCount = size(zeroSeed,2);
zeroIndices = selectMostCaptured(vectors,zeroSeed,H,zeroCount,internalIndices);
internalBasis = vectors(:,internalIndices);
zeroBasis = vectors(:,zeroIndices);
targetFrequency = sort(real(eig(target'*A*target,target'*H*target)));
internalFrequency = sort(frequency(internalIndices));
frequencyDefect = norm(internalFrequency-targetFrequency)/max(norm(targetFrequency),realmin);
projectorDefect = maximumPrincipalSine(internalBasis,target,H);
apv = trustedAPVDefect(c,direction,internalBasis);
bottom = modalBottomDefect(direction,internalBasis,frequency(internalIndices));
strong = strongModeDiagnostics(c,direction,internalBasis,frequency(internalIndices));
maximumResidual = max(relativeResidual(internalIndices),[],"all");
if isempty(zeroIndices)
    zeroResidual = 0;
    zeroAPV = 0;
    zeroBottom = 0;
    zeroStrong = 0;
    zeroResolution = Inf;
    zeroParticipation = 0;
else
    zeroResidual = max(relativeResidual(zeroIndices),[],"all");
    zeroAPV = trustedAPVDefect(c,direction,zeroBasis);
    zeroBottom = modalBottomDefect(direction,zeroBasis,frequency(zeroIndices));
    zeroStrongDiagnostics = strongModeDiagnostics(c,direction,zeroBasis,frequency(zeroIndices));
    zeroStrong = zeroStrongDiagnostics.maximumResidual;
    zeroResolution = min(abs(frequency(zeroIndices))./max(absoluteUncertainty(zeroIndices),realmin));
    zeroParticipation = maximumBottomParticipation(direction,zeroBasis);
end

iteration = emptyIteration;
iteration.frequency = frequency;
iteration.modeVectors = vectors;
iteration.residual = residual;
iteration.absoluteFrequencyUncertainty = absoluteUncertainty;
iteration.relativeResidual = relativeResidual;
iteration.internalModeIndices = internalIndices;
iteration.zeroCandidateModeIndices = zeroIndices;
iteration.internalModeBasis = internalBasis;
iteration.internalProjectorDefect = projectorDefect;
iteration.internalFrequencyDefect = frequencyDefect;
iteration.maximumInternalResidual = maximumResidual;
iteration.maximumInternalAPVDefect = apv;
iteration.maximumInternalBottomDefect = bottom;
iteration.maximumInternalStrongResidual = strong.maximumResidual;
iteration.maximumZeroCandidateResidual = zeroResidual;
iteration.maximumZeroCandidateAPVDefect = zeroAPV;
iteration.maximumZeroCandidateBottomDefect = zeroBottom;
iteration.maximumZeroCandidateStrongResidual = zeroStrong;
iteration.minimumZeroCandidateFrequencyResolutionRatio = zeroResolution;
iteration.maximumZeroCandidateBottomParticipation = zeroParticipation;
iteration.zeroCandidateStatus = "unresolved";
iteration.trialDimension = size(X,2);
iteration.passesMandatoryGates = maximumResidual <= 1e-10 ...
    && projectorDefect <= 1e-9 && frequencyDefect <= 1e-9 ...
    && apv <= 1e-8 && bottom <= 1e-8 && strong.maximumResidual <= 1e-5;
iteration.stationaryOverlapDefect = norm(G'*H*internalBasis,"fro")/max(norm(G,"fro")*norm(H*internalBasis,"fro"),realmin);
end

function correction = residualCorrection(vectors,frequency,residual,uncertainty,zeroCandidateIndices,selectedIndices,flat,zero,flatDirection,terrainScale)
correction = zeros(size(vectors));
flatZero = find(abs(flat.frequency) <= 1e3*flat.absoluteFrequencyUncertainty);
flatNonzero = setdiff((1:numel(flat.frequency))',flatZero,"stable");
zeroBasis = zero.baseBasis;
zeroFrequency = terrainScale*zero.firstOrderFrequency;
zeroUncertainty = terrainScale*zero.absoluteFrequencyUncertainty;
zeroPreconditionerBasis = [zeroBasis flat.vectors(:,flatNonzero)];
zeroPreconditionerFrequency = [zeroFrequency;flat.frequency(flatNonzero)];
zeroPreconditionerUncertainty = [zeroUncertainty;flat.absoluteFrequencyUncertainty(flatNonzero)];
zeroSelected = ismember(selectedIndices,zeroCandidateIndices);
for iColumn = 1:size(vectors,2)
    if zeroSelected(iColumn)
        basis = zeroPreconditionerBasis;
        referenceFrequency = zeroPreconditionerFrequency;
        referenceUncertainty = zeroPreconditionerUncertainty;
    else
        basis = flat.vectors;
        referenceFrequency = flat.frequency;
        referenceUncertainty = flat.absoluteFrequencyUncertainty;
    end
    denominator = referenceFrequency-frequency(iColumn);
    inseparable = abs(denominator) <= 1e3*(referenceUncertainty+uncertainty(iColumn));
    coefficient = -(basis'*residual(:,iColumn));
    coefficient(inseparable) = 0;
    coefficient(~inseparable) = coefficient(~inseparable)./denominator(~inseparable);
    correction(:,iColumn) = basis*coefficient;
end
if size(correction,1) ~= size(flatDirection.E,1)
    error("WVTerrainEnergyGalerkin:ResidualEnrichmentPreconditionerDimensionMismatch", ...
        "The residual preconditioner must retain the complete primitive ambient dimension.")
end
end

function indices = expandInseparable(indices,frequency,uncertainty)
indices = indices(:);
changed = true;
while changed
    changed = false;
    candidates = setdiff((1:numel(frequency))',indices,"stable");
    for candidate = candidates.'
        if any(abs(frequency(candidate)-frequency(indices)) <= 1e3*(uncertainty(candidate)+uncertainty(indices)))
            indices(end+1,1) = candidate; %#ok<AGROW>
            changed = true;
        end
    end
end
indices = sort(unique(indices));
end

function columns = internalColumns(internal,scale)
parts = cell(numel(internal),1);
for iBlock = 1:numel(internal)
    parts{iBlock} = internal(iBlock).baseBasis+scale*internal(iBlock).correction;
end
columns = horzcat(parts{:});
end

function columns = zeroColumns(zero,scale)
if isempty(zero.activeIndices)
    columns = zeros(size(zero.baseBasis,1),0);
else
    columns = zero.baseBasis(:,zero.activeIndices)+scale*zero.correction(:,zero.activeIndices);
end
end

function conjugate = conjugateColumns(c,columns)
n = size(columns,1);
map = sparse((1:n)',c.coordinateConjugateIndex,1,n,n);
conjugate = map*conj(columns);
end

function X = projectOut(X,basis,H)
if ~isempty(basis) && ~isempty(X)
    X = X-basis*((basis'*H*basis)\(basis'*H*X));
end
end

function X = independentColumns(X,H)
if isempty(X)
    return
end
R = chol((H+H')/2);
[~,T,p] = qr(R*X,"econ","vector");
tolerance = max(size(T))*eps(max(norm(T,2),1));
rankX = nnz(abs(diag(T)) > tolerance);
if rankX == 0
    X = zeros(size(X,1),0);
else
    X = X(:,p(1:rankX));
end
end

function X = energyOrthonormalize(X,H)
if isempty(X)
    return
end
R = chol((H+H')/2);
Q = orth(R*X);
X = R\Q;
end

function flat = flatEigensystem(direction)
H = direction.E;
A = 1i*direction.J;
R = chol((H+H')/2);
K = R'\(A/R);
[U,T] = schur(K,"complex");
frequency = real(diag(T));
[frequency,order] = sort(frequency);
vectors = R\U(:,order);
for iMode = 1:size(vectors,2)
    vectors(:,iMode) = vectors(:,iMode)/sqrt(real(vectors(:,iMode)'*H*vectors(:,iMode)));
end
residual = K*(R*vectors)-(R*vectors).*frequency.';
absoluteUncertainty = max(vecnorm(residual,2,1).'./max(vecnorm(R*vectors,2,1).',realmin),10*eps(max(norm(K,2),realmin)));
flat = struct("frequency",frequency,"vectors",vectors, ...
    "absoluteFrequencyUncertainty",absoluteUncertainty);
end

function dense = exactEigensystem(direction)
H = direction.E;
A = 1i*direction.J;
R = chol((H+H')/2);
K = R'\(A/R);
[U,T] = schur(K,"complex");
frequency = real(diag(T));
[frequency,order] = sort(frequency);
vectors = R\U(:,order);
for iMode = 1:size(vectors,2)
    vectors(:,iMode) = vectors(:,iMode)/sqrt(real(vectors(:,iMode)'*H*vectors(:,iMode)));
end
dense = struct("frequency",frequency,"vectors",vectors);
end

function indices = selectMostCaptured(vectors,target,H,count,excluded)
if count == 0
    indices = zeros(0,1);
    return
end
projector = energyProjector(target,H);
participation = real(sum(conj(vectors).*(H*projector*vectors),1)).';
participation(excluded) = -Inf;
[~,order] = maxk(participation,min(count,numel(participation)-numel(excluded)));
indices = sort(order);
end

function projector = energyProjector(basis,H)
if isempty(basis)
    projector = zeros(size(H));
else
    projector = basis*((basis'*H*basis)\(basis'*H));
end
end

function defect = maximumPrincipalSine(first,second,H)
if isempty(first) && isempty(second)
    defect = 0;
    return
elseif isempty(first) || isempty(second) || size(first,2) ~= size(second,2)
    defect = Inf;
    return
end
R = chol((H+H')/2);
firstQ = orth(R*first);
secondQ = orth(R*second);
singularValues = svd(firstQ'*secondQ);
defect = sqrt(max(0,1-min(singularValues,[],"all")^2));
end

function defect = trustedAPVDefect(c,direction,vectors)
trusted = c.layout.trustedHorizontalModes;
rows = find(repmat(trusted,c.nZ,1));
Q = direction.QProjected(rows,:);
weight = c.apvVerticalWeight(rows);
energyFactor = chol(direction.E);
scale = max(norm(sqrt(weight).*Q/energyFactor,2),realmin);
defect = max(vecnorm(sqrt(weight).*(Q*vectors),2,1).'/scale,[],"all");
end

function defect = modalBottomDefect(direction,vectors,frequency)
residual = -1i*(direction.B*vectors).*frequency.'-direction.R*vectors;
energyFactor = chol(direction.E);
scale = abs(frequency)*norm(direction.B/energyFactor,2)+norm(direction.R/energyFactor,2);
values = vecnorm(residual,2,1).'./max(scale,realmin);
defect = max(values,[],"all");
end

function participation = maximumBottomParticipation(direction,vectors)
if isempty(vectors)
    participation = 0;
    return
end
energyFactor = chol(direction.E);
scale = max(norm(direction.B/energyFactor,2),realmin);
participation = max(vecnorm(direction.B*vectors,2,1).'/scale,[],"all");
end

function diagnostics = strongModeDiagnostics(c,direction,C,frequency)
if isempty(C)
    diagnostics = struct("maximumResidual",0,"maximumByMode",zeros(0,1));
    return
end
rawState = c.N*C;
pressure = direction.pressureMap*C;
stateTendency = -1i*rawState.*frequency.';
descriptorResidual = direction.matrices.S*[stateTendency;pressure]-direction.matrices.F*C;
primitiveResidual = direction.L*C+1i*C.*frequency.';
descriptorDefect = columnRelativeResidual(descriptorResidual, ...
    {direction.matrices.S*[stateTendency;pressure],-direction.matrices.F*C});
primitiveDefect = columnRelativeResidual(primitiveResidual, ...
    {direction.L*C,1i*C.*frequency.'});
gaugeResidual = c.pressureGauge'*[zeros(c.nX,size(C,2));pressure];
gaugeDefect = vecnorm(gaugeResidual,2,1).'./max(vecnorm(pressure,2,1).',realmin);
pointwise = pointwiseStrongResiduals(c,direction,rawState,pressure,frequency);
maximum = max([descriptorDefect primitiveDefect gaugeDefect pointwise.maximumByMode],[],2);
diagnostics = struct("maximumResidual",max(maximum,[],"all"), ...
    "maximumByMode",maximum);
end

function diagnostics = pointwiseStrongResiduals(c,direction,rawState,pressureCoefficients,frequency)
H3 = repmat(c.H,c.nZ,1);
HX3 = repmat(c.HX,c.nZ,1);
HY3 = repmat(c.HY,c.nZ,1);
gamma = 1-direction.scale*H3;
gradLnGammaX = -direction.scale*HX3./gamma;
gradLnGammaY = -direction.scale*HY3./gamma;
uHat = c.Ru*rawState;
vHat = c.Rv*rawState;
wHat = c.Rwh*rawState;
eta = c.Reta*rawState;
u = uHat./gamma;
v = vHat./gamma;
w = wHat+c.xiGrid.*(gradLnGammaX.*uHat+gradLnGammaY.*vHat);
pX = c.RPX*pressureCoefficients;
pY = c.RPY*pressureCoefficients;
pXi = c.RPXi*pressureCoefficients;
Dxp = pX-c.xiGrid.*gradLnGammaX.*pXi;
Dyp = pY-c.xiGrid.*gradLnGammaY.*pXi;
N2 = c.wvt.N2Function(gamma.*c.xiGrid);
if isscalar(N2)
    N2 = repmat(N2,size(c.xiGrid));
end
N2 = N2(:);
omega = frequency.';
uDefect = weightedColumnResidual({-1i*u.*omega,-c.wvt.f*v,Dxp/c.wvt.rho0},c.volumeWeight);
vDefect = weightedColumnResidual({-1i*v.*omega,c.wvt.f*u,Dyp/c.wvt.rho0},c.volumeWeight);
wDefect = weightedColumnResidual({-1i*w.*omega,N2.*eta,pXi./(c.wvt.rho0*gamma)},c.volumeWeight);
etaDefect = weightedColumnResidual({-1i*eta.*omega,-w},c.volumeWeight);
divergence = (c.RuX+c.RvY+c.RwhXi)*rawState;
continuityDefect = weightedColumnNorm(divergence,c.volumeWeight)./max(weightedColumnNorm(c.RuX*rawState,c.volumeWeight)+weightedColumnNorm(c.RvY*rawState,c.volumeWeight)+weightedColumnNorm(c.RwhXi*rawState,c.volumeWeight),realmin);
maximumByMode = max([uDefect vDefect wDefect etaDefect continuityDefect],[],2);
diagnostics = struct("maximumByMode",maximumByMode, ...
    "maximumResidual",max(maximumByMode,[],"all"));
end

function defect = weightedColumnResidual(terms,weight)
residual = zeros(size(terms{1}));
scale = zeros(size(terms{1},2),1);
for iTerm = 1:numel(terms)
    residual = residual+terms{iTerm};
    scale = scale+weightedColumnNorm(terms{iTerm},weight);
end
defect = weightedColumnNorm(residual,weight)./max(scale,realmin);
end

function values = weightedColumnNorm(array,weight)
values = sqrt(real(sum(conj(array).*(weight.*array),1))).';
end

function defect = columnRelativeResidual(residual,terms)
scale = zeros(size(residual,2),1);
for iTerm = 1:numel(terms)
    scale = scale+vecnorm(terms{iTerm},2,1).';
end
defect = vecnorm(residual,2,1).'./max(scale,realmin);
end

function diagnostics = decompositionDiagnostics(stationary,physical,remainder,H,conjugateIndex)
orthogonality = [norm(stationary*physical,"fro"), ...
    norm(stationary*remainder,"fro"),norm(physical*remainder,"fro")];
completeness = norm(stationary+physical+remainder-eye(size(H)),"fro")/sqrt(size(H,1));
n = size(H,1);
conjugacy = sparse((1:n)',conjugateIndex,1,n,n);
projectors = {stationary,physical,remainder};
conjugacyDefect = 0;
for iProjector = 1:numel(projectors)
    residual = projectors{iProjector}*conjugacy-conjugacy*conj(projectors{iProjector});
    conjugacyDefect = max(conjugacyDefect,norm(residual,"fro")/max(2*norm(projectors{iProjector},"fro"),realmin));
end
diagnostics = struct("maximumOrthogonalityDefect",max(orthogonality), ...
    "completenessDefect",completeness,"conjugacyDefect",conjugacyDefect);
end

function iteration = emptyIteration
iteration = struct("frequency",zeros(0,1),"modeVectors",zeros(0), ...
    "residual",zeros(0),"absoluteFrequencyUncertainty",zeros(0,1), ...
    "relativeResidual",zeros(0,1),"internalModeIndices",zeros(0,1), ...
    "zeroCandidateModeIndices",zeros(0,1),"internalModeBasis",zeros(0), ...
    "internalProjectorDefect",Inf,"internalFrequencyDefect",Inf, ...
    "maximumInternalResidual",Inf,"maximumInternalAPVDefect",Inf, ...
    "maximumInternalBottomDefect",Inf,"maximumInternalStrongResidual",Inf, ...
    "maximumZeroCandidateResidual",Inf,"maximumZeroCandidateAPVDefect",Inf, ...
    "maximumZeroCandidateBottomDefect",Inf,"maximumZeroCandidateStrongResidual",Inf, ...
    "minimumZeroCandidateFrequencyResolutionRatio",0, ...
    "maximumZeroCandidateBottomParticipation",Inf, ...
    "zeroCandidateStatus","unresolved", ...
    "trialDimension",0,"passesMandatoryGates",false, ...
    "stationaryOverlapDefect",Inf);
end

function value = maximumCellField(values,field)
value = max(cellfun(@(item)item.(field),values),[],"all");
end

function diagnosis = blockerDiagnosis(structurePasses,stationaryPasses,physicalPasses,classificationPasses,primary,paddingDefect,nestedDefect,tolerance)
if ~structurePasses
    diagnosis = "The exact reduced finite-terrain forms lost their required Hermitian or skew-Hermitian structure.";
elseif ~stationaryPasses
    diagnosis = "Residual enrichment failed to preserve the complete finite-terrain stationary space.";
elseif ~classificationPasses
    diagnosis = "The stationary, enriched physical, and unresolved projectors do not form a complete conjugate-closed energy decomposition.";
elseif ~physicalPasses
    diagnosis = sprintf("Residual enrichment stopped without passing the physical gates: residual %.3g, projector %.3g, frequency %.3g, APV %.3g, bottom %.3g, strong %.3g, padding %.3g, nested %.3g.", ...
        primary.maximumInternalResidual,primary.internalProjectorDefect, ...
        primary.internalFrequencyDefect,primary.maximumInternalAPVDefect, ...
        primary.maximumInternalBottomDefect,primary.maximumInternalStrongResidual, ...
        paddingDefect,max(nestedDefect));
else
    diagnosis = sprintf("Residual enrichment did not meet the declared tolerance %.3g.",tolerance.residual);
end
end
