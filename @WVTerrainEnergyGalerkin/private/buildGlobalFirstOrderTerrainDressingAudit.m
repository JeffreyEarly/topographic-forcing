function [audit,setup] = buildGlobalFirstOrderTerrainDressingAudit(problem,trustedBounds,supportBounds,stationaryDegree,degrees,comparisonDegree,paddingFactors,terrainScales,tangentStep,quadratureOrder)
% Build the Milestone-9.4 globally coupled terrain-dressing oracle.

reference = problem.auditConvergedPhysicalSubspaces( ...
    trustedModeBounds=trustedBounds, ...
    supportModeBounds=supportBounds, ...
    stationaryPolynomialDegree=stationaryDegree, ...
    primitivePolynomialDegrees=degrees, ...
    paddingFactors=paddingFactors, ...
    terrainScales=terrainScales, ...
    quadratureOrder=quadratureOrder);
validatedTracks = find([reference.internal.tracks.isValidated]);
if isempty(validatedTracks)
    error("WVTerrainEnergyGalerkin:GlobalDressingMissingValidatedInternalBlock", ...
        "Milestone 9.4 requires at least one Milestone-9.1 validated internal-wave block.")
end

nDegree = numel(degrees);
nPadding = numel(paddingFactors);
details = cell(nDegree,nPadding);
detailSetup = cell(nDegree,nPadding);
for iDegree = 1:nDegree
    for iPadding = 1:nPadding
        [details{iDegree,iPadding},detailSetup{iDegree,iPadding}] = refinementAudit(problem,reference, ...
            validatedTracks,trustedBounds,supportBounds(iDegree,:), ...
            stationaryDegree,degrees(iDegree),paddingFactors(iPadding), ...
            terrainScales,tangentStep,quadratureOrder,iDegree);
    end
end

comparison = cell(nPadding,1);
comparisonSetup = cell(nPadding,1);
for iPadding = 1:nPadding
    [comparison{iPadding},comparisonSetup{iPadding}] = comparisonAudit(problem,reference, ...
        validatedTracks,trustedBounds,supportBounds(end,:), ...
        stationaryDegree,comparisonDegree,paddingFactors(iPadding), ...
        terrainScales,tangentStep,quadratureOrder);
end

weakSequence = cell(nPadding,1);
for iPadding = 1:nPadding
    weakSequence{iPadding} = buildBoundaryCompleteWeakEigenproblemAudit( ...
        problem,trustedBounds,supportBounds(end,:),degrees, ...
        paddingFactors(iPadding),quadratureOrder);
end

allResults = [details(:);comparison(:)];
primary = comparison{1};
convergence = convergenceDiagnostics(details,comparison);
padding = paddingDiagnostics(comparison);
weak = weakSequenceDiagnostics(weakSequence);
tolerance = struct( ...
    "tangent",1e-9, ...
    "orthogonality",1e-11, ...
    "conjugacy",1e-11, ...
    "fourierSelection",1e-11, ...
    "correction",1e-10, ...
    "bottomFirstOrder",1e-10, ...
    "stationary",1e-10, ...
    "scalingOrder",1.8, ...
    "improvement",4, ...
    "structure",1e-12, ...
    "projector",1e-8, ...
    "apv",1e-8, ...
    "bottom",1e-8, ...
    "strong",1e-5, ...
    "padding",1e-8);

coefficientPasses = maximumCellField(allResults,"tangentAgreementDefect") <= tolerance.tangent ...
    && maximumCellField(allResults,"flatEnergyOrthogonalityDefect") <= tolerance.orthogonality ...
    && maximumCellField(allResults,"conjugacyDefect") <= tolerance.conjugacy ...
    && maximumCellField(allResults,"fourierSelectionDefect") <= tolerance.fourierSelection ...
    && maximumCellField(allResults,"correctionResidual") <= tolerance.correction ...
    && maximumCellField(allResults,"firstOrderBottomDefect") <= tolerance.bottomFirstOrder ...
    && weak.maximumGreenIdentityDefect <= tolerance.stationary ...
    && weak.maximumStationaryRowDefect <= tolerance.stationary;
scalingPasses = primary.scaling.minimumDressedWeakOrder >= tolerance.scalingOrder ...
    && primary.scaling.minimumDressedBottomOrder >= tolerance.scalingOrder ...
    && primary.scaling.weakImprovement >= tolerance.improvement ...
    && primary.scaling.bottomImprovement >= tolerance.improvement;
structurePasses = maximumCellField(allResults,"reducedEnergyHermitianDefect") <= tolerance.structure ...
    && maximumCellField(allResults,"reducedExchangeSkewHermitianDefect") <= tolerance.structure;
finitePhysicalPasses = primary.finite.internalProjectorDefect <= tolerance.projector ...
    && primary.finite.maximumInternalAPVDefect <= tolerance.apv ...
    && primary.finite.maximumInternalBottomDefect <= tolerance.bottom ...
    && primary.finite.maximumInternalStrongResidual <= tolerance.strong ...
    && padding.projectorDefect <= tolerance.padding ...
    && convergence.projectorIsMonotone ...
    && convergence.bottomIsMonotone;
projectorImprovementPasses = primary.finite.internalProjectorImprovement >= 2;

if ~coefficientPasses || ~scalingPasses
    classification = "global-dressing-incompatible";
    diagnosis = "The analytic global correction does not satisfy the required first-order block, bottom, stationary, or O(h^2) residual identities.";
elseif finitePhysicalPasses && structurePasses && primary.dimension.compressionFactor >= 2
    classification = "global-dressing-acceleration";
    diagnosis = "The globally dressed blocks reproduce the validated physical projectors with at least a factor-two reduction in retained dimension.";
elseif finitePhysicalPasses && structurePasses
    classification = "global-dressing-equivalent";
    diagnosis = "The globally dressed blocks reproduce the validated physical projectors, but with less than a factor-two reduction in retained dimension.";
elseif structurePasses && projectorImprovementPasses
    classification = "global-dressing-seed";
    diagnosis = "The global first-order correction gives the required O(h^2) weak and bottom behavior and preserves exact reduced energy structure, but one correction does not yet pass every finite-amplitude physical gate.";
else
    classification = "global-dressing-incompatible";
    diagnosis = "The exact finite-amplitude restriction loses the physical Hermitian or skew-Hermitian structure.";
end

audit = struct;
audit.scope = "milestone-9.4-global-first-order-terrain-dressing";
audit.status = classification;
audit.classification = classification;
audit.isCompatible = classification ~= "global-dressing-incompatible";
audit.diagnosis = diagnosis;
audit.trustedModeBounds = trustedBounds;
audit.supportModeBounds = supportBounds;
audit.stationaryPolynomialDegree = stationaryDegree;
audit.primitivePolynomialDegrees = degrees;
audit.comparisonPolynomialDegree = comparisonDegree;
audit.paddingFactors = paddingFactors;
audit.terrainScales = terrainScales;
audit.tangentStep = tangentStep;
audit.centeredTangentSteps = 2*tangentStep./(2.^(0:2));
audit.quadratureOrder = quadratureOrder;
audit.reference = reference;
audit.validatedInternalTrackIndices = validatedTracks;
audit.details = details;
audit.comparison = comparison;
audit.primary = primary;
audit.weakSequence = weakSequence;
audit.weakSequenceDiagnostics = weak;
audit.convergence = convergence;
audit.padding = padding;
audit.requiredTolerance = tolerance;
audit.nextScope = "milestone-10-only-if-separately-authorized";
setup = struct("details",{detailSetup},"comparison",{comparisonSetup});
end

function value = maximumCellField(results,field)
value = max(cellfun(@(result)result.(field),results),[],"all");
end

function [result,setup] = refinementAudit(problem,reference,validatedTracks,trustedBounds,support,stationaryDegree,degree,paddingFactor,terrainScales,tangentStep,quadratureOrder,iDegree)
layout = retainedLayout(problem,support);
order = resolvedQuadratureOrder(degree,quadratureOrder);
[primitive,context] = buildGlobalSmallTerrainPrimitiveAudit(problem, ...
    degree,order,2*tangentStep,horizontalLayout=layout, ...
    paddingFactor=paddingFactor,trustedModeBounds=trustedBounds, ...
    rejectTerrainNyquist=true,evaluationScales=terrainScales);
flat = flatEigensystem(primitive.flatReference);
targetBlocks = cell(numel(validatedTracks),1);
for iTrack = 1:numel(validatedTracks)
    track = reference.internal.tracks(validatedTracks(iTrack));
    targetBlocks{iTrack} = track.blocks{iDegree};
end
result = dressingAudit(context,primitive,targetBlocks,trustedBounds, ...
    stationaryDegree,degree,terrainScales);
setup = struct("context",context,"primitive",primitive, ...
    "targetBlocks",{targetBlocks},"internal",result.internal, ...
    "zeroFrequency",result.zeroFrequency);
result.polynomialDegree = degree;
result.supportModeBounds = support;
result.paddingFactor = paddingFactor;
result.quadratureOrder = order;
result.flat = flat;
end

function [result,setup] = comparisonAudit(problem,reference,validatedTracks,trustedBounds,support,stationaryDegree,degree,paddingFactor,terrainScales,tangentStep,quadratureOrder)
layout = retainedLayout(problem,support);
order = resolvedQuadratureOrder(degree,quadratureOrder);
[primitive,context] = buildGlobalSmallTerrainPrimitiveAudit(problem, ...
    degree,order,2*tangentStep,horizontalLayout=layout, ...
    paddingFactor=paddingFactor,trustedModeBounds=trustedBounds, ...
    rejectTerrainNyquist=true,evaluationScales=terrainScales);
coarse = reference.finest.nestedRepresentation;
fine = nestedRepresentation(context,primitive.flatReference);
transfer = nestedTransfer(coarse,fine);
targetBlocks = cell(numel(validatedTracks),1);
for iTrack = 1:numel(validatedTracks)
    track = reference.internal.tracks(validatedTracks(iTrack));
    target = track.blocks{end};
    target.basis = transfer*target.basis;
    targetBlocks{iTrack} = target;
end
result = dressingAudit(context,primitive,targetBlocks,trustedBounds, ...
    stationaryDegree,degree,terrainScales);
setup = struct("context",context,"primitive",primitive, ...
    "targetBlocks",{targetBlocks},"internal",result.internal, ...
    "zeroFrequency",result.zeroFrequency);
result.polynomialDegree = degree;
result.supportModeBounds = support;
result.paddingFactor = paddingFactor;
result.quadratureOrder = order;
result.referenceTransfer = transfer;
result.referenceTransferDefect = norm(fine.admissibleBasis*transfer-coarseEmbeddedRaw(coarse,fine),"fro")/max(norm(coarseEmbeddedRaw(coarse,fine),"fro"),realmin);
end

function result = dressingAudit(c,primitive,targetBlocks,trustedBounds,stationaryDegree,degree,terrainScales)
flatDirection = primitive.flatReference;
first = primitive.analyticTangent;
flat = flatEigensystem(flatDirection);
pairBlocks = conjugatePairBlocks(flat);
internal = repmat(emptyDressedPair,0,1);
usedPairs = zeros(0,1);
for iTarget = 1:numel(targetBlocks)
    [pairIndex,matchDefect] = matchFlatPair(targetBlocks{iTarget},pairBlocks,flatDirection.E);
    if isempty(pairIndex) || ismember(pairIndex,usedPairs)
        continue
    end
    usedPairs(end+1,1) = pairIndex; %#ok<AGROW>
    internal(end+1,1) = dressPair(c,flatDirection,first,flat,pairBlocks(pairIndex),matchDefect); %#ok<AGROW>
end
if isempty(internal)
    error("WVTerrainEnergyGalerkin:GlobalDressingInternalBlockMatchFailure", ...
        "No validated internal-wave block could be matched to the flat eigensystem.")
end

zero = dressZeroBlock(c,flatDirection,first,flat);
exact = finiteAmplitudeAudits(c,primitive,internal,zero, ...
    trustedBounds,stationaryDegree,degree,terrainScales,targetBlocks);
scaling = scalingDiagnostics(primitive,internal,zero,terrainScales);

result = struct;
result.tangentAgreementDefect = requiredCenteredAgreement(primitive.tangentAgreement);
result.auxiliaryTangentAgreementDefect = bestCenteredAgreement(primitive.tangentAgreement);
result.flatEnergyOrthogonalityDefect = flat.energyOrthogonalityDefect;
result.internal = internal;
result.zeroFrequency = zero;
result.finiteByScale = exact.byScale;
result.finite = exact.fullScale;
result.scaling = scaling;
result.conjugacyDefect = max([internal.conjugacyDefect zero.conjugacyDefect]);
result.fourierSelectionDefect = max([internal.fourierSelectionDefect zero.fourierSelectionDefect]);
result.correctionResidual = max([internal.correctionResidual zero.correctionResidual]);
result.firstOrderBottomDefect = max([internal.firstOrderBottomDefect zero.firstOrderBottomDefect]);
result.reducedEnergyHermitianDefect = max([exact.byScale.energyHermitianDefect]);
result.reducedExchangeSkewHermitianDefect = max([exact.byScale.exchangeSkewHermitianDefect]);
result.dimension = exact.fullScale.dimension;
result.terrainScales = terrainScales;
result.projection = primitive.projection;
end

function pair = dressPair(c,flatDirection,first,flat,pairBlock,matchDefect)
positive = dressSignedBlock(c,flatDirection,first,flat,pairBlock.positiveIndices);
negative = dressSignedBlock(c,flatDirection,first,flat,pairBlock.negativeIndices);
base = [positive.baseBasis negative.baseBasis];
correction = [positive.correction negative.correction];
dressed = base+correction;
pair = emptyDressedPair;
pair.positive = positive;
pair.negative = negative;
pair.baseBasis = base;
pair.correction = correction;
pair.dressedBasisAtUnitScale = dressed;
pair.matchDefect = matchDefect;
pair.conjugacyDefect = conjugacyDefect(c,dressed,flatDirection.E);
pair.fourierSelectionDefect = max(positive.fourierSelectionDefect,negative.fourierSelectionDefect);
pair.correctionResidual = max(positive.correctionResidual,negative.correctionResidual);
pair.firstOrderBottomDefect = max(positive.firstOrderBottomDefect,negative.firstOrderBottomDefect);
pair.dimension = size(base,2);
end

function result = dressSignedBlock(c,flatDirection,first,flat,indices)
X = flat.vectors(:,indices);
omega = flat.frequency(indices);
omega0 = mean(omega);
mass = X'*flatDirection.E*X;
coupling = X'*(1i*first.J-omega0*first.E)*X;
[D,omega1] = eig(coupling,mass,"vector");
base = X*D;
for iColumn = 1:size(base,2)
    base(:,iColumn) = base(:,iColumn)/sqrt(real(base(:,iColumn)'*flatDirection.E*base(:,iColumn)));
end
complement = setdiff((1:numel(flat.frequency))',indices,"stable");
V = flat.vectors(:,complement);
denominator = flat.frequency(complement)-omega0;
correction = zeros(size(base));
correctionResidual = zeros(size(base,2),1);
complementResidual = zeros(size(base,2),1);
bottomResidual = zeros(size(base,2),1);
for iColumn = 1:size(base,2)
    right = -(1i*first.J-omega0*first.E-omega1(iColumn)*flatDirection.E)*base(:,iColumn);
    coefficient = (V'*right)./denominator;
    correction(:,iColumn) = V*coefficient;
    residual = (1i*flatDirection.J-omega0*flatDirection.E)*correction(:,iColumn)-right;
    correctionResidual(iColumn) = norm(residual)/max(norm(right),realmin);
    complementResidual(iColumn) = norm(X'*residual)/max(norm(X,"fro")*norm(residual)+norm(right),realmin);
    bottom = -1i*omega0*flatDirection.B*correction(:,iColumn)-1i*omega1(iColumn)*flatDirection.B*base(:,iColumn)-first.R*base(:,iColumn);
    scale = abs(omega0)*norm(flatDirection.B*correction(:,iColumn))+abs(omega1(iColumn))*norm(flatDirection.B*base(:,iColumn))+norm(first.R*base(:,iColumn));
    bottomResidual(iColumn) = norm(bottom)/max(scale,realmin);
end
result = struct;
result.indices = indices;
result.baseFrequency = omega0;
result.frequencySpread = max(abs(omega-omega0))/max(abs(omega0),realmin);
result.firstOrderFrequency = omega1;
result.baseBasis = base;
result.correction = correction;
result.correctionResidual = max(correctionResidual,[],"all");
result.complementResidual = max(complementResidual,[],"all");
result.firstOrderBottomDefectByDirection = bottomResidual;
result.firstOrderBottomDefect = max(bottomResidual,[],"all");
result.blockHermitianDefect = norm(coupling-coupling',"fro")/max(norm(coupling,"fro"),realmin);
result.fourierSelectionDefect = fourierSelectionDefect(c,base,correction);
end

function result = dressZeroBlock(c,flatDirection,first,flat)
indices = find(abs(flat.frequency) <= 1e3*flat.absoluteFrequencyUncertainty);
X = flat.vectors(:,indices);
mass = X'*flatDirection.E*X;
coupling = X'*(1i*first.J)*X;
[D,omega1] = eig(coupling,mass,"vector");
[omega1,order] = sort(real(omega1));
D = D(:,order);
base = X*D;
for iColumn = 1:size(base,2)
    base(:,iColumn) = base(:,iColumn)/sqrt(real(base(:,iColumn)'*flatDirection.E*base(:,iColumn)));
end
couplingResidual = coupling*D-mass*D.*omega1.';
absoluteUncertainty = vecnorm(couplingResidual,2,1).'./max(vecnorm(mass*D,2,1).',realmin);
absoluteUncertainty = max(absoluteUncertainty,10*eps(max(norm(coupling,2),realmin)));
active = find(abs(omega1) > 1e3*absoluteUncertainty);
stationary = setdiff((1:numel(omega1))',active,"stable");

complement = setdiff((1:numel(flat.frequency))',indices,"stable");
V = flat.vectors(:,complement);
denominator = flat.frequency(complement);
correction = zeros(size(base));
correctionResidual = zeros(size(base,2),1);
bottomResidual = zeros(size(base,2),1);
for iColumn = 1:size(base,2)
    right = -(1i*first.J-omega1(iColumn)*flatDirection.E)*base(:,iColumn);
    coefficient = (V'*right)./denominator;
    correction(:,iColumn) = V*coefficient;
    residual = 1i*flatDirection.J*correction(:,iColumn)-right;
    correctionResidual(iColumn) = norm(residual)/max(norm(right),realmin);
    bottom = -1i*omega1(iColumn)*flatDirection.B*base(:,iColumn)-first.R*base(:,iColumn);
    scale = abs(omega1(iColumn))*norm(flatDirection.B*base(:,iColumn))+norm(first.R*base(:,iColumn));
    bottomResidual(iColumn) = norm(bottom)/max(scale,realmin);
end
dressed = base+correction;
result = struct;
result.indices = indices;
result.baseBasis = base;
result.correction = correction;
result.dressedBasisAtUnitScale = dressed;
result.firstOrderFrequency = omega1;
result.absoluteFrequencyUncertainty = absoluteUncertainty;
result.frequencyResolutionRatio = abs(omega1)./max(absoluteUncertainty,realmin);
result.activeIndices = active;
result.stationaryIndices = stationary;
result.numberOfActiveDirections = numel(active);
result.numberOfStationaryDirections = numel(stationary);
result.correctionResidual = max(correctionResidual,[],"all");
result.firstOrderBottomDefectByDirection = bottomResidual;
if isempty(active)
    result.firstOrderBottomDefect = 0;
else
    result.firstOrderBottomDefect = max(bottomResidual(active),[],"all");
end
if isempty(stationary)
    result.stationaryBottomDiagnostic = 0;
else
    result.stationaryBottomDiagnostic = max(bottomResidual(stationary),[],"all");
end
result.blockHermitianDefect = norm(coupling-coupling',"fro")/max(norm(coupling,"fro"),realmin);
result.fourierSelectionDefect = fourierSelectionDefect(c,base,correction);
result.conjugacyDefect = conjugacyDefect(c,dressed,flatDirection.E);
end

function exact = finiteAmplitudeAudits(c,primitive,internal,zero,trustedBounds,stationaryDegree,degree,terrainScales,targetBlocks)
byScale = repmat(emptyFiniteAudit,numel(terrainScales),1);
for iScale = 1:numel(terrainScales)
    delta = terrainScales(iScale);
    direction = primitive.evaluatedDirections(iScale);
    stationary = constructCompleteStationarySpace(c,direction,trustedBounds,degree,stationaryDegree,delta);
    X = scaledDressedColumns(internal,zero,delta);
    X = projectOutStationary(X,stationary.trustedBasis,direction.E);
    [X,rankDefect] = energyOrthonormalize(X,direction.E);
    Hr = X'*direction.E*X;
    Jr = X'*direction.J*X;
    [C,frequency] = eig(1i*Jr,Hr,"vector");
    vectors = X*C;
    for iMode = 1:size(vectors,2)
        vectors(:,iMode) = vectors(:,iMode)/sqrt(real(vectors(:,iMode)'*direction.E*vectors(:,iMode)));
    end
    frequency = real(frequency);
    target = cellfun(@(block)block.basis,targetBlocks,"UniformOutput",false);
    target = horzcat(target{:});
    projectorDefect = subspaceContainmentDefect(target,X,direction.E);
    X0 = undressedColumns(internal,zero);
    X0 = projectOutStationary(X0,stationary.trustedBasis,direction.E);
    X0 = energyOrthonormalize(X0,direction.E);
    undressedProjectorDefect = subspaceContainmentDefect(target,X0,direction.E);
    internalIndices = selectMostCapturedModes(vectors,target,direction.E,size(target,2));
    apv = trustedAPVDefect(c,direction,vectors(:,internalIndices));
    bottom = modalBottomDefect(direction,vectors(:,internalIndices),frequency(internalIndices));
    strong = strongModeDiagnostics(c,direction,vectors(:,internalIndices),frequency(internalIndices));
    stationaryProjector = energyProjector(stationary.trustedBasis,direction.E);
    dressedProjector = energyProjector(X,direction.E);
    remainderProjector = eye(size(direction.E))-stationaryProjector-dressedProjector;
    decomposition = decompositionDiagnostics(stationaryProjector,dressedProjector,remainderProjector,direction.E,c.coordinateConjugateIndex);
    dimension = struct( ...
        "ambient",size(direction.E,1), ...
        "stationary",size(stationary.trustedBasis,2), ...
        "dressed",size(X,2), ...
        "unresolved",size(direction.E,1)-size(stationary.trustedBasis,2)-size(X,2), ...
        "compressionFactor",size(direction.E,1)/max(size(stationary.trustedBasis,2)+size(X,2),1));
    result = emptyFiniteAudit;
    result.terrainScale = delta;
    result.energyHermitianDefect = norm(Hr-Hr',"fro")/max(norm(Hr,"fro"),realmin);
    result.exchangeSkewHermitianDefect = norm(Jr+Jr',"fro")/max(norm(Jr,"fro"),realmin);
    result.rankDefect = rankDefect;
    result.internalProjectorDefect = projectorDefect;
    result.undressedInternalProjectorDefect = undressedProjectorDefect;
    result.internalProjectorImprovement = undressedProjectorDefect/max(projectorDefect,realmin);
    result.maximumInternalAPVDefect = apv;
    result.maximumInternalBottomDefect = bottom;
    result.maximumInternalStrongResidual = strong.maximumResidual;
    result.internalFrequency = frequency(internalIndices);
    result.internalModeIndices = internalIndices;
    result.reducedBasis = X;
    result.reducedModeVectors = vectors;
    result.frequency = frequency;
    result.stationaryBasis = stationary.trustedBasis;
    result.stationaryProjector = stationaryProjector;
    result.dressedProjector = dressedProjector;
    result.unresolvedProjector = remainderProjector;
    result.decomposition = decomposition;
    result.dimension = dimension;
    result.energyMatrix = direction.E;
    byScale(iScale) = result;
end
exact = struct("byScale",byScale,"fullScale",byScale(end));
end

function X = undressedColumns(internal,zero)
blocks = cell(numel(internal),1);
for iBlock = 1:numel(internal)
    blocks{iBlock} = internal(iBlock).baseBasis;
end
X = horzcat(blocks{:});
if ~isempty(zero.activeIndices)
    X = [X zero.baseBasis(:,zero.activeIndices)];
end
end

function X = scaledDressedColumns(internal,zero,delta)
blocks = cell(numel(internal),1);
for iBlock = 1:numel(internal)
    blocks{iBlock} = internal(iBlock).baseBasis+delta*internal(iBlock).correction;
end
X = horzcat(blocks{:});
if ~isempty(zero.activeIndices)
    X = [X zero.baseBasis(:,zero.activeIndices)+delta*zero.correction(:,zero.activeIndices)];
end
end

function scaling = scalingDiagnostics(primitive,internal,zero,terrainScales)
baseBlocks = cell(numel(internal),1);
correctionBlocks = cell(numel(internal),1);
omega0Blocks = cell(numel(internal),1);
omega1Blocks = cell(numel(internal),1);
for iBlock = 1:numel(internal)
    pair = internal(iBlock);
    baseBlocks{iBlock} = pair.baseBasis;
    correctionBlocks{iBlock} = pair.correction;
    omega0Blocks{iBlock} = [repmat(pair.positive.baseFrequency,1,size(pair.positive.baseBasis,2)) repmat(pair.negative.baseFrequency,1,size(pair.negative.baseBasis,2))];
    omega1Blocks{iBlock} = [pair.positive.firstOrderFrequency.' pair.negative.firstOrderFrequency.'];
end
base = horzcat(baseBlocks{:});
correction = horzcat(correctionBlocks{:});
omega0 = horzcat(omega0Blocks{:});
omega1 = horzcat(omega1Blocks{:});
if ~isempty(zero.activeIndices)
    base = [base zero.baseBasis(:,zero.activeIndices)];
    correction = [correction zero.correction(:,zero.activeIndices)];
    omega0 = [omega0 zeros(1,numel(zero.activeIndices))];
    omega1 = [omega1 zero.firstOrderFrequency(zero.activeIndices).'];
end
omega0 = omega0(:).';
omega1 = omega1(:).';
undressedWeak = zeros(size(terrainScales));
dressedWeak = zeros(size(terrainScales));
undressedBottom = zeros(size(terrainScales));
dressedBottom = zeros(size(terrainScales));
for iScale = 1:numel(terrainScales)
    delta = terrainScales(iScale);
    direction = primitive.evaluatedDirections(iScale);
    dressed = base+delta*correction;
    frequency = omega0+delta*omega1;
    undressedWeak(iScale) = norm(1i*direction.J*base-direction.E*base.*omega0,"fro");
    dressedWeak(iScale) = norm(1i*direction.J*dressed-direction.E*dressed.*frequency,"fro");
    undressedBottom(iScale) = norm(-1i*(direction.B*base).*omega0-direction.R*base,"fro");
    dressedBottom(iScale) = norm(-1i*(direction.B*dressed).*frequency-direction.R*dressed,"fro");
end
nFit = min(3,numel(terrainScales));
indices = 1:nFit;
weakOrder = fittedOrder(terrainScales(indices),dressedWeak(indices));
bottomOrder = fittedOrder(terrainScales(indices),dressedBottom(indices));
undressedWeakOrder = fittedOrder(terrainScales(indices),undressedWeak(indices));
undressedBottomOrder = fittedOrder(terrainScales(indices),undressedBottom(indices));
scaling = struct( ...
    "terrainScales",terrainScales, ...
    "undressedWeakResidual",undressedWeak, ...
    "dressedWeakResidual",dressedWeak, ...
    "undressedBottomResidual",undressedBottom, ...
    "dressedBottomResidual",dressedBottom, ...
    "undressedWeakOrder",undressedWeakOrder, ...
    "undressedBottomOrder",undressedBottomOrder, ...
    "minimumDressedWeakOrder",weakOrder, ...
    "minimumDressedBottomOrder",bottomOrder, ...
    "weakImprovement",undressedWeak(1)/max(dressedWeak(1),realmin), ...
    "bottomImprovement",undressedBottom(1)/max(dressedBottom(1),realmin));
end

function order = fittedOrder(scales,residual)
valid = isfinite(residual) & residual > 100*realmin;
if nnz(valid) < 2
    order = Inf;
else
    coefficient = polyfit(log(scales(valid)),log(residual(valid)),1);
    order = coefficient(1);
end
end

function diagnostics = convergenceDiagnostics(details,comparison)
n = size(details,1);
projector = zeros(n+1,1);
bottom = zeros(n+1,1);
for iDegree = 1:n
    projector(iDegree) = details{iDegree,1}.finite.internalProjectorDefect;
    bottom(iDegree) = details{iDegree,1}.finite.maximumInternalBottomDefect;
end
projector(end) = comparison{1}.finite.internalProjectorDefect;
bottom(end) = comparison{1}.finite.maximumInternalBottomDefect;
diagnostics = struct( ...
    "projectorDefect",projector, ...
    "bottomDefect",bottom, ...
    "projectorIsMonotone",all(diff(projector) <= 100*eps(max(projector,[],"all"))), ...
    "bottomIsMonotone",all(diff(bottom) <= 100*eps(max(bottom,[],"all"))));
end

function diagnostics = paddingDiagnostics(comparison)
reference = comparison{1}.finite;
projector = zeros(numel(comparison)-1,1);
frequency = zeros(numel(comparison)-1,1);
for iPadding = 2:numel(comparison)
    current = comparison{iPadding}.finite;
    projector(iPadding-1) = maximumPrincipalSine(reference.reducedBasis,current.reducedBasis,reference.energyMatrix);
    if numel(current.frequency) == numel(reference.frequency)
        frequency(iPadding-1) = max(abs(sort(current.frequency)-sort(reference.frequency)))/max(max(abs(reference.frequency)),realmin);
    else
        frequency(iPadding-1) = Inf;
    end
end
diagnostics = struct("projectorDefect",max(projector,[],"all"), ...
    "frequencyDefect",max(frequency,[],"all"));
end

function diagnostics = weakSequenceDiagnostics(audits)
green = cellfun(@(value)value.finest.trustedTangentGreenIdentityDefect,audits);
stationary = cellfun(@(value)value.finest.trustedTangentStationaryRowDefect,audits);
diagnostics = struct( ...
    "greenIdentityDefect",green, ...
    "stationaryRowDefect",stationary, ...
    "maximumGreenIdentityDefect",max(green,[],"all"), ...
    "maximumStationaryRowDefect",max(stationary,[],"all"));
end

function defect = bestCenteredAgreement(diagnostics)
names = fieldnames(diagnostics.relativeDefects);
defect = 0;
for iName = 1:numel(names)
    values = diagnostics.relativeDefects.(names{iName});
    defect = max(defect,min(values,[],"all"));
end
end

function defect = requiredCenteredAgreement(diagnostics)
names = ["E","J","R"];
defect = 0;
for name = names
    values = diagnostics.relativeDefects.(name);
    defect = max(defect,min(values,[],"all"));
end
end

function flat = flatEigensystem(direction)
H = direction.E;
J = direction.J;
[R,flag] = chol(H);
if flag ~= 0
    error("WVTerrainEnergyGalerkin:GlobalDressingEnergyNotPositive", ...
        "The flat primitive energy matrix is not positive definite.")
end
K = R'\(1i*J/R);
[U,T] = schur(K,"complex");
frequency = real(diag(T));
[frequency,order] = sort(frequency);
U = U(:,order);
vectors = R\U;
for iMode = 1:size(vectors,2)
    vectors(:,iMode) = vectors(:,iMode)/sqrt(real(vectors(:,iMode)'*H*vectors(:,iMode)));
end
residual = K*(R*vectors)-(R*vectors).*frequency.';
operatorNorm = norm(K,2);
absoluteUncertainty = max(vecnorm(residual,2,1).'./max(vecnorm(R*vectors,2,1).',realmin),10*eps(max(operatorNorm,realmin)));
flat = struct;
flat.frequency = frequency;
flat.vectors = vectors;
flat.absoluteFrequencyUncertainty = absoluteUncertainty;
flat.frequencyResolutionRatio = abs(frequency)./max(absoluteUncertainty,realmin);
flat.energyOrthogonalityDefect = norm(vectors'*H*vectors-eye(size(vectors,2)),"fro")/sqrt(size(vectors,2));
flat.operatorNorm = operatorNorm;
end

function blocks = conjugatePairBlocks(flat)
positive = find(flat.frequency > 1e3*flat.absoluteFrequencyUncertainty);
clusters = frequencyClusters(positive,flat.frequency,flat.absoluteFrequencyUncertainty);
blocks = repmat(struct("indices",zeros(0,1),"positiveIndices",zeros(0,1),"negativeIndices",zeros(0,1),"basis",zeros(size(flat.vectors,1),0)),numel(clusters),1);
for iCluster = 1:numel(clusters)
    positiveIndices = clusters{iCluster};
    negativeCandidates = find(flat.frequency < -1e3*flat.absoluteFrequencyUncertainty);
    target = -mean(flat.frequency(positiveIndices));
    tolerance = 1e3*(max(flat.absoluteFrequencyUncertainty(positiveIndices))+flat.absoluteFrequencyUncertainty(negativeCandidates));
    negativeIndices = negativeCandidates(abs(flat.frequency(negativeCandidates)-target) <= tolerance);
    if numel(negativeIndices) ~= numel(positiveIndices)
        [~,order] = sort(abs(flat.frequency(negativeCandidates)-target));
        negativeIndices = negativeCandidates(order(1:min(numel(positiveIndices),numel(order))));
    end
    indices = sort([positiveIndices;negativeIndices]);
    blocks(iCluster).indices = indices;
    blocks(iCluster).positiveIndices = positiveIndices;
    blocks(iCluster).negativeIndices = negativeIndices;
    blocks(iCluster).basis = flat.vectors(:,indices);
end
end

function clusters = frequencyClusters(indices,frequency,uncertainty)
unused = indices(:);
clusters = cell(0,1);
while ~isempty(unused)
    cluster = unused(1);
    changed = true;
    while changed
        changed = false;
        remaining = setdiff(unused,cluster,"stable");
        for candidate = remaining.'
            separation = abs(frequency(candidate)-frequency(cluster));
            tolerance = 1e3*(uncertainty(candidate)+uncertainty(cluster));
            if any(separation <= tolerance)
                cluster(end+1,1) = candidate; %#ok<AGROW>
                changed = true;
            end
        end
    end
    clusters{end+1,1} = sort(cluster); %#ok<AGROW>
    unused = setdiff(unused,cluster,"stable");
end
end

function [index,defect] = matchFlatPair(target,blocks,H)
defects = Inf(numel(blocks),1);
for iBlock = 1:numel(blocks)
    if size(blocks(iBlock).basis,2) == size(target.basis,2)
        defects(iBlock) = maximumPrincipalSine(blocks(iBlock).basis,target.basis,H);
    end
end
[defect,index] = min(defects);
if ~isfinite(defect)
    index = [];
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
[R,flag] = chol((H+H')/2);
if flag ~= 0
    error("WVTerrainEnergyGalerkin:GlobalDressingComparisonEnergyNotPositive", ...
        "The projector comparison requires a positive energy matrix.")
end
firstQ = orth(R*first);
secondQ = orth(R*second);
singularValues = svd(firstQ'*secondQ);
defect = sqrt(max(0,1-min(singularValues,[],"all")^2));
end

function defect = subspaceContainmentDefect(target,span,H)
if isempty(target)
    defect = 0;
    return
elseif isempty(span) || size(span,2) < size(target,2)
    defect = Inf;
    return
end
[R,flag] = chol((H+H')/2);
if flag ~= 0
    error("WVTerrainEnergyGalerkin:GlobalDressingContainmentEnergyNotPositive", ...
        "The subspace-containment comparison requires a positive energy matrix.")
end
targetQ = orth(R*target);
spanQ = orth(R*span);
defect = norm(targetQ-spanQ*(spanQ'*targetQ),2);
end

function defect = fourierSelectionDefect(c,base,correction)
source = false(c.nK,1);
for iK = 1:c.nK
    rows = c.layout.admissibleRanges{iK};
    source(iK) = norm(base(rows,:),"fro") > 1e-10*max(norm(base,"fro"),realmin);
end
allowed = source;
terrainModes = c.projectionDiagnostics.terrainModes;
for iSource = find(source).'
    destinations = [c.horizontalLayout.kMode(iSource)+terrainModes(:,1) c.horizontalLayout.lMode(iSource)+terrainModes(:,2)];
    for iDestination = 1:size(destinations,1)
        allowed = allowed | (c.horizontalLayout.kMode == destinations(iDestination,1) & c.horizontalLayout.lMode == destinations(iDestination,2));
    end
end
discarded = zeros(0,1);
for iK = find(~allowed).'
    discarded = [discarded c.layout.admissibleRanges{iK}]; %#ok<AGROW>
end
defect = norm(correction(discarded,:),"fro")/max(norm(correction,"fro"),realmin);
end

function defect = conjugacyDefect(c,basis,H)
if isempty(basis)
    defect = 0;
    return
end
n = size(H,1);
conjugacy = sparse((1:n)',c.coordinateConjugateIndex,1,n,n);
projector = energyProjector(basis,H);
residual = projector*conjugacy-conjugacy*conj(projector);
defect = norm(residual,"fro")/max(2*norm(projector,"fro"),realmin);
end

function projector = energyProjector(basis,H)
if isempty(basis)
    projector = zeros(size(H));
else
    gram = basis'*H*basis;
    projector = basis*(gram\(basis'*H));
end
end

function X = projectOutStationary(X,G,H)
if isempty(G)
    return
end
X = X-G*((G'*H*G)\(G'*H*X));
end

function [X,defect] = energyOrthonormalize(X,H)
[R,flag] = chol(H);
if flag ~= 0
    error("WVTerrainEnergyGalerkin:GlobalDressingFiniteEnergyNotPositive", ...
        "The finite-amplitude primitive energy matrix is not positive definite.")
end
[Q,T] = qr(R*X,0);
rankTolerance = max(size(T))*eps(max(norm(T,2),1));
if any(abs(diag(T)) <= rankTolerance)
    defect = nnz(abs(diag(T)) <= rankTolerance);
    return
end
X = R\Q;
defect = 0;
end

function indices = selectMostCapturedModes(vectors,target,H,count)
projector = energyProjector(target,H);
participation = real(sum(conj(vectors).*(H*projector*vectors),1)).';
[~,order] = maxk(participation,min(count,numel(participation)));
indices = sort(order,1);
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

function diagnostics = strongModeDiagnostics(c,direction,C,frequency)
if isempty(C)
    diagnostics = struct("maximumResidual",0,"maximumByMode",zeros(0,1));
    return
end
rawState = c.N*C;
pressure = direction.pressureMap*C;
stateTendency = -1i*rawState.*frequency.';
descriptorResidual = direction.matrices.S*[stateTendency;pressure]-direction.matrices.F*C;
momentumRows = 1:c.nX;
constraintRows = c.nX+(1:size(c.continuity,1));
descriptorTerms = {direction.matrices.S*[stateTendency;pressure],-direction.matrices.F*C};
descriptorDefect = columnRelativeResidual(descriptorResidual,descriptorTerms);
momentumDefect = columnRelativeResidual(descriptorResidual(momentumRows,:),cellfun(@(value)value(momentumRows,:),descriptorTerms,"UniformOutput",false));
continuityDefect = columnRelativeResidual(descriptorResidual(constraintRows,:),cellfun(@(value)value(constraintRows,:),descriptorTerms,"UniformOutput",false));
primitiveResidual = direction.L*C+1i*C.*frequency.';
primitiveDefect = columnRelativeResidual(primitiveResidual,{direction.L*C,1i*C.*frequency.'});
gaugeResidual = c.pressureGauge'*[zeros(c.nX,size(C,2));pressure];
gaugeDefect = vecnorm(gaugeResidual,2,1).'./max(vecnorm(pressure,2,1).',realmin);
pointwise = pointwiseStrongResiduals(c,direction,rawState,pressure,frequency);
maximum = max([descriptorDefect momentumDefect continuityDefect primitiveDefect gaugeDefect pointwise.maximumByMode],[],2);
diagnostics = struct("maximumResidual",max(maximum,[],"all"),"maximumByMode",maximum);
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
diagnostics = struct("maximumByMode",maximumByMode,"maximumResidual",max(maximumByMode,[],"all"));
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

function diagnostics = decompositionDiagnostics(stationary,dressed,remainder,H,conjugateIndex)
identity = eye(size(H));
orthogonality = [norm(stationary*dressed,"fro") norm(stationary*remainder,"fro") norm(dressed*remainder,"fro")];
completeness = norm(stationary+dressed+remainder-identity,"fro")/sqrt(size(identity,1));
n = size(H,1);
conjugacy = sparse((1:n)',conjugateIndex,1,n,n);
projectors = {stationary,dressed,remainder};
conjugacyDefect = 0;
for iProjector = 1:numel(projectors)
    residual = projectors{iProjector}*conjugacy-conjugacy*conj(projectors{iProjector});
    conjugacyDefect = max(conjugacyDefect,norm(residual,"fro")/max(2*norm(projectors{iProjector},"fro"),realmin));
end
diagnostics = struct("maximumOrthogonalityDefect",max(orthogonality),"completenessDefect",completeness,"conjugacyDefect",conjugacyDefect);
end

function representation = nestedRepresentation(c,direction)
representation = struct( ...
    "admissibleBasis",c.N, ...
    "admissibleRanges",{c.layout.admissibleRanges}, ...
    "numberOfRawStateCoefficients",c.nX, ...
    "rawStateBlockSize",c.nXBlock, ...
    "numberOfFCoordinates",c.nF, ...
    "numberOfGCoordinates",c.nG, ...
    "numberOfHCoordinates",c.nH, ...
    "horizontalLayout",c.horizontalLayout, ...
    "bottomMap",direction.B);
end

function transfer = nestedTransfer(coarse,fine)
rawMap = rawNestedMap(coarse,fine);
embeddedRaw = rawMap*coarse.admissibleBasis;
transfer = fine.admissibleBasis\embeddedRaw;
end

function embeddedRaw = coarseEmbeddedRaw(coarse,fine)
embeddedRaw = rawNestedMap(coarse,fine)*coarse.admissibleBasis;
end

function rawMap = rawNestedMap(coarse,fine)
rawMap = sparse(fine.numberOfRawStateCoefficients,coarse.numberOfRawStateCoefficients);
for iCoarse = 1:height(coarse.horizontalLayout)
    mode = [coarse.horizontalLayout.kMode(iCoarse) coarse.horizontalLayout.lMode(iCoarse)];
    iFine = find(fine.horizontalLayout.kMode == mode(1) & fine.horizontalLayout.lMode == mode(2),1);
    if isempty(iFine)
        error("WVTerrainEnergyGalerkin:NonNestedGlobalDressingSupport", ...
            "Every coarse signed Fourier mode must be retained by the comparison support.")
    end
    coarseBlock = (iCoarse-1)*coarse.rawStateBlockSize;
    fineBlock = (iFine-1)*fine.rawStateBlockSize;
    rawMap = insertIdentity(rawMap,fineBlock+(1:coarse.numberOfFCoordinates),coarseBlock+(1:coarse.numberOfFCoordinates));
    rawMap = insertIdentity(rawMap,fineBlock+fine.numberOfFCoordinates+(1:coarse.numberOfFCoordinates),coarseBlock+coarse.numberOfFCoordinates+(1:coarse.numberOfFCoordinates));
    rawMap = insertIdentity(rawMap,fineBlock+2*fine.numberOfFCoordinates+(1:coarse.numberOfGCoordinates),coarseBlock+2*coarse.numberOfFCoordinates+(1:coarse.numberOfGCoordinates));
    fineEta = fineBlock+2*fine.numberOfFCoordinates+fine.numberOfGCoordinates;
    coarseEta = coarseBlock+2*coarse.numberOfFCoordinates+coarse.numberOfGCoordinates;
    rawMap = insertIdentity(rawMap,fineEta+(1:coarse.numberOfGCoordinates),coarseEta+(1:coarse.numberOfGCoordinates));
    rawMap(fineEta+fine.numberOfHCoordinates,coarseEta+coarse.numberOfHCoordinates) = 1;
end
end

function values = insertIdentity(values,rows,columns)
values(rows,columns) = speye(numel(rows),numel(columns));
end

function order = resolvedQuadratureOrder(degree,requested)
if isempty(requested)
    order = max(2*degree+7,18);
else
    order = requested;
end
end

function layout = retainedLayout(problem,bounds)
source = problem.horizontalLayout;
mask = abs(source.kMode) <= bounds(1) & abs(source.lMode) <= bounds(2);
layout = source(mask,:);
expected = (2*bounds(1)+1)*(2*bounds(2)+1);
if height(layout) ~= expected
    error("WVTerrainEnergyGalerkin:UnavailableGlobalDressingSupport", ...
        "The originating transform does not retain the complete signed rectangle [%d %d].",bounds(1),bounds(2))
end
end

function pair = emptyDressedPair
pair = struct("positive",struct,"negative",struct, ...
    "baseBasis",zeros(0),"correction",zeros(0), ...
    "dressedBasisAtUnitScale",zeros(0), ...
    "matchDefect",NaN,"conjugacyDefect",NaN, ...
    "fourierSelectionDefect",NaN,"correctionResidual",NaN, ...
    "firstOrderBottomDefect",NaN,"dimension",0);
end

function result = emptyFiniteAudit
result = struct( ...
    "terrainScale",NaN, ...
    "energyHermitianDefect",NaN, ...
    "exchangeSkewHermitianDefect",NaN, ...
    "rankDefect",NaN, ...
    "internalProjectorDefect",NaN, ...
    "undressedInternalProjectorDefect",NaN, ...
    "internalProjectorImprovement",NaN, ...
    "maximumInternalAPVDefect",NaN, ...
    "maximumInternalBottomDefect",NaN, ...
    "maximumInternalStrongResidual",NaN, ...
    "internalFrequency",zeros(0,1), ...
    "internalModeIndices",zeros(0,1), ...
    "reducedBasis",zeros(0), ...
    "reducedModeVectors",zeros(0), ...
    "frequency",zeros(0,1), ...
    "energyMatrix",zeros(0), ...
    "stationaryBasis",zeros(0), ...
    "stationaryProjector",zeros(0), ...
    "dressedProjector",zeros(0), ...
    "unresolvedProjector",zeros(0), ...
    "decomposition",struct, ...
    "dimension",struct);
end
