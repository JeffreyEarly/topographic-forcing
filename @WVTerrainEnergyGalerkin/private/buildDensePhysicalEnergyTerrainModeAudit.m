function audit = buildDensePhysicalEnergyTerrainModeAudit(problem, ...
    trustedBounds,supportBounds,stationaryDegree,degrees,paddingFactors, ...
    terrainScales,quadratureOrder,usesFixedStationaryDegree)
% Build the Milestone-9 dense physical-energy terrain-mode oracle.

if nargin < 9
    usesFixedStationaryDegree = false;
end

nRefinement = numel(degrees);
nPadding = numel(paddingFactors);
nScale = numel(terrainScales);
details = cell(nRefinement,nPadding,nScale);
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
        [primitive,context] = buildGlobalSmallTerrainPrimitiveAudit( ...
            problem,degree,order,1e-3,horizontalLayout=layout, ...
            paddingFactor=paddingFactors(iPadding), ...
            trustedModeBounds=trustedBounds,rejectTerrainNyquist=true, ...
            evaluationScales=terrainScales);
        for iScale = 1:nScale
            direction = primitive.evaluatedDirections(iScale);
            trustedStationaryDegree = degree;
            if usesFixedStationaryDegree
                trustedStationaryDegree = stationaryDegree;
            end
            stationary = constructCompleteStationarySpace(context, ...
                direction,trustedBounds,degree,trustedStationaryDegree, ...
                terrainScales(iScale));
            stationary.minimumVerifiedDegree = stationaryDegree;
            result = denseModeAudit(context,direction,stationary, ...
                terrainScales(iScale));
            result.projection = primitive.projection;
            details{iRefinement,iPadding,iScale} = result;
            refinement(end+1,1) = summarize(result,degree,order, ...
                supportBounds(iRefinement,:),paddingFactors(iPadding), ...
                terrainScales(iScale)); %#ok<AGROW>
        end
    end
end

iFullScale = find(terrainScales == 1,1);
fullScale = cell(nRefinement,nPadding);
for iRefinement = 1:nRefinement
    for iPadding = 1:nPadding
        fullScale{iRefinement,iPadding} = ...
            details{iRefinement,iPadding,iFullScale};
    end
end
finest = fullScale{end,1};
padding = paddingDiagnostics(fullScale(end,:));
convergence = convergenceDiagnostics(fullScale(:,1));
amplitude = amplitudeDiagnostics(squeeze(details(end,1,:)));
tolerance = struct( ...
    "projection",1e-12,"structure",1e-12, ...
    "stationary",1e-10,"decomposition",1e-10, ...
    "eigenproblem",1e-10,"energyOrthogonality",1e-10, ...
    "frequencyImaginary",1e-10,"conjugacy",1e-10, ...
    "apv",1e-8,"bottom",1e-8,"strong",1e-8, ...
    "padding",1e-10,"stationaryProjector",1e-9, ...
    "minimumFrequencyResolutionRatio",1e3);

allResults = [details{:}];
structurePasses = max([allResults.projectionAdjointDefect]) ...
    <= tolerance.projection ...
    && max([allResults.exactConvolutionDefect]) <= tolerance.projection ...
    && max([allResults.energyHermitianDefect]) <= tolerance.structure ...
    && max([allResults.exchangeSkewHermitianDefect]) ...
    <= tolerance.structure ...
    && max([allResults.energyAssemblyDefect]) ...
    <= tolerance.structure ...
    && max([allResults.exchangeAssemblyDefect]) ...
    <= tolerance.structure ...
    && min([allResults.energyRcond]) > 1e-12;
modalPasses = finest.stationaryRowDefect <= tolerance.stationary ...
    && finest.stationaryEnergyOrthogonalityDefect ...
    <= tolerance.decomposition ...
    && finest.nonTangentBottomRetentionDefect ...
    <= tolerance.decomposition ...
    && finest.dynamicEnergyOrthogonalityDefect ...
    <= tolerance.energyOrthogonality ...
    && finest.stationaryDynamicOrthogonalityDefect ...
    <= tolerance.decomposition ...
    && finest.completenessDefect <= tolerance.decomposition ...
    && finest.maximumEigenResidual <= tolerance.eigenproblem ...
    && finest.maximumFrequencyImaginaryDefect ...
    <= tolerance.frequencyImaginary ...
    && finest.maximumConjugacyDefect <= tolerance.conjugacy ...
    && padding.maximumFrequencyDefect <= tolerance.padding;
stationarySeparationPasses = ...
    finest.zeroFrequencyStationaryProjectorDefect ...
    <= tolerance.stationaryProjector;
boundaryPasses = finest.numberOfClaimedBoundaryModes ...
    == finest.numberOfBoundaryCandidates ...
    && finest.maximumBoundaryCandidateAPVDefect <= tolerance.apv ...
    && finest.maximumBoundaryCandidateBottomDefect <= tolerance.bottom ...
    && finest.maximumBoundaryCandidateStrongResidual <= tolerance.strong;
frequencyPasses = finest.minimumBoundaryFrequencyResolutionRatio ...
    >= tolerance.minimumFrequencyResolutionRatio;

if ~structurePasses
    status = "implementation-unresolved";
    diagnosis = "The common projection or raw physical-energy structure failed before modal construction.";
elseif ~modalPasses
    status = "dense-terrain-mode-blocker";
    diagnosis = "The tangency-defined stationary complement did not produce a numerically valid dense physical-energy eigensystem.";
elseif ~stationarySeparationPasses && ~boundaryPasses
    status = "dense-modal-classification-blocker";
    diagnosis = "The energy eigensystem is numerically valid, but its backward-error zero-frequency projector is larger than the Milestone-8 tangency sector and its trusted non-tangent bottom candidates fail the independent APV, bottom, or strong-equation gates. Every direction remains retained and unclassified.";
elseif ~stationarySeparationPasses
    status = "stationary-space-modal-blocker";
    diagnosis = "The tangency-defined Milestone-8 stationary sector is numerically valid, but its physical-energy complement contains additional backward-error-indistinguishable zero-frequency directions. They remain retained and unclassified; they are neither added to the stationary sector nor deleted.";
elseif ~boundaryPasses
    status = "dense-boundary-mode-classification-blocker";
    diagnosis = "The dense energy eigenproblem is well resolved, but the trusted non-tangent bottom candidates do not simultaneously satisfy projected volume APV, strong bottom evolution, and the primitive strong-equation gates.";
elseif ~frequencyPasses
    status = "low-frequency-resolution-blocker";
    diagnosis = "The active bottom-mode subspace is present, but at least one claimed nonzero frequency is less than one thousand times its numerical uncertainty.";
else
    status = "dense-physical-energy-terrain-mode-oracle";
    diagnosis = "The tangency-defined stationary sector and its complete physical-energy complement yield resolved internal-wave and active bottom-mode branches.";
end

audit = struct;
audit.scope = "milestone-9-dense-physical-energy-terrain-modes";
audit.status = status;
audit.isCompatible = status == ...
    "dense-physical-energy-terrain-mode-oracle";
audit.diagnosis = diagnosis;
audit.trustedModeBounds = trustedBounds;
audit.supportModeBounds = supportBounds;
audit.stationaryPolynomialDegree = stationaryDegree;
audit.primitivePolynomialDegrees = degrees;
audit.paddingFactors = paddingFactors;
audit.terrainScales = terrainScales;
audit.quadratureOrder = quadratureOrder;
audit.usesFixedStationaryDegree = usesFixedStationaryDegree;
audit.inertialFrequency = problem.originatingTransform.f;
audit.refinement = refinement;
audit.convergence = convergence;
audit.padding = padding;
audit.amplitude = amplitude;
audit.finest = finest;
audit.details = details;
audit.requiredTolerance = tolerance;
audit.nextScope = "milestone-10-only-if-compatible-and-separately-authorized";
end

function result = denseModeAudit(c,direction,stationary,terrainScale)
H = direction.E;
J = direction.J;
[R,flag] = chol(H);
if flag ~= 0
    error("WVTerrainEnergyGalerkin:DenseModeEnergyNotPositive", ...
        "The raw finite-terrain energy matrix is not positive definite.")
end

G = stationary.trustedBasis;
Yg = R*G;
[completeQ,~] = qr(Yg);
Yd = completeQ(:,size(Yg,2)+1:end);
W = R\Yd;
K = R'\(1i*J/R);
Kmodal = Yd'*K*Yd;
[D,schurForm] = schur(Kmodal,"complex");
frequencyComplex = diag(schurForm);
[~,order] = sort(real(frequencyComplex));
frequencyComplex = frequencyComplex(order);
D = D(:,order);
D = normalizeColumns(D);
C = W*D;

frequency = real(columnRayleighQuotients(C,H,1i*J));
[~,order] = sort(frequency);
C = C(:,order);
D = D(:,order);
frequencyComplex = frequencyComplex(order);
for iMode = 1:size(C,2)
    normalization = sqrt(real(C(:,iMode)'*H*C(:,iMode)));
    C(:,iMode) = C(:,iMode)/normalization;
    D(:,iMode) = D(:,iMode)/normalization;
end
frequency = real(columnRayleighQuotients(C,H,1i*J));

energyCoordinates = R*C;
fullEnergyResidual = K*energyCoordinates-energyCoordinates.*frequency.';
operatorNorm = norm(K,2);
absoluteUncertainty = max(vecnorm(fullEnergyResidual,2,1).' ...
    ./max(vecnorm(energyCoordinates,2,1).',realmin), ...
    10*eps(max(operatorNorm,realmin)));
frequencyResolutionRatio = abs(frequency) ...
    ./max(absoluteUncertainty,realmin);
eigenResidual = vecnorm(fullEnergyResidual,2,1).' ...
    ./max((operatorNorm+abs(frequency)) ...
    .*vecnorm(energyCoordinates,2,1).',realmin);

bottomSeed = stationary.trustedBottomSeedBasis;
trustedColumns = c.layout.trustedColumns;
trustedCoordinateBasis = speye(size(H,1));
trustedCoordinateBasis = trustedCoordinateBasis(:,trustedColumns);
trustedProjector = subspaceEnergyProjector(trustedCoordinateBasis,H);
trustedParticipation = real(sum(conj(C).*(H*trustedProjector*C),1)).';
horizontalParticipation = zeros(c.nK,numel(frequency));
for iHorizontal = 1:c.nK
    horizontalBasis = speye(size(H,1));
    horizontalBasis = horizontalBasis(:, ...
        c.layout.admissibleRanges{iHorizontal});
    horizontalProjector = subspaceEnergyProjector(horizontalBasis,H);
    horizontalParticipation(iHorizontal,:) = real(sum( ...
        conj(C).*(H*horizontalProjector*C),1));
end
if isempty(bottomSeed)
    bottomParticipation = zeros(numel(frequency),1);
    bottomSeedParticipation = zeros(numel(frequency),1);
    boundaryCandidateIndices = zeros(0,1);
else
    bottomOperatorNorm = norm(direction.B/R,2);
    bottomParticipation = (vecnorm(direction.B*C,2,1).' ...
        /max(bottomOperatorNorm,realmin)).^2;
    bottomSeedProjector = subspaceEnergyProjector(bottomSeed,H);
    bottomSeedParticipation = real(sum( ...
        conj(C).*(H*bottomSeedProjector*C),1)).';
    boundaryCandidateIndices = conjugateClosedCandidates( ...
        bottomSeedParticipation,frequency,absoluteUncertainty, ...
        size(bottomSeed,2));
end

trustedAPV = trustedAPVDiagnostics(c,direction,C,G);
bottomResidual = -1i*(direction.B*C).*frequency.'-direction.R*C;
bottomScale = abs(frequency)*norm(direction.B/R,2) ...
    +norm(direction.R/R,2);
bottomDefect = vecnorm(bottomResidual,2,1).' ...
    ./max(bottomScale,realmin);
conjugacy = conjugacyDiagnostics(c,direction,C,frequency, ...
    absoluteUncertainty);
strong = strongModeDiagnostics(c,direction,C,frequency);
modePasses = eigenResidual <= 1e-10 ...
    & frequencyResolutionRatio >= 1e3 ...
    & trustedAPV.directDefect <= 1e-8 ...
    & bottomDefect <= 1e-8 ...
    & strong.maximumByMode <= 1e-8;
claimedBoundaryIndices = boundaryCandidateIndices( ...
    modePasses(boundaryCandidateIndices));
internalIndices = find(modePasses);
internalIndices = setdiff(internalIndices,claimedBoundaryIndices,"stable");
unclassifiedIndices = setdiff((1:numel(frequency))', ...
    [internalIndices;claimedBoundaryIndices],"stable");
unresolvedZeroIndices = find(abs(frequency) <= absoluteUncertainty);
if isempty(boundaryCandidateIndices)
    minimumBoundaryResolution = Inf;
    maximumBoundaryAPV = 0;
    maximumBoundaryBottom = 0;
    maximumBoundaryStrong = 0;
else
    minimumBoundaryResolution = min( ...
        frequencyResolutionRatio(boundaryCandidateIndices));
    maximumBoundaryAPV = max( ...
        trustedAPV.directDefect(boundaryCandidateIndices));
    maximumBoundaryBottom = max(bottomDefect(boundaryCandidateIndices));
    maximumBoundaryStrong = max( ...
        strong.maximumByMode(boundaryCandidateIndices));
end

stationaryProjector = G*(G'*H);
dynamicProjector = W*(W'*H);
identity = eye(size(H));
energyStationaryProjector = (R*G)*(R*G)';
if isempty(unresolvedZeroIndices)
    zeroFrequencyStationaryProjectorDefect = 0;
else
    energyZeroProjector = energyStationaryProjector ...
        +(R*C(:,unresolvedZeroIndices)) ...
        *(R*C(:,unresolvedZeroIndices))';
    zeroFrequencyStationaryProjectorDefect = ...
        norm(energyZeroProjector-energyStationaryProjector,"fro") ...
        /max(norm(energyZeroProjector,"fro"),realmin);
end
fullBottomSeed = stationary.bottomSeedBasis;
if isempty(fullBottomSeed)
    nonTangentBottomRetentionDefect = 0;
else
    discardedBottom = fullBottomSeed-dynamicProjector*fullBottomSeed;
    nonTangentBottomRetentionDefect = ...
        sqrt(max(real(trace(discardedBottom'*H*discardedBottom)),0)) ...
        /max(sqrt(max(real(trace( ...
        fullBottomSeed'*H*fullBottomSeed)),0)),realmin);
end
result = struct;
result.terrainScale = terrainScale;
result.energyMatrix = H;
result.exchangeMatrix = J;
result.stationaryBasis = G;
result.dynamicBasis = W;
result.modeVectors = C;
result.frequency = frequency;
result.frequencyComplex = frequencyComplex;
result.absoluteFrequencyUncertainty = absoluteUncertainty;
result.frequencyResolutionRatio = frequencyResolutionRatio;
result.bottomParticipation = bottomParticipation;
result.bottomSeedParticipation = bottomSeedParticipation;
result.trustedParticipation = trustedParticipation;
result.horizontalEnergyParticipation = horizontalParticipation;
result.boundaryCandidateIndices = boundaryCandidateIndices;
result.boundaryModeIndices = claimedBoundaryIndices;
result.internalWaveIndices = internalIndices;
result.unclassifiedModeIndices = unclassifiedIndices;
result.unresolvedZeroIndices = unresolvedZeroIndices;
result.numberOfStationaryModes = size(G,2);
result.numberOfDynamicModes = size(W,2);
result.numberOfRetainedNonTangentBottomDirections = ...
    size(stationary.bottomSeedBasis,2);
result.numberOfBoundaryCandidates = numel(boundaryCandidateIndices);
result.numberOfClaimedBoundaryModes = numel(claimedBoundaryIndices);
result.numberOfBoundaryModes = numel(claimedBoundaryIndices);
result.numberOfInternalWaveModes = numel(internalIndices);
result.numberOfUnclassifiedWeakNullDirections = ...
    numel(unresolvedZeroIndices);
result.numberOfNumericallyStationaryModes = ...
    size(G,2)+numel(unresolvedZeroIndices);
result.zeroFrequencyStationaryExcessDimension = ...
    numel(unresolvedZeroIndices);
result.zeroFrequencyStationaryProjectorDefect = ...
    zeroFrequencyStationaryProjectorDefect;
result.maximumUnresolvedZeroTrustedParticipation = maximumOrZero( ...
    trustedParticipation(unresolvedZeroIndices));
result.minimumBoundaryFrequencyResolutionRatio = ...
    minimumBoundaryResolution;
result.minimumAbsoluteBoundaryFrequency = minimumOrNaN( ...
    abs(frequency(boundaryCandidateIndices)));
result.maximumBoundaryCandidateAPVDefect = maximumBoundaryAPV;
result.maximumBoundaryCandidateBottomDefect = maximumBoundaryBottom;
result.maximumBoundaryCandidateStrongResidual = maximumBoundaryStrong;
result.maximumEigenResidual = max(eigenResidual,[],"all");
result.eigenResidual = eigenResidual;
result.maximumFrequencyImaginaryDefect = max(abs(imag(frequencyComplex))) ...
    /max(operatorNorm,realmin);
result.energyHermitianDefect = norm(H-H',"fro") ...
    /max(norm(H,"fro"),realmin);
result.exchangeSkewHermitianDefect = norm(J+J',"fro") ...
    /max(norm(J,"fro"),realmin);
result.reducedHermitianDefect = norm(Kmodal-Kmodal',"fro") ...
    /max(norm(Kmodal,"fro"),realmin);
result.energyRcond = rcond(H);
result.energyAssemblyDefect = ...
    direction.diagnostics.energyAssemblyDefect;
result.exchangeAssemblyDefect = ...
    direction.diagnostics.exchangeAssemblyDefect;
result.stationaryRowDefect = normalizedProduct(J*G,J,G);
stationaryEnergyResidual = K*(R*G);
result.stationaryEigenResidual = vecnorm( ...
    stationaryEnergyResidual,2,1).' ...
    /max(operatorNorm,realmin);
result.stationaryEnergyOrthogonalityDefect = ...
    norm(G'*H*G-eye(size(G,2)),"fro")/sqrt(max(size(G,2),1));
result.dynamicEnergyOrthogonalityDefect = ...
    norm(C'*H*C-eye(size(C,2)),"fro")/sqrt(max(size(C,2),1));
result.stationaryDynamicOrthogonalityDefect = ...
    norm(G'*H*W,"fro") ...
    /max(norm(G,"fro")*norm(H,"fro")*norm(W,"fro"),realmin);
result.completenessDefect = norm(stationaryProjector+dynamicProjector ...
    -identity,"fro")/sqrt(size(identity,1));
result.nonTangentBottomRetentionDefect = ...
    nonTangentBottomRetentionDefect;
result.maximumTrustedAPVDefect = trustedAPV.maximumDefect;
result.trustedAPV = trustedAPV;
result.maximumBottomDefect = max(bottomDefect,[],"all");
result.bottomDefect = bottomDefect;
result.maximumConjugacyDefect = conjugacy.maximumDefect;
result.conjugacy = conjugacy;
result.strong = strong;
result.modePassesIndependentPhysicalGates = modePasses;
result.stationary = stationary;
result.stationaryEnergyProjector = stationaryProjector;
result.dynamicEnergyProjector = dynamicProjector;
result.boundaryEnergyProjector = subspaceEnergyProjector( ...
    C(:,claimedBoundaryIndices),H);
result.projectionAdjointDefect = c.projectionDiagnostics.adjointDefect;
result.exactConvolutionDefect = ...
    c.projectionDiagnostics.maximumExactConvolutionDefect;
result.horizontalLayout = c.horizontalLayout;
result.coordinateConjugateIndex = c.coordinateConjugateIndex;
result.nestedRepresentation = struct( ...
    "admissibleBasis",c.N, ...
    "admissibleRanges",{c.layout.admissibleRanges}, ...
    "numberOfRawStateCoefficients",c.nX, ...
    "rawStateBlockSize",c.nXBlock, ...
    "numberOfFCoordinates",c.nF, ...
    "numberOfGCoordinates",c.nG, ...
    "numberOfHCoordinates",c.nH, ...
    "numberOfVerticalQuadraturePoints",c.nZ, ...
    "horizontalLayout",c.horizontalLayout, ...
    "bottomMap",direction.B, ...
    "projectedAPVMap",direction.QProjected, ...
    "projectedAPVWeight",c.apvVerticalWeight);
end

function diagnostics = trustedAPVDiagnostics(c,direction,C,G)
trustedHorizontal = c.layout.trustedHorizontalModes;
trustedRows = find(repmat(trustedHorizontal,c.nZ,1));
Q = direction.QProjected(trustedRows,:);
weights = c.apvVerticalWeight(trustedRows);
energyFactor = chol(direction.E);
weightedOperator = sqrt(weights).*Q/energyFactor;
values = Q*C;
weightedValues = sqrt(weights).*values;
scale = max(norm(weightedOperator,2),realmin);
directDefect = vecnorm(weightedValues,2,1).'/scale;
weakDefect = vecnorm(G'*direction.E*C,2,1).' ...
    /max(norm(G'*direction.E/energyFactor,2),realmin);
diagnostics = struct("directDefect",directDefect, ...
    "weakDefect",weakDefect, ...
    "maximumDirectDefect",max(directDefect,[],"all"), ...
    "maximumWeakDefect",max(weakDefect,[],"all"), ...
    "maximumDefect",max([directDefect;weakDefect],[],"all"));
end

function diagnostics = conjugacyDiagnostics(c,direction,C,frequency,uncertainty)
n = size(C,1);
conjugacyMap = sparse((1:n)',c.coordinateConjugateIndex,1,n,n);
frequencyDefect = zeros(numel(frequency),1);
vectorDefect = zeros(numel(frequency),1);
partnerIndex = zeros(numel(frequency),1);
operatorScale = max(max(abs(frequency)),realmin);
for iMode = 1:numel(frequency)
    [~,partner] = min(abs(frequency+frequency(iMode)));
    partnerIndex(iMode) = partner;
    expected = conjugacyMap*conj(C(:,iMode));
    frequencyFloor = uncertainty(iMode)+uncertainty ...
        +10*eps(max(max(abs(frequency)),realmin));
    partnerCluster = find(abs(frequency+frequency(iMode)) ...
        <= frequencyFloor);
    if isempty(partnerCluster)
        partnerCluster = partner;
    end
    overlap = C(:,partnerCluster)'*direction.E*expected;
    vectorDefect(iMode) = sqrt(max(0,1-sum(abs(overlap).^2)));
    frequencyDefect(iMode) = abs(frequency(partner)+frequency(iMode)) ...
        /operatorScale;
end
diagnostics = struct("partnerIndex",partnerIndex, ...
    "frequencyDefect",frequencyDefect,"vectorDefect",vectorDefect, ...
    "maximumFrequencyDefect",max(frequencyDefect,[],"all"), ...
    "maximumVectorDefect",max(vectorDefect,[],"all"), ...
    "operatorDefect",operatorConjugacyDefect(c,direction), ...
    "maximumDefect",max([frequencyDefect; ...
    operatorConjugacyDefect(c,direction)],[],"all"));
end

function defect = operatorConjugacyDefect(c,direction)
n = size(direction.E,1);
conjugacyMap = sparse((1:n)',c.coordinateConjugateIndex,1,n,n);
generator = direction.E\direction.J;
residual = generator*conjugacyMap ...
    -conjugacyMap*conj(generator);
defect = norm(residual,"fro")/max(2*norm(generator,"fro"),realmin);
end

function indices = conjugateClosedCandidates( ...
    participation,frequency,uncertainty,nSeed)
if nSeed == 0
    indices = zeros(0,1);
    return
end
nPositive = floor(nSeed/2);
nNegative = floor(nSeed/2);
positive = find(frequency > 0);
negative = find(frequency < 0);
[~,positiveOrder] = sort(participation(positive),"descend");
[~,negativeOrder] = sort(participation(negative),"descend");
indices = [positive(positiveOrder(1:min(nPositive,numel(positive)))); ...
    negative(negativeOrder(1:min(nNegative,numel(negative))))];
remaining = nSeed-numel(indices);
if remaining > 0
    available = setdiff((1:numel(frequency))',indices,"stable");
    [~,order] = sort(participation(available),"descend");
    indices = [indices;available(order(1:min(remaining,numel(order))))];
end
partners = zeros(numel(indices),1);
for iIndex = 1:numel(indices)
    [~,partners(iIndex)] = min(abs(frequency+frequency(indices(iIndex))));
end
indices = unique([indices;partners]);
clustered = indices;
frequencyFloor = 10*eps(max(max(abs(frequency)),realmin));
for iIndex = 1:numel(indices)
    clustered = [clustered;find(abs(frequency-frequency(indices(iIndex))) ...
        <= uncertainty+uncertainty(indices(iIndex))+frequencyFloor)]; %#ok<AGROW>
    clustered = [clustered;find(abs(frequency+frequency(indices(iIndex))) ...
        <= uncertainty+uncertainty(indices(iIndex))+frequencyFloor)]; %#ok<AGROW>
end
indices = sort(unique(clustered));
end

function diagnostics = strongModeDiagnostics(c,direction,C,frequency)
rawState = c.N*C;
pressure = direction.pressureMap*C;
stateTendency = -1i*rawState.*frequency.';
descriptorResidual = direction.matrices.S*[stateTendency;pressure] ...
    -direction.matrices.F*C;
momentumRows = 1:c.nX;
constraintRows = c.nX+(1:size(c.continuity,1));
descriptorTerms = {direction.matrices.S*[stateTendency;pressure], ...
    -direction.matrices.F*C};
descriptorDefect = columnRelativeResidual( ...
    descriptorResidual,descriptorTerms);
momentumDefect = columnRelativeResidual( ...
    descriptorResidual(momentumRows,:), ...
    cellfun(@(value)value(momentumRows,:),descriptorTerms, ...
    "UniformOutput",false));
continuityDefect = columnRelativeResidual( ...
    descriptorResidual(constraintRows,:), ...
    cellfun(@(value)value(constraintRows,:),descriptorTerms, ...
    "UniformOutput",false));
primitiveResidual = direction.L*C+1i*C.*frequency.';
primitiveDefect = columnRelativeResidual(primitiveResidual, ...
    {direction.L*C,1i*C.*frequency.'});
gaugeResidual = c.pressureGauge'*[zeros(c.nX,size(C,2));pressure];
gaugeDefect = vecnorm(gaugeResidual,2,1).' ...
    ./max(vecnorm(pressure,2,1).',realmin);

pointwise = pointwiseStrongResiduals(c,direction,rawState,pressure,frequency);
maximum = max([descriptorDefect momentumDefect continuityDefect ...
    primitiveDefect gaugeDefect pointwise.maximumByMode],[],2);
diagnostics = struct( ...
    "descriptorDefect",descriptorDefect, ...
    "momentumWeakDefect",momentumDefect, ...
    "continuityWeakDefect",continuityDefect, ...
    "primitiveGeneratorDefect",primitiveDefect, ...
    "pressureGaugeDefect",gaugeDefect, ...
    "pointwise",pointwise, ...
    "maximumByMode",maximum, ...
    "maximumResidual",max(maximum,[],"all"));
end

function diagnostics = pointwiseStrongResiduals(c,direction,rawState, ...
    pressureCoefficients,frequency)
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
uTerms = {-1i*u.*omega,-c.wvt.f*v,Dxp/c.wvt.rho0};
vTerms = {-1i*v.*omega,c.wvt.f*u,Dyp/c.wvt.rho0};
wTerms = {-1i*w.*omega,N2.*eta,pXi./(c.wvt.rho0*gamma)};
etaTerms = {-1i*eta.*omega,-w};
divergence = (c.RuX+c.RvY+c.RwhXi)*rawState;

uDefect = weightedColumnResidual(uTerms,c.volumeWeight);
vDefect = weightedColumnResidual(vTerms,c.volumeWeight);
wDefect = weightedColumnResidual(wTerms,c.volumeWeight);
etaDefect = weightedColumnResidual(etaTerms,c.volumeWeight);
continuityDefect = weightedColumnNorm(divergence,c.volumeWeight) ...
    ./max(weightedColumnNorm(c.RuX*rawState,c.volumeWeight) ...
    +weightedColumnNorm(c.RvY*rawState,c.volumeWeight) ...
    +weightedColumnNorm(c.RwhXi*rawState,c.volumeWeight),realmin);
surfaceEta = surfaceDisplacement(c,rawState);
surfaceDefect = vecnorm(surfaceEta,2,1).' ...
    /max(vecnorm(rawState,2,1).',realmin);
maximumByMode = max([uDefect vDefect wDefect etaDefect ...
    continuityDefect surfaceDefect],[],2);
diagnostics = struct("horizontalUDefect",uDefect, ...
    "horizontalVDefect",vDefect,"verticalDefect",wDefect, ...
    "displacementDefect",etaDefect, ...
    "continuityDefect",continuityDefect, ...
    "surfaceDefect",surfaceDefect, ...
    "maximumByMode",maximumByMode, ...
    "maximumResidual",max(maximumByMode,[],"all"));
end

function values = surfaceDisplacement(c,rawState)
values = zeros(c.nK,size(rawState,2));
for iK = 1:c.nK
    block = (iK-1)*c.nXBlock;
    etaColumns = block+2*c.nF+c.nG+(1:c.nH);
    values(iK,:) = c.spaces.Hendpoint(2,:)*rawState(etaColumns,:);
end
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

function values = columnRayleighQuotients(C,H,A)
values = zeros(size(C,2),1);
for iColumn = 1:size(C,2)
    values(iColumn) = (C(:,iColumn)'*A*C(:,iColumn)) ...
        /(C(:,iColumn)'*H*C(:,iColumn));
end
end

function values = normalizeColumns(values)
values = values./max(vecnorm(values,2,1),realmin);
end

function defect = columnRelativeResidual(residual,terms)
scale = zeros(size(residual,2),1);
for iTerm = 1:numel(terms)
    scale = scale+vecnorm(terms{iTerm},2,1).';
end
defect = vecnorm(residual,2,1).'./max(scale,realmin);
end

function diagnostics = convergenceDiagnostics(results)
fields = ["maximumEigenResidual","maximumTrustedAPVDefect", ...
    "maximumBottomDefect","maximumConjugacyDefect", ...
    "maximumBoundaryCandidateAPVDefect", ...
    "maximumBoundaryCandidateBottomDefect", ...
    "maximumBoundaryCandidateStrongResidual", ...
    "minimumAbsoluteBoundaryFrequency", ...
    "minimumBoundaryFrequencyResolutionRatio"];
diagnostics = struct;
for field = fields
    diagnostics.(field) = reshape(cellfun( ...
        @(result)result.(field),results),1,[]);
end
diagnostics.numberOfStationaryModes = cellfun( ...
    @(result)result.numberOfStationaryModes,results);
diagnostics.numberOfDynamicModes = cellfun( ...
    @(result)result.numberOfDynamicModes,results);
diagnostics.numberOfBoundaryModes = cellfun( ...
    @(result)result.numberOfBoundaryModes,results);
diagnostics.numberOfBoundaryCandidates = cellfun( ...
    @(result)result.numberOfBoundaryCandidates,results);
diagnostics.numberOfClaimedBoundaryModes = cellfun( ...
    @(result)result.numberOfClaimedBoundaryModes,results);
diagnostics.numberOfRetainedNonTangentBottomDirections = cellfun( ...
    @(result)result.numberOfRetainedNonTangentBottomDirections,results);
diagnostics.numberOfUnclassifiedWeakNullDirections = cellfun( ...
    @(result)result.numberOfUnclassifiedWeakNullDirections,results);
diagnostics.numberOfNumericallyStationaryModes = cellfun( ...
    @(result)result.numberOfNumericallyStationaryModes,results);
diagnostics.nonTangentBottomRetentionDefect = cellfun( ...
    @(result)result.nonTangentBottomRetentionDefect,results);
diagnostics.zeroFrequencyStationaryProjectorDefect = cellfun( ...
    @(result)result.zeroFrequencyStationaryProjectorDefect,results);
end

function diagnostics = paddingDiagnostics(results)
reference = results{1};
frequencyDefect = zeros(numel(results),1);
for iResult = 1:numel(results)
    current = results{iResult};
    if numel(current.frequency) ~= numel(reference.frequency)
        frequencyDefect(iResult) = Inf;
    else
        scale = max(max(abs(reference.frequency)),realmin);
        frequencyDefect(iResult) = ...
            max(abs(current.frequency-reference.frequency))/scale;
    end
end
diagnostics = struct("frequencyDefect",frequencyDefect, ...
    "maximumFrequencyDefect",max(frequencyDefect));
end

function diagnostics = amplitudeDiagnostics(results)
minimumFrequency = cellfun( ...
    @(result)result.minimumAbsoluteBoundaryFrequency,results);
resolutionRatio = cellfun( ...
    @(result)result.minimumBoundaryFrequencyResolutionRatio,results);
diagnostics = struct("terrainScale",cellfun( ...
    @(result)result.terrainScale,results), ...
    "minimumAbsoluteBoundaryFrequency",minimumFrequency, ...
    "minimumBoundaryFrequencyResolutionRatio",resolutionRatio);
end

function summary = summarize(result,degree,order,support,padding,scale)
summary = emptySummary;
summary.degree = degree;
summary.quadratureOrder = order;
summary.supportKMax = support(1);
summary.supportLMax = support(2);
summary.paddingFactor = padding;
summary.terrainScale = scale;
summary.numberOfStationaryModes = result.numberOfStationaryModes;
summary.numberOfDynamicModes = result.numberOfDynamicModes;
summary.numberOfBoundaryModes = result.numberOfBoundaryModes;
summary.numberOfBoundaryCandidates = result.numberOfBoundaryCandidates;
summary.numberOfClaimedBoundaryModes = result.numberOfClaimedBoundaryModes;
summary.numberOfRetainedNonTangentBottomDirections = ...
    result.numberOfRetainedNonTangentBottomDirections;
summary.numberOfUnclassifiedWeakNullDirections = ...
    result.numberOfUnclassifiedWeakNullDirections;
summary.numberOfNumericallyStationaryModes = ...
    result.numberOfNumericallyStationaryModes;
summary.nonTangentBottomRetentionDefect = ...
    result.nonTangentBottomRetentionDefect;
summary.zeroFrequencyStationaryProjectorDefect = ...
    result.zeroFrequencyStationaryProjectorDefect;
summary.maximumUnresolvedZeroTrustedParticipation = ...
    result.maximumUnresolvedZeroTrustedParticipation;
summary.minimumAbsoluteBoundaryFrequency = ...
    result.minimumAbsoluteBoundaryFrequency;
summary.minimumBoundaryFrequencyResolutionRatio = ...
    result.minimumBoundaryFrequencyResolutionRatio;
summary.maximumEigenResidual = result.maximumEigenResidual;
summary.maximumTrustedAPVDefect = result.maximumTrustedAPVDefect;
summary.maximumBottomDefect = result.maximumBottomDefect;
summary.maximumConjugacyDefect = result.maximumConjugacyDefect;
summary.maximumStrongResidual = result.strong.maximumResidual;
summary.maximumBoundaryCandidateAPVDefect = ...
    result.maximumBoundaryCandidateAPVDefect;
summary.maximumBoundaryCandidateBottomDefect = ...
    result.maximumBoundaryCandidateBottomDefect;
summary.maximumBoundaryCandidateStrongResidual = ...
    result.maximumBoundaryCandidateStrongResidual;
end

function value = subspaceEnergyProjector(C,H)
if isempty(C)
    value = zeros(size(H));
else
    gram = C'*H*C;
    value = C*(gram\(C'*H));
end
end

function value = normalizedProduct(residual,left,right)
value = norm(residual,"fro") ...
    /max(norm(left,"fro")*norm(right,"fro"),realmin);
end

function value = minimumOrNaN(values)
if isempty(values)
    value = NaN;
else
    value = min(values);
end
end

function value = maximumOrZero(values)
if isempty(values)
    value = 0;
else
    value = max(values);
end
end

function layout = retainedLayout(problem,bounds)
source = problem.horizontalLayout;
mask = abs(source.kMode) <= bounds(1) & abs(source.lMode) <= bounds(2);
layout = source(mask,:);
expected = (2*bounds(1)+1)*(2*bounds(2)+1);
if height(layout) ~= expected
    error("WVTerrainEnergyGalerkin:UnavailableDenseModeSupport", ...
        "The originating transform does not retain the complete signed rectangle [%d %d].", ...
        bounds(1),bounds(2))
end
end

function terrain = terrainSupport(problem)
wvt = problem.originatingTransform;
spectrum = fft2(problem.topographicHeight)/(wvt.Nx*wvt.Ny);
nyquist = logical(WVGeometryDoublyPeriodic.maskForNyquistModes( ...
    wvt.Nx,wvt.Ny));
scale = max(abs(spectrum),[],"all");
tolerance = 100*eps*max(scale,1);
if max(abs(spectrum(nyquist)),[],"all") > tolerance
    error("WVTerrainEnergyGalerkin:DenseModeTerrainNyquist", ...
        "The dense terrain-mode oracle requires terrain with no material Nyquist coefficient.")
end
[kMode,lMode] = ndgrid(wvt.kMode_dft,wvt.lMode_dft);
active = abs(spectrum) > tolerance & ~nyquist;
terrain = struct("kMode",kMode(active),"lMode",lMode(active), ...
    "coefficient",spectrum(active),"spectralTolerance",tolerance);
end

function validateGuard(layout,trustedBounds,terrain)
for iK = find(abs(layout.kMode) <= trustedBounds(1) ...
        & abs(layout.lMode) <= trustedBounds(2)).'
    destinations = [layout.kMode(iK)+terrain.kMode ...
        layout.lMode(iK)+terrain.lMode];
    for iDestination = 1:size(destinations,1)
        if ~any(layout.kMode == destinations(iDestination,1) ...
                & layout.lMode == destinations(iDestination,2))
            error("WVTerrainEnergyGalerkin:InsufficientDenseModeGuard", ...
                "supportModeBounds must retain every first terrain sideband of the trusted band.")
        end
    end
end
end

function summary = emptySummary
summary = struct( ...
    "degree",0,"quadratureOrder",0,"supportKMax",0,"supportLMax",0, ...
    "paddingFactor",0,"terrainScale",0, ...
    "numberOfStationaryModes",0,"numberOfDynamicModes",0, ...
    "numberOfBoundaryModes",0,"numberOfBoundaryCandidates",0, ...
    "numberOfClaimedBoundaryModes",0, ...
    "numberOfRetainedNonTangentBottomDirections",0, ...
    "numberOfUnclassifiedWeakNullDirections",0, ...
    "numberOfNumericallyStationaryModes",0, ...
    "nonTangentBottomRetentionDefect",NaN, ...
    "zeroFrequencyStationaryProjectorDefect",NaN, ...
    "maximumUnresolvedZeroTrustedParticipation",NaN, ...
    "minimumAbsoluteBoundaryFrequency",NaN, ...
    "minimumBoundaryFrequencyResolutionRatio",NaN, ...
    "maximumEigenResidual",NaN,"maximumTrustedAPVDefect",NaN, ...
    "maximumBottomDefect",NaN,"maximumConjugacyDefect",NaN, ...
    "maximumStrongResidual",NaN, ...
    "maximumBoundaryCandidateAPVDefect",NaN, ...
    "maximumBoundaryCandidateBottomDefect",NaN, ...
    "maximumBoundaryCandidateStrongResidual",NaN);
end
