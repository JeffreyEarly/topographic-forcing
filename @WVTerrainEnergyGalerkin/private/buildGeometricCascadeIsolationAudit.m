function audit = buildGeometricCascadeIsolationAudit(problem,trustedBounds, ...
    waveIndices,apvIndices,stationaryDegree,orders,degrees, ...
    comparisonDegree,paddingFactors,terrainScales, ...
    shouldRunTwoDimensionalControl,cacheDirectory,quadratureOrder)
% Build the Milestone-10.2.2 geometric-cascade isolation oracle.

cache = emptyCacheDiagnostics(cacheDirectory);
[contract,record] = cachedAuditValue(problem,cacheDirectory, ...
    "production-contract",struct("trusted",trustedBounds, ...
    "waves",waveIndices,"apv",apvIndices,"padding",paddingFactors), ...
    @() productionContract(problem,trustedBounds,waveIndices, ...
    apvIndices,paddingFactors));
cache = appendCacheRecord(cache,record);
if ~contract.isCompatible
    error("WVTerrainEnergyGalerkin:IncompleteGeometricCascadeProductionContract", ...
        "Milestone 10.2.2 requires a passing production physical-state contract.")
end

[zonal,zonalCache] = cascadeStudy(problem,contract,trustedBounds, ...
    waveIndices,stationaryDegree,orders,degrees,comparisonDegree, ...
    paddingFactors,terrainScales,cacheDirectory,quadratureOrder,"zonal");
cache = mergeCacheDiagnostics(cache,zonalCache);

twoDimensional = struct("wasRequested",shouldRunTwoDimensionalControl, ...
    "wasAttempted",false,"status","skipped-zonal-not-resolved");
classification = zonal.classification;
diagnosis = zonal.diagnosis;
if zonal.classification == "cascade-resolved" ...
        && shouldRunTwoDimensionalControl
    [controlContract,record] = cachedAuditValue(problem,cacheDirectory, ...
        "production-contract",struct("trusted",[1 1],"waves",1, ...
        "apv",apvIndices,"padding",paddingFactors), ...
        @() productionContract(problem,[1 1],1,apvIndices,paddingFactors));
    cache = appendCacheRecord(cache,record);
    if ~controlContract.isCompatible
        error("WVTerrainEnergyGalerkin:IncompleteGeometricCascadeControlContract", ...
            "The two-dimensional control requires a passing production state contract.")
    end
    [control,controlCache] = cascadeStudy(problem,controlContract,[1 1], ...
        1,stationaryDegree,orders,degrees,comparisonDegree, ...
        paddingFactors,terrainScales,cacheDirectory,quadratureOrder, ...
        "two-dimensional");
    cache = mergeCacheDiagnostics(cache,controlCache);
    twoDimensional = control;
    twoDimensional.wasRequested = true;
    twoDimensional.wasAttempted = true;
    if control.classification ~= "cascade-resolved"
        classification = control.classification;
        diagnosis = "The zonal cascade is resolved, but the required two-dimensional control returns " ...
            +control.classification+". "+control.diagnosis;
    end
elseif zonal.classification == "cascade-resolved"
    twoDimensional.status = "not-requested";
end

audit = struct;
audit.scope = "milestone-10.2.2-geometric-cascade-isolation";
audit.status = classification;
audit.classification = classification;
audit.isCompatible = classification == "cascade-resolved";
audit.diagnosis = diagnosis;
audit.zonal = zonal;
audit.twoDimensionalControl = twoDimensional;
audit.productionContract = contract;
audit.cache = finalizeCacheDiagnostics(cache);
audit.trustedModeBounds = trustedBounds;
audit.waveModeIndices = waveIndices;
audit.apvModeIndices = apvIndices;
audit.stationaryPolynomialDegree = stationaryDegree;
audit.scatteringOrders = orders;
audit.primitivePolynomialDegrees = degrees;
audit.comparisonPolynomialDegree = comparisonDegree;
audit.paddingFactors = paddingFactors;
audit.terrainScales = terrainScales;
audit.quadratureOrder = quadratureOrder;
audit.requiredTolerance = zonal.requiredTolerance;
if classification == "cascade-resolved"
    audit.nextScope = "milestone-10.3-only-if-separately-authorized";
else
    audit.nextScope = "milestone-10.2.2-stop-for-analysis";
end
end

function contract = productionContract(problem,trusted,waves,apv,padding)
support = trusted+[0 0;0 1;0 2];
contract = problem.auditProductionPhysicalStateContract( ...
    trustedModeBounds=trusted,supportModeBounds=support, ...
    waveModeIndices=waves,apvModeIndices=apv, ...
    primitivePolynomialDegrees=[8;16;24],paddingFactors=padding);
end

function [study,cache] = cascadeStudy(problem,contract,trustedBounds, ...
    waveIndices,stationaryDegree,orders,degrees,comparisonDegree, ...
    paddingFactors,terrainScales,cacheDirectory,quadratureOrder,name)
cache = emptyCacheDiagnostics(cacheDirectory);
terrain = terrainSupport(problem);
shells = scatteringShells(problem,contract,terrain,orders,trustedBounds);
maximumLayout = shells.layouts{end};

primaryConfiguration = primitiveConfiguration(name,"comparison-primary", ...
    comparisonDegree,maximumLayout,paddingFactors(1),terrainScales, ...
    trustedBounds,quadratureOrder);
[primaryBundle,record] = cachedAuditValue(problem,cacheDirectory, ...
    "primitive",primaryConfiguration,@() primitiveBundle(problem, ...
    comparisonDegree,maximumLayout,paddingFactors(1),terrainScales, ...
    trustedBounds,quadratureOrder));
cache = appendCacheRecord(cache,record);
targets = productionWaveTargetBlocks(problem,contract, ...
    primaryBundle.context,primaryBundle.primitive.flatReference);
homotopy = continuePhysicalBlock(primaryBundle.context, ...
    primaryBundle.primitive,targets,terrainScales,trustedBounds, ...
    stationaryDegree,comparisonDegree);
primary = homotopy(end);

tail = tailDiagnostics(primaryBundle.context,homotopy,shells,degrees, ...
    comparisonDegree);
[horizontal,horizontalCache] = horizontalRefinement(problem, ...
    primaryBundle,primary,shells,orders,comparisonDegree, ...
    stationaryDegree,trustedBounds,paddingFactors(1), ...
    cacheDirectory,quadratureOrder,name);
cache = mergeCacheDiagnostics(cache,horizontalCache);
[vertical,verticalCache] = verticalRefinement(problem, ...
    primaryBundle,primary,shells,degrees,comparisonDegree, ...
    stationaryDegree,trustedBounds,paddingFactors(1), ...
    cacheDirectory,quadratureOrder,name);
cache = mergeCacheDiagnostics(cache,verticalCache);

secondaryConfiguration = primitiveConfiguration(name, ...
    "comparison-padding",comparisonDegree,maximumLayout, ...
    paddingFactors(2),[0;1],trustedBounds,quadratureOrder);
[secondaryBundle,record] = cachedAuditValue(problem,cacheDirectory, ...
    "primitive",secondaryConfiguration,@() primitiveBundle(problem, ...
    comparisonDegree,maximumLayout,paddingFactors(2),[0;1], ...
    trustedBounds,quadratureOrder));
cache = appendCacheRecord(cache,record);
padding = paddingDiagnostics(primaryBundle,primary,secondaryBundle, ...
    shells,degrees,comparisonDegree,trustedBounds,stationaryDegree);

drift = driftDecomposition(primaryBundle,primary,horizontal,vertical, ...
    trustedBounds,stationaryDegree,comparisonDegree,targets);
physical = physicalDiagnostics(primaryBundle.context, ...
    primary.direction,primary.basis,primary.frequency);
isolation = isolationDiagnostics(primary);
dimension = struct("original",targets.numberOfCoordinates, ...
    "continued",size(primary.basis,2), ...
    "stationary",size(primary.stationaryBasis,2), ...
    "ambient",size(primary.direction.E,1), ...
    "compressionFactor",size(primary.direction.E,1) ...
    /max(size(primary.stationaryBasis,2)+size(primary.basis,2),1));

tolerance = struct("transfer",1e-11,"structure",1e-12, ...
    "ritz",1e-10,"apv",1e-8,"bottom",1e-10,"strong",1e-5, ...
    "tail",1e-8,"projector",1e-8,"padding",1e-10, ...
    "isolation",1e3,"compression",2,"accounting",1e-11);
structurePasses = primary.energyHermitianDefect <= tolerance.structure ...
    && primary.exchangeSkewHermitianDefect <= tolerance.structure ...
    && tail.maximumStructuralDefect <= tolerance.transfer ...
    && drift.accountingDefect <= tolerance.accounting;
physicalPasses = physical.maximumRitzResidual <= tolerance.ritz ...
    && physical.maximumAPVDefect <= tolerance.apv ...
    && physical.maximumBottomDefect <= tolerance.bottom ...
    && physical.maximumStrongResidual <= tolerance.strong;
horizontalPasses = finalTwoPass(horizontal.projectorDefect,tolerance.projector) ...
    && finalTwoPass(horizontal.tailAmplitude,tolerance.tail);
verticalPasses = finalTwoPass(vertical.projectorDefect,tolerance.projector) ...
    && finalTwoPass(vertical.tailAmplitude,tolerance.tail);
paddingPasses = padding.projectorDefect <= tolerance.padding ...
    && padding.maximumTailDifference <= tolerance.padding;
isolationPasses = isolation.minimumResolutionRatio >= tolerance.isolation;
economyPasses = dimension.compressionFactor >= tolerance.compression;
converged = structurePasses && physicalPasses && horizontalPasses ...
    && verticalPasses && paddingPasses && isolationPasses && economyPasses;
isExpanded = dimension.continued > dimension.original;
systematic = systematicDecrease(horizontal.projectorDefect) ...
    && systematicDecrease(horizontal.tailAmplitude) ...
    && systematicDecrease(vertical.projectorDefect) ...
    && systematicDecrease(vertical.tailAmplitude);

if converged && ~isExpanded
    classification = "cascade-resolved";
    diagnosis = "The original production-wave block is spectrally isolated and its horizontal and vertical geometric tails converge.";
elseif converged
    classification = "resonant-block-required";
    diagnosis = "A larger backward-error-inseparable conjugate-closed physical block is required and converges.";
elseif structurePasses && physicalPasses && systematic
    classification = "cascade-slow";
    diagnosis = "The geometric tails decrease systematically, but the available support or vertical degree does not reach every accuracy or economy gate.";
else
    classification = "cascade-nonconvergent";
    diagnosis = cascadeDiagnosis(structurePasses,physicalPasses, ...
        horizontalPasses,verticalPasses,paddingPasses,isolationPasses, ...
        economyPasses,systematic);
end

study = struct;
study.scope = "geometric-cascade-"+name;
study.status = classification;
study.classification = classification;
study.isCompatible = classification == "cascade-resolved";
study.diagnosis = diagnosis;
study.wasAttempted = true;
study.terrainSupport = terrain;
study.shells = shells;
study.homotopy = homotopy;
study.primary = primary;
study.tail = tail;
study.horizontalRefinement = horizontal;
study.verticalRefinement = vertical;
study.padding = padding;
study.isolation = isolation;
study.driftDecomposition = drift;
study.physical = physical;
study.dimension = dimension;
study.requiredTolerance = tolerance;
study.waveModeIndices = waveIndices;
study.trustedModeBounds = trustedBounds;
study.cache = finalizeCacheDiagnostics(cache);
end

function value = primitiveConfiguration(name,role,degree,layout,padding, ...
    scales,trusted,quadrature)
value = struct("study",name,"role",role,"degree",degree, ...
    "kMode",layout.kMode,"lMode",layout.lMode,"padding",padding, ...
    "scales",scales,"trusted",trusted,"quadrature",quadrature);
end

function bundle = primitiveBundle(problem,degree,layout,padding,scales, ...
    trusted,quadrature)
[primitive,context] = buildGlobalSmallTerrainPrimitiveAudit(problem, ...
    degree,quadrature,1e-3,horizontalLayout=layout, ...
    paddingFactor=padding,trustedModeBounds=trusted, ...
    rejectTerrainNyquist=true,evaluationScales=scales, ...
    shouldAuditTangent=false);
bundle = struct("primitive",primitive,"context",context);
end

function terrain = terrainSupport(problem)
wvt = problem.originatingTransform;
spectrum = fft2(problem.topographicHeight)/(wvt.Nx*wvt.Ny);
nyquist = logical(WVGeometryDoublyPeriodic.maskForNyquistModes(wvt.Nx,wvt.Ny));
scale = max(abs(spectrum),[],"all");
tolerance = 100*eps*max(scale,1);
if max(abs(spectrum(nyquist)),[],"all") > tolerance
    error("WVTerrainEnergyGalerkin:GeometricCascadeTerrainNyquist", ...
        "The geometric-cascade oracle requires terrain with no material Nyquist coefficient.")
end
[kMode,lMode] = ndgrid(wvt.kMode_dft,wvt.lMode_dft);
active = abs(spectrum) > tolerance & ~nyquist ...
    & ~(kMode == 0 & lMode == 0);
terrain = struct("kMode",kMode(active),"lMode",lMode(active), ...
    "coefficient",spectrum(active),"spectralTolerance",tolerance);
if isempty(terrain.kMode)
    error("WVTerrainEnergyGalerkin:MissingGeometricCascadeTerrainSupport", ...
        "The terrain must contain at least one nonzero Fourier coefficient.")
end
end

function shells = scatteringShells(problem,contract,terrain,orders,trusted)
source = problem.horizontalLayout;
productionHorizontal = contract.productionHorizontalLayout;
seed = unique([productionHorizontal.kMode productionHorizontal.lMode],"rows");
seed = seed(any(seed ~= 0,2),:);
terrainModes = unique([terrain.kMode terrain.lMode],"rows");
reachable = cell(max(orders)+1,1);
shell = cell(max(orders)+1,1);
reachable{1} = seed;
shell{1} = seed;
for iOrder = 1:max(orders)
    sums = zeros(0,2);
    for iTerrain = 1:size(terrainModes,1)
        sums = [sums;reachable{iOrder}+terrainModes(iTerrain,:)]; %#ok<AGROW>
    end
    reachable{iOrder+1} = unique([reachable{iOrder};sums],"rows");
    shell{iOrder+1} = setdiff(reachable{iOrder+1},reachable{iOrder},"rows","stable");
end
trustedMask = abs(source.kMode) <= trusted(1) ...
    & abs(source.lMode) <= trusted(2);
infrastructure = [source.kMode(trustedMask) source.lMode(trustedMask)];
infrastructure = unique([infrastructure;0 0],"rows");
layouts = cell(numel(reachable),1);
modeSets = cell(numel(reachable),1);
for iOrder = 1:numel(reachable)
    requested = unique([infrastructure;reachable{iOrder}],"rows");
    [isPresent,index] = ismember(requested, ...
        [source.kMode source.lMode],"rows");
    if ~all(isPresent)
        missing = requested(find(~isPresent,1),:);
        error("WVTerrainEnergyGalerkin:UnavailableGeometricCascadeSupport", ...
            "The originating Fourier layout does not retain required scattering mode (%d,%d).",missing(1),missing(2))
    end
    keep = false(height(source),1);
    keep(index) = true;
    layout = source(keep,:);
    modes = [layout.kMode layout.lMode];
    if any(~ismember(-modes,modes,"rows"))
        error("WVTerrainEnergyGalerkin:IncompleteGeometricCascadeConjugacy", ...
            "Every geometric-scattering support must retain complete Fourier conjugates.")
    end
    layouts{iOrder} = layout;
    modeSets{iOrder} = modes;
end
shells = struct("orders",(0:max(orders)).', ...
    "terrainModes",terrainModes,"seedModes",seed, ...
    "reachableModes",{reachable},"shellModes",{shell}, ...
    "modeSets",{modeSets},"layouts",{layouts});
end

function history = continuePhysicalBlock(c,primitive,targets,scales, ...
    trusted,stationaryDegree,fullDegree)
nScale = numel(scales);
history = repmat(emptyContinuedBlock,nScale,1);
previous = targets.basis;
targetDimension = targets.numberOfCoordinates;
originalDimension = targetDimension;
for iScale = 1:nScale
    direction = primitive.evaluatedDirections(iScale);
    stationary = constructCompleteStationarySpace(c,direction, ...
        trusted,fullDegree,stationaryDegree,scales(iScale));
    G = stationary.trustedBasis;
    previous = projectOut(previous,G,direction.E);
    previous = energyOrthonormalize(previous,direction.E);
    dense = physicalEigensystem(direction,G);
    selected = selectByOverlap(dense.vectors,previous,direction.E, ...
        targetDimension);
    selected = closePhysicalBlock(selected,dense,direction.E, ...
        c.coordinateConjugateIndex);
    basis = dense.vectors(:,selected);
    basis = energyOrthonormalize(basis,direction.E);
    frequency = columnRayleighQuotients(basis,direction.E,1i*direction.J);
    block = continuedBlockDiagnostics(c,direction,G,basis,frequency, ...
        dense,selected,originalDimension,scales(iScale));
    block.stationaryBasis = G;
    block.stationary = stationary;
    history(iScale) = block;
    previous = basis;
    targetDimension = size(basis,2);
end

function indices = closePhysicalBlock(indices,dense,H,conjugateIndex)
previous = zeros(0,1);
while ~isequal(previous,indices)
    previous = indices;
    indices = expandInseparable(indices,dense.frequency, ...
        dense.absoluteFrequencyUncertainty);
    indices = closeConjugacy(indices,dense.vectors,H,conjugateIndex);
end
end
end

function dense = physicalEigensystem(direction,G)
H = direction.E;
A = 1i*direction.J;
R = chol((H+H')/2);
Yg = R*G;
if isempty(Yg)
    complement = eye(size(H));
else
    [complete,~] = qr(Yg);
    complement = complete(:,size(Yg,2)+1:end);
end
K = R'\(A/R);
reduced = complement'*K*complement;
[vectors,form] = schur((reduced+reduced')/2,"complex");
frequency = real(diag(form));
[frequency,order] = sort(frequency);
vectors = R\(complement*vectors(:,order));
vectors = energyNormalizeColumns(vectors,H);
residual = A*vectors-(H*vectors).*frequency.';
energyResidual = R'\residual;
uncertainty = vecnorm(energyResidual,2,1).' ...
    ./max(vecnorm(R*vectors,2,1).',realmin);
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

function indices = expandInseparable(indices,frequency,uncertainty)
indices = unique(indices(:));
changed = true;
while changed
    changed = false;
    candidates = setdiff((1:numel(frequency))',indices,"stable");
    for candidate = candidates.'
        if any(abs(frequency(candidate)-frequency(indices)) ...
                <=1e3*(uncertainty(candidate)+uncertainty(indices)))
            indices(end+1,1) = candidate; %#ok<AGROW>
            changed = true;
        end
    end
end
indices = sort(unique(indices));
end

function indices = closeConjugacy(indices,vectors,H,conjugateIndex)
n = size(H,1);
conjugacy = sparse((1:n)',conjugateIndex,1,n,n);
indices = unique(indices(:));
for index = reshape(indices,1,[])
    partner = conjugacy*conj(vectors(:,index));
    overlap = abs(vectors'*(H*partner));
    [~,iPartner] = max(overlap);
    indices(end+1,1) = iPartner; %#ok<AGROW>
end
indices = sort(unique(indices));
end

function result = continuedBlockDiagnostics(c,direction,G,basis,frequency, ...
    dense,selected,originalDimension,scale)
H = direction.E;
A = 1i*direction.J;
R = chol((H+H')/2);
residual = A*basis-(H*basis).*frequency.';
operatorNorm = norm(R'\(A/R),2);
ritz = vecnorm(R'\residual,2,1).' ...
    ./max(operatorNorm+abs(frequency),realmin);
unselected = setdiff((1:numel(dense.frequency))',selected,"stable");
gap = Inf;
resolution = Inf;
if ~isempty(unselected)
    for index = reshape(selected,1,[])
        [currentGap,iNear] = min(abs(dense.frequency(index) ...
            -dense.frequency(unselected)));
        denominator = dense.absoluteFrequencyUncertainty(index) ...
            +dense.absoluteFrequencyUncertainty(unselected(iNear));
        gap = min(gap,currentGap);
        resolution = min(resolution,currentGap/max(denominator,realmin));
    end
end
projector = energyProjector(basis,H);
conjugacy = sparse((1:size(H,1))',c.coordinateConjugateIndex, ...
    1,size(H,1),size(H,1));
conjugacyResidual = projector*conjugacy-conjugacy*conj(projector);
result = emptyContinuedBlock;
result.scale = scale;
result.direction = direction;
result.basis = basis;
result.frequency = frequency;
result.selectedIndices = selected;
result.maximumRitzResidual = max(ritz,[],"all");
result.minimumSpectralGap = gap;
result.minimumResolutionRatio = resolution;
result.maximumBackwardUncertainty = max( ...
    dense.absoluteFrequencyUncertainty(selected),[],"all");
result.wasExpanded = size(basis,2) > originalDimension;
result.energyHermitianDefect = norm(H-H',"fro")/max(norm(H,"fro"),realmin);
result.exchangeSkewHermitianDefect = norm(direction.J+direction.J',"fro") ...
    /max(norm(direction.J,"fro"),realmin);
result.stationaryOrthogonalityDefect = norm(G'*H*basis,"fro") ...
    /max(norm(G,"fro")*norm(H*basis,"fro"),realmin);
result.conjugacyDefect = norm(conjugacyResidual,"fro") ...
    /max(2*norm(projector,"fro"),realmin);
end

function value = emptyContinuedBlock
value = struct("scale",0,"direction",struct,"basis",zeros(0), ...
    "frequency",zeros(0,1),"selectedIndices",zeros(0,1), ...
    "maximumRitzResidual",Inf,"minimumSpectralGap",0, ...
    "minimumResolutionRatio",0,"maximumBackwardUncertainty",Inf, ...
    "wasExpanded",false,"energyHermitianDefect",Inf, ...
    "exchangeSkewHermitianDefect",Inf, ...
    "stationaryOrthogonalityDefect",Inf,"conjugacyDefect",Inf, ...
    "stationaryBasis",zeros(0),"stationary",struct);
end

function diagnostics = tailDiagnostics(c,history,shells, ...
    degrees,comparisonDegree)
allDegrees = [degrees(:);comparisonDegree];
nOrder = numel(shells.orders);
nDegree = numel(allDegrees);
nScale = numel(history);
amplitude = zeros(nOrder,nDegree,nScale);
energy = zeros(nOrder,nDegree,nScale);
projectorDrift = zeros(nOrder,nDegree,nScale);
capturedEnergy = zeros(nOrder,nDegree,nScale);
rankValues = zeros(nOrder,nDegree);
reconstruction = zeros(nOrder,nDegree);
idempotence = zeros(nOrder,nDegree);
selfAdjointness = zeros(nOrder,nDegree);
conjugacy = zeros(nOrder,nDegree);
for iOrder = 1:nOrder
    for iDegree = 1:nDegree
        subspace = primitiveSubspace(c,allDegrees(iDegree), ...
            shells.modeSets{iOrder});
        rankValues(iOrder,iDegree) = size(subspace.basis,2);
        reconstruction(iOrder,iDegree) = subspace.reconstructionDefect;
        conjugacy(iOrder,iDegree) = subspace.conjugacyDefect;
        for iScale = 1:nScale
            H = history(iScale).direction.E;
            Y = history(iScale).basis;
            projected = projectInto(subspace.basis,Y,H);
            R = chol((H+H')/2);
            omitted = Y-projected;
            amplitude(iOrder,iDegree,iScale) = norm(R*omitted,2);
            energy(iOrder,iDegree,iScale) = norm(R*omitted,"fro")^2 ...
                /max(size(Y,2),1);
            capturedEnergy(iOrder,iDegree,iScale) = norm(R*projected,"fro")^2 ...
                /max(size(Y,2),1);
            projectedBasis = energyOrthonormalize(projected,H);
            projectorDrift(iOrder,iDegree,iScale) = ...
                maximumPrincipalSine(projectedBasis,Y,H);
            if iScale == nScale
                [idempotence(iOrder,iDegree), ...
                    selfAdjointness(iOrder,iDegree)] = ...
                    projectionStructureDefects(subspace.basis,H);
            end
        end
    end
end

horizontalShellAmplitude = zeros(nOrder,nScale);
horizontalShellEnergy = zeros(nOrder,nScale);
for iScale = 1:nScale
    H = history(iScale).direction.E;
    Y = history(iScale).basis;
    R = chol((H+H')/2);
    previous = zeros(size(Y));
    for iOrder = 1:nOrder
        subspace = primitiveSubspace(c,comparisonDegree, ...
            shells.modeSets{iOrder});
        current = projectInto(subspace.basis,Y,H);
        shellPart = current-previous;
        horizontalShellAmplitude(iOrder,iScale) = norm(R*shellPart,2);
        horizontalShellEnergy(iOrder,iScale) = norm(R*shellPart,"fro")^2 ...
            /max(size(Y,2),1);
        previous = current;
    end
end
accounting = max(abs(capturedEnergy+energy-1),[],"all");
scaling = shellScaling(history,horizontalShellAmplitude, ...
    horizontalShellEnergy,shells.orders);
maximumStructural = max([reconstruction(:);idempotence(:); ...
    selfAdjointness(:);conjugacy(:);accounting]);
diagnostics = struct;
diagnostics.orders = shells.orders;
diagnostics.degrees = allDegrees;
diagnostics.amplitude = amplitude;
diagnostics.energy = energy;
diagnostics.projectorDrift = projectorDrift;
diagnostics.capturedEnergy = capturedEnergy;
diagnostics.horizontalShellAmplitude = horizontalShellAmplitude;
diagnostics.horizontalShellEnergy = horizontalShellEnergy;
diagnostics.rank = rankValues;
diagnostics.reconstructionDefect = reconstruction;
diagnostics.idempotenceDefect = idempotence;
diagnostics.energySelfAdjointnessDefect = selfAdjointness;
diagnostics.conjugacyDefect = conjugacy;
diagnostics.energyAccountingDefect = accounting;
diagnostics.maximumStructuralDefect = maximumStructural;
diagnostics.scaling = scaling;
diagnostics.fullScaleAmplitude = amplitude(:,:,end);
diagnostics.fullScaleEnergy = energy(:,:,end);
diagnostics.fullScaleProjectorDrift = projectorDrift(:,:,end);
end

function subspace = primitiveSubspace(c,degree,modes)
if degree > c.nG
    error("WVTerrainEnergyGalerkin:GeometricCascadeSubspaceDegree", ...
        "A nested primitive degree cannot exceed the comparison degree.")
end
[isPresent,horizontalIndex] = ismember(modes, ...
    [c.horizontalLayout.kMode c.horizontalLayout.lMode],"rows");
if ~all(isPresent)
    error("WVTerrainEnergyGalerkin:GeometricCascadeSubspaceSupport", ...
        "Every nested primitive mode must lie in the comparison support.")
end
horizontalIndex = sort(horizontalIndex);
nColumns = 0;
for iHorizontal = reshape(horizontalIndex,1,[])
    if hypot(c.horizontalLayout.k(iHorizontal), ...
            c.horizontalLayout.l(iHorizontal)) > 0
        nColumns = nColumns+3*degree+2;
    else
        nColumns = nColumns+3*(degree+1);
    end
end
raw = zeros(c.nX,nColumns);
first = 1;
weight = c.quadrature.weight;
mass = c.spaces.F'*(weight.*c.spaces.F);
for iHorizontal = reshape(horizontalIndex,1,[])
    k = c.horizontalLayout.k(iHorizontal);
    l = c.horizontalLayout.l(iHorizontal);
    kappa = hypot(k,l);
    rawRows = (iHorizontal-1)*c.nXBlock+(1:c.nXBlock);
    if kappa > 0
        local = zeros(c.nXBlock,3*degree+2);
        nF = degree+1;
        local(1:nF,1:nF) = -(l/kappa)*eye(nF);
        local(c.nF+(1:nF),1:nF) = (k/kappa)*eye(nF);
        divergence = c.spaces.F'*(weight.*c.spaces.Gxi(:,1:degree));
        longitudinal = (1i/kappa)*(mass\divergence);
        wColumns = nF+(1:degree);
        local(1:c.nF,wColumns) = (k/kappa)*longitudinal;
        local(c.nF+(1:c.nF),wColumns) = (l/kappa)*longitudinal;
        local(2*c.nF+(1:degree),wColumns) = eye(degree);
        etaColumns = nF+degree+(1:degree+1);
        etaRows = 2*c.nF+c.nG+[1:degree c.nH];
        local(etaRows,etaColumns) = eye(degree+1);
    else
        nF = degree+1;
        local = zeros(c.nXBlock,3*nF);
        local(1:nF,1:nF) = eye(nF);
        local(c.nF+(1:nF),nF+(1:nF)) = eye(nF);
        etaColumns = 2*nF+(1:nF);
        etaRows = 2*c.nF+c.nG+[1:degree c.nH];
        local(etaRows,etaColumns) = eye(nF);
    end
    columns = first:first+size(local,2)-1;
    raw(rawRows,columns) = local;
    first = columns(end)+1;
end
T = c.N'*raw;
represented = c.N*T;
reconstruction = norm(represented-raw,"fro")/max(norm(raw,"fro"),realmin);
T = independentColumns(T,eye(size(T,1)));
n = size(T,1);
coordinateConjugacy = sparse((1:n)',c.coordinateConjugateIndex,1,n,n);
P = T*((T'*T)\T');
conjugacyResidual = P*coordinateConjugacy-coordinateConjugacy*conj(P);
subspace = struct("basis",T,"reconstructionDefect",reconstruction, ...
    "conjugacyDefect",norm(conjugacyResidual,"fro") ...
    /max(2*norm(P,"fro"),realmin));
end

function projected = projectInto(basis,values,H)
if isempty(basis)
    projected = zeros(size(values));
else
    gram = basis'*H*basis;
    projected = basis*(gram\(basis'*H*values));
end
end

function [idempotence,selfAdjointness] = projectionStructureDefects(basis,H)
if isempty(basis)
    idempotence = 0;
    selfAdjointness = 0;
    return
end
P = basis*((basis'*H*basis)\(basis'*H));
idempotence = norm(P*P-P,"fro")/max(norm(P,"fro"),realmin);
selfAdjointness = norm(P'*H-H*P,"fro") ...
    /max(norm(H*P,"fro"),realmin);
end

function scaling = shellScaling(history,amplitude,energy,orders)
scales = [history.scale].';
usable = find(scales > 0,3,"first");
amplitudeOrder = NaN(numel(orders),1);
energyOrder = NaN(numel(orders),1);
if numel(usable) >= 3
    for iOrder = 2:numel(orders)
        amplitudeOrder(iOrder) = fittedOrder(scales(usable), ...
            amplitude(iOrder,usable).');
        energyOrder(iOrder) = fittedOrder(scales(usable), ...
            energy(iOrder,usable).');
    end
end
scaling = struct("amplitudeOrder",amplitudeOrder, ...
    "energyOrder",energyOrder,"expectedAmplitudeOrder",orders, ...
    "expectedEnergyOrder",2*orders);
end

function order = fittedOrder(scale,value)
valid = value > 100*eps*max(max(value),1);
if nnz(valid) < 2
    order = NaN;
else
    coefficients = polyfit(log(scale(valid)),log(value(valid)),1);
    order = coefficients(1);
end
end

function [diagnostics,cache] = horizontalRefinement(problem,comparison, ...
    target,shells,orders,degree,stationaryDegree,trusted,padding, ...
    cacheDirectory,quadrature,name)
cache = emptyCacheDiagnostics(cacheDirectory);
values = repmat(emptyRefinement,numel(orders),1);
for iOrder = 1:numel(orders)
    layout = shells.layouts{orders(iOrder)+1};
    if orders(iOrder) == orders(end)
        values(iOrder) = comparisonRefinement(target,comparison.context, ...
            degree,orders(iOrder));
        continue
    end
    configuration = primitiveConfiguration(name, ...
        "horizontal-refinement",degree,layout,padding,1,trusted,quadrature);
    [bundle,record] = cachedAuditValue(problem,cacheDirectory, ...
        "primitive",configuration,@() primitiveBundle(problem,degree, ...
        layout,padding,1,trusted,quadrature));
    cache = appendCacheRecord(cache,record);
    values(iOrder) = refinementCase(bundle,comparison,target, ...
        trusted,stationaryDegree,degree,orders(iOrder));
end
diagnostics = summarizeRefinement(values,"scatteringOrder",orders);
end

function [diagnostics,cache] = verticalRefinement(problem,comparison, ...
    target,shells,degrees,comparisonDegree,stationaryDegree,trusted, ...
    padding,cacheDirectory,quadrature,name)
cache = emptyCacheDiagnostics(cacheDirectory);
allDegrees = [degrees(:);comparisonDegree];
values = repmat(emptyRefinement,numel(allDegrees),1);
layout = shells.layouts{end};
for iDegree = 1:numel(allDegrees)
    degree = allDegrees(iDegree);
    if degree == comparisonDegree
        values(iDegree) = comparisonRefinement(target,comparison.context, ...
            degree,shells.orders(end));
        continue
    end
    configuration = primitiveConfiguration(name,"vertical-refinement", ...
        degree,layout,padding,1,trusted,quadrature);
    [bundle,record] = cachedAuditValue(problem,cacheDirectory, ...
        "primitive",configuration,@() primitiveBundle(problem,degree, ...
        layout,padding,1,trusted,quadrature));
    cache = appendCacheRecord(cache,record);
    values(iDegree) = refinementCase(bundle,comparison,target, ...
        trusted,stationaryDegree,degree,shells.orders(end));
end
diagnostics = summarizeRefinement(values,"polynomialDegree",allDegrees);
end

function result = refinementCase(bundle,comparison,target,trusted, ...
    stationaryDegree,degree,order)
c = bundle.context;
direction = bundle.primitive.evaluatedDirections(end);
stationary = constructCompleteStationarySpace(c,direction,trusted, ...
    degree,min(stationaryDegree,degree),1);
dense = physicalEigensystem(direction,stationary.trustedBasis);
transfer = nestedTransfer(c,comparison.context,direction, ...
    target.direction);
restrictedTarget = restrictFromComparison(transfer,target.basis, ...
    target.direction.E);
restrictedTarget = energyOrthonormalize(restrictedTarget,direction.E);
selected = selectByOverlap(dense.vectors,restrictedTarget, ...
    direction.E,size(target.basis,2));
selected = expandInseparable(selected,dense.frequency, ...
    dense.absoluteFrequencyUncertainty);
selected = closeConjugacy(selected,dense.vectors,direction.E, ...
    c.coordinateConjugateIndex);
basis = energyOrthonormalize(dense.vectors(:,selected),direction.E);
frequency = columnRayleighQuotients(basis,direction.E,1i*direction.J);
commonBasis = energyOrthonormalize(transfer*basis,target.direction.E);
tailBasis = primitiveSubspace(comparison.context,degree, ...
    [c.horizontalLayout.kMode c.horizontalLayout.lMode]);
projectedTarget = projectInto(tailBasis.basis,target.basis, ...
    target.direction.E);
physical = physicalDiagnostics(c,direction,basis,frequency);
formEnergy = transfer'*target.direction.E*transfer;
formExchange = transfer'*target.direction.J*transfer;
result = emptyRefinement;
result.scatteringOrder = order;
result.polynomialDegree = degree;
result.basis = basis;
result.commonBasis = commonBasis;
result.frequency = frequency;
result.projectorDefect = maximumPrincipalSine(commonBasis, ...
    target.basis,target.direction.E);
R = chol((target.direction.E+target.direction.E')/2);
result.tailAmplitude = norm(R*(target.basis-projectedTarget),2);
result.tailEnergy = norm(R*(target.basis-projectedTarget),"fro")^2 ...
    /max(size(target.basis,2),1);
result.energyRestrictionDefect = norm(direction.E-formEnergy,"fro") ...
    /max(norm(direction.E,"fro"),realmin);
result.exchangeRestrictionDefect = norm(direction.J-formExchange,"fro") ...
    /max(norm(direction.J,"fro"),realmin);
result.transferRawDefect = transferDiagnostics(c,comparison.context, ...
    transfer,direction,target.direction);
result.maximumRitzResidual = physical.maximumRitzResidual;
result.maximumAPVDefect = physical.maximumAPVDefect;
result.maximumBottomDefect = physical.maximumBottomDefect;
result.maximumStrongResidual = physical.maximumStrongResidual;
result.blockDimension = size(basis,2);
result.ambientDimension = size(direction.E,1);
result.wasExpanded = size(basis,2) > size(target.basis,2);
end

function result = comparisonRefinement(target,c,degree,order)
result = emptyRefinement;
result.scatteringOrder = order;
result.polynomialDegree = degree;
result.basis = target.basis;
result.commonBasis = target.basis;
result.frequency = target.frequency;
result.projectorDefect = 0;
result.tailAmplitude = 0;
result.tailEnergy = 0;
result.energyRestrictionDefect = 0;
result.exchangeRestrictionDefect = 0;
result.transferRawDefect = 0;
result.maximumRitzResidual = target.maximumRitzResidual;
physical = physicalDiagnostics(c,target.direction,target.basis, ...
    target.frequency);
result.maximumAPVDefect = physical.maximumAPVDefect;
result.maximumBottomDefect = physical.maximumBottomDefect;
result.maximumStrongResidual = physical.maximumStrongResidual;
result.blockDimension = size(target.basis,2);
result.ambientDimension = size(target.direction.E,1);
result.wasExpanded = target.wasExpanded;
end

function value = emptyRefinement
value = struct("scatteringOrder",0,"polynomialDegree",0, ...
    "basis",zeros(0),"commonBasis",zeros(0),"frequency",zeros(0,1), ...
    "projectorDefect",Inf,"tailAmplitude",Inf,"tailEnergy",Inf, ...
    "energyRestrictionDefect",Inf,"exchangeRestrictionDefect",Inf, ...
    "transferRawDefect",Inf,"maximumRitzResidual",Inf, ...
    "maximumAPVDefect",Inf,"maximumBottomDefect",Inf, ...
    "maximumStrongResidual",Inf,"blockDimension",0, ...
    "ambientDimension",0,"wasExpanded",false);
end

function diagnostics = summarizeRefinement(values,coordinateName,coordinate)
diagnostics = struct;
diagnostics.values = values;
diagnostics.(coordinateName) = coordinate;
diagnostics.projectorDefect = [values.projectorDefect].';
diagnostics.tailAmplitude = [values.tailAmplitude].';
diagnostics.tailEnergy = [values.tailEnergy].';
diagnostics.energyRestrictionDefect = [values.energyRestrictionDefect].';
diagnostics.exchangeRestrictionDefect = [values.exchangeRestrictionDefect].';
diagnostics.transferRawDefect = [values.transferRawDefect].';
diagnostics.blockDimension = [values.blockDimension].';
end

function transfer = nestedTransfer(coarse,fine,coarseDirection,fineDirection)
coarseRepresentation = nestedRepresentation(coarse,coarseDirection);
fineRepresentation = nestedRepresentation(fine,fineDirection);
rawMap = rawNestedMap(coarseRepresentation,fineRepresentation);
embedded = rawMap*coarseRepresentation.admissibleBasis;
transfer = fineRepresentation.admissibleBasis\embedded;
end

function representation = nestedRepresentation(c,direction)
representation = struct("admissibleBasis",c.N, ...
    "numberOfRawStateCoefficients",c.nX, ...
    "rawStateBlockSize",c.nXBlock,"numberOfFCoordinates",c.nF, ...
    "numberOfGCoordinates",c.nG,"numberOfHCoordinates",c.nH, ...
    "horizontalLayout",c.horizontalLayout,"bottomMap",direction.B);
end

function rawMap = rawNestedMap(coarse,fine)
rawMap = sparse(fine.numberOfRawStateCoefficients, ...
    coarse.numberOfRawStateCoefficients);
for iCoarse = 1:height(coarse.horizontalLayout)
    mode = [coarse.horizontalLayout.kMode(iCoarse) ...
        coarse.horizontalLayout.lMode(iCoarse)];
    iFine = find(fine.horizontalLayout.kMode == mode(1) ...
        & fine.horizontalLayout.lMode == mode(2),1);
    if isempty(iFine)
        error("WVTerrainEnergyGalerkin:NonNestedGeometricCascadeSupport", ...
            "Every coarse horizontal mode must lie in the comparison support.")
    end
    coarseBlock = (iCoarse-1)*coarse.rawStateBlockSize;
    fineBlock = (iFine-1)*fine.rawStateBlockSize;
    rawMap = insertIdentity(rawMap, ...
        fineBlock+(1:coarse.numberOfFCoordinates), ...
        coarseBlock+(1:coarse.numberOfFCoordinates));
    rawMap = insertIdentity(rawMap, ...
        fineBlock+fine.numberOfFCoordinates+(1:coarse.numberOfFCoordinates), ...
        coarseBlock+coarse.numberOfFCoordinates+(1:coarse.numberOfFCoordinates));
    rawMap = insertIdentity(rawMap, ...
        fineBlock+2*fine.numberOfFCoordinates+(1:coarse.numberOfGCoordinates), ...
        coarseBlock+2*coarse.numberOfFCoordinates+(1:coarse.numberOfGCoordinates));
    fineEta = fineBlock+2*fine.numberOfFCoordinates+fine.numberOfGCoordinates;
    coarseEta = coarseBlock+2*coarse.numberOfFCoordinates+coarse.numberOfGCoordinates;
    rawMap = insertIdentity(rawMap,fineEta+(1:coarse.numberOfGCoordinates), ...
        coarseEta+(1:coarse.numberOfGCoordinates));
    rawMap(fineEta+fine.numberOfHCoordinates, ...
        coarseEta+coarse.numberOfHCoordinates) = 1;
end
end

function values = insertIdentity(values,rows,columns)
values(rows,columns) = speye(numel(rows),numel(columns));
end

function target = restrictFromComparison(transfer,comparison,H)
gram = transfer'*H*transfer;
target = gram\(transfer'*H*comparison);
end

function defect = transferDiagnostics(coarse,fine,transfer, ...
    coarseDirection,fineDirection)
coarseRepresentation = nestedRepresentation(coarse,coarseDirection);
fineRepresentation = nestedRepresentation(fine,fineDirection);
raw = rawNestedMap(coarseRepresentation,fineRepresentation) ...
    *coarseRepresentation.admissibleBasis;
represented = fineRepresentation.admissibleBasis*transfer;
defect = norm(represented-raw,"fro")/max(norm(raw,"fro"),realmin);
end

function diagnostics = paddingDiagnostics(primaryBundle,primary, ...
    secondaryBundle,shells,degrees,comparisonDegree,trusted,stationaryDegree)
cSecondary = secondaryBundle.context;
directionSecondary = secondaryBundle.primitive.evaluatedDirections(end);
stationary = constructCompleteStationarySpace(cSecondary, ...
    directionSecondary,trusted,comparisonDegree,stationaryDegree,1);
dense = physicalEigensystem(directionSecondary,stationary.trustedBasis);
toSecondary = nestedTransfer(primaryBundle.context,cSecondary, ...
    primary.direction,directionSecondary);
target = energyOrthonormalize(toSecondary*primary.basis, ...
    directionSecondary.E);
selected = selectByOverlap(dense.vectors,target,directionSecondary.E, ...
    size(target,2));
selected = expandInseparable(selected,dense.frequency, ...
    dense.absoluteFrequencyUncertainty);
selected = closeConjugacy(selected,dense.vectors,directionSecondary.E, ...
    cSecondary.coordinateConjugateIndex);
secondaryBasis = energyOrthonormalize(dense.vectors(:,selected), ...
    directionSecondary.E);
toPrimary = nestedTransfer(cSecondary,primaryBundle.context, ...
    directionSecondary,primary.direction);
secondaryInPrimary = energyOrthonormalize(toPrimary*secondaryBasis, ...
    primary.direction.E);
projectorDefect = maximumPrincipalSine(primary.basis,secondaryInPrimary, ...
    primary.direction.E);

secondaryBlock = emptyContinuedBlock;
secondaryBlock.scale = 1;
secondaryBlock.direction = directionSecondary;
secondaryBlock.basis = secondaryBasis;
secondaryBlock.frequency = columnRayleighQuotients(secondaryBasis, ...
    directionSecondary.E,1i*directionSecondary.J);
secondaryTail = tailDiagnostics(cSecondary,secondaryBlock,shells,degrees, ...
    comparisonDegree);
primaryTail = tailDiagnostics(primaryBundle.context,primary,shells,degrees, ...
    comparisonDegree);
maximumTailDifference = max([ ...
    abs(primaryTail.fullScaleAmplitude-secondaryTail.fullScaleAmplitude); ...
    abs(primaryTail.fullScaleEnergy-secondaryTail.fullScaleEnergy)],[],"all");
diagnostics = struct("projectorDefect",projectorDefect, ...
    "maximumTailDifference",maximumTailDifference, ...
    "secondaryBasis",secondaryBasis,"secondaryTail",secondaryTail, ...
    "transferDefect",transferDiagnostics(primaryBundle.context, ...
    cSecondary,toSecondary,primary.direction,directionSecondary));
end

function diagnostics = isolationDiagnostics(block)
diagnostics = struct("minimumSpectralGap",block.minimumSpectralGap, ...
    "maximumBackwardUncertainty",block.maximumBackwardUncertainty, ...
    "minimumResolutionRatio",block.minimumResolutionRatio, ...
    "continuedDimension",size(block.basis,2), ...
    "wasExpanded",block.wasExpanded);
end

function diagnostics = driftDecomposition(bundle,target,horizontal,vertical, ...
    trusted,stationaryDegree,degree,productionTargets)
allValues = [horizontal.values;vertical.values];
[~,iWorst] = max([allValues.projectorDefect]);
worst = allValues(iWorst);
H = target.direction.E;
if isempty(worst.commonBasis) || size(worst.commonBasis,2) ...
        ~= size(target.basis,2)
    vector = target.basis(:,1);
else
    vector = worstPrincipalVector(worst.commonBasis,target.basis,H);
end
G = target.stationaryBasis;
declared = target.basis;
bottom = target.stationary.trustedBottomSeedBasis;
flatStationary = constructCompleteStationarySpace(bundle.context, ...
    bundle.primitive.flatReference,trusted,degree,stationaryDegree,0);
flat = physicalEigensystem(bundle.primitive.flatReference, ...
    flatStationary.trustedBasis);
resolved = abs(flat.frequency) ...
    >1e3*flat.absoluteFrequencyUncertainty;
higher = flat.vectors(:,resolved);
higher = projectOut(higher,productionTargets.basis, ...
    bundle.primitive.flatReference.E);

[stationaryBasis,declaredBasis,bottomBasis,higherBasis] = ...
    orthogonalFamilyHierarchy(G,declared,bottom,higher,H);
fractions = zeros(5,1);
remainder = vector;
[fractions(1),remainder] = extractEnergy(remainder,stationaryBasis,H);
[fractions(2),remainder] = extractEnergy(remainder,declaredBasis,H);
[fractions(3),remainder] = extractEnergy(remainder,bottomBasis,H);
[fractions(4),remainder] = extractEnergy(remainder,higherBasis,H);
fractions(5) = real(remainder'*H*remainder);
fractions = fractions/max(sum(fractions),realmin);
diagnostics = struct("family",["stationary";"declared-internal"; ...
    "non-tangent-bottom";"higher-flat-wave";"unresolved"], ...
    "energyFraction",fractions,"accountingDefect", ...
    abs(sum(fractions)-1),"worstProjectorDefect", ...
    worst.projectorDefect,"worstScatteringOrder", ...
    worst.scatteringOrder,"worstPolynomialDegree", ...
    worst.polynomialDegree);
end

function vector = worstPrincipalVector(first,second,H)
R = chol((H+H')/2);
firstQ = orth(R*first);
secondQ = orth(R*second);
[U,~,~] = svd(firstQ'*secondQ,"econ");
vector = R\(firstQ*U(:,end));
vector = vector/sqrt(real(vector'*H*vector));
end

function [G,Y,B,W] = orthogonalFamilyHierarchy(G,Y,B,W,H)
G = energyOrthonormalize(G,H);
Y = energyOrthonormalize(projectOut(Y,G,H),H);
B = energyOrthonormalize(projectOut(projectOut(B,G,H),Y,H),H);
W = projectOut(W,G,H);
W = projectOut(W,Y,H);
W = projectOut(W,B,H);
W = energyOrthonormalize(W,H);
end

function [fraction,remainder] = extractEnergy(vector,basis,H)
if isempty(basis)
    fraction = 0;
    remainder = vector;
    return
end
component = projectInto(basis,vector,H);
fraction = real(component'*H*component);
remainder = vector-component;
end

function diagnostics = physicalDiagnostics(c,direction,basis,frequency)
H = direction.E;
A = 1i*direction.J;
R = chol((H+H')/2);
operatorNorm = norm(R'\(A/R),2);
residual = A*basis-(H*basis).*frequency.';
ritz = vecnorm(R'\residual,2,1).' ...
    ./max(operatorNorm+abs(frequency),realmin);
diagnostics = struct;
diagnostics.maximumRitzResidual = max(ritz,[],"all");
diagnostics.maximumAPVDefect = trustedAPVDefect(c,direction,basis);
diagnostics.maximumBottomDefect = modalBottomDefect(direction,basis,frequency);
strong = strongModeDiagnostics(c,direction,basis,frequency);
diagnostics.maximumStrongResidual = strong.maximumResidual;
end

function defect = trustedAPVDefect(c,direction,vectors)
trusted = c.layout.trustedHorizontalModes;
rows = find(repmat(trusted,c.nZ,1));
Q = direction.QProjected(rows,:);
weight = c.apvVerticalWeight(rows);
R = chol(direction.E);
scale = max(norm(sqrt(weight).*Q/R,2),realmin);
defect = max(vecnorm(sqrt(weight).*(Q*vectors),2,1).'/scale,[],"all");
end

function defect = modalBottomDefect(direction,vectors,frequency)
residual = -1i*(direction.B*vectors).*frequency.'-direction.R*vectors;
R = chol(direction.E);
scale = abs(frequency)*norm(direction.B/R,2)+norm(direction.R/R,2);
defect = max(vecnorm(residual,2,1).'./max(scale,realmin),[],"all");
end

function diagnostics = strongModeDiagnostics(c,direction,C,frequency)
if isempty(C)
    diagnostics = struct("maximumResidual",0);
    return
end
rawState = c.N*C;
pressure = direction.pressureMap*C;
stateTendency = -1i*rawState.*frequency.';
descriptorResidual = direction.matrices.S*[stateTendency;pressure] ...
    -direction.matrices.F*C;
primitiveResidual = direction.L*C+1i*C.*frequency.';
descriptorDefect = columnRelativeResidual(descriptorResidual, ...
    {direction.matrices.S*[stateTendency;pressure],-direction.matrices.F*C});
primitiveDefect = columnRelativeResidual(primitiveResidual, ...
    {direction.L*C,1i*C.*frequency.'});
gaugeResidual = c.pressureGauge'*[zeros(c.nX,size(C,2));pressure];
gaugeDefect = vecnorm(gaugeResidual,2,1).' ...
    ./max(vecnorm(pressure,2,1).',realmin);
pointwise = pointwiseStrongResiduals(c,direction,rawState,pressure,frequency);
maximum = max([descriptorDefect primitiveDefect gaugeDefect ...
    pointwise.maximumByMode],[],2);
diagnostics = struct("maximumResidual",max(maximum,[],"all"), ...
    "maximumByMode",maximum);
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
uDefect = weightedColumnResidual( ...
    {-1i*u.*omega,-c.wvt.f*v,Dxp/c.wvt.rho0},c.volumeWeight);
vDefect = weightedColumnResidual( ...
    {-1i*v.*omega,c.wvt.f*u,Dyp/c.wvt.rho0},c.volumeWeight);
wDefect = weightedColumnResidual( ...
    {-1i*w.*omega,N2.*eta,pXi./(c.wvt.rho0*gamma)},c.volumeWeight);
etaDefect = weightedColumnResidual({-1i*eta.*omega,-w},c.volumeWeight);
divergence = (c.RuX+c.RvY+c.RwhXi)*rawState;
continuityDefect = weightedColumnNorm(divergence,c.volumeWeight) ...
    ./max(weightedColumnNorm(c.RuX*rawState,c.volumeWeight) ...
    +weightedColumnNorm(c.RvY*rawState,c.volumeWeight) ...
    +weightedColumnNorm(c.RwhXi*rawState,c.volumeWeight),realmin);
maximum = max([uDefect vDefect wDefect etaDefect continuityDefect],[],2);
diagnostics = struct("maximumByMode",maximum, ...
    "maximumResidual",max(maximum,[],"all"));
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

function values = columnRayleighQuotients(C,H,A)
values = zeros(size(C,2),1);
for iColumn = 1:size(C,2)
    values(iColumn) = real((C(:,iColumn)'*A*C(:,iColumn)) ...
        /(C(:,iColumn)'*H*C(:,iColumn)));
end
end

function X = projectOut(X,basis,H)
if ~isempty(X) && ~isempty(basis)
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
[Q,~] = qr(R*X,0);
singularValues = svd(R*X,"econ");
tolerance = max(size(X))*eps(max(singularValues,[],'all'));
rankX = nnz(singularValues > tolerance);
X = R\Q(:,1:rankX);
end

function X = energyNormalizeColumns(X,H)
for iColumn = 1:size(X,2)
    X(:,iColumn) = X(:,iColumn) ...
        /sqrt(real(X(:,iColumn)'*H*X(:,iColumn)));
end
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

function passes = finalTwoPass(values,tolerance)
values = values(:);
passes = numel(values) >= 2 && all(values(end-1:end) <= tolerance);
end

function passes = systematicDecrease(values)
values = values(:);
if numel(values) < 3
    passes = false;
    return
end
tail = values(end-2:end);
slack = 100*eps(max(max(tail),1));
passes = all(diff(tail) <= slack) ...
    && tail(1) >= 4*max(tail(end),realmin);
end

function diagnosis = cascadeDiagnosis(structure,physical,horizontal, ...
    vertical,padding,isolation,economy,systematic)
if ~structure
    diagnosis = "The common primitive transfer, energy structure, conjugacy, or tail accounting failed.";
elseif ~physical
    diagnosis = "The continued comparison block fails a mandatory Ritz, APV, bottom, or strong primitive residual.";
elseif ~horizontal
    diagnosis = "The horizontal geometric-scattering tail does not reach the required convergence gate.";
elseif ~vertical
    diagnosis = "The vertical geometric-scattering tail does not reach the required convergence gate.";
elseif ~padding
    diagnosis = "Padding factors two and three do not reproduce the same physical projector and tail measures.";
elseif ~isolation
    diagnosis = "The continued block is not separated from its complement by the required backward-error margin.";
elseif ~economy
    diagnosis = "The accepted stationary-plus-wave representation does not retain factor-two compression.";
elseif ~systematic
    diagnosis = "At least one independent refinement direction lacks systematic tail decay.";
else
    diagnosis = "The available refinements do not establish a converged geometric-cascade classification.";
end
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

function first = mergeCacheDiagnostics(first,second)
first.enabled = first.enabled || second.enabled;
if first.directory == ""
    first.directory = second.directory;
end
first.records = [first.records;second.records];
end

function cache = finalizeCacheDiagnostics(cache)
if isempty(cache.records)
    cache.numberOfHits = 0;
    cache.numberOfMisses = 0;
    cache.numberOfInvalidations = 0;
    cache.bytes = 0;
    cache.secondsSaved = 0;
    return
end
status = string({cache.records.status});
cache.numberOfHits = nnz(status == "hit");
cache.numberOfMisses = nnz(status == "miss" | status == "disabled");
cache.numberOfInvalidations = nnz(status == "invalidated");
cache.bytes = sum([cache.records.bytes]);
cache.secondsSaved = sum([cache.records.secondsSaved]);
end

function [value,record] = cachedAuditValue(problem,directory,kind, ...
    configuration,builder)
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
schema = "geometric-cascade-cache-v2";
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
if isstruct(value) && isfield(value,"context")
    value.context = rmfield(value.context,["problem","wvt"]);
end
end

function value = restoreCachedValue(value,problem)
if isstruct(value) && isfield(value,"context")
    value.context.problem = problem;
    value.context.wvt = problem.originatingTransform;
end
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
updateImplementationDigest(engine);
wvt = problem.originatingTransform;
constants = [wvt.Lx wvt.Ly wvt.Lz wvt.Nx wvt.Ny wvt.Nz ...
    wvt.f wvt.rho0 wvt.g double(wvt.shouldAntialias)];
updateDigest(engine,typecast(double(constants(:)),"uint8"));
updateDigest(engine,typecast(double(problem.topographicHeight(:)),"uint8"));
updateDigest(engine,typecast(double(problem.verticalModeIndices(:)),"uint8"));
layout = problem.horizontalLayout;
updateDigest(engine,typecast(double([layout.kMode;layout.lMode]),"uint8"));
z = linspace(-wvt.Lz,0,257).';
N2 = wvt.N2Function(z);
if isscalar(N2)
    N2 = repmat(N2,size(z));
end
updateDigest(engine,typecast(double(N2(:)),"uint8"));
digest = typecast(int8(engine.digest()),"uint8");
key = string(lower(reshape(dec2hex(digest,2).',1,[])));
end

function updateImplementationDigest(engine)
privateRoot = string(fileparts(mfilename("fullpath")));
classRoot = string(fileparts(privateRoot));
files = [ ...
    fullfile(classRoot,"auditGeometricCascadeIsolation.m"); ...
    fullfile(privateRoot,"buildGeometricCascadeIsolationAudit.m"); ...
    fullfile(privateRoot,"buildGlobalSmallTerrainPrimitiveAudit.m"); ...
    fullfile(classRoot,"auditProductionPhysicalStateContract.m")];
for iFile = 1:numel(files)
    updateDigest(engine,uint8(fileread(files(iFile))));
end
end

function updateDigest(engine,bytes)
engine.update(typecast(uint8(bytes(:)),"int8"));
end
