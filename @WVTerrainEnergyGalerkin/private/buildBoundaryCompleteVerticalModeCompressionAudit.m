function audit = buildBoundaryCompleteVerticalModeCompressionAudit( ...
    problem,trustedBounds,supportBounds,stationaryDegree,degrees, ...
    comparisonDegree,modalCounts,robinLengthRatios,paddingFactors,terrainScales, ...
    solverOrders,quadratureOrder,waveCoordinateType,referenceSlope,sharedSetup)
% Build the Milestone-9.2 boundary-complete modal-compression oracle.

if nargin < 13
    waveCoordinateType = "nonhydrostatic-dirichlet";
end
if nargin < 14
    referenceSlope = [0 0];
end
if nargin < 15
    sharedSetup = [];
end
waveCoordinateType = string(waveCoordinateType);
referenceSlope = reshape(referenceSlope,1,2);

degree = comparisonDegree;
support = supportBounds(end,:);
layout = retainedLayout(problem,support);
if isempty(quadratureOrder)
    order = max(2*degree+7,18);
else
    order = quadratureOrder;
end

nCount = numel(modalCounts);
nRobin = numel(robinLengthRatios);
nPadding = numel(paddingFactors);
details = cell(nCount,nRobin,nPadding);
comparisonReferences = cell(nPadding,1);
if isempty(sharedSetup)
    reference = problem.auditConvergedPhysicalSubspaces( ...
        trustedModeBounds=trustedBounds, ...
        supportModeBounds=supportBounds, ...
        stationaryPolynomialDegree=stationaryDegree, ...
        primitivePolynomialDegrees=degrees, ...
        paddingFactors=paddingFactors,terrainScales=terrainScales, ...
        quadratureOrder=quadratureOrder);
    sharedSetup = repmat(struct("primitive",[], ...
        "context",[],"direction",[],"stationary",[], ...
        "comparisonReference",[],"providerCatalogs",[], ...
        "reference",[]),nPadding,1);
    for iPadding = 1:nPadding
        [primitive,context] = buildGlobalSmallTerrainPrimitiveAudit( ...
            problem,degree,order,1e-3,horizontalLayout=layout, ...
            paddingFactor=paddingFactors(iPadding), ...
            trustedModeBounds=trustedBounds,rejectTerrainNyquist=true, ...
            evaluationScales=terrainScales);
        direction = primitive.evaluatedDirections(end);
        stationary = constructCompleteStationarySpace(context,direction, ...
            trustedBounds,degree,stationaryDegree,terrainScales(end));
        comparisonReference = buildComparisonPrimitiveReference( ...
            context,direction,stationary,reference);
        providerCatalogs = buildVerticalCatalog( ...
            problem,context,primitive.flatReference,max(modalCounts), ...
            robinLengthRatios,solverOrders, ...
            "nonhydrostatic-dirichlet",[0 0],[]);
        sharedSetup(iPadding).primitive = primitive;
        sharedSetup(iPadding).context = context;
        sharedSetup(iPadding).direction = direction;
        sharedSetup(iPadding).stationary = stationary;
        sharedSetup(iPadding).comparisonReference = comparisonReference;
        sharedSetup(iPadding).providerCatalogs = providerCatalogs;
    end
else
    reference = sharedSetup(1).reference;
end
for iPadding = 1:nPadding
    primitive = sharedSetup(iPadding).primitive;
    context = sharedSetup(iPadding).context;
    direction = sharedSetup(iPadding).direction;
    stationary = sharedSetup(iPadding).stationary;
    comparisonReference = sharedSetup(iPadding).comparisonReference;
    comparisonReferences{iPadding} = comparisonReference;
    catalogs = buildVerticalCatalog(problem,context,primitive.flatReference, ...
        max(modalCounts),robinLengthRatios,solverOrders, ...
        waveCoordinateType,referenceSlope, ...
        sharedSetup(iPadding).providerCatalogs);
    for iCount = 1:nCount
        for iRobin = 1:nRobin
            basis = buildModalCoordinates(context,direction,stationary, ...
                catalogs(iRobin),modalCounts(iCount));
            result = reducedAudit(context,direction,stationary,basis, ...
                comparisonReference);
            result.modalCount = modalCounts(iCount);
            result.robinLengthRatio = robinLengthRatios(iRobin);
            result.robinLength = catalogs(iRobin).robinLength;
            result.paddingFactor = paddingFactors(iPadding);
            result.projection = primitive.projection;
            details{iCount,iRobin,iPadding} = result;
        end
    end
end

primaryRobin = find(robinLengthRatios == -1/4,1);
if isempty(primaryRobin)
    [~,primaryRobin] = min(abs(robinLengthRatios+1/4));
end
primary = details{end,primaryRobin,1};
convergence = countConvergence(details(:,primaryRobin,1));
padding = paddingConvergence(details(end,primaryRobin,:));
robin = robinConvergence(details(end,:,1),primaryRobin);

tolerance = struct( ...
    "provider",1e-11,"endpoint",1e-11, ...
    "structure",1e-12,"stationary",1e-10,"bottom",1e-10, ...
    "projector",1e-8,"frequency",1e-8,"padding",1e-8, ...
    "robin",1e-8,"apv",1e-8,"strong",1e-5);
allResults = [details{:}];
providerPasses = max([allResults.maximumProviderDefect]) ...
    <= tolerance.provider ...
    && max([allResults.maximumEndpointDefect]) <= tolerance.endpoint ...
    && max([allResults.negativeRobinModeCountDefect]) == 0;
structurePasses = max([allResults.energyHermitianDefect]) ...
    <= tolerance.structure ...
    && max([allResults.exchangeSkewHermitianDefect]) ...
    <= tolerance.structure ...
    && max([allResults.modalEnergyOrthogonalityDefect]) ...
    <= tolerance.structure;
compatibilityPasses = primary.stationaryRepresentationDefect ...
    <= tolerance.stationary ...
    && primary.stationaryRowDefect <= tolerance.stationary ...
    && primary.bottomEvolutionDefect <= tolerance.bottom;
physicalPasses = primary.internalProjectorDefect <= tolerance.projector ...
    && primary.internalFrequencyDefect <= tolerance.frequency ...
    && primary.maximumInternalAPVDefect <= tolerance.apv ...
    && primary.maximumInternalStrongResidual <= tolerance.strong ...
    && padding.maximumInternalProjectorDefect <= tolerance.padding;
converges = convergence.projectorDefect(end) ...
    <= tolerance.projector ...
    && convergence.frequencyDefect(end) <= tolerance.frequency;

if ~structurePasses
    classification = "modal-incompatible";
    diagnosis = "The reduced physical energy or exchange form fails before a scientific compression comparison can be made.";
elseif ~compatibilityPasses || ~converges
    classification = "modal-incompatible";
    diagnosis = "The boundary-complete modal span converges toward the internal-wave projector but does not close the finite-terrain bottom evolution or reach the required physical projector tolerance before losing useful compression.";
elseif ~providerPasses
    classification = "modal-incompatible";
    diagnosis = "The isolated InternalModesEVP provider does not meet the required order and endpoint convergence tolerance for this modal construction.";
elseif physicalPasses && primary.actualCompressionFactor >= 2
    classification = "modal-acceleration";
    diagnosis = "The boundary-complete modal family reproduces the validated physical subspaces with at least a factor-two reduction in admissible vertical degrees of freedom.";
elseif physicalPasses
    classification = "modal-equivalent";
    diagnosis = "The boundary-complete modal family reproduces the validated physical subspaces, but the stationary completion leaves less than a factor-two reduction in admissible vertical degrees of freedom.";
else
    classification = "modal-incompatible";
    diagnosis = "The one-dimensional mode families are individually valid, but their reduced terrain eigensystem does not reproduce the primitive physical projectors and independent physical residuals.";
end

audit = struct;
audit.scope = "milestone-9.2-boundary-complete-vertical-mode-compression";
audit.status = classification;
audit.classification = classification;
audit.isCompatible = classification ~= "modal-incompatible";
audit.diagnosis = diagnosis;
audit.internalModesEVPCommit = ...
    "df86687e91faa31bf65941299062d125a96904b1";
audit.providerQualification = struct( ...
    "isolatedCheckoutRequired",true, ...
    "queryOrderingAdapter","descending-z-evaluate-and-restore", ...
    "qualifiedOrders",solverOrders, ...
    "rejectedOrders",[96;192], ...
    "rejectionReason", ...
    "spurious-large-negative-Robin-eigenvalue");
audit.trustedModeBounds = trustedBounds;
audit.supportModeBounds = supportBounds;
audit.stationaryPolynomialDegree = stationaryDegree;
audit.primitivePolynomialDegrees = degrees;
audit.comparisonPolynomialDegree = comparisonDegree;
audit.modalCounts = modalCounts;
audit.robinLengthRatios = robinLengthRatios;
audit.paddingFactors = paddingFactors;
audit.terrainScales = terrainScales;
audit.internalModesEVPOrders = solverOrders;
audit.waveCoordinateType = waveCoordinateType;
audit.referenceSlope = referenceSlope;
audit.quadratureOrder = order;
audit.reference = reference;
for iPadding = 1:numel(sharedSetup)
    sharedSetup(iPadding).reference = reference;
end
audit.sharedSetup = sharedSetup;
audit.comparisonReference = comparisonReferences{1};
audit.comparisonReferences = comparisonReferences;
audit.details = details;
audit.primary = primary;
audit.convergence = convergence;
audit.padding = padding;
audit.robin = robin;
audit.robinInterpretation = struct( ...
    "role","finite-order-convergence-parameter", ...
    "isHardPhysicalGate",false, ...
    "selectionRule","train-then-freeze-for-independent-validation");
audit.requiredTolerance = tolerance;
audit.nextScope = "milestone-9.3-slope-compatible-wave-coordinates";
end

function reference = buildComparisonPrimitiveReference( ...
    c,direction,stationary,classified)
coarse = classified.finest.nestedRepresentation;
fine = struct( ...
    "admissibleBasis",c.N, ...
    "admissibleRanges",{c.layout.admissibleRanges}, ...
    "numberOfRawStateCoefficients",c.nX, ...
    "rawStateBlockSize",c.nXBlock, ...
    "numberOfFCoordinates",c.nF, ...
    "numberOfGCoordinates",c.nG, ...
    "numberOfHCoordinates",c.nH, ...
    "horizontalLayout",c.horizontalLayout);
rawMap = rawNestedMap(coarse,fine);
embeddedRaw = rawMap*coarse.admissibleBasis;
transfer = fine.admissibleBasis\embeddedRaw;
transferDefect = norm(fine.admissibleBasis*transfer-embeddedRaw,"fro") ...
    /max(norm(embeddedRaw,"fro"),realmin);

coarseInternalIndices = classified.classification.internalWaveIndices;
embeddedInternal = transfer ...
    *classified.finest.modeVectors(:,coarseInternalIndices);
H = direction.E;
J = direction.J;
energyFactor = chol(H);
Yg = energyFactor*stationary.fullBasis;
[completeQ,~] = qr(Yg);
Yd = completeQ(:,size(Yg,2)+1:end);
W = energyFactor\Yd;
K = energyFactor'\(1i*J/energyFactor);
Kdynamic = Yd'*K*Yd;
[D,T] = schur(Kdynamic,"complex");
frequencyComplex = diag(T);
[~,order] = sort(real(frequencyComplex));
D = D(:,order);
frequencyComplex = frequencyComplex(order);
modeVectors = W*D;
for iMode = 1:size(modeVectors,2)
    modeVectors(:,iMode) = modeVectors(:,iMode) ...
        /sqrt(real(modeVectors(:,iMode)'*H*modeVectors(:,iMode)));
end
embeddedProjector = energyProjector(embeddedInternal,H);
participation = real(sum(conj(modeVectors) ...
    .*(H*embeddedProjector*modeVectors),1)).';
numberOfInternal = size(embeddedInternal,2);
[~,participationOrder] = sort(participation,"descend");
internalIndices = sort(participationOrder(1:numberOfInternal));

reference = struct;
reference.scope = "primitive-modal-comparison";
reference.classification = struct( ...
    "internalWaveIndices",internalIndices);
reference.finest = struct( ...
    "modeVectors",modeVectors, ...
    "frequency",real(frequencyComplex), ...
    "energyMatrix",H);
reference.embeddedClassifiedInternalBasis = embeddedInternal;
reference.embeddedClassifiedInternalProjector = embeddedProjector;
reference.transfer = transfer;
reference.transferDefect = transferDefect;
reference.completePrimitiveDynamicDimension = size(modeVectors,2);
reference.completePrimitiveStateDimension = size(H,1);
reference.internalParticipation = participation;
end

function rawMap = rawNestedMap(coarse,fine)
rawMap = sparse(fine.numberOfRawStateCoefficients, ...
    coarse.numberOfRawStateCoefficients);
for iCoarse = 1:height(coarse.horizontalLayout)
    mode = [coarse.horizontalLayout.kMode(iCoarse) ...
        coarse.horizontalLayout.lMode(iCoarse)];
    iFine = find(fine.horizontalLayout.kMode == mode(1) ...
        &fine.horizontalLayout.lMode == mode(2),1);
    if isempty(iFine)
        error("WVTerrainEnergyGalerkin:NonNestedModalComparisonSupport", ...
            "Every classified primitive mode must be retained by the comparison support.")
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

function catalogs = buildVerticalCatalog(problem,c,flat,maximumCount,ratios,orders, ...
    waveCoordinateType,referenceSlope,baseCatalogs)
kappaValues = hypot(c.horizontalLayout.k,c.horizontalLayout.l);
uniqueKappa = unique(kappaValues(kappaValues > 0),"sorted");
nRobin = numel(ratios);
if isempty(baseCatalogs)
    catalogs = repmat(struct, nRobin,1);
else
    catalogs = baseCatalogs;
end
for iRobin = 1:nRobin
    robinLength = ratios(iRobin)*c.wvt.Lz;
    if isinf(ratios(iRobin))
        robinLength = Inf;
    end
    if isempty(baseCatalogs)
        high = cell(numel(uniqueKappa),1);
        low = cell(numel(uniqueKappa),1);
        providerDefect = zeros(numel(uniqueKappa),1);
        for iKappa = 1:numel(uniqueKappa)
            high{iKappa} = buildBoundaryCompleteVerticalReferenceFamily( ...
                problem,uniqueKappa(iKappa),maximumCount,robinLength, ...
                orders(end),c.quadrature.xi);
            low{iKappa} = buildBoundaryCompleteVerticalReferenceFamily( ...
                problem,uniqueKappa(iKappa),maximumCount,robinLength, ...
                orders(end-1),c.quadrature.xi);
            providerDefect(iKappa) = familyDifference( ...
                low{iKappa},high{iKappa},c.quadrature.weight);
        end
        diagnostics = [high{:}];
        diagnostics = [diagnostics.diagnostics];
        providerMaximumEndpointDefect = max([ ...
            diagnostics.waveEndpointDefect, ...
            diagnostics.surfaceRobinDefect, ...
            diagnostics.bottomRobinDefect, ...
            diagnostics.bottomAPVDefect, ...
            diagnostics.bottomValueDefect, ...
            diagnostics.bottomSurfaceValueDefect],[],"all");
        catalogs(iRobin).ratio = ratios(iRobin);
        catalogs(iRobin).robinLength = robinLength;
        catalogs(iRobin).kappa = uniqueKappa;
        catalogs(iRobin).families = high;
        catalogs(iRobin).providerDefect = max(providerDefect);
        catalogs(iRobin).providerMaximumEndpointDefect = ...
            providerMaximumEndpointDefect;
        catalogs(iRobin).negativeRobinModeCountDefect = ...
            max([diagnostics.negativeRobinModeCountDefect]);
    else
        high = catalogs(iRobin).families;
        providerMaximumEndpointDefect = ...
            catalogs(iRobin).providerMaximumEndpointDefect;
    end
    waveByHorizontal = cell(c.nK,1);
    waveDiagnostics = repmat(emptySlopeWaveDiagnostics,0,1);
    if contains(waveCoordinateType,"slope")
        isHydrostatic = startsWith(waveCoordinateType,"hydrostatic");
        for iHorizontal = 1:c.nK
            if hypot(c.horizontalLayout.k(iHorizontal), ...
                    c.horizontalLayout.l(iHorizontal)) == 0
                continue
            end
            waveByHorizontal{iHorizontal} = ...
                buildSlopeCompatibleWaveFamily(c,flat,iHorizontal, ...
                maximumCount,referenceSlope,isHydrostatic, ...
                high{find(abs(uniqueKappa-kappaValues(iHorizontal)) ...
                <=100*eps*max(kappaValues(iHorizontal),1),1)}.wave);
            waveDiagnostics(end+1,1) = ...
                waveByHorizontal{iHorizontal}.diagnostics; %#ok<AGROW>
        end
    end
    catalogs(iRobin).waveCoordinateType = waveCoordinateType;
    catalogs(iRobin).referenceSlope = referenceSlope;
    catalogs(iRobin).waveByHorizontal = waveByHorizontal;
    catalogs(iRobin).waveDiagnostics = waveDiagnostics;
    catalogs(iRobin).maximumEndpointDefect = ...
        providerMaximumEndpointDefect;
    if ~isempty(waveDiagnostics)
        catalogs(iRobin).maximumEndpointDefect = max( ...
            catalogs(iRobin).maximumEndpointDefect, ...
            max([waveDiagnostics.bottomResidual]));
    end
end
end

function diagnostics = emptySlopeWaveDiagnostics
diagnostics = struct("descriptorResidual",0, ...
    "momentumResidual",0,"bottomResidual",0, ...
    "backwardError",0,"flatFrequencyDefect",0, ...
    "maximumBottomDisplacement",0, ...
    "numberOfFiniteModes",0,"numberOfInfiniteModes",0, ...
    "homogeneousFiniteCount",0);
end

function family = buildSlopeCompatibleWaveFamily( ...
    c,flat,iHorizontal,numberOfModes,referenceSlope,isHydrostatic,flatWave)
% Construct boundary-only slope-compatible wave coordinates in the local
% primitive polynomial space. The periodic terrain forms are not changed.

nF = c.nF;
nH = c.nH;
nP = c.nP;
nXBlock = c.nXBlock;
nK = c.nK;
rawColumns = (iHorizontal-1)*nXBlock+(1:nXBlock);
admissible = c.layout.admissibleRanges{iHorizontal};
pressureColumns = c.nX+(iHorizontal-1)*nP+(1:nP);
uRows = (iHorizontal-1)*nF+(1:nF);
vRows = nK*nF+(iHorizontal-1)*nF+(1:nF);
wRows = 2*nK*nF+(iHorizontal-1)*nH+(1:nH);
etaRows = 2*nK*nF+nK*nH+(iHorizontal-1)*nH+(1:nH);
equationRows = [uRows vRows wRows etaRows];

Slocal = flat.matrices.S(equationRows, ...
    [rawColumns pressureColumns]);
Flocal = flat.matrices.F(equationRows,admissible);
Nlocal = c.N(rawColumns,admissible);
nA = numel(admissible);
M = [Slocal(:,1:nXBlock)*Nlocal zeros(numel(equationRows),nP)];
K = [Flocal -Slocal(:,nXBlock+(1:nP))];

k = c.horizontalLayout.k(iHorizontal);
l = c.horizontalLayout.l(iHorizontal);
sx = referenceSlope(1);
sy = referenceSlope(2);
etaBottomRaw = zeros(1,nXBlock);
etaBottomRaw(2*nF+c.nG+(1:nH)) = c.spaces.Hendpoint(1,:);
uBottomRaw = zeros(1,nXBlock);
vBottomRaw = zeros(1,nXBlock);
uBottomRaw(1:nF) = c.spaces.Fendpoint(1,:);
vBottomRaw(nF+(1:nF)) = c.spaces.Fendpoint(1,:);
B = etaBottomRaw*Nlocal;
R = (sx*uBottomRaw+sy*vBottomRaw)*Nlocal;
bottomRow = numel(equationRows);
M(bottomRow,:) = [B zeros(1,nP)];
K(bottomRow,:) = [R zeros(1,nP)];

localWRows = 2*nF+(1:nH);
if isHydrostatic
    M(localWRows,:) = 0;
end

referenceFrequency = max(abs(flatWave.omega),[],"all");
[vectors,lambda,homogeneous,scaling] = ...
    finiteHomogeneousDescriptorModes(K,M,referenceFrequency);
xCoordinates = vectors(1:nA,:);
H = flat.E(admissible,admissible);
xCoordinates = energyNormalizeColumns(xCoordinates,H);

flatCoordinates = flatWaveCoordinates(c,iHorizontal,flatWave,numberOfModes);
flatCoordinates = energyNormalizeColumns(flatCoordinates,H);
flatProjector = energyProjector(flatCoordinates,H);
participation = real(sum(conj(xCoordinates) ...
    .*(H*flatProjector*xCoordinates),1)).';
numberToRetain = min(2*numberOfModes,size(xCoordinates,2));
[~,order] = sort(participation,"descend");
selected = order(1:numberToRetain);
if numberToRetain < 2*numberOfModes
    error("WVTerrainEnergyGalerkin:InsufficientSlopeWaveModes", ...
        "The local slope descriptor returned %d finite wave candidates; %d are required.", ...
        numberToRetain,2*numberOfModes)
end
[~,frequencyOrder] = sort(real(1i*lambda(selected)));
selected = selected(frequencyOrder);
coordinates = xCoordinates(:,selected);
selectedVectors = vectors(:,selected);
selectedLambda = lambda(selected);
frequency = real(1i*selectedLambda);

descriptorResidual = columnRelativeResidual( ...
    K*selectedVectors-M*(selectedVectors.*selectedLambda.'), ...
    {K*selectedVectors,M*(selectedVectors.*selectedLambda.')});
momentumRows = 1:(numel(equationRows)-1);
momentumResidual = columnRelativeResidual( ...
    K(momentumRows,:)*selectedVectors ...
    -M(momentumRows,:)*(selectedVectors.*selectedLambda.'), ...
    {K(momentumRows,:)*selectedVectors, ...
    M(momentumRows,:)*(selectedVectors.*selectedLambda.')});
bottomResidual = abs(B*coordinates.*selectedLambda.'-R*coordinates).' ...
    ./max((abs(B*coordinates).*abs(selectedLambda.')).' ...
    +abs(R*coordinates).',realmin);
backwardError = descriptorBackwardError(K,M,selectedVectors,selectedLambda);

flatFrequencyDefect = 0;
if hypot(sx,sy) <=100*eps
    if isHydrostatic
        N2 = c.wvt.N2Function(0);
        if isscalar(N2)
            verticalWavenumber = (1:numberOfModes)'*pi/c.wvt.Lz;
            expectedMagnitude = sqrt(c.wvt.f^2 ...
                +N2*hypot(k,l)^2./verticalWavenumber.^2);
        else
            expectedMagnitude = [];
        end
    else
        expectedMagnitude = flatWave.omega(1:numberOfModes);
    end
    if ~isempty(expectedMagnitude)
        expected = sort([-expectedMagnitude(:);expectedMagnitude(:)]);
        flatFrequencyDefect = max(abs(sort(frequency)-expected)) ...
            /max(abs(expected),[],"all");
    end
end

family = struct;
family.coordinates = coordinates;
family.frequency = frequency;
family.lambda = selectedLambda;
family.embeddingDefect = zeros(2*numberOfModes,1);
family.participation = participation(selected);
family.homogeneous = homogeneous;
family.scaling = scaling;
family.diagnostics = struct( ...
    "descriptorResidual",max(descriptorResidual), ...
    "momentumResidual",max(momentumResidual), ...
    "bottomResidual",max(bottomResidual), ...
    "backwardError",max(backwardError), ...
    "flatFrequencyDefect",flatFrequencyDefect, ...
    "maximumBottomDisplacement",max(abs(B*coordinates),[],"all"), ...
    "numberOfFiniteModes",numel(lambda), ...
    "numberOfInfiniteModes",size(K,1)-numel(lambda), ...
    "homogeneousFiniteCount",nnz(homogeneous.isFinite));
end

function coordinates = flatWaveCoordinates(c,iHorizontal,wave,numberOfModes)
k = c.horizontalLayout.k(iHorizontal);
l = c.horizontalLayout.l(iHorizontal);
kappa = hypot(k,l);
admissible = c.layout.admissibleRanges{iHorizontal};
coordinates = zeros(numel(admissible),2*numberOfModes);
iColumn = 0;
for iMode = 1:numberOfModes
    for sigma = [-1 1]
        iColumn = iColumn+1;
        omega = wave.omega(iMode);
        h = wave.h(iMode);
        F = wave.F(:,iMode);
        G = wave.G(:,iMode);
        state = struct( ...
            "u",(k*omega-sigma*1i*c.wvt.f*l)/(omega*kappa)*F, ...
            "v",(l*omega+sigma*1i*c.wvt.f*k)/(omega*kappa)*F, ...
            "w",-1i*kappa*h*G, ...
            "eta",-sigma*kappa*h/omega*G);
        [globalCoordinate,~] = embedState(c,iHorizontal,state);
        coordinates(:,iColumn) = globalCoordinate(admissible);
    end
end
end

function values = energyNormalizeColumns(values,H)
energy = real(sum(conj(values).*(H*values),1));
valid = isfinite(energy) & energy > max(size(H))*eps(norm(H,2));
if ~all(valid)
    error("WVTerrainEnergyGalerkin:InvalidSlopeWaveEnergy", ...
        "Every finite local descriptor mode must have positive finite physical energy.")
end
values = values./sqrt(energy);
end

function [vectors,lambda,homogeneous,scaling] = ...
    finiteHomogeneousDescriptorModes(K,M,referenceFrequency)
frequencyScale = max(abs(referenceFrequency),eps);
Ks = K/frequencyScale;
Ms = M;
columnNorm = sqrt(sum(abs(Ks).^2+abs(Ms).^2,1));
columnScale = 1./max(columnNorm,sqrt(realmin));
Ks = Ks.*columnScale;
Ms = Ms.*columnScale;
rowNorm = sqrt(sum(abs(Ks).^2+abs(Ms).^2,2));
rowScale = 1./max(rowNorm,sqrt(realmin));
Ks = rowScale.*Ks;
Ms = rowScale.*Ms;

[AA,BB] = qz(Ks,Ms,"complex");
alpha = diag(AA);
beta = diag(BB);
homogeneousTolerance = max(size(K))*eps*100;
isFinite = abs(beta) > homogeneousTolerance ...
    .*max(abs(alpha)+abs(beta),realmin);

[scaledVectors,scaledLambda] = eig(Ks,Ms,"vector");
finiteEigenvector = isfinite(scaledLambda) & ~isnan(scaledLambda);
vectors = columnScale.'.*scaledVectors(:,finiteEigenvector);
lambda = frequencyScale*scaledLambda(finiteEigenvector);
vectors = vectors./max(vecnorm(vectors),sqrt(realmin));
homogeneous = struct("alpha",alpha,"beta",beta, ...
    "isFinite",isFinite,"tolerance",homogeneousTolerance);
scaling = struct("frequency",frequencyScale, ...
    "row",rowScale,"column",columnScale.');
if nnz(isFinite) ~= numel(lambda)
    error("WVTerrainEnergyGalerkin:HomogeneousDescriptorCountMismatch", ...
        "Generalized Schur and eigenvector finite counts disagree (%d versus %d).", ...
        nnz(isFinite),numel(lambda))
end
end

function errorValue = descriptorBackwardError(K,M,vectors,lambda)
residual = K*vectors-M*(vectors.*lambda.');
scale = (norm(K,2)+abs(lambda).'*norm(M,2)) ...
    .*max(vecnorm(vectors,2,1).',realmin);
errorValue = vecnorm(residual,2,1).'./max(scale,realmin);
end

function value = familyDifference(first,second,weight)
fields = {first.wave.F,second.wave.F; ...
    first.wave.G,second.wave.G; ...
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
if size(first,2) ~= size(second,2)
    value = Inf;
    return
end
value = norm(first-second*(second'*(weight.*first)),2);
end

function values = weightedOrthonormalize(values,weight)
gram = values'*(weight.*values);
[R,flag] = chol(gram);
if flag ~= 0
    [left,singularMatrix,~] = svd(sqrt(weight).*values,"econ");
    singularValues = diag(singularMatrix);
    rankValue = nnz(singularValues ...
        >max(size(values))*eps(max(singularValues))*100);
    values = left(:,1:rankValue)./sqrt(weight);
else
    values = values/R;
end
end

function basis = buildModalCoordinates(c,direction,stationary,catalog,count)
nState = size(direction.E,1);
columns = zeros(nState,0);
labels = strings(0,1);
embeddingDefect = zeros(0,1);
for iHorizontal = 1:c.nK
    k = c.horizontalLayout.k(iHorizontal);
    l = c.horizontalLayout.l(iHorizontal);
    kappa = hypot(k,l);
    admissible = c.layout.admissibleRanges{iHorizontal};
    if kappa == 0
        meanColumns = speye(nState);
        meanColumns = meanColumns(:,admissible);
        columns = [columns full(meanColumns)]; %#ok<AGROW>
        labels = [labels;repmat("mean-sector",numel(admissible),1)]; %#ok<AGROW>
        embeddingDefect = [embeddingDefect;zeros(numel(admissible),1)]; %#ok<AGROW>
        continue
    end
    familyIndex = find(abs(catalog.kappa-kappa) ...
        <=100*eps*max(kappa,1),1);
    family = catalog.families{familyIndex};
    if contains(catalog.waveCoordinateType,"slope")
        slopeWave = catalog.waveByHorizontal{iHorizontal};
        for iMode = 1:(2*count)
            coordinate = zeros(nState,1);
            coordinate(admissible) = slopeWave.coordinates(:,iMode);
            columns(:,end+1) = coordinate; %#ok<AGROW>
            labels(end+1,1) = catalog.waveCoordinateType; %#ok<AGROW>
            embeddingDefect(end+1,1) = ...
                slopeWave.embeddingDefect(iMode); %#ok<AGROW>
        end
    else
        for iMode = 1:count
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
                [coordinate,defect] = embedState(c,iHorizontal,state);
                columns(:,end+1) = coordinate; %#ok<AGROW>
                labels(end+1,1) = "wave-"+string(sigma); %#ok<AGROW>
                embeddingDefect(end+1,1) = defect; %#ok<AGROW>
            end
        end
    end
    for iMode = 1:count
        F = family.robin.F(:,iMode);
        state = struct("u",-1i*l*F,"v",1i*k*F, ...
            "w",zeros(size(F)),"eta",family.robin.eta(:,iMode));
        [coordinate,defect] = embedState(c,iHorizontal,state);
        columns(:,end+1) = coordinate; %#ok<AGROW>
        labels(end+1,1) = "robin-geostrophic"; %#ok<AGROW>
        embeddingDefect(end+1,1) = defect; %#ok<AGROW>
    end
    state = struct("u",-1i*l*family.bottom.psi, ...
        "v",1i*k*family.bottom.psi, ...
        "w",zeros(size(family.bottom.psi)), ...
        "eta",family.bottom.eta);
    [coordinate,defect] = embedState(c,iHorizontal,state);
    columns(:,end+1) = coordinate; %#ok<AGROW>
    labels(end+1,1) = "zero-apv-bottom"; %#ok<AGROW>
    embeddingDefect(end+1,1) = defect; %#ok<AGROW>
end

H = direction.E;
energyFactor = chol(H);
columnScale = vecnorm(energyFactor*columns,2,1).';
if any(~isfinite(columnScale)) || any(columnScale <= 0)
    error("WVTerrainEnergyGalerkin:ModalReferenceEnergyFailure", ...
        "Every boundary-complete reference column must have positive finite physical energy.")
end
columns = columns./columnScale.';
[~,coordinateFactor] = qr(energyFactor*columns,0);
normalizedDiagonal = abs(diag(coordinateFactor)) ...
    /max(abs(diag(coordinateFactor)));
if min(normalizedDiagonal) ...
        <=max(size(columns))*eps*100
    error("WVTerrainEnergyGalerkin:ModalReferenceRankFailure", ...
        "The count-%d, ell_b/D=%g reference has minimum normalized energy QR diagonal %.3g.", ...
        count,catalog.ratio,min(normalizedDiagonal))
end
columns = columns/coordinateFactor;
stationaryResidual = stationary.fullBasis ...
    -columns*(columns'*H*stationary.fullBasis);
[stationaryCompletion,completionDiagnostics] = ...
    independentEnergyColumns(stationaryResidual,H);
allColumns = [columns stationaryCompletion];
[~,allFactor] = qr(energyFactor*allColumns,0);
allDiagonal = abs(diag(allFactor))/max(abs(diag(allFactor)));
if min(allDiagonal) <=max(size(allColumns))*eps*100
    error("WVTerrainEnergyGalerkin:ModalCompletionRankFailure", ...
        "The stationary completion is not independent of the modal reference columns.")
end
allColumns = allColumns/allFactor;

basis = struct;
basis.columns = allColumns;
basis.referenceColumns = columns;
basis.stationaryCompletion = stationaryCompletion;
basis.labels = labels;
basis.maximumEmbeddingDefect = max(embeddingDefect,[],"all");
basis.embeddingDefect = embeddingDefect;
basis.numberOfNominalModalCoordinates = size(columns,2);
basis.numberOfStationaryCompletionCoordinates = ...
    size(stationaryCompletion,2);
basis.numberOfCoordinates = size(allColumns,2);
basis.primitiveDimension = nState;
basis.actualCompressionFactor = nState/size(allColumns,2);
basis.referenceGramConditionNumber = cond(coordinateFactor)^2;
basis.completedGramConditionNumber = cond(allFactor)^2;
basis.minimumReferenceQRDiagonal = min(normalizedDiagonal);
basis.minimumCompletedQRDiagonal = min(allDiagonal);
basis.stationaryCompletionDiagnostics = completionDiagnostics;
basis.providerDefect = catalog.providerDefect;
basis.maximumEndpointDefect = catalog.maximumEndpointDefect;
basis.negativeRobinModeCountDefect = ...
    catalog.negativeRobinModeCountDefect;
basis.waveCoordinateType = catalog.waveCoordinateType;
basis.referenceSlope = catalog.referenceSlope;
if isempty(catalog.waveDiagnostics)
    basis.maximumWaveDescriptorResidual = 0;
    basis.maximumWaveMomentumResidual = 0;
    basis.maximumWaveBottomResidual = 0;
    basis.maximumWaveBackwardError = 0;
    basis.maximumWaveFlatFrequencyDefect = 0;
    basis.maximumWaveBottomDisplacement = 0;
else
    basis.maximumWaveDescriptorResidual = max( ...
        [catalog.waveDiagnostics.descriptorResidual]);
    basis.maximumWaveMomentumResidual = max( ...
        [catalog.waveDiagnostics.momentumResidual]);
    basis.maximumWaveBottomResidual = max( ...
        [catalog.waveDiagnostics.bottomResidual]);
    basis.maximumWaveBackwardError = max( ...
        [catalog.waveDiagnostics.backwardError]);
    basis.maximumWaveFlatFrequencyDefect = max( ...
        [catalog.waveDiagnostics.flatFrequencyDefect]);
    basis.maximumWaveBottomDisplacement = max( ...
        [catalog.waveDiagnostics.maximumBottomDisplacement]);
end
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
localCoordinate = localBasis\rawBlock;
represented = localBasis*localCoordinate;
coordinate = zeros(size(c.N,2),1);
coordinate(admissible) = localCoordinate;
defect = norm(represented-rawBlock) ...
    /max(norm(rawBlock),realmin);
end

function coefficient = weightedFit(reconstruction,values,weight)
coefficient = (reconstruction'*(weight.*reconstruction)) ...
    \(reconstruction'*(weight.*values));
end

function [values,diagnostics] = independentEnergyColumns(values,H)
if isempty(values)
    diagnostics = struct("rank",0,"singularValues",zeros(0,1), ...
        "tolerance",0);
    return
end
R = chol(H);
[~,singularMatrix,right] = svd(R*values,"econ");
singularValues = diag(singularMatrix);
tolerance = max(size(values))*eps(max(singularValues))*100;
rankValue = nnz(singularValues > tolerance);
values = values*right(:,1:rankValue) ...
    /singularMatrix(1:rankValue,1:rankValue);
diagnostics = struct("rank",rankValue, ...
    "singularValues",singularValues,"tolerance",tolerance);
end

function result = reducedAudit(c,direction,stationary,basis,reference)
V = basis.columns;
primitiveEnergyFactor = chol(direction.E);
energyCoordinates = primitiveEnergyFactor*V;
H = energyCoordinates'*energyCoordinates;
J = V'*direction.J*V;
[R,flag] = chol(H);
if flag ~= 0
    error("WVTerrainEnergyGalerkin:ReducedModalEnergyNotPositive", ...
        "The reduced physical-energy Gram matrix is not positive definite.")
end
K = R'\(1i*J/R);
[D,T] = schur(K,"complex");
frequencyComplex = diag(T);
[~,order] = sort(real(frequencyComplex));
D = D(:,order);
frequencyComplex = frequencyComplex(order);
C = V*(R\D);
frequency = real(frequencyComplex);
for iMode = 1:size(C,2)
    C(:,iMode) = C(:,iMode) ...
        /sqrt(real(C(:,iMode)'*direction.E*C(:,iMode)));
end

G = stationary.trustedBasis;
Gcoordinate = H\(V'*direction.E*G);
stationaryRepresentation = V*Gcoordinate-G;
stationaryRepresentationDefect = energyNorm( ...
    stationaryRepresentation,direction.E) ...
    /max(energyNorm(G,direction.E),realmin);
stationaryRowDefect = norm(Gcoordinate'*J,"fro") ...
    /max(norm(Gcoordinate,"fro")*norm(J,"fro"),realmin);
L = H\J;
bottomResidual = direction.B*V*L-direction.R*V;
bottomEvolutionDefect = norm(bottomResidual,"fro") ...
    /max(norm(direction.B*V*L,"fro") ...
    +norm(direction.R*V,"fro"),realmin);

referenceIndices = reference.classification.internalWaveIndices;
referenceModes = reference.finest.modeVectors(:,referenceIndices);
referenceFrequency = reference.finest.frequency(referenceIndices);
numberOfInternal = size(referenceModes,2);
referenceProjector = energyProjector( ...
    referenceModes,direction.E);
participation = real(sum(conj(C).*(direction.E* ...
    referenceProjector*C),1)).';
[~,participationOrder] = sort(participation,"descend");
internalIndices = sort(participationOrder(1:min( ...
    numberOfInternal,numel(participationOrder))));
internalModes = C(:,internalIndices);
if size(internalModes,2) == size(referenceModes,2)
    internalProjectorDefect = maximumPrincipalSine( ...
        internalModes,referenceModes,direction.E);
    internalFrequencyDefect = sortedFrequencyDefect( ...
        frequency(internalIndices),referenceFrequency);
else
    internalProjectorDefect = Inf;
    internalFrequencyDefect = Inf;
end

trustedRows = find(repmat(c.layout.trustedHorizontalModes,c.nZ,1));
weightedQ = sqrt(c.apvVerticalWeight(trustedRows)) ...
    .*direction.QProjected(trustedRows,:);
weightedQScale = max(norm(weightedQ/chol(direction.E),2),realmin);
apvDefect = vecnorm(weightedQ*internalModes,2,1).' ...
    /weightedQScale;
strong = strongModeDiagnostics(c,direction,internalModes, ...
    frequency(internalIndices));

result = struct;
result.modalBasis = basis;
result.modeVectors = C;
result.frequency = frequency;
result.frequencyComplex = frequencyComplex;
result.internalModeIndices = internalIndices;
result.internalWaveProjector = energyProjector( ...
    internalModes,direction.E);
result.referenceInternalWaveProjector = referenceProjector;
result.internalProjectorDefect = internalProjectorDefect;
result.internalFrequencyDefect = internalFrequencyDefect;
result.internalAPVDefect = apvDefect;
if isempty(apvDefect)
    result.maximumInternalAPVDefect = Inf;
else
    result.maximumInternalAPVDefect = max(apvDefect,[],"all");
end
result.internalStrong = strong;
result.maximumInternalStrongResidual = strong.maximumResidual;
result.stationaryRepresentationDefect = ...
    stationaryRepresentationDefect;
result.stationaryRowDefect = stationaryRowDefect;
result.bottomEvolutionDefect = bottomEvolutionDefect;
result.energyHermitianDefect = norm(H-H',"fro") ...
    /max(norm(H,"fro"),realmin);
result.exchangeSkewHermitianDefect = norm(J+J',"fro") ...
    /max(norm(J,"fro"),realmin);
result.modalEnergyOrthogonalityDefect = ...
    norm(C'*direction.E*C-eye(size(C,2)),"fro") ...
    /sqrt(size(C,2));
result.maximumEmbeddingDefect = basis.maximumEmbeddingDefect;
result.maximumProviderDefect = basis.providerDefect;
result.maximumEndpointDefect = basis.maximumEndpointDefect;
result.negativeRobinModeCountDefect = ...
    basis.negativeRobinModeCountDefect;
result.waveCoordinateType = basis.waveCoordinateType;
result.referenceSlope = basis.referenceSlope;
result.maximumWaveDescriptorResidual = ...
    basis.maximumWaveDescriptorResidual;
result.maximumWaveMomentumResidual = ...
    basis.maximumWaveMomentumResidual;
result.maximumWaveBottomResidual = ...
    basis.maximumWaveBottomResidual;
result.maximumWaveBackwardError = ...
    basis.maximumWaveBackwardError;
result.maximumWaveFlatFrequencyDefect = ...
    basis.maximumWaveFlatFrequencyDefect;
result.maximumWaveBottomDisplacement = ...
    basis.maximumWaveBottomDisplacement;
result.primitiveDimension = basis.primitiveDimension;
result.nominalModalDimension = ...
    basis.numberOfNominalModalCoordinates;
result.stationaryCompletionDimension = ...
    basis.numberOfStationaryCompletionCoordinates;
result.actualModalDimension = basis.numberOfCoordinates;
result.actualCompressionFactor = basis.actualCompressionFactor;
result.referenceGramConditionNumber = ...
    basis.referenceGramConditionNumber;
result.completedGramConditionNumber = ...
    basis.completedGramConditionNumber;
result.minimumReferenceQRDiagonal = ...
    basis.minimumReferenceQRDiagonal;
result.minimumCompletedQRDiagonal = ...
    basis.minimumCompletedQRDiagonal;
result.allModeParticipationInReferenceInternalSpace = participation;
end

function value = energyProjector(basis,H)
if isempty(basis)
    value = zeros(size(H));
else
    value = basis*((basis'*H*basis)\(basis'*H));
end
end

function value = maximumPrincipalSine(first,second,H)
if size(first,2) ~= size(second,2)
    value = Inf;
    return
end
first = orthonormalize(first,H);
second = orthonormalize(second,H);
R = chol(H);
value = norm(R*(first-second*(second'*H*first)),2);
end

function values = orthonormalize(values,H)
R = chol(H);
[~,singularMatrix,right] = svd(R*values,"econ");
values = values*right/singularMatrix;
end

function value = sortedFrequencyDefect(first,second)
first = sort(first(:));
second = sort(second(:));
if isempty(first) && isempty(second)
    value = 0;
    return
end
if numel(first) ~= numel(second)
    value = Inf;
    return
end
value = max(abs(first-second)) ...
    /max([abs(first);abs(second);realmin]);
end

function value = energyNorm(columns,H)
value = sqrt(max(real(trace(columns'*H*columns)),0));
end

function diagnostics = countConvergence(results)
results = results(:);
diagnostics.projectorDefect = cellfun( ...
    @(result)result.internalProjectorDefect,results);
diagnostics.frequencyDefect = cellfun( ...
    @(result)result.internalFrequencyDefect,results);
diagnostics.actualDimension = cellfun( ...
    @(result)result.actualModalDimension,results);
diagnostics.compressionFactor = cellfun( ...
    @(result)result.actualCompressionFactor,results);
end

function diagnostics = paddingConvergence(results)
results = results(:);
reference = results{1};
defect = zeros(numel(results),1);
for iResult = 2:numel(results)
    defect(iResult) = projectorDifference( ...
        reference.internalWaveProjector, ...
        results{iResult}.internalWaveProjector);
end
diagnostics.internalProjectorDefect = defect;
diagnostics.maximumInternalProjectorDefect = max(defect);
end

function diagnostics = robinConvergence(results,referenceIndex)
results = results(:);
reference = results{referenceIndex};
defect = zeros(numel(results),1);
for iResult = 1:numel(results)
    defect(iResult) = projectorDifference( ...
        reference.internalWaveProjector, ...
        results{iResult}.internalWaveProjector);
end
diagnostics.internalProjectorDefect = defect;
diagnostics.maximumInternalProjectorDefect = max(defect);
end

function value = projectorDifference(first,second)
value = norm(first-second,"fro") ...
    /max(norm(first,"fro"),realmin);
end

function diagnostics = strongModeDiagnostics(c,direction,C,frequency)
if isempty(C)
    diagnostics = struct("maximumByMode",zeros(0,1), ...
        "maximumResidual",0);
    return
end
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
maximum = max([descriptorDefect momentumDefect continuityDefect ...
    primitiveDefect gaugeDefect],[],2);
diagnostics = struct( ...
    "descriptorDefect",descriptorDefect, ...
    "momentumWeakDefect",momentumDefect, ...
    "continuityWeakDefect",continuityDefect, ...
    "primitiveGeneratorDefect",primitiveDefect, ...
    "pressureGaugeDefect",gaugeDefect, ...
    "maximumByMode",maximum, ...
    "maximumResidual",max(maximum,[],"all"));
end

function defect = columnRelativeResidual(residual,terms)
scale = zeros(size(residual,2),1);
for iTerm = 1:numel(terms)
    scale = scale+vecnorm(terms{iTerm},2,1).';
end
defect = vecnorm(residual,2,1).'./max(scale,realmin);
end

function layout = retainedLayout(problem,bounds)
source = problem.horizontalLayout;
mask = abs(source.kMode) <= bounds(1) ...
    & abs(source.lMode) <= bounds(2);
layout = source(mask,:);
expected = (2*bounds(1)+1)*(2*bounds(2)+1);
if height(layout) ~= expected
    error("WVTerrainEnergyGalerkin:UnavailableModalCompressionSupport", ...
        "The originating transform does not retain the complete signed rectangle [%d %d].", ...
        bounds(1),bounds(2))
end
end
