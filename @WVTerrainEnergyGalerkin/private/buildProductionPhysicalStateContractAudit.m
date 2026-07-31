function audit = buildProductionPhysicalStateContractAudit(problem,trustedBounds,supportBounds,waveIndices,apvIndices,degrees,paddingFactors,quadratureOrder)
% Build the Milestone-10.1 production physical-state contract audit.

production = productionLayout(problem,trustedBounds,waveIndices,apvIndices);
configurations = refinementConfigurations(degrees,supportBounds,paddingFactors);
details = cell(numel(configurations),1);
[nativeCache,quadratureOrders] = nativeCacheForDegrees( ...
    problem,production,degrees,quadratureOrder);
for iConfiguration = 1:numel(configurations)
    configuration = configurations(iConfiguration);
    degree = configuration.degree;
    support = configuration.support;
    padding = configuration.padding;
    layout = retainedLayout(problem,support);
    iDegree = find(degrees == degree,1);
    order = quadratureOrders(iDegree);
    [primitive,context] = buildGlobalSmallTerrainPrimitiveAudit(problem,degree,order,1e-3, ...
        horizontalLayout=layout,paddingFactor=padding, ...
        trustedModeBounds=trustedBounds,rejectTerrainNyquist=true);
    detail = contractAtResolution(problem,production,context,primitive.flatReference,nativeCache{iDegree});
    detail.degree = degree;
    detail.supportModeBounds = support;
    detail.paddingFactor = padding;
    detail.quadratureOrder = order;
    details{iConfiguration} = detail;
end

details = vertcat(details{:});
vertical = details([configurations.kind] == "vertical");
support = details([configurations.kind] == "support");
padding = details([configurations.kind] == "padding");
primary = vertical(end);
tolerance = struct("roundTrip",1e-12,"energy",1e-12, ...
    "endpoint",1e-12,"zeroAPV",5e-12,"embedding",1e-8, ...
    "support",1e-10,"padding",1e-10,"conjugacy",1e-11, ...
    "nativeMatch",5e-10);
verticalDefect = [vertical.maximumPhysicalEmbeddingDefect].';
verticalConverges = verticalDefect(end) <= tolerance.embedding ...
    && verticalDefect(end) <= verticalDefect(1)/4 ...
    && all(diff(verticalDefect) <= 100*eps(max(verticalDefect)));
supportDefect = maximumReferenceGramDifference(support,support(end));
paddingDefect = maximumReferenceGramDifference(padding,padding(end));
passes = primary.roundTripDefect <= tolerance.roundTrip ...
    && primary.energyNormalizationDefect <= tolerance.energy ...
    && primary.maximumEndpointDefect <= tolerance.endpoint ...
    && primary.maximumZeroAPVDefect <= tolerance.zeroAPV ...
    && primary.conjugacyDefect <= tolerance.conjugacy ...
    && primary.maximumNativeMatchDefect <= tolerance.nativeMatch ...
    && verticalConverges && supportDefect <= tolerance.support ...
    && paddingDefect <= tolerance.padding;

if passes
    status = "complete-production-state-contract";
    diagnosis = "Every declared production coordinate has a native physical precursor, a complete conjugacy assignment, and a converged primitive-oracle image.";
else
    status = "production-state-blocker";
    diagnosis = blockerDiagnosis(primary,verticalConverges,supportDefect,paddingDefect,tolerance);
end

audit = struct;
audit.scope = "milestone-10.1-production-physical-state-contract";
audit.status = status;
audit.isCompatible = passes;
audit.diagnosis = diagnosis;
audit.trustedModeBounds = trustedBounds;
audit.supportModeBounds = supportBounds;
audit.waveModeIndices = waveIndices;
audit.apvModeIndices = apvIndices;
audit.primitivePolynomialDegrees = degrees;
audit.paddingFactors = paddingFactors;
audit.productionHorizontalLayout = production.horizontalLayout;
audit.productionLayout = production.layout;
audit.conjugateCoordinateIndex = production.conjugateIndex;
audit.degreeOfFreedomBudget = production.degreeOfFreedomBudget;
audit.numberOfProductionCoordinates = height(production.layout);
audit.numberOfIndependentRealDegreesOfFreedom = production.numberOfIndependentRealDegreesOfFreedom;
audit.guardLayouts = guardLayouts(problem,trustedBounds,supportBounds);
audit.details = details;
audit.verticalConvergence = struct("degrees",degrees, ...
    "maximumPhysicalEmbeddingDefect",verticalDefect, ...
    "converges",verticalConverges);
audit.supportConvergence = struct("bounds",supportBounds, ...
    "maximumReferenceGramDefect",supportDefect);
audit.paddingConvergence = struct("factors",paddingFactors, ...
    "maximumReferenceGramDefect",paddingDefect);
audit.primary = primary;
audit.requiredTolerance = tolerance;
audit.nextScope = "milestone-10.2-complete-internal-wave-coverage";
end

function [cache,orders] = nativeCacheForDegrees(problem,production,degrees,requestedOrder)
cache = cell(numel(degrees),1);
zByDegree = cell(numel(degrees),1);
orders = zeros(numel(degrees),1);
for iDegree = 1:numel(degrees)
    if isempty(requestedOrder)
        orders(iDegree) = max(2*degrees(iDegree)+7,18);
    else
        orders(iDegree) = requestedOrder;
    end
    zByDegree{iDegree} = gaussLegendrePoints(orders(iDegree), ...
        problem.originatingTransform.Lz);
end
zCommon = unique([vertcat(zByDegree{:}); ...
    problem.originatingTransform.z(:)],"sorted");
nativeCommon = nativePhysicalColumns(problem,production,zCommon);
names = ["u","v","w","eta"];
for iDegree = 1:numel(degrees)
    [isPresent,index] = ismember(zByDegree{iDegree},zCommon);
    if ~all(isPresent)
        error("WVTerrainEnergyGalerkin:ProductionQuadratureMismatch", ...
            "The common native-mode evaluation grid does not contain every oracle quadrature point.")
    end
    for name = names
        cache{iDegree}.(name) = nativeCommon.(name)(index,:);
    end
    cache{iDegree}.waveVortexMatchDefect = ...
        nativeCommon.waveVortexMatchDefect;
    cache{iDegree}.bottomInversionMatchDefect = ...
        nativeCommon.bottomInversionMatchDefect;
end
end

function z = gaussLegendrePoints(order,D)
index = (1:order-1)';
offDiagonal = index./sqrt(4*index.^2-1);
[~,values] = eig(diag(offDiagonal,1)+diag(offDiagonal,-1),"vector");
r = sort(values);
z = D*(r-1)/2;
end

function configurations = refinementConfigurations(degrees,supportBounds,paddingFactors)
prototype = struct("degree",0,"support",[0 0],"padding",0,"kind","");
configurations = repmat(prototype,0,1);
for degree = reshape(degrees,1,[])
    configurations(end+1,1) = struct("degree",degree, ...
        "support",supportBounds(1,:),"padding",paddingFactors(1), ...
        "kind","vertical"); %#ok<AGROW>
end
for iSupport = 1:size(supportBounds,1)
    if isequal(supportBounds(iSupport,:),supportBounds(end,:))
        continue
    end
    configurations(end+1,1) = struct("degree",degrees(1), ...
        "support",supportBounds(iSupport,:),"padding",paddingFactors(1), ...
        "kind","support"); %#ok<AGROW>
end
configurations(end+1,1) = struct("degree",degrees(1), ...
    "support",supportBounds(end,:),"padding",paddingFactors(1), ...
    "kind","support");
for iPadding = 2:numel(paddingFactors)
    configurations(end+1,1) = struct("degree",degrees(1), ...
        "support",supportBounds(end,:),"padding",paddingFactors(iPadding), ...
        "kind","padding"); %#ok<AGROW>
end
configurations(end+1,1) = struct("degree",degrees(1), ...
    "support",supportBounds(end,:),"padding",paddingFactors(1), ...
    "kind","padding");
end

function production = productionLayout(problem,bounds,waveIndices,apvIndices)
horizontal = retainedLayout(problem,bounds);
horizontalIndex = zeros(0,1);
family = strings(0,1);
component = strings(0,1);
j = zeros(0,1);
branchSign = zeros(0,1);
units = strings(0,1);
nativeHorizontalIndex = zeros(0,1);
nativeCoefficientIndex = NaN(0,1);
isPrimary = false(0,1);
wvt = problem.originatingTransform;
for iK = 1:height(horizontal)
    kappa = hypot(horizontal.k(iK),horizontal.l(iK));
    if kappa > 0
        families = [repmat("internal-wave",2*numel(waveIndices),1); ...
            repmat("apv-balanced",numel(apvIndices),1);"bottom-inversion"];
        components = [repmat(["Ap";"Am"],numel(waveIndices),1); ...
            repmat("A0",numel(apvIndices),1);"etaB"];
        jValues = [repelem(waveIndices(:),2,1);apvIndices;NaN];
        signs = [repmat([1;-1],numel(waveIndices),1);zeros(numel(apvIndices)+1,1)];
        unitValues = [repmat("m/s",2*numel(waveIndices),1); ...
            repmat("m^2/s",numel(apvIndices),1);"m"];
    else
        inertialIndices = [0;waveIndices];
        positiveMDA = apvIndices(apvIndices > 0);
        families = [repmat("inertial",2*numel(inertialIndices),1); ...
            repmat("mda",numel(positiveMDA),1);"mean-bottom-compatible"];
        components = [repmat(["Ap";"Am"],numel(inertialIndices),1); ...
            repmat("A0",numel(positiveMDA),1);"etaB"];
        jValues = [repelem(inertialIndices(:),2,1);positiveMDA;0];
        signs = [repmat([1;-1],numel(inertialIndices),1);zeros(numel(positiveMDA)+1,1)];
        unitValues = [repmat("m/s",2*numel(inertialIndices),1); ...
            repmat("m",numel(positiveMDA)+1,1)];
    end
    n = numel(families);
    horizontalIndex(end+(1:n),1) = iK;
    family(end+(1:n),1) = families;
    component(end+(1:n),1) = components;
    j(end+(1:n),1) = jValues;
    branchSign(end+(1:n),1) = signs;
    units(end+(1:n),1) = unitValues;
    nativeHorizontalIndex(end+(1:n),1) = horizontal.nativeIndex(iK);
    isPrimary(end+(1:n),1) = horizontal.isPrimary(iK);
    for iLocal = 1:n
        if components(iLocal) ~= "etaB"
            iJ = find(wvt.j == jValues(iLocal),1);
            nativeCoefficientIndex(end+1,1) = sub2ind(wvt.spectralMatrixSize,iJ,horizontal.nativeIndex(iK)); %#ok<AGROW>
        else
            nativeCoefficientIndex(end+1,1) = NaN; %#ok<AGROW>
        end
    end
end
index = (1:numel(family)).';
layout = table(index,horizontalIndex,family,component,j,branchSign,units, ...
    nativeHorizontalIndex,nativeCoefficientIndex,isPrimary);
conjugateIndex = productionConjugateMap(layout,horizontal);
layout.conjugateIndex = conjugateIndex;
counts = groupsummary(layout,"family");
counts.Properties.VariableNames{2} = 'numberOfComplexCoordinates';
counts.numberOfIndependentRealDegreesOfFreedom = counts.numberOfComplexCoordinates;
production = struct("horizontalLayout",horizontal,"layout",layout, ...
    "conjugateIndex",conjugateIndex, ...
    "degreeOfFreedomBudget",counts, ...
    "numberOfIndependentRealDegreesOfFreedom",height(layout));
end

function map = productionConjugateMap(layout,horizontal)
map = zeros(height(layout),1);
for row = 1:height(layout)
    iK = layout.horizontalIndex(row);
    partnerK = find(horizontal.kMode == -horizontal.kMode(iK) ...
        & horizontal.lMode == -horizontal.lMode(iK),1);
    partnerComponent = layout.component(row);
    if partnerComponent == "Ap"
        partnerComponent = "Am";
    elseif partnerComponent == "Am"
        partnerComponent = "Ap";
    end
    candidates = find(layout.horizontalIndex == partnerK ...
        & layout.component == partnerComponent ...
        & layout.family == layout.family(row));
    partner = candidates(layout.j(candidates) == layout.j(row) ...
        | (isnan(layout.j(candidates)) & isnan(layout.j(row))));
    if numel(partner) ~= 1
        error("WVTerrainEnergyGalerkin:IncompleteProductionConjugacy", ...
            "Production coordinate %d has %d conjugate candidates.",row,numel(partner))
    end
    map(row) = partner;
end
if any(map(map) ~= (1:height(layout)).')
    error("WVTerrainEnergyGalerkin:InvalidProductionConjugacy", ...
        "The production conjugate map is not an involution.")
end
end

function native = nativePhysicalColumns(problem,production,z)
wvt = problem.originatingTransform;
modes = nativeWVTModeShapes(wvt,production.horizontalLayout,z);
nCoordinate = height(production.layout);
nZ = numel(z);
names = ["u","v","w","eta"];
for name = names
    native.(name) = zeros(nZ,nCoordinate);
end
[bottom,bottomMatchDefect] = ...
    bottomProfiles(problem,production.horizontalLayout,z);
for row = find(production.layout.isPrimary).'
    state = nativeState(wvt,modes,production,row,bottom);
    for name = names
        native.(name)(:,row) = state.(name);
    end
end
for row = find(~production.layout.isPrimary).'
    partner = production.conjugateIndex(row);
    for name = names
        native.(name)(:,row) = conj(native.(name)(:,partner));
    end
end
native.waveVortexMatchDefect = modes.maximumMatchDefect;
native.bottomInversionMatchDefect = bottomMatchDefect;
end

function modes = nativeWVTModeShapes(wvt,horizontal,z)
sourceRows = sourceGridRows(z,wvt.z);
sampleZ = linspace(-wvt.Lz,0,max(17,wvt.Nz)).';
sampleN2 = wvt.N2Function(sampleZ);
isConstant = max(abs(sampleN2-mean(sampleN2))) ...
    <= 100*eps*max(abs(mean(sampleN2)),1);
if isConstant
    verticalModes = InternalModesConstantStratification( ...
        N0=sqrt(mean(sampleN2)),zIn=[-wvt.Lz 0],zOut=z, ...
        latitude=wvt.latitude,rho0=wvt.rho0,nModes=wvt.Nj, ...
        rotationRate=wvt.rotationRate,g=wvt.g);
    verticalModes.normalization = Normalization.kConstant;
else
    verticalModes = InternalModesWKBSpectral(N2=wvt.N2Function, ...
        zIn=[-wvt.Lz 0],zOut=z,latitude=wvt.latitude,rho0=wvt.rho0, ...
        nModes=wvt.Nj,nEVP=max(256,floor(2.1*wvt.Nz)), ...
        rotationRate=wvt.rotationRate,g=wvt.g);
    verticalModes.normalization = Normalization.geostrophic;
end
verticalModes.upperBoundary = UpperBoundary.rigidLid;
[F0,G0] = verticalModes.modesAtFrequency(0);
modes = struct;
modes.z = z;
[modes.F0,modes.G0,modes.maximumMatchDefect] = matchNativeNormalization( ...
    [ones(numel(z),1) F0(:,1:wvt.Nj-1)], ...
    [zeros(numel(z),1) G0(:,1:wvt.Nj-1)], ...
    wvt.FinvMatrix,wvt.GinvMatrix,sourceRows);
modes.wave = cell(0,3);
kappa = unique(hypot(horizontal.k,horizontal.l),"sorted");
for currentKappa = reshape(kappa,1,[])
    iHorizontal = find(abs(hypot(horizontal.k,horizontal.l)-currentKappa) ...
        <= 100*eps*max(currentKappa,1),1);
    verticalModes.normalization = Normalization.kConstant;
    [Fw,Gw] = verticalModes.modesAtWavenumber(currentKappa);
    [Fw,Gw,matchDefect] = matchNativeNormalization( ...
        [ones(numel(z),1) Fw(:,1:wvt.Nj-1)], ...
        [zeros(numel(z),1) Gw(:,1:wvt.Nj-1)], ...
        wvt.FwInvMatrix(horizontal.kMode(iHorizontal), ...
        horizontal.lMode(iHorizontal)), ...
        wvt.GwInvMatrix(horizontal.kMode(iHorizontal), ...
        horizontal.lMode(iHorizontal)),sourceRows);
    modes.wave(end+1,:) = {currentKappa, ...
        Fw,Gw};
    modes.maximumMatchDefect = max(modes.maximumMatchDefect,matchDefect);
end
end

function [F,G,defect] = matchNativeNormalization(F,G,Fsource,Gsource,sourceRows)
for iMode = 1:size(F,2)
    values = [F(sourceRows,iMode);G(sourceRows,iMode)];
    reference = [Fsource(:,iMode);Gsource(:,iMode)];
    scaling = values\reference;
    F(:,iMode) = scaling*F(:,iMode);
    G(:,iMode) = scaling*G(:,iMode);
end
represented = [F(sourceRows,:);G(sourceRows,:)];
reference = [Fsource;Gsource];
defect = norm(represented-reference,"fro") ...
    /max(norm(reference,"fro"),realmin);
end

function rows = sourceGridRows(targetZ,sourceZ)
[isPresent,rows] = ismember(sourceZ(:),targetZ(:));
if ~all(isPresent)
    error("WVTerrainEnergyGalerkin:ProductionNativeGridMismatch", ...
        "The native evaluation grid does not contain every originating WaveVortexModel depth.")
end
end

function state = nativeState(wvt,modes,production,row,bottom)
layout = production.layout;
iK = layout.horizontalIndex(row);
horizontal = production.horizontalLayout(iK,:);
iNative = horizontal.nativeIndex;
iJ = find(wvt.j == layout.j(row),1);
component = layout.component(row);
kappa = hypot(horizontal.k,horizontal.l);
nZ = numel(modes.z);
state = struct("u",zeros(nZ,1),"v",zeros(nZ,1), ...
    "w",zeros(nZ,1),"eta",zeros(nZ,1));
if component == "Ap" || component == "Am"
    waveKappa = cell2mat(modes.wave(:,1));
    iWave = find(abs(waveKappa-kappa) <= 100*eps*max(kappa,1),1);
    F = modes.wave{iWave,2};
    G = modes.wave{iWave,3};
    if component == "Ap"
        state.u = F(:,iJ)*wvt.UAp(iJ,iNative);
        state.v = F(:,iJ)*wvt.VAp(iJ,iNative);
        state.w = G(:,iJ)*wvt.WAp(iJ,iNative);
        state.eta = G(:,iJ)*wvt.NAp(iJ,iNative);
    else
        state.u = F(:,iJ)*wvt.UAm(iJ,iNative);
        state.v = F(:,iJ)*wvt.VAm(iJ,iNative);
        state.w = G(:,iJ)*wvt.WAm(iJ,iNative);
        state.eta = G(:,iJ)*wvt.NAm(iJ,iNative);
    end
elseif component == "A0"
    state.u = modes.F0(:,iJ)*wvt.UA0(iJ,iNative);
    state.v = modes.F0(:,iJ)*wvt.VA0(iJ,iNative);
    state.eta = modes.G0(:,iJ)*wvt.NA0(iJ,iNative);
elseif kappa > 0
    bottomKappa = cell2mat(bottom(:,1));
    iBottom = find(abs(bottomKappa-kappa) <= 100*eps*max(kappa,1),1);
    profile = bottom{iBottom,2};
    state.u = -1i*horizontal.l*profile.psi;
    state.v = 1i*horizontal.k*profile.psi;
    state.eta = profile.eta;
else
    state.eta = -modes.z/wvt.Lz;
end
end

function [values,matchDefect] = bottomProfiles(problem,horizontal,z)
kappa = unique(hypot(horizontal.k,horizontal.l),"sorted");
kappa = kappa(kappa > 0);
values = cell(numel(kappa),2);
if isempty(kappa)
    matchDefect = 0;
    return
end
wvt = problem.originatingTransform;
N2Bottom = wvt.N2Function(-wvt.Lz);
modeProblem = IMSurfaceGeostrophicModes.atWavenumber( ...
    N2=wvt.N2Function,zDomain=[-wvt.Lz 0],f0=wvt.f,g=wvt.g, ...
    k=kappa,g0=Inf,gd=N2Bottom,surfaceAnomaly="noFreeSurface");
basis = IMSolverSpectral(nEVP=max(128,4*wvt.Nz)). ...
    solveSurfaceGeostrophicModes(modeProblem);
[zDescending,order] = sort(z,"descend");
inverseOrder = zeros(size(order));
inverseOrder(order) = 1:numel(order);
psi = basis.F(zDescending);
eta = (wvt.f/wvt.g)*basis.G(zDescending);
psi = psi(inverseOrder,:);
eta = eta(inverseOrder,:);
endpoint = (wvt.f/wvt.g)*basis.G([-wvt.Lz;0]);
normalization = endpoint(1,:);
psi = psi./normalization;
eta = eta./normalization;
for iK = 1:numel(kappa)
    values{iK,1} = kappa(iK);
    values{iK,2} = struct("psi",psi(:,iK),"eta",eta(:,iK));
end
sourceRows = sourceGridRows(z,wvt.z);
fullHorizontal = problem.horizontalLayout;
matchDefect = 0;
for iK = 1:numel(kappa)
    iHorizontal = find(abs(hypot(fullHorizontal.k,fullHorizontal.l)-kappa(iK)) ...
        <= 100*eps*max(kappa(iK),1),1);
    reference = problem.bottomInversionProfiles{iHorizontal};
    represented = [values{iK,2}.psi(sourceRows); ...
        values{iK,2}.eta(sourceRows)];
    expected = [reference.psi;reference.eta];
    matchDefect = max(matchDefect,norm(represented-expected) ...
        /max(norm(expected),realmin));
end
end

function result = contractAtResolution(~,production,c,flat,native)
nProduction = height(production.layout);
X = zeros(size(c.N,2),nProduction);
embedding = zeros(nProduction,1);
bottomValue = zeros(nProduction,1);
surfaceValue = zeros(nProduction,1);
expectedBottom = double(ismember(production.layout.family, ...
    ["bottom-inversion","mean-bottom-compatible"]));
for row = 1:nProduction
    iProductionK = production.layout.horizontalIndex(row);
    kMode = production.horizontalLayout.kMode(iProductionK);
    lMode = production.horizontalLayout.lMode(iProductionK);
    iContextK = find(c.horizontalLayout.kMode == kMode ...
        & c.horizontalLayout.lMode == lMode,1);
    state = struct("u",native.u(:,row),"v",native.v(:,row), ...
        "w",native.w(:,row),"eta",native.eta(:,row));
    [X(:,row),embedding(row),bottomValue(row),surfaceValue(row)] = ...
        embedState(c,iContextK,state,expectedBottom(row));
end

H = X'*flat.E*X;
energy = real(diag(H));
if any(~isfinite(energy)) || any(energy <= 0)
    error("WVTerrainEnergyGalerkin:NonpositiveProductionEnergy", ...
        "Every production coordinate must have positive finite flat physical energy.")
end
energyScale = 1./sqrt(energy);
Xnormalized = X.*energyScale.';
Hnormalized = Xnormalized'*flat.E*Xnormalized;
leftInverseNormalized = Hnormalized\(Xnormalized'*flat.E);
leftInverse = diag(energyScale)*leftInverseNormalized;
roundTrip = leftInverseNormalized*Xnormalized;
QX = flat.QProjected*Xnormalized;
zeroAPV = ismember(production.layout.family, ...
    ["internal-wave","inertial","bottom-inversion"]);
zeroAPVDefect = vecnorm(QX(:,zeroAPV),2,1).' ...
    /max(norm(flat.QProjected/chol(flat.E),2),realmin);
conjugacy = conjugacyDefect(c,X,production.conjugateIndex);
referenceGram = Hnormalized;
result = struct;
result.productionToPrimitive = X;
result.primitiveToProduction = leftInverse;
result.energyNormalizedProductionToPrimitive = Xnormalized;
result.primitiveToEnergyNormalizedProduction = leftInverseNormalized;
result.physicalEnergyGram = H;
result.normalizedPhysicalEnergyGram = Hnormalized;
result.familyProjectors = familyProjectors(production.layout,X,flat.E);
result.nativeColumnEnergy = energy;
result.maximumPhysicalEmbeddingDefect = max(embedding);
result.physicalEmbeddingDefect = embedding;
result.roundTripDefect = norm(roundTrip-eye(nProduction),"fro")/sqrt(nProduction);
result.energyNormalizationDefect = max(abs(diag(Hnormalized)-1));
result.scaledEnergyReciprocalConditionNumber = rcond(Hnormalized);
result.maximumZeroAPVDefect = max(zeroAPVDefect,[],"all");
result.zeroAPVDefect = zeroAPVDefect;
result.maximumEndpointDefect = max([abs(bottomValue-expectedBottom);abs(surfaceValue)]);
result.bottomValue = bottomValue;
result.surfaceValue = surfaceValue;
result.conjugacyDefect = conjugacy;
result.waveVortexMatchDefect = native.waveVortexMatchDefect;
result.bottomInversionMatchDefect = native.bottomInversionMatchDefect;
result.maximumNativeMatchDefect = max( ...
    native.waveVortexMatchDefect,native.bottomInversionMatchDefect);
result.referenceGram = referenceGram;
result.numberOfProductionCoordinates = nProduction;
result.numberOfPrimitiveCoordinates = size(c.N,2);
result.numberOfDiscardedPrimitiveCoordinates = size(c.N,2)-rank(X);
end

function [coordinate,defect,bottomValue,surfaceValue] = embedState(c,iHorizontal,state,expectedBottom)
weight = c.quadrature.weight;
uCoefficient = weightedFit(c.spaces.F,state.u,weight);
vCoefficient = weightedFit(c.spaces.F,state.v,weight);
wCoefficient = weightedFit(c.spaces.G,state.w,weight);
chi = c.spaces.H(:,end);
etaInterior = weightedFit(c.spaces.H(:,1:end-1), ...
    state.eta-expectedBottom*chi,weight);
etaCoefficient = [etaInterior;expectedBottom];
rawBlock = [uCoefficient;vCoefficient;wCoefficient;etaCoefficient];
rawRows = (iHorizontal-1)*c.nXBlock+(1:c.nXBlock);
admissible = c.layout.admissibleRanges{iHorizontal};
localBasis = c.N(rawRows,admissible);
localCoordinate = localBasis\rawBlock;
represented = localBasis*localCoordinate;
coordinate = zeros(size(c.N,2),1);
coordinate(admissible) = localCoordinate;
defect = physicalEmbeddingDefect(c,state,represented);
etaRows = 2*c.nF+c.nG+(1:c.nH);
bottomValue = c.spaces.Hendpoint(1,:)*represented(etaRows);
surfaceValue = c.spaces.Hendpoint(2,:)*represented(etaRows);
end

function defect = physicalEmbeddingDefect(c,state,represented)
u = c.spaces.F*represented(1:c.nF);
v = c.spaces.F*represented(c.nF+(1:c.nF));
w = c.spaces.G*represented(2*c.nF+(1:c.nG));
eta = c.spaces.H*represented(2*c.nF+c.nG+(1:c.nH));
weight = c.quadrature.weight;
N2 = c.wvt.N2Function(c.quadrature.xi);
errorSquared = sum(weight.*(abs(u-state.u).^2+abs(v-state.v).^2 ...
    +abs(w-state.w).^2+N2.*abs(eta-state.eta).^2));
scaleSquared = sum(weight.*(abs(state.u).^2+abs(state.v).^2 ...
    +abs(state.w).^2+N2.*abs(state.eta).^2));
defect = sqrt(max(real(errorSquared),0)/max(real(scaleSquared),realmin));
end

function coefficient = weightedFit(reconstruction,values,weight)
coefficient = (reconstruction'*(weight.*reconstruction)) ...
    \(reconstruction'*(weight.*values));
end

function defect = conjugacyDefect(c,X,productionConjugate)
nPrimitive = size(X,1);
nProduction = size(X,2);
primitiveConjugate = sparse((1:nPrimitive)',c.coordinateConjugateIndex,1,nPrimitive,nPrimitive);
productionMap = sparse((1:nProduction)',productionConjugate,1,nProduction,nProduction);
defect = norm(primitiveConjugate*conj(X)-X*productionMap,"fro") ...
    /max(norm(X,"fro"),realmin);
end

function projectors = familyProjectors(layout,X,H)
families = unique(layout.family,"stable");
projectors = struct;
for family = reshape(families,1,[])
    columns = X(:,layout.family == family);
    gram = columns'*H*columns;
    field = matlab.lang.makeValidName(family);
    projectors.(field) = columns*(gram\(columns'*H));
end
end

function value = maximumReferenceGramDifference(results,reference)
value = 0;
for iResult = 1:numel(results)
    value = max(value,norm(results(iResult).referenceGram-reference.referenceGram,"fro") ...
        /max(norm(reference.referenceGram,"fro"),realmin));
end
end

function layouts = guardLayouts(problem,trustedBounds,supportBounds)
layouts = cell(size(supportBounds,1),1);
for iSupport = 1:size(supportBounds,1)
    support = retainedLayout(problem,supportBounds(iSupport,:));
    isTrusted = abs(support.kMode) <= trustedBounds(1) ...
        & abs(support.lMode) <= trustedBounds(2);
    layouts{iSupport} = support(~isTrusted,:);
end
end

function diagnosis = blockerDiagnosis(primary,verticalConverges,supportDefect,paddingDefect,tolerance)
if primary.roundTripDefect > tolerance.roundTrip || primary.energyNormalizationDefect > tolerance.energy
    diagnosis = "The declared production coordinates do not form a stable positive-energy coordinate system in the primitive oracle.";
elseif primary.maximumEndpointDefect > tolerance.endpoint
    diagnosis = "At least one production family fails its required flat displacement endpoint.";
elseif primary.maximumZeroAPVDefect > tolerance.zeroAPV
    diagnosis = "A wave, inertial, or explicit bottom coordinate fails the flat zero-APV classification.";
elseif primary.conjugacyDefect > tolerance.conjugacy
    diagnosis = "The native production columns do not commute with the declared Fourier and branch conjugacy map.";
elseif primary.maximumNativeMatchDefect > tolerance.nativeMatch
    diagnosis = "The continuous production profiles do not reproduce the originating WaveVortexModel modes or stored bottom inversion.";
elseif ~verticalConverges
    diagnosis = "The native physical columns do not converge under independent primitive vertical refinement.";
elseif supportDefect > tolerance.support || paddingDefect > tolerance.padding
    diagnosis = "The trusted production contract depends materially on guard support or padding.";
else
    diagnosis = "The production physical-state contract fails an unclassified acceptance condition.";
end
end

function layout = retainedLayout(problem,bounds)
source = problem.horizontalLayout;
mask = abs(source.kMode) <= bounds(1) ...
    & abs(source.lMode) <= bounds(2);
layout = source(mask,:);
expected = (2*bounds(1)+1)*(2*bounds(2)+1);
if height(layout) ~= expected
    error("WVTerrainEnergyGalerkin:UnavailableProductionSupport", ...
        "The originating transform does not retain the complete signed rectangle [%d %d].", ...
        bounds(1),bounds(2))
end
end
