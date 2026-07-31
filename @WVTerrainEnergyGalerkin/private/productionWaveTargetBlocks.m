function result = productionWaveTargetBlocks(problem,productionContract,c,flat)
% Map the declared native production waves into one primitive context.

layout = productionContract.productionLayout;
isWave = layout.family == "internal-wave";
waveRows = find(isWave);
if isempty(waveRows)
    error("WVTerrainEnergyGalerkin:MissingProductionInternalWaves", ...
        "The production contract must contain at least one internal-wave coordinate.")
end
productionHorizontal = productionContract.productionHorizontalLayout;
z = c.quadrature.xi;
modes = nativeWaveShapes(problem.originatingTransform,productionHorizontal,z);
X = zeros(size(c.N,2),numel(waveRows));
embedding = zeros(numel(waveRows),1);
for iWave = find(layout.isPrimary(waveRows)).'
    row = waveRows(iWave);
    state = nativeWaveState(problem.originatingTransform,modes, ...
        productionHorizontal,layout(row,:));
    [X(:,iWave),embedding(iWave)] = embedWaveState(c,state, ...
        productionHorizontal.kMode(layout.horizontalIndex(row)), ...
        productionHorizontal.lMode(layout.horizontalIndex(row)));
end
primitiveConjugate = sparse((1:size(X,1))',c.coordinateConjugateIndex, ...
    1,size(X,1),size(X,1));
for iWave = find(~layout.isPrimary(waveRows)).'
    partnerRow = productionContract.conjugateCoordinateIndex(waveRows(iWave));
    partner = find(waveRows == partnerRow,1);
    if isempty(partner)
        error("WVTerrainEnergyGalerkin:IncompleteProductionWaveConjugacy", ...
            "Every secondary production wave must have its primary partner in the retained wave set.")
    end
    X(:,iWave) = primitiveConjugate*conj(X(:,partner));
    embedding(iWave) = embedding(partner);
end

wvt = problem.originatingTransform;
nativeFrequency = abs(wvt.Omega(layout.nativeCoefficientIndex(waveRows)));
frequencyScale = max(nativeFrequency,[],"all");
frequencyTolerance = 1e3*eps(max(frequencyScale,realmin));
unused = (1:numel(waveRows)).';
blocks = cell(0,1);
blockIndex = zeros(numel(waveRows),1);
H = flat.E;
while ~isempty(unused)
    seed = unused(1);
    members = unused(abs(nativeFrequency(unused)-nativeFrequency(seed)) ...
        <= frequencyTolerance);
    gram = X(:,members)'*H*X(:,members);
    [R,flag] = chol((gram+gram')/2);
    if flag ~= 0
        error("WVTerrainEnergyGalerkin:ProductionWaveEnergyNotPositive", ...
            "Each mapped production-wave block must have a positive physical-energy Gram matrix.")
    end
    X(:,members) = X(:,members)/R;
    block = struct;
    block.basis = X(:,members);
    block.productionCoordinateIndices = waveRows(members);
    block.productionWaveColumnIndices = members;
    block.nativeFrequency = mean(nativeFrequency(members));
    block.dimension = numel(members);
    blocks{end+1,1} = block; %#ok<AGROW>
    blockIndex(members) = numel(blocks);
    unused = setdiff(unused,members,"stable");
end
roundTrip = (X'*H*X)\(X'*H*X);

assignment = table(waveRows,blockIndex,nativeFrequency, ...
    layout.component(waveRows),layout.j(waveRows), ...
    VariableNames=["productionCoordinateIndex","blockIndex", ...
    "nativeFrequency","component","j"]);
result = struct;
result.basis = X;
result.blocks = blocks;
result.assignment = assignment;
result.productionCoordinateIndices = waveRows;
result.numberOfCoordinates = numel(waveRows);
result.maximumEmbeddingDefect = max(embedding,[],"all");
result.embeddingDefect = embedding;
result.maximumNativeMatchDefect = modes.maximumMatchDefect;
result.energyOrthogonalityDefect = norm(X'*H*X-eye(size(X,2)),"fro") ...
    /sqrt(size(X,2));
result.roundTripDefect = norm(roundTrip-eye(size(roundTrip)),"fro") ...
    /sqrt(size(roundTrip,1));
end

function modes = nativeWaveShapes(wvt,horizontal,z)
zCommon = unique([z(:);wvt.z(:)],"sorted");
[isTarget,targetRows] = ismember(z(:),zCommon);
if ~all(isTarget)
    error("WVTerrainEnergyGalerkin:ProductionWaveQuadratureMismatch", ...
        "The common wave evaluation grid must contain every primitive quadrature point.")
end
sourceRows = sourceGridRows(zCommon,wvt.z);
sampleZ = linspace(-wvt.Lz,0,max(17,wvt.Nz)).';
sampleN2 = wvt.N2Function(sampleZ);
isConstant = max(abs(sampleN2-mean(sampleN2))) ...
    <= 100*eps*max(abs(mean(sampleN2)),1);
if isConstant
    verticalModes = InternalModesConstantStratification( ...
        N0=sqrt(mean(sampleN2)),zIn=[-wvt.Lz 0],zOut=zCommon, ...
        latitude=wvt.latitude,rho0=wvt.rho0,nModes=wvt.Nj, ...
        rotationRate=wvt.rotationRate,g=wvt.g);
else
    verticalModes = InternalModesWKBSpectral(N2=wvt.N2Function, ...
        zIn=[-wvt.Lz 0],zOut=zCommon,latitude=wvt.latitude,rho0=wvt.rho0, ...
        nModes=wvt.Nj,nEVP=max(256,floor(2.1*wvt.Nz)), ...
        rotationRate=wvt.rotationRate,g=wvt.g);
end
verticalModes.upperBoundary = UpperBoundary.rigidLid;
verticalModes.normalization = Normalization.kConstant;
modes = struct;
modes.wave = cell(0,3);
modes.maximumMatchDefect = 0;
kappa = unique(hypot(horizontal.k,horizontal.l),"sorted");
for currentKappa = reshape(kappa,1,[])
    iHorizontal = find(abs(hypot(horizontal.k,horizontal.l)-currentKappa) ...
        <= 100*eps*max(currentKappa,1),1);
    [Fw,Gw] = verticalModes.modesAtWavenumber(currentKappa);
    [Fw,Gw,matchDefect] = matchNativeNormalization( ...
        [ones(numel(zCommon),1) Fw(:,1:wvt.Nj-1)], ...
        [zeros(numel(zCommon),1) Gw(:,1:wvt.Nj-1)], ...
        wvt.FwInvMatrix(horizontal.kMode(iHorizontal), ...
        horizontal.lMode(iHorizontal)), ...
        wvt.GwInvMatrix(horizontal.kMode(iHorizontal), ...
        horizontal.lMode(iHorizontal)),sourceRows);
    modes.wave(end+1,:) = {currentKappa,Fw(targetRows,:),Gw(targetRows,:)};
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
    error("WVTerrainEnergyGalerkin:ProductionWaveNativeGridMismatch", ...
        "The native wave grid does not contain every originating WaveVortexModel depth.")
end
end

function state = nativeWaveState(wvt,modes,horizontal,layout)
iK = layout.horizontalIndex;
iNative = horizontal.nativeIndex(iK);
iJ = find(wvt.j == layout.j,1);
kappa = hypot(horizontal.k(iK),horizontal.l(iK));
waveKappa = cell2mat(modes.wave(:,1));
iWave = find(abs(waveKappa-kappa) <= 100*eps*max(kappa,1),1);
F = modes.wave{iWave,2};
G = modes.wave{iWave,3};
if layout.component == "Ap"
    state = struct("u",F(:,iJ)*wvt.UAp(iJ,iNative), ...
        "v",F(:,iJ)*wvt.VAp(iJ,iNative), ...
        "w",G(:,iJ)*wvt.WAp(iJ,iNative), ...
        "eta",G(:,iJ)*wvt.NAp(iJ,iNative));
else
    state = struct("u",F(:,iJ)*wvt.UAm(iJ,iNative), ...
        "v",F(:,iJ)*wvt.VAm(iJ,iNative), ...
        "w",G(:,iJ)*wvt.WAm(iJ,iNative), ...
        "eta",G(:,iJ)*wvt.NAm(iJ,iNative));
end
end

function [coordinate,defect] = embedWaveState(c,state,kMode,lMode)
iHorizontal = find(c.horizontalLayout.kMode == kMode ...
    & c.horizontalLayout.lMode == lMode,1);
if isempty(iHorizontal)
    error("WVTerrainEnergyGalerkin:ProductionWaveOutsidePrimitiveSupport", ...
        "Every production wave must lie inside the primitive support.")
end
weight = c.quadrature.weight;
uCoefficient = weightedFit(c.spaces.F,state.u,weight);
vCoefficient = weightedFit(c.spaces.F,state.v,weight);
wCoefficient = weightedFit(c.spaces.G,state.w,weight);
etaCoefficient = [weightedFit(c.spaces.H(:,1:end-1),state.eta,weight);0];
rawBlock = [uCoefficient;vCoefficient;wCoefficient;etaCoefficient];
rawRows = (iHorizontal-1)*c.nXBlock+(1:c.nXBlock);
admissible = c.layout.admissibleRanges{iHorizontal};
localBasis = c.N(rawRows,admissible);
localCoordinate = localBasis\rawBlock;
represented = localBasis*localCoordinate;
coordinate = zeros(size(c.N,2),1);
coordinate(admissible) = localCoordinate;
defect = physicalEmbeddingDefect(c,state,represented);
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
