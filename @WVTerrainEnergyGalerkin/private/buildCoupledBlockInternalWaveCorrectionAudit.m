function audit = buildCoupledBlockInternalWaveCorrectionAudit(problem,trustedBounds,supportBounds,stationaryDegree,degrees,paddingFactors,terrainScales,tangentStep,maximumOuterIterations,maximumKrylovRestarts,shouldCompareSlopeSeed,quadratureOrder,productionContract)
% Build the Milestone-10.2.1 coupled-block correction oracle.

comparisonDegree = degrees(end);
[baseline,setup] = buildResidualEnrichedTerrainModeAudit(problem, ...
    trustedBounds,supportBounds,stationaryDegree,degrees, ...
    comparisonDegree,paddingFactors,terrainScales,tangentStep, ...
    maximumOuterIterations,quadratureOrder,productionContract);

nDegree = numel(degrees);
nPadding = numel(paddingFactors);
exactDetails = cell(nDegree,nPadding);
for iDegree = 1:nDegree
    for iPadding = 1:nPadding
        exactDetails{iDegree,iPadding} = exactCoupledCase( ...
            setup.details{iDegree,iPadding},baseline.details{iDegree,iPadding}, ...
            maximumOuterIterations);
    end
end
primary = exactDetails{end,1};

paddingDefect = maximumPrincipalSine( ...
    exactDetails{end,1}.internalModeBasis, ...
    exactDetails{end,2}.internalModeBasis,primary.energyMatrix);
nestedDefect = zeros(nDegree-1,1);
for iDegree = 1:nDegree-1
    transfer = baseline.seed.reference.transfers{iDegree}.matrix;
    nestedDefect(iDegree) = maximumPrincipalSine( ...
        transfer*exactDetails{iDegree,1}.internalModeBasis, ...
        exactDetails{iDegree+1,1}.internalModeBasis, ...
        exactDetails{iDegree+1,1}.energyMatrix);
end

tolerance = struct("structure",1e-12,"stationary",1e-10, ...
    "residual",1e-10,"projector",1e-9,"apv",1e-8, ...
    "bottom",1e-10,"strong",1e-5,"padding",1e-8, ...
    "nested",1e-8,"conjugacy",1e-10,"correction",1e-10);
allCases = [exactDetails{:}];
isolationPasses = ~any([allCases.isolationFailure]) ...
    && primary.internalProjectorDefect <= tolerance.projector ...
    && paddingDefect <= tolerance.padding ...
    && max(nestedDefect,[],"all") <= tolerance.nested;
structurePasses = max([allCases.energyHermitianDefect]) <= tolerance.structure ...
    && max([allCases.exchangeSkewHermitianDefect]) <= tolerance.structure;
stationaryPasses = max([allCases.stationaryOrthogonalityDefect]) <= tolerance.stationary;
correctionPasses = max([allCases.maximumCorrectionResidual]) <= tolerance.correction ...
    && max([allCases.maximumCorrectionOrthogonalityDefect]) <= tolerance.correction;
physicalPasses = primary.maximumInternalResidual <= tolerance.residual ...
    && primary.maximumInternalAPVDefect <= tolerance.apv ...
    && primary.maximumInternalBottomDefect <= tolerance.bottom ...
    && primary.maximumInternalStrongResidual <= tolerance.strong ...
    && primary.conjugacyDefect <= tolerance.conjugacy;
formulationPasses = structurePasses && stationaryPasses ...
    && correctionPasses && physicalPasses;
exactPasses = isolationPasses && formulationPasses;

if ~isolationPasses
    classification = "coupled-block-isolation-blocker";
    diagnosis = exactIsolationDiagnosis(primary,paddingDefect,nestedDefect,tolerance);
elseif ~formulationPasses
    classification = "coupled-block-formulation-blocker";
    diagnosis = exactBlockerDiagnosis(primary,paddingDefect,nestedDefect,tolerance);
else
    classification = "exact-coupled-pass";
    diagnosis = "The exact complementary Sylvester oracle passes; the iterative coupled realization is required before a final Milestone-10.2.1 classification.";
end

audit = struct;
audit.scope = "milestone-10.2.1-coupled-block-internal-waves";
audit.status = classification;
audit.classification = classification;
audit.isCompatible = false;
audit.diagnosis = diagnosis;
audit.baseline = baseline;
audit.exact = struct("details",{exactDetails},"primary",primary, ...
    "passes",exactPasses,"paddingProjectorDefect",paddingDefect, ...
    "nestedProjectorDefect",nestedDefect, ...
    "formulationPasses",formulationPasses, ...
    "isolationPasses",isolationPasses);
audit.iterative = struct("wasAttempted",false,"passes",false);
audit.slopeControl = struct("wasAttempted",false, ...
    "wasRequested",shouldCompareSlopeSeed);
audit.requiredTolerance = tolerance;
audit.economy = struct( ...
    "NJoint",primary.dimension.maximumJointTrial, ...
    "NIndependent",baseline.economy.independentTrialDimension, ...
    "NAmbient",primary.dimension.ambient);
audit.maximumOuterIterations = maximumOuterIterations;
audit.maximumKrylovRestarts = maximumKrylovRestarts;
audit.productionContract = productionContract;
audit.trustedModeBounds = trustedBounds;
audit.supportModeBounds = supportBounds;
audit.stationaryPolynomialDegree = stationaryDegree;
audit.primitivePolynomialDegrees = degrees;
audit.paddingFactors = paddingFactors;
audit.terrainScales = terrainScales;
audit.tangentStep = tangentStep;
audit.quadratureOrder = quadratureOrder;
if exactPasses
    audit.nextScope = "milestone-10.2.1-iterative-coupled-realization";
else
    audit.nextScope = "milestone-10.2.1-stop-for-analysis";
end
end

function result = exactCoupledCase(setup,baseline,maximumIterations)
c = setup.context;
primitive = setup.primitive;
direction = primitive.evaluatedDirections(end);
H = direction.E;
G = baseline.stationaryBasis;
target = baseline.denseTargetBasis;
targetFrequency = baseline.denseTargetFrequency;
targetDimension = setup.productionTargets.numberOfCoordinates;
X = internalColumns(setup.internal,1);
X = projectOut(X,G,H);
X = energyOrthonormalize(X,H);
previous = X;
history = repmat(emptyExactIteration,maximumIterations+1,1);
maximumCorrectionResidual = 0;
maximumCorrectionOrthogonality = 0;
isolationFailure = false;
for iIteration = 0:maximumIterations
    ritz = trackedRitz(direction,X,previous,targetDimension);
    history(iIteration+1) = ritz;
    if ritz.maximumResidual <= 1e-10 || iIteration == maximumIterations
        isolationFailure = ritz.isolationFailure;
        break
    end
    correction = exactComplementaryCorrection(direction,G,ritz.basis, ...
        ritz.theta,ritz.residual);
    maximumCorrectionResidual = max(maximumCorrectionResidual, ...
        correction.residualDefect);
    maximumCorrectionOrthogonality = max( ...
        maximumCorrectionOrthogonality,correction.orthogonalityDefect);
    X = energyOrthonormalize(ritz.basis+correction.values,H);
end
history = history(1:iIteration+1);
final = history(end);
internal = final.basis;
projectorDefect = maximumPrincipalSine(internal,target,H);
internalFrequency = sort(final.frequency);
targetFrequency = sort(targetFrequency);
if numel(internalFrequency) == numel(targetFrequency)
    frequencyDefect = norm(internalFrequency-targetFrequency) ...
        /max(norm(targetFrequency),realmin);
else
    frequencyDefect = Inf;
end
apv = trustedAPVDefect(c,direction,internal);
bottom = modalBottomDefect(direction,internal,final.frequency);
strong = strongModeDiagnostics(c,direction,internal,final.frequency);
stationaryOrthogonality = norm(G'*H*internal,"fro") ...
    /max(norm(G,"fro")*norm(H*internal,"fro"),realmin);
projector = energyProjector(internal,H);
conjugacy = sparse((1:size(H,1))',c.coordinateConjugateIndex, ...
    1,size(H,1),size(H,1));
conjugacyResidual = projector*conjugacy-conjugacy*conj(projector);

result = struct;
result.internalModeBasis = internal;
result.internalFrequency = final.frequency;
result.internalProjectorDefect = projectorDefect;
result.internalFrequencyDefect = frequencyDefect;
result.maximumInternalResidual = final.maximumResidual;
result.maximumInternalAPVDefect = apv;
result.maximumInternalBottomDefect = bottom;
result.maximumInternalStrongResidual = strong.maximumResidual;
result.stationaryOrthogonalityDefect = stationaryOrthogonality;
result.conjugacyDefect = norm(conjugacyResidual,"fro") ...
    /max(2*norm(projector,"fro"),realmin);
result.energyHermitianDefect = norm(H-H',"fro")/max(norm(H,"fro"),realmin);
result.exchangeSkewHermitianDefect = norm(direction.J+direction.J',"fro") ...
    /max(norm(direction.J,"fro"),realmin);
result.maximumCorrectionResidual = maximumCorrectionResidual;
result.maximumCorrectionOrthogonalityDefect = maximumCorrectionOrthogonality;
result.isolationFailure = isolationFailure;
result.numberOfIterations = numel(history)-1;
result.history = history;
result.energyMatrix = H;
result.denseTargetBasis = target;
result.denseTargetFrequency = targetFrequency;
result.dimension = struct("ambient",size(H,1),"stationary",size(G,2), ...
    "validatedInternal",size(internal,2), ...
    "maximumJointTrial",max([history.trialDimension]), ...
    "compressionFactor",size(H,1) ...
    /max(size(G,2)+max([history.trialDimension]),1));
end

function iteration = trackedRitz(direction,X,previous,targetDimension)
H = direction.E;
A = 1i*direction.J;
Hr = X'*H*X;
Ar = X'*A*X;
[R,flag] = chol((Hr+Hr')/2);
if flag ~= 0
    error("WVTerrainEnergyGalerkin:CoupledBlockReducedEnergyNotPositive", ...
        "The coupled reduced energy matrix must remain positive definite.")
end
K = R'\(Ar/R);
[U,T] = schur(K,"complex");
frequency = real(diag(T));
[frequency,order] = sort(frequency);
vectors = X*(R\U(:,order));
vectors = energyNormalizeColumns(vectors,H);
residual = A*vectors-(H*vectors).*frequency.';
energyFactor = chol((H+H')/2);
energyResidual = energyFactor'\residual;
absoluteUncertainty = vecnorm(energyResidual,2,1).' ...
    ./max(vecnorm(energyFactor*vectors,2,1).',realmin);
absoluteUncertainty = max(absoluteUncertainty, ...
    10*eps(max(norm(energyFactor'\(A/energyFactor),2),realmin)));
projector = energyProjector(previous,H);
participation = real(sum(conj(vectors).*(H*projector*vectors),1)).';
[~,selected] = maxk(participation,min(targetDimension,numel(participation)));
selected = expandInseparable(sort(selected),frequency,absoluteUncertainty);
isolationFailure = numel(selected) ~= targetDimension;
basis = vectors(:,selected);
theta = basis'*A*basis;
selectedResidual = A*basis-H*basis*theta;
energySelectedResidual = energyFactor'\selectedResidual;
operatorNorm = norm(energyFactor'\(A/energyFactor),2);
relativeResidual = vecnorm(energySelectedResidual,2,1).' ...
    ./max(operatorNorm+vecnorm(theta,2,1).',realmin);
iteration = emptyExactIteration;
iteration.frequency = real(diag(theta));
iteration.basis = basis;
iteration.theta = theta;
iteration.residual = selectedResidual;
iteration.maximumResidual = max(relativeResidual,[],"all");
iteration.absoluteFrequencyUncertainty = absoluteUncertainty(selected);
iteration.isolationFailure = isolationFailure;
iteration.trialDimension = size(X,2);
end

function correction = exactComplementaryCorrection(direction,G,Y,theta,residual)
H = direction.E;
A = 1i*direction.J;
energyFactor = chol((H+H')/2);
U = energyOrthonormalize([G Y],H);
Ubar = energyFactor*U;
[complete,~] = qr(Ubar);
Q = complete(:,size(Ubar,2)+1:end);
Abar = energyFactor'\(A/energyFactor);
residualBar = energyFactor'\residual;
reducedA = Q'*Abar*Q;
right = -Q'*residualBar;
warningState = warning("off","MATLAB:sylvester:EquationNotSolved");
cleanup = onCleanup(@()warning(warningState));
Z = sylvester(reducedA,-theta,right);
clear cleanup
deltaBar = Q*Z;
values = energyFactor\deltaBar;
lambda = -Ubar'*(residualBar+Abar*deltaBar-deltaBar*theta);
equationResidual = Abar*deltaBar-deltaBar*theta+Ubar*lambda+residualBar;
scale = norm(Abar*deltaBar,"fro")+norm(deltaBar*theta,"fro") ...
    +norm(Ubar*lambda,"fro")+norm(residualBar,"fro");
correction = struct;
correction.values = values;
correction.lambda = lambda;
correction.residualDefect = norm(equationResidual,"fro")/max(scale,realmin);
correction.orthogonalityDefect = norm(Ubar'*deltaBar,"fro") ...
    /max(norm(deltaBar,"fro"),realmin);
end

function columns = internalColumns(internal,scale)
parts = cell(numel(internal),1);
for iBlock = 1:numel(internal)
    parts{iBlock} = internal(iBlock).baseBasis ...
        +scale*internal(iBlock).correction;
end
columns = horzcat(parts{:});
end

function indices = expandInseparable(indices,frequency,uncertainty)
indices = indices(:);
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

function X = projectOut(X,basis,H)
if ~isempty(basis) && ~isempty(X)
    X = X-basis*((basis'*H*basis)\(basis'*H*X));
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
firstResidual = firstQ-secondQ*(secondQ'*firstQ);
secondResidual = secondQ-firstQ*(firstQ'*secondQ);
defect = max(norm(firstResidual,2),norm(secondResidual,2));
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
scale = abs(frequency)*norm(direction.B/energyFactor,2) ...
    +norm(direction.R/energyFactor,2);
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
pointwise = pointwiseStrongResiduals( ...
    c,direction,rawState,pressure,frequency);
maximum = max([descriptorDefect primitiveDefect gaugeDefect ...
    pointwise.maximumByMode],[],2);
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
maximumByMode = max( ...
    [uDefect vDefect wDefect etaDefect continuityDefect],[],2);
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

function iteration = emptyExactIteration
iteration = struct("frequency",zeros(0,1),"basis",zeros(0), ...
    "theta",zeros(0),"residual",zeros(0),"maximumResidual",Inf, ...
    "absoluteFrequencyUncertainty",zeros(0,1), ...
    "isolationFailure",false,"trialDimension",0);
end

function diagnosis = exactBlockerDiagnosis(primary,padding,nested,tolerance)
diagnosis = sprintf(strjoin([ ...
    "The exact complementary oracle failed a physical gate: correction %.3g, " ...
    "orthogonality %.3g, Ritz %.3g, projector %.3g, APV %.3g, bottom %.3g, " ...
    "strong %.3g, padding %.3g, nested %.3g (required Ritz %.3g)."],""), ...
    primary.maximumCorrectionResidual, ...
    primary.maximumCorrectionOrthogonalityDefect, ...
    primary.maximumInternalResidual,primary.internalProjectorDefect, ...
    primary.maximumInternalAPVDefect,primary.maximumInternalBottomDefect, ...
    primary.maximumInternalStrongResidual,padding,max(nested,[],"all"), ...
    tolerance.residual);
end

function diagnosis = exactIsolationDiagnosis(primary,padding,nested,tolerance)
diagnosis = sprintf(strjoin([ ...
    "The exact coupled correction reaches the requested physical subspace at each fixed discretization, but that subspace is not isolated under refinement: " ...
    "projector %.3g, padding %.3g, nested-degree %.3g (required %.3g). " ...
    "The iterative coupled realization is therefore not attempted."],""), ...
    primary.internalProjectorDefect,padding,max(nested,[],"all"), ...
    tolerance.nested);
end
