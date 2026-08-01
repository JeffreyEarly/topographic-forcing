function audit = buildWKBPrimitiveForwardStateAudit(problem,trustedBounds, ...
    supportBounds,wkbDegrees,legendreDegrees,paddingFactors,terrainScales, ...
    cacheDirectory)
% Build the Milestone-10.2.6 complete spectral representation oracle.

configuration = struct("trustedBounds",trustedBounds, ...
    "supportBounds",supportBounds,"wkbDegrees",wkbDegrees, ...
    "legendreDegrees",legendreDegrees, ...
    "paddingFactors",paddingFactors,"terrainScales",terrainScales);
cache = emptyCacheDiagnostics(cacheDirectory);
comparisonOrder = max(65,2*max([wkbDegrees;legendreDegrees])+11);
baselineSupport = supportBounds(1,:);
baselinePadding = paddingFactors(1);

[legendre,cache] = representationSweep(problem,"legendre", ...
    legendreDegrees,baselineSupport,baselinePadding,trustedBounds, ...
    terrainScales,comparisonOrder,cacheDirectory,cache,"degree");
[wkb,cache] = representationSweep(problem,"wkb",wkbDegrees, ...
    baselineSupport,baselinePadding,trustedBounds,terrainScales, ...
    comparisonOrder,cacheDirectory,cache,"degree");

tolerance = requiredTolerance;
legendre = convergenceDiagnostics(legendre,legendre(end),tolerance);
wkb = convergenceDiagnostics(wkb,legendre(end),tolerance);

allControlProblems = standardizedControlProblems(problem);
controlProblems = allControlProblems(1:end-1);
twoDimensionalProblem = allControlProblems(end);
[legendreStandardized,cache] = standardizedDegreeSweep( ...
    controlProblems,"legendre",legendreDegrees,baselineSupport, ...
    baselinePadding,trustedBounds,comparisonOrder,cacheDirectory,cache,[]);
[wkbStandardized,cache] = standardizedDegreeSweep(controlProblems,"wkb", ...
    wkbDegrees,baselineSupport,baselinePadding,trustedBounds, ...
    comparisonOrder,cacheDirectory,cache,legendreStandardized(end));
[legendreSelected,controls,cache] = selectedCaseWithControls(problem, ...
    "legendre",legendreDegrees,legendre,legendreStandardized, ...
    supportBounds,paddingFactors,baselineSupport,baselinePadding, ...
    trustedBounds,terrainScales,comparisonOrder,cacheDirectory,cache, ...
    tolerance,twoDimensionalProblem);
[wkbSelected,wkbControls,cache] = selectedCaseWithControls(problem, ...
    "wkb",wkbDegrees,wkb,wkbStandardized,supportBounds,paddingFactors, ...
    baselineSupport,baselinePadding,trustedBounds,terrainScales, ...
    comparisonOrder,cacheDirectory,cache,tolerance,twoDimensionalProblem);
constantNEquivalence = constantNEquivalenceDiagnostics( ...
    legendreStandardized,wkbStandardized);
if wkbSelected > 0 && ~constantNEquivalence.passes
    wkbSelected = 0;
end
legendrePasses = legendreSelected > 0;
wkbPasses = wkbSelected > 0;
standardized = selectedStandardizedControls(controlProblems, ...
    legendreStandardized,legendreSelected,wkbStandardized,wkbSelected);
if wkbPasses && legendrePasses
    coordinateRatio = legendre(legendreSelected).numberOfAdmissibleCoordinates ...
        /wkb(wkbSelected).numberOfAdmissibleCoordinates;
elseif wkbPasses
    coordinateRatio = Inf;
else
    coordinateRatio = 0;
end
matched = matchedAccuracyDiagnostics(legendre,legendreSelected,wkb, ...
    wkbSelected,coordinateRatio);
if wkbPasses && coordinateRatio >= 2
    classification = "wkb-spectral-acceleration";
    selectedCoordinate = "wkb";
    diagnosis = "The complete WKB--Chebyshev state passes with at least a factor-two coordinate reduction.";
elseif wkbPasses
    classification = "wkb-spectral-equivalent";
    if legendrePasses && matched.shouldPreferLegendre
        selectedCoordinate = "legendre";
    else
        selectedCoordinate = "wkb";
    end
    diagnosis = "The complete WKB--Chebyshev state passes, but not with factor-two compression.";
elseif legendrePasses
    classification = "wkb-spectral-fallback";
    selectedCoordinate = "legendre";
    diagnosis = "The WKB coordinate is unqualified or does not pass the physical gates; the complete Legendre state passes.";
else
    classification = "complete-spectral-blocker";
    selectedCoordinate = "";
    diagnosis = "Neither complete spectral representation passes the required physical and convergence gates.";
end

audit = struct;
audit.scope = "milestone-10.2.6-complete-wkb-primitive-oracle";
audit.status = classification;
audit.classification = classification;
audit.isCompatible = classification ~= "complete-spectral-blocker";
audit.diagnosis = diagnosis;
audit.selectedVerticalCoordinate = selectedCoordinate;
audit.selectedWKBIndex = wkbSelected;
audit.selectedLegendreIndex = legendreSelected;
audit.coordinateReductionFactor = coordinateRatio;
audit.matchedAccuracy = matched;
audit.constantNEquivalence = constantNEquivalence;
audit.wkbDegreeConvergence = wkb;
audit.legendreDegreeConvergence = legendre;
audit.wkbControls = wkbControls;
audit.legendreControls = controls;
audit.wkbControlDiagnostics = summarizeControls(wkbControls);
audit.legendreControlDiagnostics = summarizeControls(controls);
audit.standardizedControls = standardized;
audit.twoDimensionalControls = struct( ...
    "legendre",controls.twoDimensional, ...
    "wkb",wkbControls.twoDimensional);
audit.standardizedLegendreConvergence = legendreStandardized;
audit.standardizedWKBConvergence = wkbStandardized;
audit.controlCases = standardizedControlDefinitions(problem);
audit.requiredTolerance = tolerance;
audit.trustedModeBounds = trustedBounds;
audit.supportModeBounds = supportBounds;
audit.verticalDegrees = wkbDegrees;
audit.legendreReferenceDegrees = legendreDegrees;
audit.paddingFactors = paddingFactors;
audit.terrainScales = terrainScales;
audit.comparisonQuadratureOrder = comparisonOrder;
audit.projectorProbeDimension = 12;
audit.projectorProbeDescription = ...
    "energy-weighted mean-horizontal Legendre vertical degrees zero through two for each physical field";
audit.cache = finalizeCacheDiagnostics(cache);
audit.configuration = configuration;
audit.nextScope = "milestone-10.3-only-after-review";
end

function [sweep,cache] = standardizedDegreeSweep(controls,coordinate,degrees, ...
    support,padding,trusted,comparisonOrder,cacheDirectory,cache,reference)
sweep = repmat(emptyStandardizedDegree,numel(degrees),1);
for iDegree = 1:numel(degrees)
    results = repmat(emptyCase,numel(controls),1);
    for iControl = 1:numel(controls)
        control = controls(iControl);
        controlSupport = max(support,control.minimumSupportBounds);
        layout = retainedLayout(control.problem,controlSupport);
        controlTrusted = min(trusted,control.trustedBounds);
        scales = [0;1];
        configuration = struct("name","standardized-"+control.name, ...
            "coordinate",coordinate,"degree",degrees(iDegree), ...
            "support",controlSupport,"padding",padding, ...
            "trusted",controlTrusted, ...
            "scales",scales,"comparisonOrder",comparisonOrder);
        [results(iControl),record] = cachedCase(control.problem,cacheDirectory, ...
            configuration,@()representationCase(control.problem,layout, ...
            controlTrusted,coordinate,degrees(iDegree),padding,scales, ...
            comparisonOrder));
        cache = appendCacheRecord(cache,record);
    end
    sweep(iDegree) = struct("degree",degrees(iDegree), ...
        "controls",results,"passes",false);
end
for iControl = 1:numel(controls)
    if isempty(reference)
        referenceCase = sweep(end).controls(iControl);
    else
        referenceCase = reference.controls(iControl);
    end
    values = repmat(emptyCase,numel(degrees),1);
    for iDegree = 1:numel(degrees)
        values(iDegree) = sweep(iDegree).controls(iControl);
    end
    values = convergenceDiagnostics(values,referenceCase,requiredTolerance);
    for iDegree = 1:numel(degrees)
        sweep(iDegree).controls(iControl) = values(iDegree);
    end
end
for iDegree = 1:numel(degrees)
    sweep(iDegree).passes = all([sweep(iDegree).controls.passes]);
end
end

function [index,controls,cache] = selectedCaseWithControls(problem, ...
    coordinate,degrees,primary,standardized,supportBounds,paddingFactors, ...
    baselineSupport,baselinePadding,trusted,scales,comparisonOrder, ...
    cacheDirectory,cache,tolerance,twoDimensionalProblem)
index = 0;
controls = emptyControls;
for iCase = 1:numel(degrees)
    if ~verticalCandidatePasses(primary(iCase),tolerance) ...
            || ~standardized(iCase).passes
        continue
    end
    [candidate.support,cache] = controlSweep(problem,coordinate, ...
        degrees(iCase),supportBounds,baselinePadding,trusted,scales, ...
        comparisonOrder,cacheDirectory,cache,"support");
    [candidate.padding,cache] = controlSweep(problem,coordinate, ...
        degrees(1),paddingFactors,supportBounds(end,:),trusted,scales, ...
        comparisonOrder,cacheDirectory,cache,"padding");
    [candidate.twoDimensional,record] = standardizedCase( ...
        twoDimensionalProblem,coordinate,min(degrees(1),8),baselineSupport, ...
        baselinePadding,trusted,comparisonOrder,cacheDirectory);
    cache = appendCacheRecord(cache,record);
    controls = candidate;
    if controlsPass(candidate,tolerance)
        index = iCase;
        return
    end
end
end

function [result,record] = standardizedCase(control,coordinate,degree, ...
    support,padding,trusted,comparisonOrder,cacheDirectory)
controlSupport = max(support,control.minimumSupportBounds);
layout = retainedLayout(control.problem,controlSupport);
controlTrusted = min(trusted,control.trustedBounds);
scales = [0;1];
configuration = struct("name","standardized-"+control.name, ...
    "coordinate",coordinate,"degree",degree,"support",controlSupport, ...
    "padding",padding,"trusted",controlTrusted,"scales",scales, ...
    "comparisonOrder",comparisonOrder);
[result,record] = cachedCase(control.problem,cacheDirectory,configuration, ...
    @()representationCase(control.problem,layout,controlTrusted,coordinate, ...
    degree,padding,scales,comparisonOrder));
if result.isQualified
    result = convergenceDiagnostics(result,result,requiredTolerance);
end
end

function results = selectedStandardizedControls(controls,legendreSweep, ...
    legendreIndex,wkbSweep,wkbIndex)
results = repmat(emptyStandardizedControl,numel(controls),1);
for iControl = 1:numel(controls)
    results(iControl).name = controls(iControl).name;
    if legendreIndex > 0
        results(iControl).legendre = ...
            legendreSweep(legendreIndex).controls(iControl);
        results(iControl).legendrePasses = ...
            results(iControl).legendre.passes;
    end
    if wkbIndex > 0
        results(iControl).wkb = wkbSweep(wkbIndex).controls(iControl);
        results(iControl).wkbPasses = results(iControl).wkb.passes;
    end
end
end

function [results,cache] = representationSweep(problem,coordinate,degrees, ...
    support,padding,trusted,scales,comparisonOrder,cacheDirectory,cache,name)
layout = retainedLayout(problem,support);
results = repmat(emptyCase,numel(degrees),1);
for iDegree = 1:numel(degrees)
    caseConfiguration = struct("name",name,"coordinate",coordinate, ...
        "degree",degrees(iDegree),"support",support,"padding",padding, ...
        "trusted",trusted,"scales",scales, ...
        "comparisonOrder",comparisonOrder);
    [results(iDegree),record] = cachedCase(problem,cacheDirectory, ...
        caseConfiguration,@()representationCase(problem,layout,trusted, ...
        coordinate,degrees(iDegree),padding,scales,comparisonOrder));
    cache = appendCacheRecord(cache,record);
end
end

function [results,cache] = controlSweep(problem,coordinate,degree,control, ...
    fixed,trusted,scales,comparisonOrder,cacheDirectory,cache,name)
if name == "support"
    count = size(control,1);
else
    count = numel(control);
end
results = repmat(emptyCase,count,1);
for iControl = 1:count
    if name == "support"
        support = control(iControl,:);
        padding = fixed;
    else
        support = fixed;
        padding = control(iControl);
    end
    layout = retainedLayout(problem,support);
    caseConfiguration = struct("name",name,"coordinate",coordinate, ...
        "degree",degree,"support",support,"padding",padding, ...
        "trusted",trusted,"scales",scales, ...
        "comparisonOrder",comparisonOrder);
    [results(iControl),record] = cachedCase(problem,cacheDirectory, ...
        caseConfiguration,@()representationCase(problem,layout,trusted, ...
        coordinate,degree,padding,scales,comparisonOrder));
    cache = appendCacheRecord(cache,record);
end
if ~isempty(results)
    results = convergenceDiagnostics(results,results(end),requiredTolerance);
end
end

function result = representationCase(problem,layout,trusted,coordinate, ...
    degree,padding,scales,comparisonOrder)
result = emptyCase;
result.verticalCoordinate = coordinate;
result.degree = degree;
result.supportModeBounds = [max(abs(layout.kMode)) ...
    max(abs(layout.lMode))];
result.paddingFactor = padding;
result.terrainScales = scales;
quadratureOrder = max(2*degree+7,18);
result.quadratureOrder = quadratureOrder;
timer = tic;
try
    [primitive,c] = buildGlobalSmallTerrainPrimitiveAudit(problem,degree, ...
        quadratureOrder,1e-3,horizontalLayout=layout, ...
        paddingFactor=padding,trustedModeBounds=trusted, ...
        rejectTerrainNyquist=true,evaluationScales=scales, ...
        shouldAuditTangent=false,verticalCoordinate=coordinate);
catch exception
    if coordinate == "wkb" && ismember(string(exception.identifier), ...
            ["WVTerrainEnergyGalerkin:WKBPrimitiveStratificationNotPositive", ...
            "WVTerrainEnergyGalerkin:InvalidWKBPrimitiveIntegral", ...
            "WVTerrainEnergyGalerkin:WKBPrimitiveCoordinateFailure"])
        result.qualificationReason = string(exception.message);
        result.buildSeconds = toc(timer);
        return
    end
    rethrow(exception)
end
result.isQualified = true;
result.verticalDiagnostics = c.verticalDiagnostics;
result.numberOfRawCoordinates = c.nX;
result.numberOfAdmissibleCoordinates = size(c.N,2);
result.numberOfPressureCoordinates = c.nPressure;
result.continuityDefect = c.continuityDiagnostics.nullspaceDefect;
result.roundTripDefect = coefficientRoundTrip(c);
result.quadratureAdjointDefect = quadratureAdjointDefect(c);
result.sequence = repmat(emptySequence,numel(scales),1);
result.physicalBases = cell(numel(scales),1);
comparison = trustedPhysicalContext(problem,c,coordinate,degree, ...
    comparisonOrder);
for iScale = 1:numel(scales)
    direction = primitive.evaluatedDirections(iScale);
    scalarDegree = min(4,degree-1);
    sequence = auditPrimitiveSequenceDirection(c,direction,trusted, ...
        scalarDegree,scales(iScale));
    result.sequence(iScale) = sequenceSummary(sequence);
    [basis,energyDefect] = trustedPhysicalBasis(c,direction, ...
        scales(iScale),comparison);
    result.physicalBases{iScale} = basis;
    result.sequence(iScale).comparisonEnergyDefect = energyDefect;
end
result.conjugacyDefect = max([result.sequence.conjugacyDefect]);
result.maximumStructuralDefect = max([[result.sequence.energyHermitianDefect] ...
    [result.sequence.exchangeSkewHermitianDefect]]);
result.maximumStateRepresentationDefect = ...
    max([result.sequence.stateRepresentationDefect]);
result.maximumEnergyDefect = max([result.sequence.energyDefect]);
result.maximumAPVDefect = max([result.sequence.strongAPVDefect]);
result.maximumBottomDefect = max([result.sequence.bottomDefect]);
result.maximumGreenDefect = max([result.sequence.greenIdentityDefect]);
result.maximumStationaryDefect = max([result.sequence.stationaryRowDefect]);
result.maximumProjectionDefect = max([[result.sequence.projectionAdjointDefect] ...
    [result.sequence.exactConvolutionDefect]]);
result.energyConditionNumber = 1/min([result.sequence.energyRcond]);
result.buildSeconds = toc(timer);
end

function summary = sequenceSummary(value)
fields = ["terrainScale","stateRepresentationDefect", ...
    "greenIdentityDefect","stationaryRowDefect","strongAPVDefect", ...
    "weakStrongAgreementDefect","weakEvolutionDefect", ...
    "primitiveGeneratorAgreementDefect","energyDefect","bottomDefect", ...
    "bottomResidualNorm","conjugacyDefect","enstrophyDefect", ...
    "energyHermitianDefect","exchangeSkewHermitianDefect", ...
    "energyRcond", ...
    "continuityTangencyDefect","projectionAdjointDefect", ...
    "exactConvolutionDefect"];
summary = emptySequence;
for field = fields
    summary.(field) = value.(field);
end
summary.comparisonEnergyDefect = NaN;
end

function comparison = trustedPhysicalContext(problem,c,coordinate,degree, ...
    comparisonOrder)
[xi,verticalWeight] = physicalQuadrature(comparisonOrder,c.wvt.Lz);
spaces = evaluatePrimitiveVerticalSpaces(c.wvt,degree,coordinate,xi);
Nx = c.layout.oversampledSize(1);
Ny = c.layout.oversampledSize(2);
[x,y] = ndgrid((0:Nx-1)'*c.wvt.Lx/Nx, ...
    (0:Ny-1)'*c.wvt.Ly/Ny);
phase = exp(1i*(x(:)*c.horizontalLayout.k.' ...
    +y(:)*c.horizontalLayout.l.'));
trusted = c.layout.trustedColumns;
raw = c.N(:,trusted);
nColumn = numel(trusted);
nXY = Nx*Ny;
nZ = numel(xi);
uHat = zeros(nXY*nZ,nColumn);
vHat = zeros(nXY*nZ,nColumn);
wHat = zeros(nXY*nZ,nColumn);
eta = zeros(nXY*nZ,nColumn);
for iK = 1:c.nK
    rows = (iK-1)*c.nXBlock+(1:c.nXBlock);
    block = raw(rows,:);
    iu = 1:c.nF;
    iv = c.nF+(1:c.nF);
    iw = 2*c.nF+(1:c.nG);
    ieta = 2*c.nF+c.nG+(1:c.nH);
    uHat = uHat+kron(spaces.F,phase(:,iK))*block(iu,:);
    vHat = vHat+kron(spaces.F,phase(:,iK))*block(iv,:);
    wHat = wHat+kron(spaces.G,phase(:,iK))*block(iw,:);
    eta = eta+kron(spaces.H,phase(:,iK))*block(ieta,:);
end
h = interpft(interpft(problem.topographicHeight,Nx,1),Ny,2);
[hX,hY] = horizontalDerivatives(real(h),c.wvt.Lx,c.wvt.Ly);
xiGrid = kron(xi,ones(nXY,1));
weight = kron(verticalWeight,ones(nXY,1)/nXY);
comparison = struct("xiGrid",xiGrid,"weight",weight, ...
    "h",real(h(:)),"hX",hX(:),"hY",hY(:),"nZ",nZ, ...
    "nXY",nXY,"depth",c.wvt.Lz, ...
    "uHat",uHat,"vHat",vHat,"wHat",wHat,"eta",eta);
[comparison.probeUHat,comparison.probeVHat,comparison.probeWHat, ...
    comparison.probeEta] = admissiblePhysicalProbes(c,xi,phase);
end

function [basis,energyDefect] = trustedPhysicalBasis(c,direction, ...
    scale,comparison)
gamma2 = 1-scale*comparison.h/c.wvt.Lz;
gamma = repmat(gamma2,comparison.nZ,1);
gammaX = repmat(-scale*comparison.hX/c.wvt.Lz,comparison.nZ,1);
gammaY = repmat(-scale*comparison.hY/c.wvt.Lz,comparison.nZ,1);
w = comparison.wHat+comparison.xiGrid.*((gammaX./gamma) ...
    .*comparison.uHat+(gammaY./gamma).*comparison.vHat);
N2 = c.wvt.N2Function(gamma.*comparison.xiGrid);
if isscalar(N2)
    N2 = repmat(N2,size(gamma));
end
weight = comparison.weight;
embedded = sqrt(c.wvt.rho0)*[sqrt(weight./gamma).*comparison.uHat; ...
    sqrt(weight./gamma).*comparison.vHat; ...
    sqrt(weight.*gamma).*w; ...
    sqrt(weight.*gamma.*N2(:)).*comparison.eta];
gram = embedded'*embedded;
trusted = c.layout.trustedColumns;
expected = direction.E(trusted,trusted);
energyDefect = norm(gram-expected,"fro") ...
    /max(norm(gram,"fro")+norm(expected,"fro"),realmin);
stateBasis = orthonormalColumns(embedded);
probes = smoothPhysicalProbes(comparison,gamma,N2(:),scale);
basis = probes'*(stateBasis*(stateBasis'*probes));
end

function probes = smoothPhysicalProbes(comparison,gamma,N2,scale)
weight = comparison.weight;
nGrid = numel(comparison.xiGrid);
probeW = comparison.probeWHat+comparison.xiGrid.*( ...
    (repmat(-scale*comparison.hX/comparison.depth,comparison.nZ,1)./gamma) ...
    .*comparison.probeUHat ...
    +(repmat(-scale*comparison.hY/comparison.depth,comparison.nZ,1)./gamma) ...
    .*comparison.probeVHat);
probes = [sqrt(weight./gamma).*comparison.probeUHat; ...
    sqrt(weight./gamma).*comparison.probeVHat; ...
    sqrt(weight.*gamma).*probeW; ...
    sqrt(weight.*gamma.*N2).*comparison.probeEta];
if size(probes,1) ~= 4*nGrid
    error("WVTerrainEnergyGalerkin:InvalidWKBPrimitiveProbeSize", ...
        "The trusted physical probes have an inconsistent grid size.")
end
probes = orthonormalColumns(probes);
end

function [uHat,vHat,wHat,eta] = admissiblePhysicalProbes(c,xi,phase)
% Fixed mapped-divergence-free probes compare physical subspaces, not raw fields.
t = 2*xi/c.wvt.Lz+1;
F = [ones(size(t)) t (3*t.^2-1)/2 (5*t.^3-3*t)/2];
dF = [zeros(size(t)) 2/c.wvt.Lz+zeros(size(t)) ...
    6*t/c.wvt.Lz (15*t.^2-3)/c.wvt.Lz];
G = (1-t.^2).*F;
dG = (-4*t/c.wvt.Lz).*F+(1-t.^2).*dF;
H = [G(:,1:3) (1-t)/2];
trustedModes = find(c.layout.trustedHorizontalModes);
kappa = hypot(c.horizontalLayout.k(trustedModes), ...
    c.horizontalLayout.l(trustedModes));
iNonzero = find(kappa > 0,1);
if isempty(iNonzero)
    iMode = trustedModes(1);
    horizontalPhase = phase(:,iMode);
    zero = zeros(numel(horizontalPhase)*numel(xi),4);
    uHat = [kron(F,horizontalPhase) zero zero];
    vHat = [zero kron(F,horizontalPhase) zero];
    wHat = [zero zero zero];
    eta = [zero zero kron(H,horizontalPhase)];
    return
end
iMode = trustedModes(iNonzero);
k = c.horizontalLayout.k(iMode);
l = c.horizontalLayout.l(iMode);
kappa2 = k^2+l^2;
horizontalPhase = phase(:,iMode);
zero = zeros(numel(horizontalPhase)*numel(xi),4);
rotationalU = (-l/sqrt(kappa2))*kron(F,horizontalPhase);
rotationalV = ( k/sqrt(kappa2))*kron(F,horizontalPhase);
divergentU = (1i*k/kappa2)*kron(dG,horizontalPhase);
divergentV = (1i*l/kappa2)*kron(dG,horizontalPhase);
uHat = [rotationalU divergentU zero];
vHat = [rotationalV divergentV zero];
wHat = [zero kron(G,horizontalPhase) zero];
eta = [zero zero kron(H,horizontalPhase)];
end

function values = orthonormalColumns(values)
gram = values'*values;
gram = (gram+gram')/2;
factor = chol(gram);
values = values/factor;
end

function value = coefficientRoundTrip(c)
trusted = c.layout.trustedColumns;
R = [c.Ru*c.N(:,trusted);c.Rv*c.N(:,trusted); ...
    c.Rwh*c.N(:,trusted);c.Reta*c.N(:,trusted)];
weight = repmat(c.volumeWeight,4,1);
gram = R'*(weight.*R);
analysis = gram\(R'.*weight.');
coefficient = sin((1:numel(trusted))')+1i*cos((1:numel(trusted))');
recovered = analysis*(R*coefficient);
value = norm(recovered-coefficient)/max(norm(coefficient),realmin);
end

function value = quadratureAdjointDefect(c)
trusted = c.layout.trustedColumns;
R = [c.Ru*c.N(:,trusted);c.Rv*c.N(:,trusted); ...
    c.Rwh*c.N(:,trusted);c.Reta*c.N(:,trusted)];
weight = repmat(c.volumeWeight,4,1);
coefficient = sin((1:numel(trusted))')+1i*cos((1:numel(trusted))');
field = cos((1:size(R,1))'/max(size(R,1),1)) ...
    +1i*sin((1:size(R,1))'/max(size(R,1),1));
left = (R*coefficient)'*(weight.*field);
right = coefficient'*(R'*(weight.*field));
value = abs(left-right)/max(abs(left)+abs(right),realmin);
end

function results = convergenceDiagnostics(results,reference,tolerance)
for iResult = 1:numel(results)
    if ~results(iResult).isQualified
        continue
    end
    defects = zeros(numel(reference.physicalBases),1);
    for iScale = 1:numel(defects)
        defects(iScale) = projectorActionDefect( ...
            results(iResult).physicalBases{iScale}, ...
            reference.physicalBases{iScale});
    end
    results(iResult).maximumProjectorDefect = max(defects);
    results(iResult).projectorDefects = defects;
    results(iResult).passes = casePasses(results(iResult),tolerance);
end
end

function value = projectorActionDefect(candidate,reference)
if isempty(candidate) || isempty(reference) ...
        || ~isequal(size(candidate),size(reference))
    value = Inf;
    return
end
value = norm(candidate-reference,"fro") ...
    /max(norm(candidate,"fro")+norm(reference,"fro"),realmin);
end

function tf = casePasses(value,tolerance)
tf = value.isQualified ...
    && value.verticalDiagnostics.forwardMapDefect <= tolerance.coordinate ...
    && value.verticalDiagnostics.inverseMapDefect <= tolerance.coordinate ...
    && value.verticalDiagnostics.jacobianDefect <= tolerance.coordinate ...
    && value.verticalDiagnostics.endpointDefect <= tolerance.coordinate ...
    && value.verticalDiagnostics.quadratureDefect <= tolerance.coordinate ...
    && value.quadratureAdjointDefect <= tolerance.coordinate ...
    && value.continuityDefect <= tolerance.continuity ...
    && value.conjugacyDefect <= tolerance.conjugacy ...
    && value.maximumStructuralDefect <= tolerance.structure ...
    && value.maximumProjectionDefect <= tolerance.structure ...
    && value.maximumEnergyDefect <= tolerance.energy ...
    && value.maximumGreenDefect <= tolerance.stationary ...
    && value.maximumStationaryDefect <= tolerance.stationary ...
    && value.maximumAPVDefect <= tolerance.apv ...
    && value.maximumBottomDefect <= tolerance.bottom ...
    && value.roundTripDefect <= tolerance.roundTrip ...
    && max([value.sequence.comparisonEnergyDefect]) <= tolerance.projector ...
    && value.maximumProjectorDefect <= tolerance.projector;
end

function tf = controlsPass(controls,tolerance)
if isempty(controls.support) || isempty(controls.padding)
    tf = false;
    return
end
support = controls.support;
padding = controls.padding;
supportDrift = 0;
if numel(support) > 1
    for iScale = 1:numel(support(end).physicalBases)
        supportDrift = max(supportDrift,projectorActionDefect( ...
            support(end-1).physicalBases{iScale}, ...
            support(end).physicalBases{iScale}));
    end
end
tf = support(end).passes && all([padding.passes]) ...
    && controls.twoDimensional.passes ...
    && supportDrift <= tolerance.drift;
end

function tf = verticalCandidatePasses(value,tolerance)
tf = value.isQualified ...
    && value.verticalDiagnostics.forwardMapDefect <= tolerance.coordinate ...
    && value.verticalDiagnostics.inverseMapDefect <= tolerance.coordinate ...
    && value.verticalDiagnostics.jacobianDefect <= tolerance.coordinate ...
    && value.verticalDiagnostics.endpointDefect <= tolerance.coordinate ...
    && value.verticalDiagnostics.quadratureDefect <= tolerance.coordinate ...
    && value.quadratureAdjointDefect <= tolerance.coordinate ...
    && value.continuityDefect <= tolerance.continuity ...
    && value.conjugacyDefect <= tolerance.conjugacy ...
    && value.maximumStructuralDefect <= tolerance.structure ...
    && value.maximumProjectionDefect <= tolerance.structure ...
    && value.maximumEnergyDefect <= tolerance.energy ...
    && value.maximumGreenDefect <= tolerance.stationary ...
    && value.maximumStationaryDefect <= tolerance.stationary ...
    && value.roundTripDefect <= tolerance.roundTrip ...
    && max([value.sequence.comparisonEnergyDefect]) <= tolerance.projector ...
    && value.maximumProjectorDefect <= tolerance.projector;
end

function diagnostics = summarizeControls(controls)
diagnostics = struct("supportProjectorDrift",Inf, ...
    "paddingProjectorDrift",Inf,"maximumSupportDefect",Inf, ...
    "maximumPaddingDefect",Inf,"twoDimensionalPasses",false, ...
    "passes",false);
if isempty(controls.support) || isempty(controls.padding)
    return
end
diagnostics.supportProjectorDrift = maximumAdjacentProjectorDrift( ...
    controls.support);
diagnostics.paddingProjectorDrift = maximumAdjacentProjectorDrift( ...
    controls.padding);
diagnostics.maximumSupportDefect = ...
    max([controls.support.maximumProjectorDefect]);
diagnostics.maximumPaddingDefect = ...
    max([controls.padding.maximumProjectorDefect]);
diagnostics.twoDimensionalPasses = controls.twoDimensional.passes;
diagnostics.passes = controlsPass(controls,requiredTolerance);
end

function value = maximumAdjacentProjectorDrift(cases)
value = 0;
if numel(cases) < 2
    return
end
for iCase = 2:numel(cases)
    if ~cases(iCase-1).isQualified || ~cases(iCase).isQualified ...
            || numel(cases(iCase-1).physicalBases) ...
            ~= numel(cases(iCase).physicalBases)
        value = Inf;
        return
    end
    for iScale = 1:numel(cases(iCase).physicalBases)
        value = max(value,projectorActionDefect( ...
            cases(iCase-1).physicalBases{iScale}, ...
            cases(iCase).physicalBases{iScale}));
    end
end
end

function diagnostics = constantNEquivalenceDiagnostics(legendre,wkb)
commonDegrees = intersect([legendre.degree],[wkb.degree]);
entries = repmat(struct("degree",0,"projectorDefects",zeros(0,1), ...
    "maximumProjectorDefect",Inf),numel(commonDegrees),1);
for iDegree = 1:numel(commonDegrees)
    iLegendre = find([legendre.degree] == commonDegrees(iDegree),1);
    iWKB = find([wkb.degree] == commonDegrees(iDegree),1);
    legendreCase = legendre(iLegendre).controls(1);
    wkbCase = wkb(iWKB).controls(1);
    defects = Inf(max(numel(legendreCase.physicalBases), ...
        numel(wkbCase.physicalBases)),1);
    if legendreCase.isQualified && wkbCase.isQualified ...
            && numel(legendreCase.physicalBases) ...
            == numel(wkbCase.physicalBases)
        for iScale = 1:numel(defects)
            defects(iScale) = max( ...
                projectorActionDefect(legendreCase.physicalBases{iScale}, ...
                wkbCase.physicalBases{iScale}), ...
                projectorActionDefect(wkbCase.physicalBases{iScale}, ...
                legendreCase.physicalBases{iScale}));
        end
    end
    entries(iDegree) = struct("degree",commonDegrees(iDegree), ...
        "projectorDefects",defects, ...
        "maximumProjectorDefect",max(defects));
end
if isempty(entries)
    maximum = Inf;
else
    maximum = max([entries.maximumProjectorDefect]);
end
diagnostics = struct("commonDegrees",commonDegrees(:), ...
    "entries",entries,"maximumProjectorDefect",maximum, ...
    "passes",maximum <= 1e-10);
end

function diagnostics = matchedAccuracyDiagnostics(legendre,legendreIndex, ...
    wkb,wkbIndex,coordinateRatio)
diagnostics = struct("legendrePasses",legendreIndex > 0, ...
    "wkbPasses",wkbIndex > 0,"legendreDegree",0,"wkbDegree",0, ...
    "legendreCoordinates",0,"wkbCoordinates",0, ...
    "coordinateReductionFactor",coordinateRatio, ...
    "legendreBuildSeconds",Inf,"wkbBuildSeconds",Inf, ...
    "wkbToLegendreTimeRatio",Inf,"relativeCoordinateDifference",Inf, ...
    "relativeTimingDifference",Inf,"shouldPreferLegendre",false);
if legendreIndex > 0
    diagnostics.legendreDegree = legendre(legendreIndex).degree;
    diagnostics.legendreCoordinates = ...
        legendre(legendreIndex).numberOfAdmissibleCoordinates;
    diagnostics.legendreBuildSeconds = legendre(legendreIndex).buildSeconds;
end
if wkbIndex > 0
    diagnostics.wkbDegree = wkb(wkbIndex).degree;
    diagnostics.wkbCoordinates = ...
        wkb(wkbIndex).numberOfAdmissibleCoordinates;
    diagnostics.wkbBuildSeconds = wkb(wkbIndex).buildSeconds;
end
if legendreIndex == 0
    return
elseif wkbIndex == 0
    diagnostics.shouldPreferLegendre = true;
    return
end
diagnostics.wkbToLegendreTimeRatio = diagnostics.wkbBuildSeconds ...
    /max(diagnostics.legendreBuildSeconds,realmin);
diagnostics.relativeCoordinateDifference = abs( ...
    diagnostics.legendreCoordinates-diagnostics.wkbCoordinates) ...
    /max(diagnostics.legendreCoordinates,diagnostics.wkbCoordinates);
diagnostics.relativeTimingDifference = abs( ...
    diagnostics.legendreBuildSeconds-diagnostics.wkbBuildSeconds) ...
    /max(diagnostics.legendreBuildSeconds,diagnostics.wkbBuildSeconds);
if diagnostics.wkbCoordinates < 0.9*diagnostics.legendreCoordinates
    diagnostics.shouldPreferLegendre = false;
elseif diagnostics.legendreCoordinates < 0.9*diagnostics.wkbCoordinates
    diagnostics.shouldPreferLegendre = true;
elseif diagnostics.relativeTimingDifference < 0.1
    diagnostics.shouldPreferLegendre = true;
else
    diagnostics.shouldPreferLegendre = ...
        diagnostics.legendreBuildSeconds <= diagnostics.wkbBuildSeconds;
end
end

function tolerance = requiredTolerance
tolerance = struct("coordinate",1e-12,"continuity",1e-11, ...
    "conjugacy",1e-11,"structure",1e-12,"energy",1e-10, ...
    "stationary",1e-10,"apv",1e-10,"bottom",1e-10, ...
    "roundTrip",1e-10,"projector",1e-8,"drift",1e-8);
end

function definitions = standardizedControlDefinitions(problem)
wvt = problem.originatingTransform;
z = linspace(-wvt.Lz,0,257).';
N2 = wvt.N2Function(z);
if isscalar(N2)
    N2 = repmat(N2,size(z));
end
meanN2 = trapz(z,N2)/wvt.Lz;
amplitude = max(abs(problem.topographicHeight),[],"all");
standardizedAmplitude = min(amplitude,wvt.Lz/5000);
if standardizedAmplitude == 0
    standardizedAmplitude = wvt.Lz/5000;
end
definitions = struct("meanN2",meanN2,"terrainAmplitude",amplitude, ...
    "standardizedTerrainAmplitude",standardizedAmplitude, ...
    "constantN2Description","depth-mean input N2", ...
    "exponentialN2Description","mean-normalized exp(2z/D)", ...
    "thermoclineDescription", ...
    "depth-mean-normalized 0.5+1.5 sech^2((z+0.2D)/(0.2D))", ...
    "terrain2DDescription", ...
    "a[cos(2 pi x/Lx)+cos(2 pi y/Ly)]/2");
end

function controls = standardizedControlProblems(problem)
wvt = problem.originatingTransform;
D = wvt.Lz;
zReference = linspace(-D,0,513).';
inputN2 = wvt.N2Function(zReference);
if isscalar(inputN2)
    inputN2 = repmat(inputN2,size(zReference));
end
positive = inputN2(isfinite(inputN2) & inputN2 > 0);
if isempty(positive)
    referenceN2 = max(wvt.f^2,1e-8);
else
    referenceN2 = trapz(zReference,max(inputN2,0))/D;
    if referenceN2 <= 0
        referenceN2 = mean(positive);
    end
end
constantN2 = @(z)referenceN2+0*z;
exponentialNormalization = 2/(1-exp(-2));
exponentialN2 = @(z)referenceN2*exponentialNormalization.*exp(2*z/D);
rawThermocline = @(z)0.5+1.5./cosh((z+0.2*D)/(0.2*D)).^2;
thermoclineMean = trapz(zReference,rawThermocline(zReference))/D;
thermoclineN2 = @(z)referenceN2*rawThermocline(z)/thermoclineMean;

amplitude = min(max(abs(problem.topographicHeight),[],"all"),D/5000);
if amplitude == 0
    amplitude = D/5000;
end
names = ["constant-sinusoidal","exponential-sinusoidal", ...
    "thermocline-sinusoidal","uniform-depth","two-dimensional-terrain"];
profiles = {constantN2,exponentialN2,thermoclineN2,constantN2,constantN2};
terrainKinds = ["sinusoidal","sinusoidal","sinusoidal","uniform","two-dimensional"];
controls = repmat(struct("name","","problem",problem, ...
    "trustedBounds",[Inf Inf],"minimumSupportBounds",[0 0]), ...
    numel(names),1);
for iControl = 1:numel(names)
    if terrainKinds(iControl) == "two-dimensional"
        horizontalResolution = [max(wvt.Nx,8) max(wvt.Ny,8)];
    else
        horizontalResolution = [wvt.Nx wvt.Ny];
    end
    transform = WVTransformBoussinesq( ...
        [wvt.Lx wvt.Ly wvt.Lz],[horizontalResolution wvt.Nz], ...
        N2Function=profiles{iControl},latitude=wvt.latitude, ...
        rho0=wvt.rho0,g=wvt.g,shouldAntialias=wvt.shouldAntialias, ...
        z=wvt.z);
    [x,y] = ndgrid((0:transform.Nx-1)'*transform.Lx/transform.Nx, ...
        (0:transform.Ny-1)'*transform.Ly/transform.Ny);
    if terrainKinds(iControl) == "uniform"
        terrain = amplitude*ones(transform.Nx,transform.Ny);
    elseif terrainKinds(iControl) == "two-dimensional"
        terrain = (amplitude/2)*(cos(2*pi*x/transform.Lx) ...
            +cos(2*pi*y/transform.Ly));
    else
        terrain = amplitude*cos(2*pi*y/transform.Ly);
    end
    controls(iControl).name = names(iControl);
    controls(iControl).problem = WVTerrainEnergyGalerkin.fromTopography( ...
        transform,topographicHeight=terrain, ...
        verticalModeIndices=[0;min(transform.j(transform.j > 0))]);
end
controls(end).trustedBounds = [0 0];
controls(end).minimumSupportBounds = [2 2];
end

function layout = retainedLayout(problem,bounds)
source = problem.horizontalLayout;
mask = abs(source.kMode) <= bounds(1) ...
    & abs(source.lMode) <= bounds(2);
layout = source(mask,:);
expected = (2*bounds(1)+1)*(2*bounds(2)+1);
if height(layout) ~= expected
    error("WVTerrainEnergyGalerkin:UnavailableWKBPrimitiveSupport", ...
        "The originating transform does not retain the complete signed rectangle [%d %d].", ...
        bounds(1),bounds(2))
end
end

function [xi,weight] = physicalQuadrature(order,D)
index = (1:order-1)';
offDiagonal = index./sqrt(4*index.^2-1);
[vectors,values] = eig(diag(offDiagonal,1)+diag(offDiagonal,-1),"vector");
[r,permutation] = sort(values);
vectors = vectors(:,permutation);
weightR = 2*(vectors(1,:)').^2;
xi = D*(r-1)/2;
weight = D*weightR/2;
end

function [xDerivative,yDerivative] = horizontalDerivatives(values,Lx,Ly)
[Nx,Ny] = size(values);
k = 2*pi*[0:floor((Nx-1)/2) -floor(Nx/2):-1]'/Lx;
l = 2*pi*[0:floor((Ny-1)/2) -floor(Ny/2):-1]/Ly;
spectrum = fft2(values);
xDerivative = real(ifft2((1i*k).*spectrum));
yDerivative = real(ifft2(spectrum.*(1i*l)));
end

function value = emptyCase
value = struct("verticalCoordinate","","degree",0, ...
    "supportModeBounds",[0 0],"paddingFactor",0, ...
    "terrainScales",zeros(0,1),"quadratureOrder",0, ...
    "isQualified",false,"qualificationReason","", ...
    "verticalDiagnostics",emptyVerticalDiagnostics, ...
    "numberOfRawCoordinates",0,"numberOfAdmissibleCoordinates",0, ...
    "numberOfPressureCoordinates",0,"continuityDefect",Inf, ...
    "conjugacyDefect",Inf,"roundTripDefect",Inf, ...
    "quadratureAdjointDefect",Inf,"energyConditionNumber",Inf, ...
    "sequence",repmat(emptySequence,0,1),"physicalBases",{{}}, ...
    "maximumStructuralDefect",Inf,"maximumEnergyDefect",Inf, ...
    "maximumStateRepresentationDefect",Inf, ...
    "maximumAPVDefect",Inf,"maximumBottomDefect",Inf, ...
    "maximumGreenDefect",Inf,"maximumStationaryDefect",Inf, ...
    "maximumProjectionDefect",Inf,"maximumProjectorDefect",Inf, ...
    "projectorDefects",zeros(0,1),"passes",false,"buildSeconds",0);
end

function value = emptyVerticalDiagnostics
value = struct("name","","isQualified",false,"reason","", ...
    "conditionNumber",Inf,"forwardMapDefect",Inf, ...
    "inverseMapDefect",Inf,"jacobianDefect",Inf,"endpointDefect",Inf, ...
    "quadratureDefect",Inf,"integratedBuoyancyFrequency",NaN, ...
    "minimumN2",NaN,"maximumN2",NaN);
end

function value = emptySequence
value = struct("terrainScale",0,"stateRepresentationDefect",Inf, ...
    "greenIdentityDefect",Inf,"stationaryRowDefect",Inf, ...
    "strongAPVDefect",Inf,"weakStrongAgreementDefect",Inf, ...
    "weakEvolutionDefect",Inf,"primitiveGeneratorAgreementDefect",Inf, ...
    "energyDefect",Inf,"bottomDefect",Inf,"bottomResidualNorm",Inf, ...
    "conjugacyDefect",Inf,"enstrophyDefect",Inf, ...
    "energyHermitianDefect",Inf,"exchangeSkewHermitianDefect",Inf, ...
    "energyRcond",0, ...
    "continuityTangencyDefect",Inf,"projectionAdjointDefect",Inf, ...
    "exactConvolutionDefect",Inf,"comparisonEnergyDefect",Inf);
end

function value = emptyStandardizedControl
value = struct("name","","legendre",emptyCase,"wkb",emptyCase, ...
    "legendrePasses",false,"wkbPasses",false);
end

function value = emptyStandardizedDegree
value = struct("degree",0,"controls",repmat(emptyCase,0,1), ...
    "passes",false);
end

function value = emptyControls
value = struct("support",repmat(emptyCase,0,1), ...
    "padding",repmat(emptyCase,0,1),"twoDimensional",emptyCase);
end

function [value,record] = cachedCase(problem,directory,configuration,builder)
schema = "wkb-primitive-forward-case-v2";
key = cacheKey(problem,configuration,schema);
path = "";
wasInvalidated = false;
if directory ~= ""
    if ~isfolder(directory)
        mkdir(directory);
    end
    path = string(fullfile(directory,"case-"+key+".mat"));
    if isfile(path)
        try
            loaded = load(path,"entry");
            if loaded.entry.schema == schema && loaded.entry.key == key
                value = loaded.entry.value;
                information = dir(path);
                record = cacheRecord("hit",key,path,information.bytes, ...
                    loaded.entry.buildSeconds,loaded.entry.buildSeconds);
                return
            end
            wasInvalidated = true;
        catch
            % Rebuild an invalid or incomplete disposable cache entry.
            wasInvalidated = true;
        end
    end
end
timer = tic;
value = builder();
seconds = toc(timer);
if directory == ""
    record = cacheRecord("disabled",key,"",0,seconds,0);
    return
end
entry = struct("schema",schema,"key",key,"configuration",configuration, ...
    "buildSeconds",seconds,"value",value);
temporary = string(tempname(directory))+".mat";
cleanup = onCleanup(@()deleteIfPresent(temporary));
save(temporary,"entry","-v7.3");
movefile(temporary,path,"f");
clear cleanup
information = dir(path);
if wasInvalidated
    status = "invalidated";
else
    status = "miss";
end
record = cacheRecord(status,key,path,information.bytes,seconds,0);
end

function key = cacheKey(problem,configuration,schema)
engine = java.security.MessageDigest.getInstance("SHA-256");
updateDigest(engine,uint8(char(schema)));
updateDigest(engine,uint8(char(jsonencode(configuration))));
files = [string(mfilename("fullpath"))+".m"; ...
    fullfile(fileparts(fileparts(mfilename("fullpath"))), ...
    "auditWKBPrimitiveForwardState.m"); ...
    fullfile(fileparts(mfilename("fullpath")), ...
    "buildGlobalSmallTerrainPrimitiveAudit.m"); ...
    fullfile(fileparts(mfilename("fullpath")), ...
    "buildPrimitiveVerticalDiscretization.m"); ...
    fullfile(fileparts(mfilename("fullpath")), ...
    "evaluatePrimitiveVerticalSpaces.m"); ...
    fullfile(fileparts(mfilename("fullpath")), ...
    "legendreValuesForPrimitiveAudit.m"); ...
    fullfile(fileparts(mfilename("fullpath")), ...
    "auditPrimitiveSequenceDirection.m")];
for iFile = 1:numel(files)
    updateDigest(engine,uint8(fileread(files(iFile))));
end
wvt = problem.originatingTransform;
updateDigest(engine,typecast(double([wvt.Lx wvt.Ly wvt.Lz ...
    wvt.Nx wvt.Ny wvt.Nz wvt.f wvt.g wvt.rho0]),"uint8"));
updateDigest(engine,typecast(double(problem.topographicHeight(:)),"uint8"));
z = linspace(-wvt.Lz,0,max(65,4*wvt.Nz+1)).';
N2 = wvt.N2Function(z);
if isscalar(N2)
    N2 = repmat(N2,size(z));
end
updateDigest(engine,typecast(double(N2(:)),"uint8"));
digest = typecast(engine.digest(),"uint8");
key = lower(join(compose("%02x",digest),""));
end

function updateDigest(engine,bytes)
engine.update(bytes);
end

function cache = emptyCacheDiagnostics(directory)
cache = struct("enabled",directory ~= "","directory",directory, ...
    "records",repmat(cacheRecord("","","",0,0,0),0,1));
end

function cache = appendCacheRecord(cache,record)
cache.records(end+1,1) = record;
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
status = [cache.records.status];
cache.numberOfHits = nnz(status == "hit");
cache.numberOfMisses = nnz(status == "miss" | status == "invalidated" ...
    | status == "disabled");
cache.numberOfInvalidations = nnz(status == "invalidated");
cache.bytes = sum([cache.records.bytes]);
cache.secondsSaved = sum([cache.records.secondsSaved]);
end

function record = cacheRecord(status,key,path,bytes,buildSeconds,secondsSaved)
record = struct("status",status,"key",key,"path",path, ...
    "bytes",bytes,"buildSeconds",buildSeconds, ...
    "secondsSaved",secondsSaved);
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
