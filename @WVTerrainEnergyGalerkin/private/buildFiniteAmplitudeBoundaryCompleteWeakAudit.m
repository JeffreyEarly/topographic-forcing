function audit = buildFiniteAmplitudeBoundaryCompleteWeakAudit(problem,trustedBounds,supportBounds,scalarDegree,degrees,paddingFactors,terrainScales,quadratureOrder)
% Build the finite-amplitude boundary-complete primitive weak oracle.

nRefinement = numel(degrees);
nPadding = numel(paddingFactors);
nScale = numel(terrainScales);
details = cell(nRefinement,nPadding,nScale);
flatDetails = cell(nRefinement,nPadding);
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
        [primitive,context] = buildGlobalSmallTerrainPrimitiveAudit(problem, ...
            degree,order,1e-3,horizontalLayout=layout, ...
            paddingFactor=paddingFactors(iPadding),trustedModeBounds=trustedBounds, ...
            rejectTerrainNyquist=true,evaluationScales=terrainScales);
        flatDetails{iRefinement,iPadding} = finiteSequenceAudit( ...
            context,primitive.flatReference,trustedBounds,scalarDegree,0);
        tangent = finiteTangent(context,primitive,trustedBounds,scalarDegree);
        for iScale = 1:nScale
            result = finiteSequenceAudit(context,primitive.evaluatedDirections(iScale), ...
                trustedBounds,scalarDegree,terrainScales(iScale));
            result.tangentRemainder = tangentRemainder( ...
                result,flatDetails{iRefinement,iPadding},tangent,terrainScales(iScale));
            result.projection = primitive.projection;
            details{iRefinement,iPadding,iScale} = result;
            refinement(end+1,1) = summarize(result,degree,order, ...
                supportBounds(iRefinement,:),paddingFactors(iPadding),terrainScales(iScale)); %#ok<AGROW>
        end
    end
end

iFullScale = find(terrainScales == 1,1);
fullScale = cell(nRefinement,nPadding);
for iRefinement = 1:nRefinement
    for iPadding = 1:nPadding
        fullScale{iRefinement,iPadding} = details{iRefinement,iPadding,iFullScale};
    end
end
finest = fullScale{end,1};
flat = flatDetails{end,1};
convergence = convergenceDiagnostics(fullScale(:,1));
padding = paddingDiagnostics(fullScale(end,:));
amplitude = amplitudeDiagnostics(squeeze(details(end,1,:)));
tolerance = struct("projection",1e-12,"structure",1e-12,"weak",1e-10, ...
    "energy",1e-11,"bottom",1e-10,"conjugacy",1e-11, ...
    "representation",1e-10,"green",1e-10,"stationary",1e-10, ...
    "apv",1e-10,"convergentAPV",1e-8,"padding",1e-10, ...
    "minimumReduction",2);

allFull = [fullScale{:}];
allFinest = [fullScale{end,:}];
bottomPasses = all(arrayfun(@(result)result.bottomDefect <= tolerance.bottom ...
    || result.bottomResidualNorm <= 1e-12,allFinest));
structurePasses = max([allFull.projectionAdjointDefect]) <= tolerance.projection ...
    && max([allFull.exactConvolutionDefect]) <= tolerance.projection ...
    && max([allFull.energyHermitianDefect]) <= tolerance.structure ...
    && max([allFull.exchangeSkewHermitianDefect]) <= tolerance.structure ...
    && min([allFull.energyRcond]) > 1e-12 ...
    && max([allFull.weakEvolutionDefect]) <= tolerance.weak ...
    && max([allFull.energyDefect]) <= tolerance.energy ...
    && max([allFull.conjugacyDefect]) <= tolerance.conjugacy ...
    && bottomPasses ...
    && padding.maximumDefect <= tolerance.padding ...
    && flat.greenIdentityDefect <= tolerance.green ...
    && flat.stationaryRowDefect <= tolerance.stationary ...
    && flat.weakEvolutionDefect <= tolerance.weak ...
    && flat.energyDefect <= tolerance.energy ...
    && flat.bottomResidualNorm <= 1e-12;
exactAPVPasses = finest.stateRepresentationDefect <= tolerance.representation ...
    && finest.greenIdentityDefect <= tolerance.green ...
    && finest.stationaryRowDefect <= tolerance.stationary ...
    && finest.weakAPVDefect <= tolerance.apv ...
    && finest.strongAPVDefect <= tolerance.apv ...
    && finest.weakStrongAgreementDefect <= tolerance.apv;
convergentAPVPasses = finest.stateRepresentationDefect <= tolerance.convergentAPV ...
    && finest.greenIdentityDefect <= tolerance.convergentAPV ...
    && finest.stationaryRowDefect <= tolerance.convergentAPV ...
    && finest.strongAPVDefect <= tolerance.convergentAPV ...
    && finest.weakStrongAgreementDefect <= tolerance.convergentAPV ...
    && all(convergence.stateRepresentationDefectReduction ...
        >= tolerance.minimumReduction);
strongDiagnosticsConverge = passesOrDecreases( ...
        convergence.primitiveGeneratorAgreementDefect,tolerance.weak) ...
    && passesOrDecreases(convergence.primitiveWeakEvolutionDefect,tolerance.weak) ...
    && passesOrDecreases(convergence.primitiveEnergyDefect,tolerance.energy);

if ~structurePasses
    status = "implementation-unresolved";
    diagnosis = "A projection, form-structure, primitive-equivalence, energy, bottom, conjugacy, padding, or flat-limit gate failed.";
elseif exactAPVPasses && strongDiagnosticsConverge
    status = "compatible-finite-amplitude-weak-oracle";
    diagnosis = "The finite-amplitude primitive weak evolution preserves energy, derives stationary APV, and converges to the independently evaluated bottom and strong primitive equations.";
elseif convergentAPVPasses && strongDiagnosticsConverge
    status = "convergent-finite-amplitude-weak-oracle";
    diagnosis = "The boundary-complete APV identities converge under independent primitive enrichment at finite terrain amplitude.";
else
    status = "finite-amplitude-weak-sequence-blocker";
    diagnosis = "The primitive energy and bottom identities close, but the finite-amplitude geostrophic inclusion or APV Green identity does not converge to the required tolerance.";
end

audit = struct;
audit.scope = "milestone-7-finite-amplitude-boundary-complete-weak-gate";
audit.status = status;
audit.isCompatible = status == "compatible-finite-amplitude-weak-oracle" ...
    || status == "convergent-finite-amplitude-weak-oracle";
audit.diagnosis = diagnosis;
audit.trustedModeBounds = trustedBounds;
audit.supportModeBounds = supportBounds;
audit.scalarPolynomialDegree = scalarDegree;
audit.primitivePolynomialDegrees = degrees;
audit.paddingFactors = paddingFactors;
audit.terrainScales = terrainScales;
audit.quadratureOrder = quadratureOrder;
audit.refinement = refinement;
audit.convergence = convergence;
audit.padding = padding;
audit.amplitude = amplitude;
audit.flat = flat;
audit.finest = finest;
audit.details = details;
audit.requiredTolerance = tolerance;
audit.nextScope = "terrain-mode-construction-only-if-compatible-and-separately-authorized";
end

function result = finiteSequenceAudit(c,direction,trustedBounds,scalarDegree,terrainScale)
scalar = scalarTestFields(c,trustedBounds,scalarDegree);
[G,representation] = projectGeostrophicState(c,scalar,terrainScale);
H3 = repmat(c.H,c.nZ,1);
gamma = 1-terrainScale*H3;
mass = c.wvt.rho0*(scalar.values'*(c.volumeWeight.*gamma.*scalar.values));
weakRow = G'*direction.E;
strongRow = -c.wvt.rho0*(scalar.values'*(c.volumeWeight.*gamma.*direction.Q));
strongAPV = -(mass\strongRow);
energyFactor = chol((direction.E+direction.E')/2);
weakEnergyGenerator = energyFactor'\(direction.J/energyFactor);
L = energyFactor\(weakEnergyGenerator*energyFactor);
energyGenerator = energyFactor*L/energyFactor;
weakEnergyRow = weakRow/energyFactor;
strongEnergyRow = strongRow/energyFactor;
strongEnergyAPV = strongAPV/energyFactor;

trustedRows = scalar.trustedRows;
trustedColumns = c.layout.trustedColumns;
greenIdentityDefect = pairDefect(weakRow,strongRow,trustedRows,trustedColumns);
stationaryRowDefect = normalizedProduct(G'*direction.J,G,direction.J, ...
    trustedRows,trustedColumns);
energyColumns = 1:size(energyGenerator,2);
weakEvolutionResidual = direction.E*L-direction.J;
stationaryEvolution = (G'*direction.J+G'*weakEvolutionResidual)/energyFactor;
greenEvolutionCorrection = ((strongRow-weakRow)/energyFactor)*energyGenerator;
weakAPVDefect = normalizedResidual(stationaryEvolution,trustedRows,energyColumns, ...
    norm(weakEnergyRow(trustedRows,:),"fro")*norm(energyGenerator,"fro"));
strongAPVEvolution = stationaryEvolution+greenEvolutionCorrection;
strongAPVDefect = normalizedResidual(strongAPVEvolution,trustedRows,energyColumns, ...
    norm(strongEnergyRow(trustedRows,:),"fro")*norm(energyGenerator,"fro"));
weakStrongAgreementDefect = normalizedResidual(greenEvolutionCorrection, ...
    trustedRows,energyColumns, ...
    (norm(weakEnergyRow(trustedRows,:),"fro") ...
    +norm(strongEnergyRow(trustedRows,:),"fro"))*norm(energyGenerator,"fro"));

n = size(L,1);
C = sparse((1:n)',c.coordinateConjugateIndex,1,n,n);
energyResidual = energyGenerator'+energyGenerator;
bottomResidual = direction.B*L-direction.R;
conjugacyResidual = L*C-C*conj(L);
strongAPVCoefficientEvolution = -(mass\strongAPVEvolution);
trustedMass = mass(trustedRows,trustedRows);
trustedEnergyAPV = strongEnergyAPV(trustedRows,:);
trustedAPVEvolution = strongAPVCoefficientEvolution(trustedRows,:);
Z = trustedEnergyAPV'*trustedMass*trustedEnergyAPV;
enstrophyResidual = trustedEnergyAPV'*trustedMass*trustedAPVEvolution ...
    +trustedAPVEvolution'*trustedMass*trustedEnergyAPV;
strongPrimitiveEnergyGenerator = energyFactor*direction.L/energyFactor;
primitiveAgreement = strongPrimitiveEnergyGenerator-energyGenerator;
primitiveRawAgreement = direction.L-L;
massHermitian = (mass+mass')/2;

result = struct;
result.terrainScale = terrainScale;
result.scalar = scalar.diagnostics;
result.G = G;
result.massMatrix = mass;
result.weakAPVRow = weakRow;
result.strongAPVRow = strongRow;
result.strongAPVMap = strongAPV;
result.generator = L;
result.strongPrimitiveGenerator = direction.L;
result.energyGenerator = energyGenerator;
result.energyMatrix = direction.E;
result.exchangeMatrix = direction.J;
result.stateRepresentationDefect = representation.trustedMaximumDefect;
result.stateRepresentation = representation;
result.greenIdentityDefect = greenIdentityDefect;
result.stationaryRowDefect = stationaryRowDefect;
result.weakAPVDefect = weakAPVDefect;
result.strongAPVDefect = strongAPVDefect;
result.weakStrongAgreementDefect = weakStrongAgreementDefect;
result.weakEvolutionDefect = productDefect(weakEvolutionResidual,{direction.E*L,direction.J});
result.primitiveGeneratorAgreementDefect = norm( ...
    primitiveRawAgreement(:,trustedColumns),"fro") ...
    /max(norm(L(:,trustedColumns),"fro") ...
    +norm(direction.L(:,trustedColumns),"fro"),realmin);
result.fullSupportPrimitiveGeneratorAgreementDefect = norm(primitiveAgreement,"fro") ...
    /max(norm(energyGenerator,"fro")+norm(strongPrimitiveEnergyGenerator,"fro"),realmin);
result.energyDefect = productDefect(energyResidual,{energyGenerator',energyGenerator});
result.bottomDefect = bottomProductDefect(direction.B,L,direction.R, ...
    bottomResidual,trustedColumns);
result.fullSupportBottomDefect = bottomProductDefect(direction.B,L, ...
    direction.R,bottomResidual,1:size(L,2));
result.bottomResidualNorm = norm(bottomResidual(:,trustedColumns),"fro");
result.bottomLeftNorm = norm(direction.B*L(:,trustedColumns),"fro");
result.bottomRightNorm = norm(direction.R(:,trustedColumns),"fro");
result.conjugacyDefect = norm(conjugacyResidual,"fro")/max(2*norm(L,"fro"),realmin);
result.enstrophyDefect = norm(enstrophyResidual,"fro") ...
    /max(2*norm(Z,"fro")*norm(energyGenerator,"fro"),realmin);
result.energyHermitianDefect = norm(direction.E-direction.E',"fro")/max(norm(direction.E,"fro"),realmin);
result.exchangeSkewHermitianDefect = norm(direction.J+direction.J',"fro")/max(norm(direction.J,"fro"),realmin);
result.energyRcond = rcond((direction.E+direction.E')/2);
result.massHermitianDefect = norm(mass-mass',"fro")/max(norm(mass,"fro"),realmin);
result.massRcond = rcond(massHermitian);
result.geostrophicBottomTangencyDefect = norm(direction.R*G(:,trustedRows),"fro") ...
    /max(norm(direction.R,"fro")*norm(G(:,trustedRows),"fro"),realmin);
result.continuityTangencyDefect = direction.diagnostics.continuityTangencyDefect;
result.saddleResidual = direction.diagnostics.gaugedSaddleResidual;
result.primitiveWeakEvolutionDefect = direction.diagnostics.weakEvolutionDefect;
result.primitiveEnergyDefect = direction.diagnostics.energyDefect;
result.primitiveAPVDefect = direction.diagnostics.apvDefect;
result.primitiveBottomDefect = direction.diagnostics.bottomDefect;
result.primitiveEnstrophyDefect = direction.diagnostics.enstrophyDefect;
result.projectionAdjointDefect = c.projectionDiagnostics.adjointDefect;
result.exactConvolutionDefect = c.projectionDiagnostics.maximumExactConvolutionDefect;
end

function tangent = finiteTangent(c,primitive,trustedBounds,scalarDegree)
scalar = scalarTestFields(c,trustedBounds,scalarDegree);
step = 1e-5;
[Gplus,~] = projectGeostrophicState(c,scalar,step);
[Gminus,~] = projectGeostrophicState(c,scalar,-step);
tangent = struct("energyMatrix",primitive.analyticTangent.E, ...
    "exchangeMatrix",primitive.analyticTangent.J, ...
    "G",(Gplus-Gminus)/(2*step));
end

function diagnostics = tangentRemainder(result,flat,tangent,scale)
diagnostics = struct( ...
    "energy",blockDefect(result.energyMatrix, ...
    flat.energyMatrix+scale*tangent.energyMatrix), ...
    "exchange",blockDefect(result.exchangeMatrix, ...
    flat.exchangeMatrix+scale*tangent.exchangeMatrix), ...
    "geostrophicInclusion",blockDefect(result.G,flat.G+scale*tangent.G));
end

function scalar = scalarTestFields(c,trustedBounds,scalarDegree)
n = (0:scalarDegree)';
bottomValue = (-1).^n;
topDerivative = n.*(n+1)/c.wvt.Lz;
constraint = [bottomValue.';topDerivative.'];
[~,singularValues,V] = svd(constraint);
tolerance = max(size(constraint))*eps(max(diag(singularValues)));
rankConstraint = nnz(diag(singularValues) > tolerance);
trustedCoefficients = canonicalColumnPhases(V(:,rankConstraint+1:end));
coefficients = zeros(c.nF,size(trustedCoefficients,2));
coefficients(1:scalarDegree+1,:) = trustedCoefficients;
S = c.spaces.F*coefficients;
Sxi = c.spaces.Fxi*coefficients;
nS = size(S,2);
nScalar = c.nK*nS;
nGrid = c.nXY*c.nZ;
values = zeros(nGrid,nScalar);
valuesXi = zeros(nGrid,nScalar);
valuesX = zeros(nGrid,nScalar);
valuesY = zeros(nGrid,nScalar);
trustedRows = zeros(0,1);
trustedHorizontal = abs(c.horizontalLayout.kMode) <= trustedBounds(1) ...
    & abs(c.horizontalLayout.lMode) <= trustedBounds(2);
for iK = 1:c.nK
    columns = (iK-1)*nS+(1:nS);
    values(:,columns) = kron(S,c.phase(:,iK));
    valuesXi(:,columns) = kron(Sxi,c.phase(:,iK));
    valuesX(:,columns) = 1i*c.horizontalLayout.k(iK)*values(:,columns);
    valuesY(:,columns) = 1i*c.horizontalLayout.l(iK)*values(:,columns);
    if trustedHorizontal(iK)
        trustedRows = [trustedRows columns]; %#ok<AGROW>
    end
end
surfaceDerivative = topDerivative.'*trustedCoefficients;
bottomScalar = bottomValue.'*trustedCoefficients;
diagnostics = struct("numberOfVerticalFunctions",nS, ...
    "numberOfScalarTests",nScalar,"numberOfTrustedTests",numel(trustedRows), ...
    "trustedScalarDegree",scalarDegree,"constraintRank",rankConstraint, ...
    "bottomValueDefect",norm(bottomScalar,"fro")/max(norm(coefficients,"fro"),realmin), ...
    "surfaceDerivativeDefect",norm(surfaceDerivative,"fro")/max(norm(coefficients,"fro"),realmin), ...
    "boundaryDefect",max(norm(bottomScalar,"fro"),norm(surfaceDerivative,"fro")) ...
        /max(norm(coefficients,"fro"),realmin));
scalar = struct("values",values,"valuesXi",valuesXi,"valuesX",valuesX, ...
    "valuesY",valuesY,"trustedRows",trustedRows,"diagnostics",diagnostics);
end

function [G,diagnostics] = projectGeostrophicState(c,scalar,terrainScale)
H3 = repmat(c.H,c.nZ,1);
HX3 = repmat(c.HX,c.nZ,1);
HY3 = repmat(c.HY,c.nZ,1);
gamma = 1-terrainScale*H3;
gammaX = -terrainScale*HX3;
gammaY = -terrainScale*HY3;
N2 = c.wvt.N2Function(gamma.*c.xiGrid);
if isscalar(N2)
    N2 = repmat(N2,size(c.xiGrid));
end
N2 = N2(:);
u = -gamma.*scalar.valuesY+c.xiGrid.*gammaY.*scalar.valuesXi;
v = gamma.*scalar.valuesX-c.xiGrid.*gammaX.*scalar.valuesXi;
w = c.xiGrid.*(gammaX.*scalar.valuesY-gammaY.*scalar.valuesX);
eta = -c.wvt.f*scalar.valuesXi./(gamma.*N2);

[uColumns,vColumns,wColumns,etaColumns] = componentColumns(c);
projectU = weightedProjector(c.Ru(:,uColumns),c.volumeWeight);
projectV = weightedProjector(c.Rv(:,vColumns),c.volumeWeight);
projectW = weightedProjector(c.Rwh(:,wColumns),c.volumeWeight);
projectEta = weightedProjector(c.Reta(:,etaColumns),c.volumeWeight);
raw = zeros(c.nX,size(u,2));
raw(uColumns,:) = projectU*u;
raw(vColumns,:) = projectV*v;
raw(wColumns,:) = projectW*w;
raw(etaColumns,:) = projectEta*eta;
G = c.N'*raw;
representedRaw = c.N*G;
represented = [c.Ru*representedRaw,c.Rv*representedRaw, ...
    c.Rwh*representedRaw,c.Reta*representedRaw];
targets = [u,v,w,eta];
trustedByField = arrayfun(@(iField)(iField-1)*size(u,2)+scalar.trustedRows, ...
    1:4,"UniformOutput",false);
trustedColumns = horzcat(trustedByField{:});
diagnostics = struct("rawCoordinateDefect",blockDefect(representedRaw,raw), ...
    "maximumDefect",blockDefect(represented,targets), ...
    "trustedMaximumDefect",blockDefect( ...
        represented(:,trustedColumns),targets(:,trustedColumns)), ...
    "continuityDefect",norm(c.continuity*representedRaw,"fro") ...
        /max(norm(c.continuity,"fro")*norm(representedRaw,"fro"),realmin));
end

function projector = weightedProjector(reconstruction,weight)
mass = reconstruction'*(weight.*reconstruction);
projector = mass\(reconstruction'.*weight.');
end

function [uColumns,vColumns,wColumns,etaColumns] = componentColumns(c)
uColumns = zeros(c.nK*c.nF,1);
vColumns = zeros(c.nK*c.nF,1);
wColumns = zeros(c.nK*c.nG,1);
etaColumns = zeros(c.nK*c.nH,1);
for iK = 1:c.nK
    block = (iK-1)*c.nXBlock;
    uColumns((iK-1)*c.nF+(1:c.nF)) = block+(1:c.nF);
    vColumns((iK-1)*c.nF+(1:c.nF)) = block+c.nF+(1:c.nF);
    wColumns((iK-1)*c.nG+(1:c.nG)) = block+2*c.nF+(1:c.nG);
    etaColumns((iK-1)*c.nH+(1:c.nH)) = block+2*c.nF+c.nG+(1:c.nH);
end
end

function diagnostics = convergenceDiagnostics(results)
fields = ["stateRepresentationDefect","greenIdentityDefect", ...
    "stationaryRowDefect","strongAPVDefect","weakStrongAgreementDefect", ...
    "primitiveGeneratorAgreementDefect","primitiveWeakEvolutionDefect", ...
    "primitiveEnergyDefect","bottomDefect"];
diagnostics = struct;
for field = fields
    values = reshape(cellfun(@(result)result.(field),results),1,[]);
    diagnostics.(field) = values;
    diagnostics.(field+"Reduction") = values(1:end-1)./max(values(2:end),realmin);
end
end

function diagnostics = paddingDiagnostics(results)
reference = results{1};
fields = ["generator","G","weakAPVRow","strongAPVRow"];
diagnostics = struct;
maximum = 0;
for field = fields
    values = zeros(numel(results)-1,1);
    for iResult = 2:numel(results)
        values(iResult-1) = blockDefect(results{iResult}.(field),reference.(field));
    end
    diagnostics.(field+"Defect") = values;
    maximum = max(maximum,max(values));
end
diagnostics.maximumDefect = maximum;
end

function diagnostics = amplitudeDiagnostics(results)
scales = cellfun(@(result)result.terrainScale,results);
fields = ["energy","exchange","geostrophicInclusion"];
diagnostics = struct("terrainScales",scales);
for field = fields
    values = cellfun(@(result)result.tangentRemainder.(field),results);
    diagnostics.(field+"Remainder") = values;
    diagnostics.(field+"Reduction") = values(2:end)./max(values(1:end-1),realmin);
end
end

function summary = summarize(result,degree,order,support,padding,scale)
summary = emptySummary;
summary.primitivePolynomialDegree = degree;
summary.quadratureOrder = order;
summary.supportModeBounds = support;
summary.paddingFactor = padding;
summary.terrainScale = scale;
summary.stateRepresentationDefect = result.stateRepresentationDefect;
summary.greenIdentityDefect = result.greenIdentityDefect;
summary.stationaryRowDefect = result.stationaryRowDefect;
summary.weakAPVDefect = result.weakAPVDefect;
summary.strongAPVDefect = result.strongAPVDefect;
summary.weakStrongAgreementDefect = result.weakStrongAgreementDefect;
summary.weakEvolutionDefect = result.weakEvolutionDefect;
summary.primitiveGeneratorAgreementDefect = result.primitiveGeneratorAgreementDefect;
summary.energyDefect = result.energyDefect;
summary.bottomDefect = result.bottomDefect;
summary.conjugacyDefect = result.conjugacyDefect;
summary.enstrophyDefect = result.enstrophyDefect;
end

function summary = emptySummary
summary = struct("primitivePolynomialDegree",0,"quadratureOrder",0, ...
    "supportModeBounds",[0 0],"paddingFactor",0,"terrainScale",0, ...
    "stateRepresentationDefect",NaN,"greenIdentityDefect",NaN, ...
    "stationaryRowDefect",NaN,"weakAPVDefect",NaN,"strongAPVDefect",NaN, ...
    "weakStrongAgreementDefect",NaN,"weakEvolutionDefect",NaN, ...
    "primitiveGeneratorAgreementDefect",NaN,"energyDefect",NaN, ...
    "bottomDefect",NaN,"conjugacyDefect",NaN,"enstrophyDefect",NaN);
end

function value = pairDefect(first,second,rows,columns)
value = norm(first(rows,columns)-second(rows,columns),"fro") ...
    /max(norm(first(rows,columns),"fro")+norm(second(rows,columns),"fro"),realmin);
end

function value = normalizedProduct(residual,left,right,rows,columns)
value = norm(residual(rows,columns),"fro") ...
    /max(norm(left(:,rows),"fro")*norm(right(:,columns),"fro"),realmin);
end

function value = normalizedResidual(residual,rows,columns,naturalScale)
value = norm(residual(rows,columns),"fro")/max(naturalScale,realmin);
end

function value = productDefect(residual,terms)
scale = 0;
for iTerm = 1:numel(terms)
    scale = scale+norm(terms{iTerm},"fro");
end
value = norm(residual,"fro")/max(scale,realmin);
end

function value = bottomProductDefect(B,L,R,residual,columns)
scale = norm(B*L(:,columns),"fro")+norm(R(:,columns),"fro");
naturalScale = norm(B,"fro")*norm(L(:,columns),"fro")+norm(R(:,columns),"fro");
value = norm(residual(:,columns),"fro")/max([scale 100*eps*naturalScale realmin]);
end

function value = blockDefect(first,second)
value = norm(first-second,"fro")/max(norm(first,"fro")+norm(second,"fro"),realmin);
end

function tf = passesOrDecreases(values,tolerance)
tf = values(end) <= tolerance || values(end) < values(1);
end

function values = canonicalColumnPhases(values)
for iColumn = 1:size(values,2)
    [~,pivot] = max(abs(values(:,iColumn)));
    if values(pivot,iColumn) ~= 0
        values(:,iColumn) = values(:,iColumn)*conj(values(pivot,iColumn)) ...
            /abs(values(pivot,iColumn));
    end
end
end

function layout = retainedLayout(problem,bounds)
source = problem.horizontalLayout;
mask = abs(source.kMode) <= bounds(1) & abs(source.lMode) <= bounds(2);
layout = source(mask,:);
expected = (2*bounds(1)+1)*(2*bounds(2)+1);
if height(layout) ~= expected
    error("WVTerrainEnergyGalerkin:UnavailableFiniteAmplitudeSupport", ...
        "The originating transform does not retain the complete signed rectangle [%d %d].", ...
        bounds(1),bounds(2))
end
end

function terrain = terrainSupport(problem)
wvt = problem.originatingTransform;
spectrum = fft2(problem.topographicHeight)/(wvt.Nx*wvt.Ny);
nyquist = logical(WVGeometryDoublyPeriodic.maskForNyquistModes(wvt.Nx,wvt.Ny));
scale = max(abs(spectrum),[],"all");
tolerance = 100*eps*max(scale,1);
if max(abs(spectrum(nyquist)),[],"all") > tolerance
    error("WVTerrainEnergyGalerkin:FiniteAmplitudeTerrainNyquist", ...
        "The finite-amplitude weak oracle requires terrain with no material Nyquist coefficient.")
end
[kMode,lMode] = ndgrid(wvt.kMode_dft,wvt.lMode_dft);
active = abs(spectrum) > tolerance & ~nyquist;
terrain = struct("kMode",kMode(active),"lMode",lMode(active), ...
    "coefficient",spectrum(active),"spectralTolerance",tolerance);
end

function validateGuard(layout,trustedBounds,terrain)
for iK = find(abs(layout.kMode) <= trustedBounds(1) ...
        & abs(layout.lMode) <= trustedBounds(2)).'
    destinations = [layout.kMode(iK)+terrain.kMode layout.lMode(iK)+terrain.lMode];
    for iDestination = 1:size(destinations,1)
        if ~any(layout.kMode == destinations(iDestination,1) ...
                & layout.lMode == destinations(iDestination,2))
            error("WVTerrainEnergyGalerkin:InsufficientFiniteAmplitudeGuard", ...
                "supportModeBounds must retain every first terrain sideband of the trusted band.")
        end
    end
end
end
