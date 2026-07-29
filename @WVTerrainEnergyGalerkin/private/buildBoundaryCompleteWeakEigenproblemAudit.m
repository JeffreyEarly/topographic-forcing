function audit = buildBoundaryCompleteWeakEigenproblemAudit(problem,trustedBounds,supportBounds,degrees,paddingFactor,quadratureOrder)
% Build the boundary-complete weak terrain eigenproblem oracle.

layout = retainedLayout(problem,supportBounds);
terrain = terrainSupport(problem);
validateGuard(layout,trustedBounds,terrain);

refinement = repmat(emptyRefinement,0,1);
details = cell(numel(degrees),1);
for iDegree = 1:numel(degrees)
    degree = degrees(iDegree);
    if isempty(quadratureOrder)
        order = max(2*degree+7,18);
    else
        order = quadratureOrder;
    end
    [primitive,context] = buildGlobalSmallTerrainPrimitiveAudit(problem, ...
        degree,order,1e-3,horizontalLayout=layout,paddingFactor=paddingFactor, ...
        trustedModeBounds=trustedBounds,rejectTerrainNyquist=true);
    details{iDegree} = weakSequenceAudit(context,primitive,trustedBounds,degrees(1));
    refinement(end+1,1) = summarize(details{iDegree},primitive); %#ok<AGROW>
end

finest = details{end};
finestSummary = refinement(end);
convergence = convergenceDiagnostics(refinement);
tolerance = struct("scalarBoundary",1e-12,"stateRepresentation",1e-10, ...
    "greenIdentity",1e-10,"stationaryRows",1e-10,"energy",1e-11, ...
    "bottom",1e-11,"conjugacy",1e-11,"eigenproblem",1e-10, ...
    "modalAPV",1e-8,"refinementFactor",2);
structurePasses = max([refinement.scalarBoundaryDefect]) <= tolerance.scalarBoundary ...
    && finestSummary.trustedStateRepresentationDefect <= tolerance.stateRepresentation ...
    && finest.flatGreenIdentityDefect <= tolerance.greenIdentity ...
    && finest.flatStationaryRowDefect <= tolerance.stationaryRows ...
    && max([refinement.energyDefect]) <= tolerance.energy ...
    && max([refinement.bottomDefect]) <= tolerance.bottom ...
    && max([refinement.conjugacyDefect]) <= tolerance.conjugacy ...
    && finest.modeDiagnostics.eigenproblemDefect <= tolerance.eigenproblem;
tangentPasses = finest.trustedTangentGreenIdentityDefect <= tolerance.greenIdentity ...
    && finest.trustedTangentStationaryRowDefect <= tolerance.stationaryRows ...
    && finest.modeDiagnostics.trustedWaveAPVDefect <= tolerance.modalAPV;
convergencePasses = all(convergence.greenIdentityReduction >= tolerance.refinementFactor) ...
    && all(convergence.stationaryRowReduction >= tolerance.refinementFactor);

if ~structurePasses
    status = "implementation-unresolved";
    diagnosis = "The derived geostrophic test sequence failed a boundary, representation, flat Green-identity, energy, bottom, conjugacy, or eigenproblem gate.";
elseif tangentPasses
    status = "compatible-boundary-complete-weak-oracle";
    diagnosis = "The terrain-dependent geostrophic test sequence makes APV a consequence of the primitive weak equations on the trusted band.";
elseif convergencePasses
    status = "convergent-boundary-complete-weak-oracle";
    diagnosis = "The finite weak sequence is not closed at the finest degree, but its terrain Green-identity and stationary-row defects converge under vertical refinement.";
else
    status = "weak-sequence-representation-blocker";
    diagnosis = "The exact continuous Green identity is established, but the current primitive polynomial spaces do not represent its terrain-dependent geostrophic test sequence with the required convergence.";
end

audit = struct;
audit.scope = "boundary-complete-weak-terrain-eigenproblem-tangent-only";
audit.status = status;
audit.isCompatible = status == "compatible-boundary-complete-weak-oracle" ...
    || status == "convergent-boundary-complete-weak-oracle";
audit.diagnosis = diagnosis;
audit.trustedModeBounds = trustedBounds;
audit.supportModeBounds = supportBounds;
audit.polynomialDegrees = degrees;
audit.paddingFactor = paddingFactor;
audit.quadratureOrder = quadratureOrder;
audit.terrain = terrain;
audit.refinement = refinement;
audit.convergence = convergence;
audit.finest = finest;
audit.requiredTolerance = tolerance;
audit.nextScope = "milestone-7-only-if-compatible";
end

function result = weakSequenceAudit(c,primitive,trustedBounds,scalarDegree)
[scalar,fields] = scalarTestFields(c,trustedBounds,scalarDegree);
[G0,G1,representation] = projectGeostrophicStates(c,fields,scalar);

flat = primitive.flatReference;
first = primitive.analyticTangent;
H3 = repmat(c.H,c.nZ,1);
direct0 = -c.wvt.rho0*(scalar.values'*(c.volumeWeight.*flat.Q));
direct1 = -c.wvt.rho0*(scalar.values'*(c.volumeWeight.*(first.Q-H3.*flat.Q)));
weak0 = G0'*flat.E;
weak1 = G1'*flat.E+G0'*first.E;

trustedRows = scalar.apvRows;
trustedColumns = c.layout.trustedColumns;
flatGreen = pairDefect(weak0,direct0,trustedRows,trustedColumns);
tangentGreen = pairDefect(weak1,direct1,trustedRows,trustedColumns);
greenTerms = struct( ...
    "flatWeakNorm",norm(weak0(trustedRows,trustedColumns),"fro"), ...
    "flatDirectNorm",norm(direct0(trustedRows,trustedColumns),"fro"), ...
    "tangentWeakNorm",norm(weak1(trustedRows,trustedColumns),"fro"), ...
    "tangentDirectNorm",norm(direct1(trustedRows,trustedColumns),"fro"), ...
    "tangentDifferenceNorm",norm( ...
        weak1(trustedRows,trustedColumns)-direct1(trustedRows,trustedColumns),"fro"));
tangentStationary0 = G1'*flat.J;
tangentStationary1 = G0'*first.J;
tangentStationaryDifference = tangentStationary0-tangentStationary1;
flatStationary = normalizedProduct(G0'*flat.J,G0,flat.J,trustedRows,trustedColumns);
tangentStationaryMatrix = tangentStationary0+tangentStationary1;
tangentStationary = normalizedSum(tangentStationaryMatrix, ...
    {tangentStationary0,tangentStationary1},trustedRows,trustedColumns);
stationaryTerms = struct( ...
    "stateDerivativeNorm",norm(tangentStationary0(trustedRows,trustedColumns),"fro"), ...
    "formDerivativeNorm",norm(tangentStationary1(trustedRows,trustedColumns),"fro"), ...
    "sumNorm",norm(tangentStationaryMatrix(trustedRows,trustedColumns),"fro"), ...
    "differenceNorm",norm(tangentStationaryDifference(trustedRows,trustedColumns),"fro"));

weakAPVResidual = weak0*first.L+weak1*flat.L;
directAPVResidual = direct0*first.L+direct1*flat.L;
weakAPVDefect = normalizedSum(weakAPVResidual, ...
    {weak0*first.L,weak1*flat.L},trustedRows,trustedColumns);
directAPVDefect = normalizedSum(directAPVResidual, ...
    {direct0*first.L,direct1*flat.L},trustedRows,trustedColumns);
weakStrongAgreement = normalizedSum(weakAPVResidual-directAPVResidual, ...
    {weak0*first.L,weak1*flat.L,direct0*first.L,direct1*flat.L}, ...
    trustedRows,trustedColumns);

modeDiagnostics = tangentModeDiagnostics(flat,first,G0,G1,direct0,direct1, ...
    trustedRows,trustedColumns);
result = struct;
result.scalar = scalar.diagnostics;
result.representation = representation;
result.G0 = G0;
result.G1 = G1;
result.weakAPVMap0 = weak0;
result.weakAPVMap1 = weak1;
result.directAPVMap0 = direct0;
result.directAPVMap1 = direct1;
result.flatGreenIdentityDefect = flatGreen;
result.trustedTangentGreenIdentityDefect = tangentGreen;
result.greenTerms = greenTerms;
result.flatStationaryRowDefect = flatStationary;
result.trustedTangentStationaryRowDefect = tangentStationary;
result.stationaryTerms = stationaryTerms;
result.trustedWeakAPVDefect = weakAPVDefect;
result.trustedDirectAPVDefect = directAPVDefect;
result.trustedWeakStrongAgreementDefect = weakStrongAgreement;
result.modeDiagnostics = modeDiagnostics;
end

function [scalar,fields] = scalarTestFields(c,trustedBounds,scalarDegree)
nF = c.nF;
n = (0:scalarDegree)';
bottomValue = (-1).^n;
topDerivative = n.*(n+1)/c.wvt.Lz;
constraint = [bottomValue.';topDerivative.'];
[~,singularValues,V] = svd(constraint);
tolerance = max(size(constraint))*eps(max(diag(singularValues)));
rankConstraint = nnz(diag(singularValues) > tolerance);
trustedCoefficients = canonicalColumnPhases(V(:,rankConstraint+1:end));
scalarCoefficients = zeros(nF,size(trustedCoefficients,2));
scalarCoefficients(1:scalarDegree+1,:) = trustedCoefficients;
S = c.spaces.F*scalarCoefficients;
Sxi = c.spaces.Fxi*scalarCoefficients;
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
apvHorizontal = trustedHorizontal;
terrainModes = c.projectionDiagnostics.terrainModes;
for iK = find(trustedHorizontal).'
    destinations = [c.horizontalLayout.kMode(iK)+terrainModes(:,1) ...
        c.horizontalLayout.lMode(iK)+terrainModes(:,2)];
    for iDestination = 1:size(destinations,1)
        apvHorizontal = apvHorizontal ...
            | (c.horizontalLayout.kMode == destinations(iDestination,1) ...
            & c.horizontalLayout.lMode == destinations(iDestination,2));
    end
end
apvRows = zeros(0,1);
for iK = 1:c.nK
    columns = (iK-1)*nS+(1:nS);
    values(:,columns) = kron(S,c.phase(:,iK));
    valuesXi(:,columns) = kron(Sxi,c.phase(:,iK));
    valuesX(:,columns) = 1i*c.horizontalLayout.k(iK)*values(:,columns);
    valuesY(:,columns) = 1i*c.horizontalLayout.l(iK)*values(:,columns);
    if trustedHorizontal(iK)
        trustedRows = [trustedRows columns]; %#ok<AGROW>
    end
    if apvHorizontal(iK)
        apvRows = [apvRows columns]; %#ok<AGROW>
    end
end
H3 = repmat(c.H,c.nZ,1);
HX3 = repmat(c.HX,c.nZ,1);
HY3 = repmat(c.HY,c.nZ,1);
N20 = kron(c.N20,ones(c.nXY,1));
N2Xi = kron(c.N2Xi,ones(c.nXY,1));
stratificationFactor = H3.*(1+c.xiGrid.*N2Xi./N20);

fields = struct;
fields.u0 = -valuesY;
fields.v0 = valuesX;
fields.w0 = zeros(size(values));
fields.eta0 = -c.wvt.f*valuesXi./N20;
fields.u1 = H3.*valuesY-c.xiGrid.*HY3.*valuesXi;
fields.v1 = -H3.*valuesX+c.xiGrid.*HX3.*valuesXi;
fields.w1 = c.xiGrid.*(-HX3.*valuesY+HY3.*valuesX);
fields.eta1 = stratificationFactor.*fields.eta0;

surfaceDerivative = topDerivative.'*trustedCoefficients;
bottomScalar = bottomValue.'*trustedCoefficients;
diagnostics = struct("numberOfVerticalFunctions",nS, ...
    "numberOfScalarTests",nScalar,"numberOfTrustedTests",numel(trustedRows), ...
    "numberOfAPVTests",numel(apvRows),"trustedScalarDegree",scalarDegree, ...
    "constraintRank",rankConstraint, ...
    "bottomValueDefect",norm(bottomScalar,"fro")/max(norm(scalarCoefficients,"fro"),realmin), ...
    "surfaceDerivativeDefect",norm(surfaceDerivative,"fro")/max(norm(scalarCoefficients,"fro"),realmin), ...
    "boundaryDefect",max(norm(bottomScalar,"fro"),norm(surfaceDerivative,"fro")) ...
        /max(norm(scalarCoefficients,"fro"),realmin));
scalar = struct("values",values,"trustedRows",trustedRows,"apvRows",apvRows, ...
    "verticalCoefficients",scalarCoefficients,"diagnostics",diagnostics);
end

function [G0,G1,diagnostics] = projectGeostrophicStates(c,fields,scalar)
[uColumns,vColumns,wColumns,etaColumns] = componentColumns(c);
projectU = weightedProjector(c.Ru(:,uColumns),c.volumeWeight);
projectV = weightedProjector(c.Rv(:,vColumns),c.volumeWeight);
projectW = weightedProjector(c.Rwh(:,wColumns),c.volumeWeight);
projectEta = weightedProjector(c.Reta(:,etaColumns),c.volumeWeight);

[G0,flat] = projectState(fields.u0,fields.v0,fields.w0,fields.eta0);
[G1,terrain] = projectState(fields.u1,fields.v1,fields.w1,fields.eta1);
diagnostics = struct("flat",flat,"terrain",terrain, ...
    "trustedMaximumDefect",max(flat.trustedMaximumDefect,terrain.trustedMaximumDefect));

    function [coordinates,result] = projectState(u,v,w,eta)
        raw = zeros(c.nX,size(u,2));
        raw(uColumns,:) = projectU*u;
        raw(vColumns,:) = projectV*v;
        raw(wColumns,:) = projectW*w;
        raw(etaColumns,:) = projectEta*eta;
        coordinates = c.N'*raw;
        representedRaw = c.N*coordinates;
        represented = [c.Ru*representedRaw,c.Rv*representedRaw, ...
            c.Rwh*representedRaw,c.Reta*representedRaw];
        targets = [u,v,w,eta];
        fieldDefects = zeros(4,1);
        trustedFieldDefects = zeros(4,1);
        for iField = 1:4
            fieldColumns = (iField-1)*size(u,2)+(1:size(u,2));
            fieldDefects(iField) = blockDefect(represented(:,fieldColumns),targets(:,fieldColumns));
            trusted = scalar.apvRows;
            trustedFieldDefects(iField) = blockDefect( ...
                represented(:,(iField-1)*size(u,2)+trusted),targets(:,(iField-1)*size(u,2)+trusted));
        end
        trustedByField = arrayfun(@(iField)(iField-1)*size(u,2)+scalar.apvRows, ...
            1:4,"UniformOutput",false);
        trustedColumns = horzcat(trustedByField{:});
        result = struct("rawCoordinateDefect",blockDefect(representedRaw,raw), ...
            "fieldDefects",fieldDefects,"maximumDefect",max(fieldDefects), ...
            "trustedFieldDefects",trustedFieldDefects, ...
            "trustedMaximumDefect",blockDefect( ...
                represented(:,trustedColumns),targets(:,trustedColumns)), ...
            "continuityDefect",norm(c.continuity*representedRaw,"fro") ...
                /max(norm(c.continuity,"fro")*norm(representedRaw,"fro"),realmin));
    end
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

function diagnostics = tangentModeDiagnostics(flat,first,G0,G1,direct0,direct1,trustedRows,trustedColumns)
scales = [0.1;0.05;0.025];
balancedResidual = zeros(size(scales));
waveAPV = zeros(size(scales));
bottomResidual = zeros(size(scales));
frequencyImaginaryDefect = zeros(size(scales));
eigenproblemDefect = zeros(size(scales));
orthogonalityDefect = zeros(size(scales));
distinctFrequencyOverlap = zeros(size(scales));
numberOfWaves = zeros(size(scales));
numberOfBoundaryActiveWaves = zeros(size(scales));
for iScale = 1:numel(scales)
    delta = scales(iScale);
    H = flat.E+delta*first.E;
    J = flat.J+delta*first.J;
    G = G0+delta*G1;
    balancedResidual(iScale) = normalizedSum(G'*J, ...
        {G0'*flat.J,delta*G1'*flat.J,delta*G0'*first.J,delta^2*G1'*first.J}, ...
        trustedRows,trustedColumns);

    [C,omega] = eig(1i*J,H,"vector");
    for iMode = 1:size(C,2)
        energyNorm = real(C(:,iMode)'*H*C(:,iMode));
        if ~isfinite(energyNorm) || energyNorm <= 0
            error("WVTerrainEnergyGalerkin:NonpositiveWeakEigenproblemMode", ...
                "A tangent terrain mode has nonpositive finite-terrain energy at scale %.3g.",delta)
        end
        C(:,iMode) = C(:,iMode)/sqrt(energyNorm);
    end
    residual = 1i*J*C-H*C.*omega.';
    eigenproblemDefect(iScale) = norm(residual,"fro") ...
        /max(norm(J,"fro")*norm(C,"fro")+norm(H,"fro")*norm(C.*omega.',"fro"),realmin);
    frequencyImaginaryDefect(iScale) = max(abs(imag(omega))) ...
        /max(max(abs(real(omega))),1);
    gram = C'*H*C;
    frequencySeparation = abs(omega-omega.');
    distinct = frequencySeparation > 1e-8*max(max(abs(omega)),1);
    distinctFrequencyOverlap(iScale) = norm(gram(distinct)) ...
        /max(sqrt(nnz(distinct)),1);
    orthogonalityResidual = (conj(omega)-omega.').*gram;
    orthogonalityDefect(iScale) = norm(orthogonalityResidual,"fro") ...
        /max(2*max(abs(omega))*norm(gram,"fro"),realmin);

    frequencyTolerance = 1e-8*max(max(abs(real(omega))),1);
    waves = abs(real(omega)) > frequencyTolerance;
    numberOfWaves(iScale) = nnz(waves);
    weakAPV = direct0+delta*direct1;
    waveAPV(iScale) = norm(weakAPV(trustedRows,:)*C(:,waves),"fro") ...
        /max(norm(weakAPV(trustedRows,:),"fro")*norm(C(:,waves),"fro"),realmin);
    modalBottom = -1i*(flat.B*C(:,waves)).*omega(waves).' ...
        -delta*first.R*C(:,waves);
    bottomResidual(iScale) = norm(modalBottom,"fro") ...
        /max(norm((flat.B*C(:,waves)).*omega(waves).',"fro") ...
            +norm(delta*first.R*C(:,waves),"fro"),realmin);
    bottomAmplitude = sqrt(sum(abs(flat.B*C(:,waves)).^2,1));
    numberOfBoundaryActiveWaves(iScale) = nnz(bottomAmplitude > 1e-8*max(max(bottomAmplitude),realmin));
end
diagnostics = struct("terrainScales",scales, ...
    "balancedResidual",balancedResidual, ...
    "balancedResidualReduction",scales(1:end-1).^2./scales(2:end).^2, ...
    "observedBalancedResidualReduction",balancedResidual(1:end-1)./max(balancedResidual(2:end),realmin), ...
    "trustedWaveAPVByScale",waveAPV, ...
    "trustedWaveAPVDefect",waveAPV(end), ...
    "bottomResidualByScale",bottomResidual, ...
    "frequencyImaginaryDefectByScale",frequencyImaginaryDefect, ...
    "eigenproblemDefectByScale",eigenproblemDefect, ...
    "eigenproblemDefect",max(eigenproblemDefect), ...
    "orthogonalityDefectByScale",orthogonalityDefect, ...
    "distinctFrequencyOverlapByScale",distinctFrequencyOverlap, ...
    "numberOfWaves",numberOfWaves, ...
    "numberOfBoundaryActiveWaves",numberOfBoundaryActiveWaves);
end

function summary = summarize(result,primitive)
summary = emptyRefinement;
summary.polynomialDegree = primitive.polynomialDegree;
summary.quadratureOrder = primitive.quadratureOrder;
summary.numberOfScalarTests = result.scalar.numberOfScalarTests;
summary.scalarBoundaryDefect = result.scalar.boundaryDefect;
summary.trustedStateRepresentationDefect = result.representation.trustedMaximumDefect;
summary.flatGreenIdentityDefect = result.flatGreenIdentityDefect;
summary.trustedTangentGreenIdentityDefect = result.trustedTangentGreenIdentityDefect;
summary.flatStationaryRowDefect = result.flatStationaryRowDefect;
summary.trustedTangentStationaryRowDefect = result.trustedTangentStationaryRowDefect;
summary.trustedWeakAPVDefect = result.trustedWeakAPVDefect;
summary.trustedDirectAPVDefect = result.trustedDirectAPVDefect;
summary.trustedWeakStrongAgreementDefect = result.trustedWeakStrongAgreementDefect;
summary.energyDefect = primitive.compatibility.energyDefect;
summary.bottomDefect = primitive.compatibility.bottomDefect;
summary.conjugacyDefect = primitive.compatibility.conjugacyDefect;
summary.modeDiagnostics = result.modeDiagnostics;
end

function summary = emptyRefinement
summary = struct("polynomialDegree",0,"quadratureOrder",0, ...
    "numberOfScalarTests",0,"scalarBoundaryDefect",NaN, ...
    "trustedStateRepresentationDefect",NaN,"flatGreenIdentityDefect",NaN, ...
    "trustedTangentGreenIdentityDefect",NaN,"flatStationaryRowDefect",NaN, ...
    "trustedTangentStationaryRowDefect",NaN,"trustedWeakAPVDefect",NaN, ...
    "trustedDirectAPVDefect",NaN,"trustedWeakStrongAgreementDefect",NaN, ...
    "energyDefect",NaN,"bottomDefect",NaN,"conjugacyDefect",NaN, ...
    "modeDiagnostics",struct);
end

function diagnostics = convergenceDiagnostics(refinement)
diagnostics = struct( ...
    "greenIdentityReduction",reduction([refinement.trustedTangentGreenIdentityDefect]), ...
    "stationaryRowReduction",reduction([refinement.trustedTangentStationaryRowDefect]), ...
    "weakAPVReduction",reduction([refinement.trustedWeakAPVDefect]), ...
    "directAPVReduction",reduction([refinement.trustedDirectAPVDefect]));
end

function values = reduction(defects)
values = defects(1:end-1)./max(defects(2:end),realmin);
end

function value = pairDefect(first,second,rows,columns)
value = norm(first(rows,columns)-second(rows,columns),"fro") ...
    /max(norm(first(rows,columns),"fro")+norm(second(rows,columns),"fro"),realmin);
end

function value = normalizedProduct(residual,left,right,rows,columns)
value = norm(residual(rows,columns),"fro") ...
    /max(norm(left(:,rows),"fro")*norm(right(:,columns),"fro"),realmin);
end

function value = normalizedSum(residual,terms,rows,columns)
scale = 0;
for iTerm = 1:numel(terms)
    scale = scale+norm(terms{iTerm}(rows,columns),"fro");
end
value = norm(residual(rows,columns),"fro")/max(scale,realmin);
end

function value = blockDefect(first,second)
value = norm(first-second,"fro")/max(norm(first,"fro")+norm(second,"fro"),realmin);
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
    error("WVTerrainEnergyGalerkin:UnavailableWeakEigenproblemSupport", ...
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
    error("WVTerrainEnergyGalerkin:WeakEigenproblemTerrainNyquist", ...
        "The boundary-complete weak oracle requires terrain with no material Nyquist coefficient.")
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
            error("WVTerrainEnergyGalerkin:InsufficientWeakEigenproblemGuard", ...
                "supportModeBounds must retain every first terrain sideband of the trusted band.")
        end
    end
end
end
