function result = auditPrimitiveSequenceDirection(c,direction, ...
    trustedBounds,scalarDegree,terrainScale)
% Audit one complete primitive direction using its own vertical coordinate.
scalar = scalarTestFields(c,trustedBounds,scalarDegree);
[G,representation] = projectGeostrophicState(c,scalar,terrainScale);
H3 = repmat(c.H,c.nZ,1);
gamma = 1-terrainScale*H3;
mass = c.wvt.rho0*(scalar.values'*(c.volumeWeight.*gamma.*scalar.values));
weakRow = G'*direction.E;
strongRow = -c.wvt.rho0 ...
    *(scalar.values'*(c.volumeWeight.*gamma.*direction.Q));
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
greenIdentityDefect = pairDefect(weakRow,strongRow, ...
    trustedRows,trustedColumns);
stationaryRowDefect = normalizedProduct(G'*direction.J,G,direction.J, ...
    trustedRows,trustedColumns);
weakEvolutionResidual = direction.E*L-direction.J;
stationaryEvolution = (G'*direction.J+G'*weakEvolutionResidual) ...
    /energyFactor;
greenEvolutionCorrection = ((strongRow-weakRow)/energyFactor) ...
    *energyGenerator;
energyColumns = 1:size(energyGenerator,2);
weakAPVDefect = normalizedResidual(stationaryEvolution, ...
    trustedRows,energyColumns, ...
    norm(weakEnergyRow(trustedRows,:),"fro")*norm(energyGenerator,"fro"));
strongAPVEvolution = stationaryEvolution+greenEvolutionCorrection;
strongAPVDefect = normalizedResidual(strongAPVEvolution, ...
    trustedRows,energyColumns, ...
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

result = struct;
result.terrainScale = terrainScale;
result.scalar = scalar.diagnostics;
result.G = G;
result.generator = L;
result.energyMatrix = direction.E;
result.exchangeMatrix = direction.J;
result.stateRepresentationDefect = representation.trustedMaximumDefect;
result.greenIdentityDefect = greenIdentityDefect;
result.stationaryRowDefect = stationaryRowDefect;
result.weakAPVDefect = weakAPVDefect;
result.strongAPVDefect = strongAPVDefect;
result.weakStrongAgreementDefect = weakStrongAgreementDefect;
result.weakEvolutionDefect = productDefect(weakEvolutionResidual, ...
    {direction.E*L,direction.J});
result.primitiveGeneratorAgreementDefect = norm( ...
    primitiveAgreement(:,trustedColumns),"fro") ...
    /max(norm(energyGenerator(:,trustedColumns),"fro") ...
    +norm(strongPrimitiveEnergyGenerator(:,trustedColumns),"fro"),realmin);
result.energyDefect = productDefect(energyResidual, ...
    {energyGenerator',energyGenerator});
result.bottomDefect = bottomProductDefect(direction.B,L,direction.R, ...
    bottomResidual,trustedColumns);
result.bottomResidualNorm = norm(bottomResidual(:,trustedColumns),"fro");
result.conjugacyDefect = norm(conjugacyResidual,"fro") ...
    /max(2*norm(L,"fro"),realmin);
result.enstrophyDefect = norm(enstrophyResidual,"fro") ...
    /max(2*norm(Z,"fro")*norm(energyGenerator,"fro"),realmin);
result.energyHermitianDefect = norm(direction.E-direction.E',"fro") ...
    /max(norm(direction.E,"fro"),realmin);
result.exchangeSkewHermitianDefect = norm(direction.J+direction.J',"fro") ...
    /max(norm(direction.J,"fro"),realmin);
result.energyRcond = rcond((direction.E+direction.E')/2);
result.geostrophicBottomTangencyDefect = norm( ...
    direction.R*G(:,trustedRows),"fro") ...
    /max(norm(direction.R,"fro")*norm(G(:,trustedRows),"fro"),realmin);
result.continuityTangencyDefect = direction.diagnostics.continuityTangencyDefect;
result.projectionAdjointDefect = c.projectionDiagnostics.adjointDefect;
result.exactConvolutionDefect = ...
    c.projectionDiagnostics.maximumExactConvolutionDefect;
result.verticalDiagnostics = c.verticalDiagnostics;
end

function scalar = scalarTestFields(c,trustedBounds,scalarDegree)
if scalarDegree >= c.nF
    error("WVTerrainEnergyGalerkin:WKBScalarDegreeOutsidePrimitiveSpace", ...
        "The scalar test degree must be smaller than the horizontal-velocity basis dimension.")
end
bottomValue = c.spaces.Fendpoint(1,1:scalarDegree+1);
topDerivative = c.spaces.FendpointXi(2,1:scalarDegree+1);
constraint = [bottomValue;topDerivative];
[~,singularValues,V] = svd(constraint);
singularValues = diag(singularValues);
tolerance = max(size(constraint))*eps(max(singularValues));
rankConstraint = nnz(singularValues > tolerance);
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
diagnostics = struct("numberOfVerticalFunctions",nS, ...
    "numberOfScalarTests",nScalar,"numberOfTrustedTests",numel(trustedRows), ...
    "trustedScalarDegree",scalarDegree,"constraintRank",rankConstraint, ...
    "bottomValueDefect",norm(bottomValue*trustedCoefficients,"fro") ...
    /max(norm(coefficients,"fro"),realmin), ...
    "surfaceDerivativeDefect",norm(topDerivative*trustedCoefficients,"fro") ...
    /max(norm(coefficients,"fro"),realmin));
scalar = struct("values",values,"valuesXi",valuesXi, ...
    "valuesX",valuesX,"valuesY",valuesY, ...
    "trustedRows",trustedRows,"diagnostics",diagnostics);
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
trustedByField = arrayfun(@(iField)(iField-1)*size(u,2) ...
    +scalar.trustedRows,1:4,"UniformOutput",false);
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
    etaColumns((iK-1)*c.nH+(1:c.nH)) = ...
        block+2*c.nF+c.nG+(1:c.nH);
end
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

function value = pairDefect(first,second,rows,columns)
value = norm(first(rows,columns)-second(rows,columns),"fro") ...
    /max(norm(first(rows,columns),"fro") ...
    +norm(second(rows,columns),"fro"),realmin);
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
naturalScale = norm(B,"fro")*norm(L(:,columns),"fro") ...
    +norm(R(:,columns),"fro");
value = norm(residual(:,columns),"fro") ...
    /max(naturalScale,realmin);
end

function value = blockDefect(first,second)
value = norm(first-second,"fro") ...
    /max(norm(first,"fro")+norm(second,"fro"),realmin);
end
