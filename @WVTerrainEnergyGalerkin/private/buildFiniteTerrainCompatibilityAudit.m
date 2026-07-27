function audit = buildFiniteTerrainCompatibilityAudit(problem)
% Audit the unmodified finite-terrain Galerkin generator.

forms = problem.finiteTerrainForms;
E = forms.energyMatrix;
J = forms.exchangeMatrix;
Q = forms.apvMatrix;
QUnweighted = forms.unweightedAPVMatrix;
Z = Q'*Q;
L = E\J;
nState = size(L,1);

[B,R,bottomDiagnostics] = buildBottomConstraintMaps(problem);
coordinateConjugacy = sparse((1:nState)',problem.conjugateCoordinateIndex,1,nState,nState);

energyResidual = L'*E+E*L;
apvResidual = Q*L;
unweightedAPVResidual = QUnweighted*L;
bottomResidual = B*L-R;
enstrophyResidual = L'*Z+Z*L;
conjugacyResidual = L*coordinateConjugacy-coordinateConjugacy*conj(L);

diagnostics = struct;
diagnostics.energyDefect = relativeProductDefect(energyResidual,norm(E,"fro")*norm(L,"fro"));
diagnostics.apvDefect = relativeProductDefect(apvResidual,norm(Q,"fro")*norm(L,"fro"));
diagnostics.unweightedAPVDefect = relativeProductDefect(unweightedAPVResidual,norm(QUnweighted,"fro")*norm(L,"fro"));
diagnostics.bottomDefect = relativeProductDefect(bottomResidual,norm(B,"fro")*norm(L,"fro")+norm(R,"fro"));
diagnostics.enstrophyDefect = relativeProductDefect(enstrophyResidual,norm(Z,"fro")*norm(L,"fro"));
diagnostics.conjugacyDefect = relativeProductDefect(conjugacyResidual,norm(L,"fro"));

rawE = forms.rawEnergyMatrix;
rawJ = forms.rawExchangeMatrix;
preRestorationGenerator = rawE\rawJ;
rawEnergyResidual = preRestorationGenerator'*rawE+rawE*preRestorationGenerator;
rawAPVResidual = Q*preRestorationGenerator;
rawBottomResidual = B*preRestorationGenerator-R;
rawEnstrophyResidual = preRestorationGenerator'*Z+Z*preRestorationGenerator;
rawConjugacyResidual = preRestorationGenerator*coordinateConjugacy-coordinateConjugacy*conj(preRestorationGenerator);
diagnostics.preRestorationEnergyDefect = relativeProductDefect(rawEnergyResidual,norm(rawE,"fro")*norm(preRestorationGenerator,"fro"));
diagnostics.preRestorationAPVDefect = relativeProductDefect(rawAPVResidual,norm(Q,"fro")*norm(preRestorationGenerator,"fro"));
diagnostics.preRestorationBottomDefect = relativeProductDefect(rawBottomResidual,norm(B,"fro")*norm(preRestorationGenerator,"fro")+norm(R,"fro"));
diagnostics.preRestorationEnstrophyDefect = relativeProductDefect(rawEnstrophyResidual,norm(Z,"fro")*norm(preRestorationGenerator,"fro"));
diagnostics.preRestorationConjugacyDefect = relativeProductDefect(rawConjugacyResidual,norm(preRestorationGenerator,"fro"));
diagnostics.rawRestoredGeneratorAgreement = norm(preRestorationGenerator-L,"fro")/max(norm(L,"fro"),realmin);

S = chol(E);
energyGenerator = S*L/S;
energyAPVMap = Q/S;
energyBottomValueMap = B/S;
energyBottomTendencyMap = R/S;
energyEnstrophyMatrix = energyAPVMap'*energyAPVMap;
apvSingularValues = svd(energyAPVMap);
if isempty(apvSingularValues)
    apvRankTolerance = 0;
else
    apvRankTolerance = max(size(energyAPVMap))*eps(max(apvSingularValues));
end
apvRank = nnz(apvSingularValues > apvRankTolerance);
diagnostics.apvSingularValues = apvSingularValues;
diagnostics.apvRankTolerance = apvRankTolerance;
diagnostics.apvRank = apvRank;
diagnostics.apvNullity = nState-apvRank;
diagnostics.generatorRank = rank(energyGenerator);
diagnostics.generatorRangeFitsAPVNullspace = diagnostics.generatorRank <= diagnostics.apvNullity;

diagnostics.randomStates = randomStateDiagnostics(E,L,Q,Z,B,R,coordinateConjugacy);
diagnostics.subspace = subspaceDiagnostics(problem,S,energyGenerator,energyAPVMap,energyBottomValueMap,energyBottomTendencyMap,energyEnstrophyMatrix);
diagnostics.apvResidualByVerticalLevel = residualByVerticalLevel(apvResidual,forms.oversampledSize);
diagnostics.unweightedAPVResidualByVerticalLevel = residualByVerticalLevel(unweightedAPVResidual,forms.oversampledSize);
diagnostics.bottomResidualByHorizontalMode = vecnorm(bottomResidual,2,2);
diagnostics.bottom = bottomDiagnostics;
diagnostics.rawHermitianDefect = problem.constructionDiagnostics.finiteTerrain.hermitianDefect;
diagnostics.rawSkewHermitianDefect = problem.constructionDiagnostics.finiteTerrain.skewHermitianDefect;
diagnostics.scaledEnergyRcond = problem.constructionDiagnostics.finiteTerrain.scaledEnergyRcond;
diagnostics.quadratureEnergyRelativeResidual = problem.constructionDiagnostics.finiteTerrain.quadratureEnergyRelativeResidual;
diagnostics.quadratureExchangeRelativeResidual = problem.constructionDiagnostics.finiteTerrain.quadratureExchangeRelativeResidual;

failedConditions = strings(0,1);
if diagnostics.rawHermitianDefect > 1e-12 || diagnostics.rawSkewHermitianDefect > 1e-12 || diagnostics.scaledEnergyRcond <= 1e-12
    failedConditions(end+1,1) = "finite-form-structure";
end
if diagnostics.energyDefect > 1e-12
    failedConditions(end+1,1) = "energy-compatibility";
end
if diagnostics.apvDefect > 1e-10
    failedConditions(end+1,1) = "apv-compatibility";
end
if diagnostics.bottomDefect > 1e-10
    failedConditions(end+1,1) = "bottom-evolution-compatibility";
end
if diagnostics.enstrophyDefect > 1e-10
    failedConditions(end+1,1) = "quadratic-enstrophy-compatibility";
end
if diagnostics.conjugacyDefect > 1e-12
    failedConditions(end+1,1) = "conjugacy";
end
isCompatible = isempty(failedConditions);
status = "compatible";
if ~isCompatible
    status = "incompatible";
end

audit = struct;
audit.status = status;
audit.isCompatible = isCompatible;
audit.failedConditions = failedConditions;
audit.rawGenerator = L;
audit.preRestorationGenerator = preRestorationGenerator;
audit.quadraticEnstrophyMatrix = Z;
audit.bottomValueMatrix = B;
audit.bottomTendencyMatrix = R;
audit.energyResidual = energyResidual;
audit.apvResidual = apvResidual;
audit.unweightedAPVResidual = unweightedAPVResidual;
audit.bottomResidual = bottomResidual;
audit.enstrophyResidual = enstrophyResidual;
audit.conjugacyResidual = conjugacyResidual;
audit.preRestorationEnergyResidual = rawEnergyResidual;
audit.preRestorationAPVResidual = rawAPVResidual;
audit.preRestorationBottomResidual = rawBottomResidual;
audit.preRestorationEnstrophyResidual = rawEnstrophyResidual;
audit.preRestorationConjugacyResidual = rawConjugacyResidual;
audit.diagnostics = diagnostics;
end

function value = relativeProductDefect(residual,scale)
value = norm(residual,"fro")/max(scale,realmin);
end

function diagnostics = randomStateDiagnostics(E,L,Q,Z,B,R,C)
n = size(L,1);
stream = RandStream("mt19937ar","Seed",481903);
complexState = randn(stream,n,1)+1i*randn(stream,n,1);
candidate = randn(stream,n,1)+1i*randn(stream,n,1);
physicalState = (candidate+C*conj(candidate))/2;
diagnostics = struct;
diagnostics.complex = oneStateDiagnostics(complexState,E,L,Q,Z,B,R);
diagnostics.physical = oneStateDiagnostics(physicalState,E,L,Q,Z,B,R);
diagnostics.physicalConjugacyDefect = norm(physicalState-C*conj(physicalState))/max(norm(physicalState),realmin);
end

function diagnostics = oneStateDiagnostics(a,E,L,Q,Z,B,R)
stateNorm = norm(a);
stateNormSquared = stateNorm^2;
tendency = L*a;
diagnostics = struct;
diagnostics.energyTendency = abs(real(a'*E*tendency))/max(norm(E,"fro")*norm(L,"fro")*stateNormSquared,realmin);
diagnostics.apvTendency = norm(Q*tendency)/max(norm(Q,"fro")*norm(L,"fro")*stateNorm,realmin);
diagnostics.bottomTendency = norm((B*L-R)*a)/max((norm(B,"fro")*norm(L,"fro")+norm(R,"fro"))*stateNorm,realmin);
diagnostics.enstrophyTendency = abs(real(a'*Z*tendency))/max(norm(Z,"fro")*norm(L,"fro")*stateNormSquared,realmin);
end

function diagnostics = subspaceDiagnostics(problem,S,K,Q,B,R,G)
nState = size(K,1);
names = ["waves","apv-bearing-geostrophic","bottom-inversion","mean-density"];
coordinates = repmat({zeros(nState,0)},size(names));
for iK = 1:numel(problem.flatModeBlocks)
    rows = problem.stateLayout.horizontalIndex == iK;
    block = problem.flatModeBlocks{iK};
    kappa = hypot(problem.horizontalLayout.k(iK),problem.horizontalLayout.l(iK));
    indices = {block.waveIndices,[],[],[]};
    if kappa > 0
        indices{2} = block.apvBearingStationaryIndices;
        indices{3} = block.zeroAPVStationaryIndices;
    else
        indices{4} = block.stationaryIndices;
    end
    for iSubspace = 1:numel(names)
        if isempty(indices{iSubspace})
            continue
        end
        vectors = zeros(nState,numel(indices{iSubspace}));
        vectors(rows,:) = block.eigenvectors(:,indices{iSubspace});
        coordinates{iSubspace} = [coordinates{iSubspace} vectors];
    end
end

apvResidual = Q*K;
bottomResidual = B*K-R;
enstrophyResidual = K'*G+G*K;
energyResidual = K'+K;
diagnostics = repmat(struct("name","","dimension",0,"energyDefect",0,"apvDefect",0,"bottomDefect",0,"enstrophyDefect",0),numel(names),1);
for iSubspace = 1:numel(names)
    [U,~] = qr(S*coordinates{iSubspace},0);
    diagnostics(iSubspace).name = names(iSubspace);
    diagnostics(iSubspace).dimension = size(U,2);
    diagnostics(iSubspace).energyDefect = norm(energyResidual*U,"fro")/max(norm(K,"fro"),realmin);
    diagnostics(iSubspace).apvDefect = norm(apvResidual*U,"fro")/max(norm(Q,"fro")*norm(K,"fro"),realmin);
    diagnostics(iSubspace).bottomDefect = norm(bottomResidual*U,"fro")/max(norm(B,"fro")*norm(K,"fro")+norm(R,"fro"),realmin);
    diagnostics(iSubspace).enstrophyDefect = norm(enstrophyResidual*U,"fro")/max(norm(G,"fro")*norm(K,"fro"),realmin);
end
end

function values = residualByVerticalLevel(residual,oversampledSize)
nXY = prod(oversampledSize(1:2));
Nz = oversampledSize(3);
values = zeros(Nz,1);
scale = max(norm(residual,"fro"),realmin);
for iz = 1:Nz
    rows = (iz-1)*nXY+(1:nXY);
    values(iz) = norm(residual(rows,:),"fro")/scale;
end
end
