function audit = buildQuadraticInvariantClosureAudit(problem)
% Audit an energy-, quadratic-enstrophy-, and bottom-compatible closure.

forms = problem.finiteTerrainForms;
E = forms.energyMatrix;
J = forms.exchangeMatrix;
Q = forms.apvMatrix;
nState = size(E,1);

[stateRealification,stateRealificationDiagnostics] = buildRealConjugacyCoordinates(problem.conjugateCoordinateIndex);
[bottomValueMatrix,bottomTendencyMatrix,bottomDiagnostics] = buildBottomConstraintMaps(problem);
bottomConjugateIndex = horizontalConjugateIndex(problem);
[bottomRealification,bottomRealificationDiagnostics] = buildRealConjugacyCoordinates(bottomConjugateIndex);

complexRealEnergy = stateRealification'*E*stateRealification;
complexRealExchange = stateRealification'*J*stateRealification;
energyImaginaryDefect = imaginaryDefect(complexRealEnergy);
exchangeImaginaryDefect = imaginaryDefect(complexRealExchange);
realEnergyMatrix = real(complexRealEnergy);
realExchangeMatrix = real(complexRealExchange);
realEnergyMatrix = (realEnergyMatrix+realEnergyMatrix')/2;
realExchangeMatrix = (realExchangeMatrix-realExchangeMatrix')/2;
realEnergyFactor = chol(realEnergyMatrix);
energyFactor = realEnergyFactor*stateRealification';

rawGenerator = (realEnergyFactor'\realExchangeMatrix)/realEnergyFactor;
rawGeneratorSkewDefect = relativeSkewDefect(rawGenerator);
rawGenerator = (rawGenerator-rawGenerator')/2;

complexAPVMap = Q*stateRealification/realEnergyFactor;
apvImaginaryDefect = imaginaryDefect(complexAPVMap);
apvMap = real(complexAPVMap);
quadraticEnstrophyMatrix = Q'*Q;
enstrophyMetric = apvMap'*apvMap;
enstrophyMetric = (enstrophyMetric+enstrophyMetric')/2;
quadraticFormDefect = norm(quadraticEnstrophyMatrix-energyFactor'*enstrophyMetric*energyFactor,"fro")/max(norm(quadraticEnstrophyMatrix,"fro"),realmin);

complexBottomValueMap = bottomRealification'*bottomValueMatrix*stateRealification/realEnergyFactor;
complexBottomTendencyMap = bottomRealification'*bottomTendencyMatrix*stateRealification/realEnergyFactor;
bottomValueImaginaryDefect = imaginaryDefect(complexBottomValueMap);
bottomTendencyImaginaryDefect = imaginaryDefect(complexBottomTendencyMap);
bottomValueMap = real(complexBottomValueMap);
bottomTendencyMap = real(complexBottomTendencyMap);

realificationTolerance = 1e-12;
maximumRealificationDefect = max([ ...
    stateRealificationDiagnostics.unitaryDefect, ...
    stateRealificationDiagnostics.conjugacyDefect, ...
    bottomRealificationDiagnostics.unitaryDefect, ...
    bottomRealificationDiagnostics.conjugacyDefect, ...
    energyImaginaryDefect,exchangeImaginaryDefect,apvImaginaryDefect, ...
    bottomValueImaginaryDefect,bottomTendencyImaginaryDefect,quadraticFormDefect]);
if maximumRealificationDefect > realificationTolerance
    error("WVTerrainEnergyGalerkin:QuadraticClosureRealificationFailure", ...
        "The independent-real coordinate forms have a maximum structural defect of %.3g.",maximumRealificationDefect)
end

[eigenvectors,eigenvalues] = rightSingularEigensystem(apvMap);
eigenResiduals = vecnorm(enstrophyMetric*eigenvectors-eigenvectors.*eigenvalues.',2,1).';
backwardFloor = 10*nState*eps(max(norm(enstrophyMetric,2),realmin));
backwardRadii = eigenResiduals+backwardFloor;
groups = eigenvalueGroups(eigenvalues,backwardRadii);

nGroups = numel(groups);
blockDiagnostics = repmat(emptyBlockDiagnostics,1,nGroups);
minimumResidualGenerator = zeros(nState);
constrainedGenerator = zeros(nState);
allBlocksFeasible = true;
constraintTolerance = 1e-12*max([norm(bottomTendencyMap,"fro"),norm(bottomValueMap,"fro")*norm(rawGenerator,"fro"),realmin]);
minimumConstraintResidualSquared = 0;
for iGroup = 1:nGroups
    indices = groups(iGroup).firstIndex:groups(iGroup).lastIndex;
    U = eigenvectors(:,indices);
    B = bottomValueMap*U;
    R = bottomTendencyMap*U;
    rawBlock = U'*rawGenerator*U;
    rawBlock = (rawBlock-rawBlock')/2;
    [constraintMatrix,rawCoordinates] = realSkewConstraintSystem(B,rawBlock);
    reconstructedRawBlock = realSkewMatrixFromCoordinates(rawCoordinates,numel(indices));
    rawCoordinateRoundTripDefect = norm(reconstructedRawBlock-rawBlock,"fro")/max(norm(rawBlock,"fro"),realmin);
    rawBlockBottomResidual = norm(B*rawBlock-R,"fro");
    target = R(:);
    [leftVectors,singularValues,rightVectors,numericalRank,rankTolerance] = singularFactors(constraintMatrix);
    minimumResidualCoordinates = applyPseudoinverse(leftVectors,singularValues,rightVectors,numericalRank,target);
    minimumResidualBlock = realSkewMatrixFromCoordinates(minimumResidualCoordinates,numel(indices));
    if numericalRank == 0
        projectedTarget = zeros(size(target));
    else
        projectedTarget = leftVectors(:,1:numericalRank)*(leftVectors(:,1:numericalRank)'*target);
    end
    leastSquaresResidual = projectedTarget-target;
    minimumResidual = norm(leastSquaresResidual);
    minimumRelativeResidual = minimumResidual/max(norm(R,"fro"),constraintTolerance);
    blockIsFeasible = minimumResidual <= constraintTolerance;
    allBlocksFeasible = allBlocksFeasible && blockIsFeasible;
    minimumConstraintResidualSquared = minimumConstraintResidualSquared+minimumResidual^2;

    if numericalRank == 0
        normalEquationRelativeResidual = 0;
    else
        normalEquationRelativeResidual = norm(leftVectors(:,1:numericalRank)'*leastSquaresResidual)/max(norm(target),realmin);
    end
    pseudoinverseReconstructionResidual = norm(constraintMatrix*minimumResidualCoordinates-projectedTarget)/max(norm(projectedTarget),realmin);
    minimumResidualGenerator = minimumResidualGenerator+U*minimumResidualBlock*U';

    correctionCoordinates = [];
    constrainedCoordinates = [];
    correctionResidual = NaN;
    if blockIsFeasible
        correctionTarget = target-constraintMatrix*rawCoordinates;
        correctionCoordinates = applyPseudoinverse(leftVectors,singularValues,rightVectors,numericalRank,correctionTarget);
        constrainedCoordinates = rawCoordinates+correctionCoordinates;
        constrainedBlock = realSkewMatrixFromCoordinates(constrainedCoordinates,numel(indices));
        correctionResidual = norm(constraintMatrix*constrainedCoordinates-target);
        constrainedGenerator = constrainedGenerator+U*constrainedBlock*U';
    end

    diagnostics = emptyBlockDiagnostics;
    diagnostics.firstIndex = groups(iGroup).firstIndex;
    diagnostics.lastIndex = groups(iGroup).lastIndex;
    diagnostics.multiplicity = numel(indices);
    diagnostics.minimumEigenvalue = eigenvalues(indices(1));
    diagnostics.maximumEigenvalue = eigenvalues(indices(end));
    diagnostics.maximumBackwardRadius = max(backwardRadii(indices));
    diagnostics.separationBelow = groups(iGroup).separationBelow;
    diagnostics.separationAbove = groups(iGroup).separationAbove;
    diagnostics.numberOfSkewCoordinates = size(constraintMatrix,2);
    diagnostics.constraintRank = numericalRank;
    diagnostics.rankTolerance = rankTolerance;
    diagnostics.singularValues = singularValues;
    diagnostics.minimumResidual = minimumResidual;
    diagnostics.minimumRelativeResidual = minimumRelativeResidual;
    diagnostics.rawCoordinateRoundTripDefect = rawCoordinateRoundTripDefect;
    diagnostics.rawBlockBottomResidual = rawBlockBottomResidual;
    diagnostics.normalEquationRelativeResidual = normalEquationRelativeResidual;
    diagnostics.pseudoinverseReconstructionResidual = pseudoinverseReconstructionResidual;
    diagnostics.isFeasible = blockIsFeasible;
    diagnostics.rawCoordinates = rawCoordinates;
    diagnostics.minimumResidualCoordinates = minimumResidualCoordinates;
    diagnostics.correctionCoordinates = correctionCoordinates;
    diagnostics.constrainedCoordinates = constrainedCoordinates;
    diagnostics.correctionResidual = correctionResidual;
    blockDiagnostics(iGroup) = diagnostics;
end
minimumResidualGenerator = (minimumResidualGenerator-minimumResidualGenerator')/2;
minimumConstraintResidual = sqrt(minimumConstraintResidualSquared);
minimumConstraintRelativeResidual = minimumConstraintResidual/max(norm(bottomTendencyMap,"fro"),constraintTolerance);
constraintSystemIsFeasible = allBlocksFeasible && minimumConstraintResidual <= constraintTolerance;
minimumResidualRealizationResidual = norm(bottomValueMap*minimumResidualGenerator-bottomTendencyMap,"fro");
rawEnstrophyDefect = relativeCommutatorDefect(enstrophyMetric,rawGenerator);
rawBottomDefect = relativeBottomDefect(bottomValueMap,rawGenerator,bottomTendencyMap);
rawSatisfiesConstraints = rawGeneratorSkewDefect <= 1e-13 && rawEnstrophyDefect <= 1e-12 && rawBottomDefect <= 1e-12;

constrainedExchangeMatrix = [];
correctionNorm = NaN;
energyConstraintResidual = NaN;
enstrophyConstraintResidual = NaN;
bottomConstraintResidual = NaN;
constrainedConjugacyDefect = NaN;
maximumFrequencyImaginaryPart = NaN;
apvRedistributionNorm = NaN;
if constraintSystemIsFeasible
    if rawSatisfiesConstraints
        constrainedGenerator = rawGenerator;
    else
        constrainedGenerator = (constrainedGenerator-constrainedGenerator')/2;
    end
    constrainedExchangeMatrix = energyFactor'*constrainedGenerator*energyFactor;
    correctionNorm = norm(constrainedGenerator-rawGenerator,"fro")/max(norm(rawGenerator,"fro"),realmin);
    energyConstraintResidual = relativeSkewDefect(constrainedGenerator);
    enstrophyConstraintResidual = relativeCommutatorDefect(enstrophyMetric,constrainedGenerator);
    bottomConstraintResidual = relativeBottomDefect(bottomValueMap,constrainedGenerator,bottomTendencyMap);
    fullGenerator = energyFactor\(constrainedGenerator*energyFactor);
    coordinateConjugacy = sparse((1:nState)',problem.conjugateCoordinateIndex,1,nState,nState);
    constrainedConjugacyDefect = norm(fullGenerator*coordinateConjugacy-coordinateConjugacy*conj(fullGenerator),"fro")/max(norm(fullGenerator,"fro"),realmin);
    frequencies = eig(1i*constrainedGenerator);
    maximumFrequencyImaginaryPart = max(abs(imag(frequencies)))/max(1,max(abs(real(frequencies))));
    apvRedistributionNorm = norm(apvMap*constrainedGenerator,"fro")/max(norm(apvMap,"fro")*norm(constrainedGenerator,"fro"),realmin);
else
    constrainedGenerator = [];
end

failedConditions = strings(0,1);
if ~constraintSystemIsFeasible
    failedConditions(end+1,1) = "incompatible-bottom-constraint";
end
if constraintSystemIsFeasible && energyConstraintResidual > 1e-13
    failedConditions(end+1,1) = "energy-constraint-residual";
end
if constraintSystemIsFeasible && enstrophyConstraintResidual > 1e-12
    failedConditions(end+1,1) = "enstrophy-commutator-residual";
end
if constraintSystemIsFeasible && bottomConstraintResidual > 1e-12
    failedConditions(end+1,1) = "bottom-constraint-residual";
end
if constraintSystemIsFeasible && constrainedConjugacyDefect > 1e-12
    failedConditions(end+1,1) = "constrained-conjugacy-defect";
end
isScientificallyAdmissible = isempty(failedConditions);
status = "admissible";
if ~isScientificallyAdmissible
    status = "incompatible";
    if constraintSystemIsFeasible
        constrainedGenerator = [];
        constrainedExchangeMatrix = [];
    end
end

diagnostics = struct;
diagnostics.stateRealification = stateRealificationDiagnostics;
diagnostics.bottomRealification = bottomRealificationDiagnostics;
diagnostics.energyImaginaryDefect = energyImaginaryDefect;
diagnostics.exchangeImaginaryDefect = exchangeImaginaryDefect;
diagnostics.apvImaginaryDefect = apvImaginaryDefect;
diagnostics.bottomValueImaginaryDefect = bottomValueImaginaryDefect;
diagnostics.bottomTendencyImaginaryDefect = bottomTendencyImaginaryDefect;
diagnostics.quadraticFormDefect = quadraticFormDefect;
diagnostics.maximumRealificationDefect = maximumRealificationDefect;
diagnostics.rawGeneratorSkewDefect = rawGeneratorSkewDefect;
diagnostics.rawEnstrophyDefect = rawEnstrophyDefect;
diagnostics.rawBottomDefect = rawBottomDefect;
diagnostics.rawSatisfiesConstraints = rawSatisfiesConstraints;
diagnostics.eigenvalues = eigenvalues;
diagnostics.eigenResiduals = eigenResiduals;
diagnostics.backwardRadii = backwardRadii;
diagnostics.eigenvalueGroups = groups;
diagnostics.block = blockDiagnostics;
diagnostics.minimumConstraintResidual = minimumConstraintResidual;
diagnostics.minimumConstraintRelativeResidual = minimumConstraintRelativeResidual;
diagnostics.minimumResidualRealizationResidual = minimumResidualRealizationResidual;
diagnostics.constraintTolerance = constraintTolerance;
diagnostics.correctionNorm = correctionNorm;
diagnostics.energyConstraintResidual = energyConstraintResidual;
diagnostics.enstrophyConstraintResidual = enstrophyConstraintResidual;
diagnostics.bottomConstraintResidual = bottomConstraintResidual;
diagnostics.constrainedConjugacyDefect = constrainedConjugacyDefect;
diagnostics.maximumFrequencyImaginaryPart = maximumFrequencyImaginaryPart;
diagnostics.apvRedistributionNorm = apvRedistributionNorm;
diagnostics.bottom = bottomDiagnostics;

audit = struct;
audit.status = status;
audit.constraintSystemIsFeasible = constraintSystemIsFeasible;
audit.isScientificallyAdmissible = isScientificallyAdmissible;
audit.failedConditions = failedConditions;
audit.stateRealification = stateRealification;
audit.bottomRealification = bottomRealification;
audit.realEnergyMatrix = realEnergyMatrix;
audit.realExchangeMatrix = realExchangeMatrix;
audit.realEnergyFactor = realEnergyFactor;
audit.energyFactor = energyFactor;
audit.rawGenerator = rawGenerator;
audit.minimumResidualGenerator = minimumResidualGenerator;
audit.constrainedGenerator = constrainedGenerator;
audit.constrainedExchangeMatrix = constrainedExchangeMatrix;
audit.quadraticEnstrophyMatrix = quadraticEnstrophyMatrix;
audit.enstrophyMetric = enstrophyMetric;
audit.apvMap = apvMap;
audit.bottomValueMatrix = bottomValueMatrix;
audit.bottomTendencyMatrix = bottomTendencyMatrix;
audit.bottomValueMap = bottomValueMap;
audit.bottomTendencyMap = bottomTendencyMap;
audit.diagnostics = diagnostics;
end

function conjugateIndex = horizontalConjugateIndex(problem)
nHorizontal = height(problem.horizontalLayout);
conjugateIndex = zeros(nHorizontal,1);
for iK = 1:nHorizontal
    rows = find(problem.stateLayout.horizontalIndex == iK);
    bottomRow = rows(problem.stateLayout.component(rows) == "etaB");
    partnerRow = problem.conjugateCoordinateIndex(bottomRow);
    conjugateIndex(iK) = problem.stateLayout.horizontalIndex(partnerRow);
end
end

function groups = eigenvalueGroups(eigenvalues,backwardRadii)
n = numel(eigenvalues);
groups = repmat(struct("firstIndex",0,"lastIndex",0,"separationBelow",Inf,"separationAbove",Inf),n,1);
first = 1;
numberOfGroups = 0;
while first <= n
    last = first;
    upper = eigenvalues(first)+backwardRadii(first);
    while last < n && eigenvalues(last+1)-backwardRadii(last+1) <= upper
        last = last+1;
        upper = max(upper,eigenvalues(last)+backwardRadii(last));
    end
    group = struct;
    group.firstIndex = first;
    group.lastIndex = last;
    group.separationBelow = Inf;
    group.separationAbove = Inf;
    if first > 1
        group.separationBelow = eigenvalues(first)-eigenvalues(first-1);
    end
    if last < n
        group.separationAbove = eigenvalues(last+1)-eigenvalues(last);
    end
    numberOfGroups = numberOfGroups+1;
    groups(numberOfGroups) = group;
    first = last+1;
end
groups = groups(1:numberOfGroups);
end

function diagnostics = emptyBlockDiagnostics
diagnostics = struct( ...
    "firstIndex",0, ...
    "lastIndex",0, ...
    "multiplicity",0, ...
    "minimumEigenvalue",NaN, ...
    "maximumEigenvalue",NaN, ...
    "maximumBackwardRadius",NaN, ...
    "separationBelow",NaN, ...
    "separationAbove",NaN, ...
    "numberOfSkewCoordinates",0, ...
    "constraintRank",0, ...
    "rankTolerance",NaN, ...
    "singularValues",[], ...
    "minimumResidual",NaN, ...
    "minimumRelativeResidual",NaN, ...
    "rawCoordinateRoundTripDefect",NaN, ...
    "rawBlockBottomResidual",NaN, ...
    "normalEquationRelativeResidual",NaN, ...
    "pseudoinverseReconstructionResidual",NaN, ...
    "isFeasible",false, ...
    "rawCoordinates",[], ...
    "minimumResidualCoordinates",[], ...
    "correctionCoordinates",[], ...
    "constrainedCoordinates",[], ...
    "correctionResidual",NaN);
end

function [A,x0] = realSkewConstraintSystem(B,H0)
r = size(H0,1);
nCoordinates = r*(r-1)/2;
A = zeros(numel(B),nCoordinates);
x0 = zeros(nCoordinates,1);
coordinate = 0;
for j = 1:r-1
    for k = j+1:r
        coordinate = coordinate+1;
        T = sparse([j k],[k j],[1 -1]/sqrt(2),r,r);
        A(:,coordinate) = reshape(B*T,[],1);
        x0(coordinate) = sum(T.*H0,"all");
    end
end
end

function H = realSkewMatrixFromCoordinates(x,r)
H = zeros(r);
coordinate = 0;
for j = 1:r-1
    for k = j+1:r
        coordinate = coordinate+1;
        H(j,k) = x(coordinate)/sqrt(2);
        H(k,j) = -x(coordinate)/sqrt(2);
    end
end
end

function [U,singularValues,V,numericalRank,tolerance] = singularFactors(A)
if isempty(A)
    U = zeros(size(A,1),0);
    singularValues = zeros(0,1);
    V = zeros(size(A,2),0);
    numericalRank = 0;
    tolerance = 0;
    return
end
[U,S,V] = svd(A,"econ");
singularValues = diag(S);
tolerance = max(size(A))*eps(max(singularValues));
numericalRank = nnz(singularValues > tolerance);
end

function [rightVectors,eigenvalues] = rightSingularEigensystem(A)
n = size(A,2);
if size(A,1) >= n
    [~,S,rightVectors] = svd(A,"econ");
    singularValues = diag(S);
else
    [~,S,rightVectors] = svd(A);
    singularValues = [diag(S);zeros(n-size(A,1),1)];
end
[singularValues,order] = sort(singularValues);
rightVectors = real(rightVectors(:,order));
eigenvalues = singularValues.^2;
end

function x = applyPseudoinverse(U,singularValues,V,numericalRank,rhs)
if numericalRank == 0
    x = zeros(size(V,1),1);
    return
end
x = V(:,1:numericalRank)*((U(:,1:numericalRank)'*rhs)./singularValues(1:numericalRank));
end

function defect = imaginaryDefect(A)
defect = norm(imag(A),"fro")/max(norm(A,"fro"),realmin);
end

function defect = relativeSkewDefect(K)
defect = norm(K+K',"fro")/max(norm(K,"fro"),realmin);
end

function defect = relativeCommutatorDefect(G,K)
defect = norm(G*K-K*G,"fro")/max(norm(G,"fro")*norm(K,"fro"),realmin);
end

function defect = relativeBottomDefect(B,K,R)
defect = norm(B*K-R,"fro")/max(norm(B,"fro")*norm(K,"fro")+norm(R,"fro"),realmin);
end
