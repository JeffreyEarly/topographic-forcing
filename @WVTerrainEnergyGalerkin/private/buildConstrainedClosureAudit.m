function audit = buildConstrainedClosureAudit(problem)
% Audit a dense energy-, APV-, and bottom-compatible closure.

forms = problem.finiteTerrainForms;
E = forms.energyMatrix;
J = forms.exchangeMatrix;
nState = size(E,1);
S = chol(E);
rawGenerator = (S'\J)/S;
rawGeneratorSkewDefect = norm(rawGenerator+rawGenerator',"fro")/max(norm(rawGenerator,"fro"),realmin);
rawGenerator = (rawGenerator-rawGenerator')/2;

apvMap = forms.apvMatrix/S;
[apvSingularMatrix,apvRightVectors] = rightSingularFactors(apvMap);
apvSingularValues = diag(apvSingularMatrix);
apvRankTolerance = numericalRankTolerance(apvSingularValues,size(apvMap));
apvRank = nnz(apvSingularValues > apvRankTolerance);
apvNullBasis = apvRightVectors(:,apvRank+1:end);
apvNullity = size(apvNullBasis,2);
apvProjector = apvNullBasis*apvNullBasis';
expectedWaveDimension = flatWaveDimension(problem);
expectedBottomInversionDimension = nnz(hypot(problem.horizontalLayout.k,problem.horizontalLayout.l) > 0);
expectedZeroAPVDimension = expectedWaveDimension+expectedBottomInversionDimension;
waveSpaceIsComplete = apvNullity >= expectedWaveDimension;

[bottomValueMatrix,bottomTendencyMatrix,bottomDiagnostics] = buildBottomConstraintMaps(problem);
bottomValueMap = bottomValueMatrix/S;
bottomTendencyMap = bottomTendencyMatrix/S;
identity = eye(nState);
rightCompatibilityMatrix = bottomTendencyMap*(identity-apvProjector);
rightCompatibilityDefect = norm(rightCompatibilityMatrix,"fro")/max(norm(bottomTendencyMap,"fro"),realmin);
skewCompatibilityMatrix = bottomTendencyMap*bottomValueMap'+bottomValueMap*bottomTendencyMap';
skewCompatibilityDefect = norm(skewCompatibilityMatrix,"fro")/max(norm(bottomTendencyMap*bottomValueMap',"fro"),realmin);

coordinateConjugacy = sparse((1:nState)',problem.conjugateCoordinateIndex,1,nState,nState);
energyCoordinateConjugacy = S*coordinateConjugacy/conj(S);
rawConjugacyDefect = conjugacyDefect(rawGenerator,energyCoordinateConjugacy);

reducedBottomValue = bottomValueMap*apvNullBasis;
reducedBottomTendency = bottomTendencyMap*apvNullBasis;
reducedRawGenerator = apvNullBasis'*rawGenerator*apvNullBasis;
reducedRawGenerator = (reducedRawGenerator-reducedRawGenerator')/2;
[constraintMatrix,rawCoordinates] = skewConstraintSystem(reducedBottomValue,reducedRawGenerator);
target = [real(reducedBottomTendency(:));imag(reducedBottomTendency(:))];
[constraintLeftVectors,constraintSingularMatrix,constraintRightVectors] = svd(constraintMatrix,"econ");
constraintSingularValues = diag(constraintSingularMatrix);
constraintRankTolerance = numericalRankTolerance(constraintSingularValues,size(constraintMatrix));
constraintRank = nnz(constraintSingularValues > constraintRankTolerance);

minimumResidualCoordinates = applyPseudoinverse(constraintLeftVectors,constraintSingularValues,constraintRightVectors,constraintRank,target);
reducedMinimumResidual = norm(constraintMatrix*minimumResidualCoordinates-target);
rightMinimumResidual = norm(rightCompatibilityMatrix,"fro");
minimumConstraintResidual = hypot(reducedMinimumResidual,rightMinimumResidual);
constraintTolerance = 1e-12*max(1,norm(bottomTendencyMap,"fro"));
constraintSystemIsFeasible = minimumConstraintResidual <= constraintTolerance;

constrainedGenerator = [];
constrainedExchangeMatrix = [];
correctionNorm = NaN;
apvConstraintResidual = NaN;
bottomConstraintResidual = NaN;
energyConstraintResidual = NaN;
constrainedConjugacyDefect = NaN;
if constraintSystemIsFeasible
    rawAPVConstraintResidual = norm(apvMap*rawGenerator,"fro");
    rawBottomConstraintResidual = norm(bottomValueMap*rawGenerator-bottomTendencyMap,"fro");
    if max(rawAPVConstraintResidual,rawBottomConstraintResidual) <= constraintTolerance
        constrainedGenerator = rawGenerator;
    else
        correctionTarget = target-constraintMatrix*rawCoordinates;
        correctionCoordinates = applyPseudoinverse(constraintLeftVectors,constraintSingularValues,constraintRightVectors,constraintRank,correctionTarget);
        constrainedCoordinates = rawCoordinates+correctionCoordinates;
        reducedConstrainedGenerator = skewMatrixFromCoordinates(constrainedCoordinates,apvNullity);
        constrainedGenerator = apvNullBasis*reducedConstrainedGenerator*apvNullBasis';
    end
    constrainedExchangeMatrix = S'*constrainedGenerator*S;
    correctionNorm = norm(constrainedGenerator-rawGenerator,"fro")/max(norm(rawGenerator,"fro"),realmin);
    apvConstraintResidual = norm(apvMap*constrainedGenerator,"fro")/max(1,norm(constrainedGenerator,"fro"));
    bottomConstraintResidual = norm(bottomValueMap*constrainedGenerator-bottomTendencyMap,"fro")/max(1,norm(bottomTendencyMap,"fro"));
    energyConstraintResidual = norm(constrainedGenerator+constrainedGenerator',"fro")/max(1,norm(constrainedGenerator,"fro"));
    constrainedConjugacyDefect = conjugacyDefect(constrainedGenerator,energyCoordinateConjugacy);
end

failedConditions = strings(0,1);
if ~waveSpaceIsComplete
    failedConditions(end+1,1) = "insufficient-apv-nullspace";
end
if ~constraintSystemIsFeasible
    failedConditions(end+1,1) = "incompatible-bottom-constraint";
end
if rawConjugacyDefect > 1e-12
    failedConditions(end+1,1) = "raw-conjugacy-defect";
end
if constraintSystemIsFeasible && constrainedConjugacyDefect > 1e-12
    failedConditions(end+1,1) = "constrained-conjugacy-defect";
end
if constraintSystemIsFeasible && apvConstraintResidual > 1e-12
    failedConditions(end+1,1) = "apv-constraint-residual";
end
if constraintSystemIsFeasible && bottomConstraintResidual > 1e-12
    failedConditions(end+1,1) = "bottom-constraint-residual";
end
if constraintSystemIsFeasible && energyConstraintResidual > 1e-13
    failedConditions(end+1,1) = "energy-constraint-residual";
end
isScientificallyAdmissible = isempty(failedConditions);
status = "admissible";
if ~isScientificallyAdmissible
    status = "incompatible";
end

diagnostics = struct;
diagnostics.rawGeneratorSkewDefect = rawGeneratorSkewDefect;
diagnostics.rawAPVDefect = norm(apvMap*rawGenerator,"fro")/max(1,norm(rawGenerator,"fro"));
diagnostics.rawBottomDefect = norm(bottomValueMap*rawGenerator-bottomTendencyMap,"fro")/max(1,norm(bottomTendencyMap,"fro"));
diagnostics.apvSingularValues = apvSingularValues;
diagnostics.apvRankTolerance = apvRankTolerance;
diagnostics.apvRank = apvRank;
diagnostics.apvNullity = apvNullity;
diagnostics.expectedWaveDimension = expectedWaveDimension;
diagnostics.expectedBottomInversionDimension = expectedBottomInversionDimension;
diagnostics.expectedZeroAPVDimension = expectedZeroAPVDimension;
diagnostics.rightCompatibilityDefect = rightCompatibilityDefect;
diagnostics.skewCompatibilityDefect = skewCompatibilityDefect;
diagnostics.minimumConstraintResidual = minimumConstraintResidual;
diagnostics.minimumConstraintRelativeResidual = minimumConstraintResidual/max(1,norm(bottomTendencyMap,"fro"));
diagnostics.constraintTolerance = constraintTolerance;
diagnostics.constraintSingularValues = constraintSingularValues;
diagnostics.constraintRankTolerance = constraintRankTolerance;
diagnostics.constraintRank = constraintRank;
diagnostics.correctionNorm = correctionNorm;
diagnostics.apvConstraintResidual = apvConstraintResidual;
diagnostics.bottomConstraintResidual = bottomConstraintResidual;
diagnostics.energyConstraintResidual = energyConstraintResidual;
diagnostics.rawConjugacyDefect = rawConjugacyDefect;
diagnostics.constrainedConjugacyDefect = constrainedConjugacyDefect;
diagnostics.bottom = bottomDiagnostics;

audit = struct;
audit.status = status;
audit.constraintSystemIsFeasible = constraintSystemIsFeasible;
audit.waveSpaceIsComplete = waveSpaceIsComplete;
audit.isScientificallyAdmissible = isScientificallyAdmissible;
audit.failedConditions = failedConditions;
audit.energyFactor = S;
audit.rawGenerator = rawGenerator;
audit.constrainedGenerator = constrainedGenerator;
audit.constrainedExchangeMatrix = constrainedExchangeMatrix;
audit.apvMap = apvMap;
audit.bottomValueMatrix = bottomValueMatrix;
audit.bottomTendencyMatrix = bottomTendencyMatrix;
audit.bottomValueMap = bottomValueMap;
audit.bottomTendencyMap = bottomTendencyMap;
audit.diagnostics = diagnostics;
end

function count = flatWaveDimension(problem)
count = 0;
for iK = 1:numel(problem.flatModeBlocks)
    block = problem.flatModeBlocks{iK};
    tolerance = 100*eps*size(block.E,1)*max([abs(problem.originatingTransform.f);abs(block.frequency)]);
    count = count+nnz(abs(block.frequency) > tolerance);
end
end

function [A,x0] = skewConstraintSystem(C,H0)
r = size(H0,1);
nCoordinates = r^2;
nEquations = 2*numel(C*H0);
A = zeros(nEquations,nCoordinates);
x0 = zeros(nCoordinates,1);
coordinate = 0;
for j = 1:r
    coordinate = coordinate+1;
    T = sparse(j,j,1i,r,r);
    [A(:,coordinate),x0(coordinate)] = constraintColumn(C,T,H0);
end
for j = 1:r-1
    for k = j+1:r
        coordinate = coordinate+1;
        T = sparse([j k],[k j],[1 -1]/sqrt(2),r,r);
        [A(:,coordinate),x0(coordinate)] = constraintColumn(C,T,H0);
        coordinate = coordinate+1;
        T = sparse([j k],[k j],[1i 1i]/sqrt(2),r,r);
        [A(:,coordinate),x0(coordinate)] = constraintColumn(C,T,H0);
    end
end
end

function [column,coordinate] = constraintColumn(C,T,H0)
product = C*T;
column = [real(product(:));imag(product(:))];
coordinate = real(sum(conj(T).*H0,"all"));
end

function H = skewMatrixFromCoordinates(x,r)
H = zeros(r);
coordinate = 0;
for j = 1:r
    coordinate = coordinate+1;
    H(j,j) = 1i*x(coordinate);
end
for j = 1:r-1
    for k = j+1:r
        coordinate = coordinate+1;
        H(j,k) = H(j,k)+x(coordinate)/sqrt(2);
        H(k,j) = H(k,j)-x(coordinate)/sqrt(2);
        coordinate = coordinate+1;
        H(j,k) = H(j,k)+1i*x(coordinate)/sqrt(2);
        H(k,j) = H(k,j)+1i*x(coordinate)/sqrt(2);
    end
end
end

function x = applyPseudoinverse(U,singularValues,V,numericalRank,rhs)
if numericalRank == 0
    x = zeros(size(V,1),1);
    return
end
x = V(:,1:numericalRank)*((U(:,1:numericalRank)'*rhs)./singularValues(1:numericalRank));
end

function tolerance = numericalRankTolerance(singularValues,matrixSize)
if isempty(singularValues) || max(singularValues) == 0
    tolerance = 0;
    return
end
tolerance = max(matrixSize)*eps(max(singularValues));
end

function [singularMatrix,rightVectors] = rightSingularFactors(A)
if size(A,1) >= size(A,2)
    [~,singularMatrix,rightVectors] = svd(A,"econ");
else
    [~,singularMatrix,rightVectors] = svd(A);
end
end

function defect = conjugacyDefect(K,C)
defect = norm(K*C-C*conj(K),"fro")/max(norm(K*C,"fro"),realmin);
end
