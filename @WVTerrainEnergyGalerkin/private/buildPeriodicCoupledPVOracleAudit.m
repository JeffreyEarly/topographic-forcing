function audit = buildPeriodicCoupledPVOracleAudit(problem,degree,quadratureOrder)
% Build and classify the projected periodic volume--boundary PV oracle.

wvt = problem.originatingTransform;
quadrature = legendreQuadrature(quadratureOrder,wvt.Lz);
layout = pvLayout(problem,degree);
[blocks,forms] = inversionForms(problem,layout,degree,quadrature);
terrain = terrainSpectrum(problem);
[exactBottomOperator,sidebands] = exactBottomJacobian(problem,layout,terrain);
pseudospectral = pseudospectralBottomJacobian(problem,layout);
pseudospectral.diagnostics.exactConvolutionDefect = matrixDefect( ...
    pseudospectral.bottomOperator-exactBottomOperator,exactBottomOperator);
dynamics = coupledDynamics(layout,forms,exactBottomOperator);
fourier = fourierDiagnostics(layout,terrain,exactBottomOperator,sidebands,dynamics);
bridge = balancedBasisBridge(problem,layout,blocks,dynamics);

requiredTolerance = struct("structure",1e-12,"dynamics",1e-12, ...
    "pseudospectral",1e-11,"frequencyImaginary",1e-11, ...
    "basisReconstruction",2e-2);
inversionPasses = forms.diagnostics.maximumInversionDefect <= requiredTolerance.structure ...
    && forms.diagnostics.maximumGreenIdentityDefect <= requiredTolerance.structure ...
    && forms.diagnostics.minimumEnergyEigenvalue > 0;
exactPasses = dynamics.diagnostics.bottomOperatorSkewDefect <= requiredTolerance.structure ...
    && dynamics.diagnostics.energyDefect <= requiredTolerance.dynamics ...
    && dynamics.diagnostics.apvDefect <= requiredTolerance.dynamics ...
    && dynamics.diagnostics.enstrophyDefect <= requiredTolerance.dynamics ...
    && dynamics.diagnostics.bottomDefect <= requiredTolerance.dynamics ...
    && dynamics.diagnostics.conjugacyDefect <= requiredTolerance.dynamics ...
    && dynamics.diagnostics.frequencyImaginaryDefect <= requiredTolerance.frequencyImaginary;
pseudospectralPasses = pseudospectral.diagnostics.exactConvolutionDefect <= requiredTolerance.pseudospectral ...
    && pseudospectral.diagnostics.meanTendencyDefect <= requiredTolerance.pseudospectral;
basisPasses = isempty(bridge.resolvedIndices) ...
    || (bridge.diagnostics.rank == bridge.diagnostics.expectedRank ...
    && bridge.diagnostics.maximumBottomValueDefect <= requiredTolerance.structure ...
    && bridge.diagnostics.resolvedReconstructionDefect <= requiredTolerance.basisReconstruction);

if ~inversionPasses
    status = "inversion-blocker";
    diagnosis = "At least one volume--boundary PV inversion block fails positivity or its Green identity.";
elseif ~exactPasses
    status = "exact-convolution-blocker";
    diagnosis = "The exact projected terrain convolution fails energy, APV, enstrophy, bottom, frequency, or conjugacy closure.";
elseif ~pseudospectralPasses
    status = "pseudospectral-blocker";
    diagnosis = "The oversampled pseudospectral bottom Jacobian does not reproduce the exact projected convolution.";
elseif ~basisPasses
    status = "basis-bridge-blocker";
    diagnosis = "The periodic coupled-PV oracle passes, but its active modes are not resolved by the existing balanced-plus-bottom coordinates.";
else
    status = "passed";
    diagnosis = "The projected periodic volume--boundary PV system closes at Fourier edges and is represented by the existing boundary-complete balanced basis.";
end

audit = struct();
audit.scope = "periodic-qg-coupled-volume-boundary-pv-oracle";
audit.status = status;
audit.isCompatible = status == "passed";
audit.diagnosis = diagnosis;
audit.polynomialDegree = degree;
audit.quadratureOrder = quadratureOrder;
audit.layout = layout;
audit.inversionBlocks = blocks;
audit.forms = forms;
audit.terrain = terrain;
audit.exactConvolution = struct("bottomOperator",exactBottomOperator, ...
    "sidebands",sidebands);
audit.pseudospectral = pseudospectral;
audit.dynamics = dynamics;
audit.fourier = fourier;
audit.bridge = bridge;
audit.requiredTolerance = requiredTolerance;
audit.nextScope = "hybrid-primitive-pv-representation-only-if-passed";
end

function layout = pvLayout(problem,degree)
activeHorizontalIndex = find(hypot(problem.horizontalLayout.k,problem.horizontalLayout.l) > 0);
nK = numel(activeHorizontalIndex);
nq = degree-1;
nBlock = nq+1;
ranges = cell(nK,1);
qRanges = cell(nK,1);
bottomRows = zeros(nK,1);
for i = 1:nK
    ranges{i} = (i-1)*nBlock+(1:nBlock);
    qRanges{i} = ranges{i}(1:nq);
    bottomRows(i) = ranges{i}(end);
end

activeLookup = zeros(height(problem.horizontalLayout),1);
activeLookup(activeHorizontalIndex) = (1:nK)';
conjugateHorizontalIndex = zeros(nK,1);
for i = 1:nK
    iK = activeHorizontalIndex(i);
    partner = find(problem.horizontalLayout.kMode == -problem.horizontalLayout.kMode(iK) ...
        & problem.horizontalLayout.lMode == -problem.horizontalLayout.lMode(iK),1);
    conjugateHorizontalIndex(i) = activeLookup(partner);
end
coordinateConjugateIndex = zeros(nK*nBlock,1);
for i = 1:nK
    coordinateConjugateIndex(ranges{i}) = ranges{conjugateHorizontalIndex(i)};
end

layout = struct("activeHorizontalIndex",activeHorizontalIndex, ...
    "numberOfHorizontalModes",nK,"numberOfVolumeCoefficientsPerMode",nq, ...
    "blockSize",nBlock,"numberOfStateCoefficients",nK*nBlock, ...
    "ranges",{ranges},"volumeRanges",{qRanges},"bottomRows",bottomRows, ...
    "conjugateHorizontalIndex",conjugateHorizontalIndex, ...
    "coordinateConjugateIndex",coordinateConjugateIndex, ...
    "kMode",problem.horizontalLayout.kMode(activeHorizontalIndex), ...
    "lMode",problem.horizontalLayout.lMode(activeHorizontalIndex), ...
    "k",problem.horizontalLayout.k(activeHorizontalIndex), ...
    "l",problem.horizontalLayout.l(activeHorizontalIndex));
end

function [blocks,forms] = inversionForms(problem,layout,degree,quadrature)
nK = layout.numberOfHorizontalModes;
nState = layout.numberOfStateCoefficients;
nq = layout.numberOfVolumeCoefficientsPerMode;
blocks = cell(nK,1);
E = zeros(nState);
Z = zeros(nState);
Q = zeros(nK*nq,nState);
B = zeros(nK,nState);
psiBottomMap = zeros(nK,nState);
inversionDefect = zeros(nK,1);
greenDefect = zeros(nK,1);
surfaceFluxDefect = zeros(nK,1);
bottomFluxDefect = zeros(nK,1);
for i = 1:nK
    block = buildCoupledPVInversionBlock(problem.originatingTransform,hypot(layout.k(i),layout.l(i)),degree,quadrature);
    blocks{i} = block;
    rows = layout.ranges{i};
    qRows = layout.volumeRanges{i};
    outputQ = (i-1)*nq+(1:nq);
    E(rows,rows) = block.energyMatrix;
    Z(qRows,qRows) = block.volumeMassMatrix;
    Q(outputQ,qRows) = eye(nq);
    B(i,layout.bottomRows(i)) = 1;
    psiBottomMap(i,rows) = block.stateToPsiBottom;
    inversionDefect(i) = block.diagnostics.inversionDefect;
    greenDefect(i) = block.diagnostics.greenIdentityDefect;
    surfaceFluxDefect(i) = block.diagnostics.surfaceFluxDefect;
    bottomFluxDefect(i) = block.diagnostics.bottomFluxDefect;
end
energyEigenvalue = eig((E+E')/2);
forms = struct("energyMatrix",E,"potentialEnstrophyMatrix",Z, ...
    "volumeAPVSelector",Q,"bottomSelector",B,"psiBottomMap",psiBottomMap);
forms.diagnostics = struct("inversionDefect",inversionDefect, ...
    "greenIdentityDefect",greenDefect, ...
    "surfaceFluxDefect",surfaceFluxDefect, ...
    "bottomFluxDefect",bottomFluxDefect, ...
    "maximumInversionDefect",max(inversionDefect), ...
    "maximumGreenIdentityDefect",max(greenDefect), ...
    "maximumSurfaceFluxDefect",max(surfaceFluxDefect), ...
    "maximumBottomFluxDefect",max(bottomFluxDefect), ...
    "energyHermitianDefect",matrixDefect(E-E',E), ...
    "enstrophyHermitianDefect",matrixDefect(Z-Z',Z), ...
    "minimumEnergyEigenvalue",min(real(energyEigenvalue)), ...
    "maximumEnergyEigenvalue",max(real(energyEigenvalue)));
end

function terrain = terrainSpectrum(problem)
wvt = problem.originatingTransform;
spectrum = fft2(problem.topographicHeight)/(wvt.Nx*wvt.Ny);
nyquistMask = logical(WVGeometryDoublyPeriodic.maskForNyquistModes(wvt.Nx,wvt.Ny));
scale = max(abs(spectrum),[],"all");
tolerance = 100*eps*max(scale,1);
nyquistAmplitude = max(abs(spectrum(nyquistMask)),[],"all");
if nyquistAmplitude > tolerance
    error("WVTerrainEnergyGalerkin:PeriodicPVOracleTerrainNyquist", ...
        "The periodic coupled-PV oracle requires terrain with no material Nyquist coefficient.")
end
[kMode,lMode] = ndgrid(wvt.kMode_dft,wvt.lMode_dft);
mask = abs(spectrum) > tolerance & ~nyquistMask;
coefficient = spectrum(mask);
kMode = kMode(mask);
lMode = lMode(mask);
isMean = kMode == 0 & lMode == 0;
terrain = struct("coefficient",coefficient(~isMean),"kMode",kMode(~isMean), ...
    "lMode",lMode(~isMean),"meanCoefficient",spectrum(1,1), ...
    "spectralTolerance",tolerance,"nyquistAmplitude",nyquistAmplitude, ...
    "numberOfActiveModes",nnz(~isMean));
end

function [C,sidebands] = exactBottomJacobian(problem,layout,terrain)
wvt = problem.originatingTransform;
nK = layout.numberOfHorizontalModes;
C = zeros(nK);
discardedNorm = zeros(nK,1);
meanCoefficient = zeros(nK,1);
discardedCount = zeros(nK,1);
for iIn = 1:nK
    pMode = [layout.kMode(iIn) layout.lMode(iIn)];
    p = [layout.k(iIn) layout.l(iIn)];
    for iTerrain = 1:terrain.numberOfActiveModes
        qMode = [terrain.kMode(iTerrain) terrain.lMode(iTerrain)];
        q = 2*pi*[qMode(1)/wvt.Lx qMode(2)/wvt.Ly];
        destination = pMode+qMode;
        factor = wvt.f*(p(1)*q(2)-p(2)*q(1))*terrain.coefficient(iTerrain);
        iOut = find(layout.kMode == destination(1) & layout.lMode == destination(2),1);
        if ~isempty(iOut)
            C(iOut,iIn) = C(iOut,iIn)+factor;
        elseif all(destination == 0)
            meanCoefficient(iIn) = meanCoefficient(iIn)+factor;
        else
            discardedNorm(iIn) = hypot(discardedNorm(iIn),abs(factor));
            discardedCount(iIn) = discardedCount(iIn)+1;
        end
    end
end
sidebands = struct("discardedOperatorNormByInput",discardedNorm, ...
    "discardedCountByInput",discardedCount, ...
    "meanCoefficientByInput",meanCoefficient, ...
    "interiorHorizontalModes",discardedCount == 0, ...
    "edgeHorizontalModes",discardedCount > 0, ...
    "numberOfInteriorHorizontalModes",nnz(discardedCount == 0), ...
    "numberOfEdgeHorizontalModes",nnz(discardedCount > 0), ...
    "maximumDiscardedOperatorNorm",max(discardedNorm), ...
    "maximumMeanCoefficient",max(abs(meanCoefficient)));
end

function pseudo = pseudospectralBottomJacobian(problem,layout)
wvt = problem.originatingTransform;
oversampling = problem.horizontalOversamplingFactor;
Nx = oversampling*wvt.Nx;
Ny = oversampling*wvt.Ny;
[x,y] = ndgrid((0:Nx-1)'*wvt.Lx/Nx,(0:Ny-1)'*wvt.Ly/Ny);
phase = exp(1i*(x(:)*layout.k.'+y(:)*layout.l.'));
h = interpft(interpft(problem.topographicHeight,Nx,1),Ny,2);
[hX,hY] = horizontalDerivatives(h,wvt.Lx,wvt.Ly);
nK = layout.numberOfHorizontalModes;
C = zeros(nK);
meanTendency = zeros(nK,1);
for iIn = 1:nK
    tendency = -(1i*layout.k(iIn)*phase(:,iIn)).*(wvt.f*hY(:)) ...
        +(1i*layout.l(iIn)*phase(:,iIn)).*(wvt.f*hX(:));
    C(:,iIn) = phase'*tendency/(Nx*Ny);
    meanTendency(iIn) = sum(tendency)/(Nx*Ny);
end
pseudo = struct("bottomOperator",C,"oversampledSize",[Nx Ny], ...
    "phaseGramDefect",norm(phase'*phase/(Nx*Ny)-eye(nK),"fro")/sqrt(nK), ...
    "meanTendency",meanTendency);
pseudo.diagnostics = struct("exactConvolutionDefect",NaN, ...
    "meanTendencyDefect",max(abs(meanTendency))/max(norm(C,"fro"),realmin), ...
    "terrainInterpolationImaginaryDefect",max(abs(imag(h)),[],"all")/max(max(abs(h),[],"all"),realmin));
end

function dynamics = coupledDynamics(layout,forms,C)
nState = layout.numberOfStateCoefficients;
L = zeros(nState);
R = C*forms.psiBottomMap;
L(layout.bottomRows,:) = R;
E = forms.energyMatrix;
Q = forms.volumeAPVSelector;
Z = forms.potentialEnstrophyMatrix;
B = forms.bottomSelector;
lambda = eig(L);
frequencyComplex = 1i*lambda;
frequencyScale = max(max(abs(real(frequencyComplex))),norm(L,2)*eps);
if frequencyScale == 0
    frequencyImaginaryDefect = 0;
else
    frequencyImaginaryDefect = max(abs(imag(frequencyComplex)))/frequencyScale;
end
ci = layout.coordinateConjugateIndex;
dynamics = struct("generator",L,"bottomTendencyMap",R, ...
    "frequency",real(frequencyComplex),"frequencyComplex",frequencyComplex);
dynamics.diagnostics = struct( ...
    "bottomOperatorSkewDefect",matrixDefect(C+C',C), ...
    "energyDefect",productDefect(L'*E+E*L,{L'*E,E*L}), ...
    "apvDefect",norm(Q*L,"fro")/max(norm(Q,"fro")*norm(L,"fro"),realmin), ...
    "enstrophyDefect",productDefect(L'*Z+Z*L,{L'*Z,Z*L}), ...
    "bottomDefect",productDefect(B*L-R,{B*L,R}), ...
    "conjugacyDefect",matrixDefect(L(ci,ci)-conj(L),L), ...
    "frequencyImaginaryDefect",frequencyImaginaryDefect);
end

function diagnostics = fourierDiagnostics(layout,terrain,C,sidebands,dynamics)
allowed = false(size(C));
for iOut = 1:layout.numberOfHorizontalModes
    for iIn = 1:layout.numberOfHorizontalModes
        difference = [layout.kMode(iOut)-layout.kMode(iIn), ...
            layout.lMode(iOut)-layout.lMode(iIn)];
        allowed(iOut,iIn) = any(terrain.kMode == difference(1) & terrain.lMode == difference(2));
    end
end
disallowed = C;
disallowed(allowed) = 0;
nonzeroFrequency = abs(dynamics.frequencyComplex) > 1e-10*max(norm(dynamics.generator,2),realmin);
if any(nonzeroFrequency)
    [V,D] = eig(dynamics.generator,"vector");
    selected = abs(D) > 1e-10*max(norm(dynamics.generator,2),realmin);
    qRows = setdiff(1:layout.numberOfStateCoefficients,layout.bottomRows);
    nonzeroModeVolumeAPVDefect = norm(V(qRows,selected),"fro")/max(norm(V(:,selected),"fro"),realmin);
else
    nonzeroModeVolumeAPVDefect = 0;
end
diagnostics = struct("allowedHorizontalCouplings",allowed, ...
    "couplingLeakage",norm(disallowed,"fro")/max(norm(C,"fro"),realmin), ...
    "terrainModes",[terrain.kMode terrain.lMode], ...
    "interiorHorizontalModes",sidebands.interiorHorizontalModes, ...
    "edgeHorizontalModes",sidebands.edgeHorizontalModes, ...
    "nonzeroModeVolumeAPVDefect",nonzeroModeVolumeAPVDefect, ...
    "numberOfNonzeroFrequencies",nnz(nonzeroFrequency), ...
    "projectedEdgeAPVDefect",dynamics.diagnostics.apvDefect);
end

function bridge = balancedBasisBridge(problem,layout,blocks,dynamics)
wvt = problem.originatingTransform;
nState = layout.numberOfStateCoefficients;
coefficientMap = zeros(height(problem.stateLayout),nState);
targetMaps = cell(layout.numberOfHorizontalModes,1);
reconstructionMaps = cell(layout.numberOfHorizontalModes,1);
rankByMode = zeros(layout.numberOfHorizontalModes,1);
expectedRankByMode = zeros(layout.numberOfHorizontalModes,1);
bottomValueDefectByMode = zeros(layout.numberOfHorizontalModes,1);
rNative = 2*wvt.z(:)/wvt.Lz+1;
for i = 1:layout.numberOfHorizontalModes
    iK = layout.activeHorizontalIndex(i);
    block = blocks{i};
    [Pnative,PrNative] = legendreValues(rNative,size(block.polynomial.psiTransformation,1)-1);
    psiMap = Pnative*block.polynomial.psiTransformation*block.inversionMap;
    psiXiMap = (2/wvt.Lz)*PrNative*block.polynomial.psiTransformation*block.inversionMap;
    uMap = -1i*layout.l(i)*psiMap;
    vMap = 1i*layout.k(i)*psiMap;
    etaMap = -(wvt.f./wvt.N2(:)).*psiXiMap;
    bottomMap = zeros(1,layout.blockSize);
    bottomMap(end) = -1/wvt.f;
    etaMap(1,:) = bottomMap;
    weight = wvt.z_int(:);
    target = [sqrt(weight).*uMap;sqrt(weight).*vMap;sqrt(weight.*wvt.N2(:)).*etaMap];

    publicRows = find(problem.stateLayout.horizontalIndex == iK);
    components = problem.stateLayout.component(publicRows);
    selectedLocal = [find(components == "A0");find(components == "etaB",1)];
    selectedRows = publicRows(selectedLocal);
    basis = problem.basisBlocks{iK};
    R = [sqrt(weight).*basis.uHat(:,selectedLocal); ...
        sqrt(weight).*basis.vHat(:,selectedLocal); ...
        sqrt(weight.*wvt.N2(:)).*basis.etaHat(:,selectedLocal)];
    balancedCoefficient = R(:,1:end-1)\(target-R(:,end).*bottomMap);
    localCoefficient = [balancedCoefficient;bottomMap];
    localReconstruction = R*localCoefficient;
    coefficientMap(selectedRows,layout.ranges{i}) = localCoefficient;
    targetMaps{i} = target;
    reconstructionMaps{i} = localReconstruction;
    rankByMode(i) = rank(R,max(size(R))*eps(max(norm(R,2),1)));
    expectedRankByMode(i) = size(R,2);
    bottomValueDefectByMode(i) = norm(localCoefficient(end,:)-bottomMap)/max(norm(bottomMap),realmin);
end

[V,lambda] = eig(dynamics.generator,"vector");
frequency = 1i*lambda;
nonzero = find(abs(frequency) > 1e-10*max(norm(dynamics.generator,2),realmin));
[~,order] = sort(abs(frequency(nonzero)),"descend");
resolved = nonzero(order(1:min(4,numel(order))));
defect = zeros(numel(resolved),1);
for iMode = 1:numel(resolved)
    targetNorm2 = 0;
    residualNorm2 = 0;
    for i = 1:layout.numberOfHorizontalModes
        local = V(layout.ranges{i},resolved(iMode));
        target = targetMaps{i}*local;
        residual = reconstructionMaps{i}*local-target;
        targetNorm2 = targetNorm2+norm(target)^2;
        residualNorm2 = residualNorm2+norm(residual)^2;
    end
    defect(iMode) = sqrt(residualNorm2/max(targetNorm2,realmin));
end
if isempty(defect)
    resolvedDefect = 0;
else
    resolvedDefect = max(defect);
end

bridge = struct("coefficientMap",coefficientMap,"resolvedIndices",resolved, ...
    "resolvedFrequency",real(frequency(resolved)), ...
    "resolvedReconstructionDefectByMode",defect);
bridge.diagnostics = struct("rankByHorizontalMode",rankByMode, ...
    "expectedRankByHorizontalMode",expectedRankByMode, ...
    "rank",sum(rankByMode),"expectedRank",sum(expectedRankByMode), ...
    "bottomValueDefectByHorizontalMode",bottomValueDefectByMode, ...
    "maximumBottomValueDefect",max(bottomValueDefectByMode), ...
    "resolvedReconstructionDefect",resolvedDefect, ...
    "geostrophicBalanceDefect",0,"hydrostaticBalanceDefect",0);
end

function quadrature = legendreQuadrature(order,D)
index = (1:order-1)';
offDiagonal = index./sqrt(4*index.^2-1);
[vectors,values] = eig(diag(offDiagonal,1)+diag(offDiagonal,-1),"vector");
[r,permutation] = sort(values);
vectors = vectors(:,permutation);
weightR = 2*(vectors(1,:)').^2;
quadrature = struct("r",r,"xi",D*(r-1)/2,"weight",D*weightR/2);
end

function [P,Pr] = legendreValues(r,degree)
P = zeros(numel(r),degree+1);
Pr = zeros(numel(r),degree+1);
P(:,1) = 1;
if degree == 0
    return
end
P(:,2) = r;
Pr(:,2) = 1;
for n = 2:degree
    P(:,n+1) = ((2*n-1)*r.*P(:,n)-(n-1)*P(:,n-1))/n;
    Pr(:,n+1) = ((2*n-1)*(P(:,n)+r.*Pr(:,n))-(n-1)*Pr(:,n-1))/n;
end
end

function [xDerivative,yDerivative] = horizontalDerivatives(values,Lx,Ly)
[Nx,Ny] = size(values);
k = 2*pi*[0:floor((Nx-1)/2) -floor(Nx/2):-1]'/Lx;
l = 2*pi*[0:floor((Ny-1)/2) -floor(Ny/2):-1]/Ly;
spectrum = fft2(values);
xDerivative = ifft2((1i*k).*spectrum);
yDerivative = ifft2(spectrum.*(1i*l));
end

function value = productDefect(residual,terms)
scale = 0;
for iTerm = 1:numel(terms)
    scale = scale+norm(terms{iTerm},"fro");
end
value = norm(residual,"fro")/max(scale,realmin);
end

function value = matrixDefect(residual,reference)
value = norm(residual,"fro")/max(norm(reference,"fro"),realmin);
end
