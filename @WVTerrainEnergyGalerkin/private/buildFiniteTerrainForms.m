function [forms,diagnostics] = buildFiniteTerrainForms(problem)
% Construct the dense finite-terrain weak forms by common quadrature.

wvt = problem.originatingTransform;
nState = height(problem.stateLayout);
oversampling = problem.horizontalOversamplingFactor;
Nx = oversampling*wvt.Nx;
Ny = oversampling*wvt.Ny;
nXY = Nx*Ny;
Nz = wvt.Nz;
nGrid = nXY*Nz;

[x,y] = ndgrid((0:Nx-1)'*wvt.Lx/Nx,(0:Ny-1)'*wvt.Ly/Ny);
phase = exp(1i*(x(:)*problem.horizontalLayout.k.'+y(:)*problem.horizontalLayout.l.'));
[h,gamma,gammaX,gammaY,terrainInterpolationImaginaryPart] = oversampledTerrain(problem,Nx,Ny);
gradLnGammaX = gammaX./gamma;
gradLnGammaY = gammaY./gamma;

xi = wvt.z(:);
xiGrid = kron(xi,ones(nXY,1));
gammaGrid = repmat(gamma(:),Nz,1);
gradLnGammaXGrid = repmat(gradLnGammaX(:),Nz,1);
gradLnGammaYGrid = repmat(gradLnGammaY(:),Nz,1);
mappedDepth = gammaGrid.*xiGrid;
N2 = mappedStratification(wvt,mappedDepth);

uHat = zeros(nGrid,nState);
vHat = zeros(nGrid,nState);
wHat = zeros(nGrid,nState);
etaHat = zeros(nGrid,nState);
uHatXi = zeros(nGrid,nState);
vHatXi = zeros(nGrid,nState);
etaHatXi = zeros(nGrid,nState);
k = zeros(nState,1);
l = zeros(nState,1);
for iK = 1:height(problem.horizontalLayout)
    rows = find(problem.stateLayout.horizontalIndex == iK);
    block = problem.basisBlocks{iK};
    horizontalPhase = phase(:,iK);
    for iLocal = 1:numel(rows)
        row = rows(iLocal);
        uHat(:,row) = kron(block.uHat(:,iLocal),horizontalPhase);
        vHat(:,row) = kron(block.vHat(:,iLocal),horizontalPhase);
        wHat(:,row) = kron(block.wHat(:,iLocal),horizontalPhase);
        etaHat(:,row) = kron(block.etaHat(:,iLocal),horizontalPhase);
        uHatXi(:,row) = kron(block.uHatXi(:,iLocal),horizontalPhase);
        vHatXi(:,row) = kron(block.vHatXi(:,iLocal),horizontalPhase);
        etaHatXi(:,row) = kron(block.etaHatXi(:,iLocal),horizontalPhase);
        k(row) = problem.horizontalLayout.k(iK);
        l(row) = problem.horizontalLayout.l(iK);
    end
end

u = uHat./gammaGrid;
v = vHat./gammaGrid;
w = wHat+xiGrid.*(gradLnGammaXGrid.*uHat+gradLnGammaYGrid.*vHat);
volumeWeight = kron(wvt.z_int(:),ones(nXY,1)/nXY);
gammaWeight = volumeWeight.*gammaGrid;

rawEnergyMatrix = wvt.rho0*(uHat'*(volumeWeight./gammaGrid.*uHat)+vHat'*(volumeWeight./gammaGrid.*vHat) ...
    +w'*(gammaWeight.*w)+etaHat'*(gammaWeight.*N2.*etaHat));
exchange = wvt.f*u'*(gammaWeight.*v)+etaHat'*(gammaWeight.*N2.*w);
rawExchangeMatrix = wvt.rho0*(exchange-exchange');
quadratureEnergyMatrix = rawEnergyMatrix;
quadratureExchangeMatrix = rawExchangeMatrix;
if max(gamma,[],"all")-min(gamma,[],"all") <= 10*eps(max(abs(gamma),[],"all"))
    [rawEnergyMatrix,rawExchangeMatrix] = constantGammaForms(problem,gamma(1),N2(1:nXY:end));
end
if max(abs(problem.topographicHeight),[],"all") == 0
    [rawEnergyMatrix,rawExchangeMatrix] = globalFlatForms(problem);
end

uY = uHat.*(1i*l.')-gradLnGammaYGrid.*uHat;
vX = vHat.*(1i*k.')-gradLnGammaXGrid.*vHat;
mappedUY = uY./gammaGrid-xiGrid.*gradLnGammaYGrid.*(uHatXi./gammaGrid);
mappedVX = vX./gammaGrid-xiGrid.*gradLnGammaXGrid.*(vHatXi./gammaGrid);
unweightedAPVMatrix = mappedVX-mappedUY-(wvt.f./gammaGrid).*etaHatXi;
apvWeight = gammaWeight;
apvMatrix = sqrt(apvWeight).*unweightedAPVMatrix;

energyScale = 1./sqrt(real(diag(rawEnergyMatrix)));
scaledEnergy = (energyScale.*rawEnergyMatrix).*energyScale.';
hermitianDefect = norm(rawEnergyMatrix-rawEnergyMatrix',"fro")/max(norm(rawEnergyMatrix,"fro"),realmin);
skewHermitianDefect = norm(rawExchangeMatrix+rawExchangeMatrix',"fro")/max(norm(rawExchangeMatrix,"fro"),realmin);
scaledEnergyRcond = rcond((scaledEnergy+scaledEnergy')/2);
if hermitianDefect > 1e-12 || skewHermitianDefect > 1e-12
    error("WVTerrainEnergyGalerkin:FiniteTerrainStructureFailure", ...
        "Raw finite-terrain forms failed their structural gate: E %.3g, J %.3g.",hermitianDefect,skewHermitianDefect)
end
if scaledEnergyRcond <= 1e-12
    error("WVTerrainEnergyGalerkin:IllConditionedTerrainEnergy", ...
        "The scaled finite-terrain energy form is ill-conditioned; reciprocal condition number is %.3g.",scaledEnergyRcond)
end

energyMatrix = (rawEnergyMatrix+rawEnergyMatrix')/2;
exchangeMatrix = (rawExchangeMatrix-rawExchangeMatrix')/2;
[flatEnergyMatrix,flatExchangeMatrix] = globalFlatForms(problem);
flatEnergyRelativeResidual = norm(energyMatrix-flatEnergyMatrix,"fro")/max(norm(flatEnergyMatrix,"fro"),realmin);
flatExchangeRelativeResidual = norm(exchangeMatrix-flatExchangeMatrix,"fro")/max(norm(flatExchangeMatrix,"fro"),realmin);
if max(abs(problem.topographicHeight),[],"all") == 0 && max(flatEnergyRelativeResidual,flatExchangeRelativeResidual) > 1e-12
    error("WVTerrainEnergyGalerkin:FlatLimitFailure", ...
        "The finite-terrain forms do not recover the flat oracle: E %.3g, J %.3g (absolute J %.3g, flat norm %.3g).", ...
        flatEnergyRelativeResidual,flatExchangeRelativeResidual,norm(exchangeMatrix-flatExchangeMatrix,"fro"),norm(flatExchangeMatrix,"fro"))
end

testVector = complex((1:nState)',mod((1:nState)',7)-3);
testField = phase*complex((1:height(problem.horizontalLayout))',mod((1:height(problem.horizontalLayout))',5)-2);
leftAdjointProduct = (phase*testVectorForHorizontalLayout(problem,testVector))'*testField/nXY;
rightAdjointProduct = testVectorForHorizontalLayout(problem,testVector)'*(phase'*testField/nXY);
oversamplingAdjointDefect = abs(leftAdjointProduct-rightAdjointProduct)/max([abs(leftAdjointProduct),abs(rightAdjointProduct),realmin]);

forms = struct();
forms.rawEnergyMatrix = rawEnergyMatrix;
forms.rawExchangeMatrix = rawExchangeMatrix;
forms.quadratureEnergyMatrix = quadratureEnergyMatrix;
forms.quadratureExchangeMatrix = quadratureExchangeMatrix;
forms.energyMatrix = energyMatrix;
forms.exchangeMatrix = exchangeMatrix;
forms.apvMatrix = apvMatrix;
forms.unweightedAPVMatrix = unweightedAPVMatrix;
forms.apvWeight = apvWeight;
forms.flatEnergyMatrix = flatEnergyMatrix;
forms.flatExchangeMatrix = flatExchangeMatrix;
forms.oversampledTopographicHeight = h;
forms.oversampledGamma = gamma;
forms.oversampledGradLnGammaX = gradLnGammaX;
forms.oversampledGradLnGammaY = gradLnGammaY;
forms.oversampledN2 = reshape(N2,[Nx Ny Nz]);
forms.oversampledSize = [Nx Ny Nz];

diagnostics = struct();
diagnostics.hermitianDefect = hermitianDefect;
diagnostics.skewHermitianDefect = skewHermitianDefect;
diagnostics.scaledEnergyRcond = scaledEnergyRcond;
diagnostics.flatEnergyRelativeResidual = flatEnergyRelativeResidual;
diagnostics.flatExchangeRelativeResidual = flatExchangeRelativeResidual;
diagnostics.quadratureEnergyRelativeResidual = norm(quadratureEnergyMatrix-energyMatrix,"fro")/max(norm(energyMatrix,"fro"),realmin);
diagnostics.quadratureExchangeRelativeResidual = norm(quadratureExchangeMatrix-exchangeMatrix,"fro")/max(norm(exchangeMatrix,"fro"),realmin);
diagnostics.oversamplingAdjointDefect = oversamplingAdjointDefect;
diagnostics.maximumTerrainInterpolationImaginaryPart = terrainInterpolationImaginaryPart;
diagnostics.numberOfDegreesOfFreedom = nState;
diagnostics.estimatedDenseFormBytes = 16*(7*nState*nState+nGrid*nState);
end

function [h,gamma,gammaX,gammaY,imaginaryPart] = oversampledTerrain(problem,Nx,Ny)
wvt = problem.originatingTransform;
h = interpft(interpft(problem.topographicHeight,Nx,1),Ny,2);
imaginaryPart = max(abs(imag(h)),[],"all");
if imaginaryPart > 1e-12*max(1,max(abs(h),[],"all"))
    error("WVTerrainEnergyGalerkin:TerrainInterpolationFailure", ...
        "Periodic terrain interpolation produced a material imaginary component.")
end
h = real(h);
gamma = 1-h/wvt.Lz;
if any(gamma <= 0,"all")
    error("WVTerrainEnergyGalerkin:NonpositiveOversampledGamma", ...
        "The periodically interpolated terrain must satisfy gamma>0 on the oversampled grid.")
end

k = 2*pi*[0:floor((Nx-1)/2) -floor(Nx/2):-1]'/wvt.Lx;
l = 2*pi*[0:floor((Ny-1)/2) -floor(Ny/2):-1]/wvt.Ly;
hDFT = fft2(h);
gammaX = real(ifft2((-1i*k).*hDFT/wvt.Lz));
gammaY = real(ifft2(hDFT.*(-1i*l)/wvt.Lz));
end

function N2 = mappedStratification(wvt,mappedDepth)
try
    N2 = wvt.N2Function(mappedDepth);
catch exception
    error("WVTerrainEnergyGalerkin:InvalidMappedStratification", ...
        "N2Function could not be evaluated at mapped terrain depths: %s",exception.message)
end
if isscalar(N2)
    N2 = repmat(N2,size(mappedDepth));
end
if ~isa(N2,"double") || ~isequal(size(N2),size(mappedDepth)) || ~isreal(N2) || any(~isfinite(N2),"all") || any(N2 <= 0,"all")
    error("WVTerrainEnergyGalerkin:InvalidMappedStratification", ...
        "N2Function must return a real, finite, positive double array matching the mapped depths.")
end
end

function [E0,J0] = globalFlatForms(problem)
nState = height(problem.stateLayout);
E0 = zeros(nState);
J0 = zeros(nState);
for iK = 1:height(problem.horizontalLayout)
    rows = find(problem.stateLayout.horizontalIndex == iK);
    E0(rows,rows) = problem.flatModeBlocks{iK}.E;
    J0(rows,rows) = problem.flatModeBlocks{iK}.J;
end
end

function [E,J] = constantGammaForms(problem,gamma,N2)
wvt = problem.originatingTransform;
zWeight = wvt.z_int(:);
N2Weight = zWeight.*N2(:);
nState = height(problem.stateLayout);
E = zeros(nState);
J = zeros(nState);
for iK = 1:height(problem.horizontalLayout)
    rows = find(problem.stateLayout.horizontalIndex == iK);
    basis = problem.basisBlocks{iK};
    U = basis.uHat;
    V = basis.vHat;
    W = basis.wHat;
    Eta = basis.etaHat;
    E(rows,rows) = wvt.rho0*((U'*(zWeight.*U)+V'*(zWeight.*V))/gamma ...
        +gamma*(W'*(zWeight.*W)+Eta'*(N2Weight.*Eta)));
    exchange = (wvt.f/gamma)*U'*(zWeight.*V)+gamma*Eta'*(N2Weight.*W);
    J(rows,rows) = wvt.rho0*(exchange-exchange');
end
end

function values = testVectorForHorizontalLayout(problem,stateVector)
values = zeros(height(problem.horizontalLayout),1);
for iK = 1:height(problem.horizontalLayout)
    values(iK) = sum(stateVector(problem.stateLayout.horizontalIndex == iK));
end
end
