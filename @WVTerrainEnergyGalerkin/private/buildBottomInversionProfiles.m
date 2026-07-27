function [profiles,diagnostics] = buildBottomInversionProfiles(problem)
% Construct complete flat balanced bottom-inversion profiles.

wvt = problem.originatingTransform;
xi = wvt.z(:);
D = wvt.Lz;
f = wvt.f;
g = wvt.g;
N2 = wvt.N2(:);
kappa = hypot(problem.horizontalLayout.k,problem.horizontalLayout.l);
positiveKappa = unique(kappa(kappa > 0),"sorted");
nEVP = max(128,4*wvt.Nz);

profiles = cell(height(problem.horizontalLayout),1);
if isempty(positiveKappa)
    diagnostics = struct("solverOrder",nEVP,"uniquePositiveKappa",positiveKappa, ...
        "maximumBottomValueError",0,"maximumSurfaceValueError",0, ...
        "maximumAPVResidual",0,"maximumInversionResidual",0);
    return
end

N2Bottom = wvt.N2Function(-D);
if ~isscalar(N2Bottom) || ~isreal(N2Bottom) || ~isfinite(N2Bottom) || N2Bottom <= 0
    error("WVTerrainEnergyGalerkin:InvalidBottomStratification", ...
        "N2Function must return one positive finite value at the bottom.")
end

modeProblem = IMSurfaceGeostrophicModes.atWavenumber( ...
    N2=wvt.N2Function,zDomain=[-D 0],f0=f,g=g,k=positiveKappa, ...
    g0=Inf,gd=N2Bottom,surfaceAnomaly="noFreeSurface");
modeBasis = IMSolverSpectral(nEVP=nEVP).solveSurfaceGeostrophicModes(modeProblem);
psi = modeBasis.F(xi);
eta = (f/g)*modeBasis.G(xi);

normalization = eta(1,:);
if any(abs(normalization) <= 100*eps)
    error("WVTerrainEnergyGalerkin:BottomInversionFailure", ...
        "The bottom-inversion solver returned a vanishing bottom displacement.")
end
psi = psi./normalization;
eta = eta./normalization;
eta(1,:) = 1;
eta(end,:) = 0;

bottomError = abs(eta(1,:)-1);
surfaceError = abs(eta(end,:));
apvResidual = zeros(size(positiveKappa));
inversionResidual = zeros(size(positiveKappa));
for iMode = 1:numel(positiveKappa)
    currentKappa = positiveKappa(iMode);
    psiXi = -(N2/f).*eta(:,iMode);
    etaXi = -(currentKappa^2/f)*psi(:,iMode);
    q = -currentKappa^2*psi(:,iMode)-f*etaXi;
    apvScale = max(currentKappa^2*norm(psi(:,iMode))+abs(f)*norm(etaXi),realmin);
    apvResidual(iMode) = norm(q)/apvScale;

    nativePsi = modeBasis.nativeModes(:,iMode)/normalization(iMode);
    solver = modeBasis.solver;
    zNative = solver.zNative;
    nativePsiValues = solver.evaluateNativeModes(nativePsi,zNative);
    nativeN2 = wvt.N2Function(zNative);
    nativeN2 = nativeN2(:);
    nativeN2Xi = solver.differentiateGridValues(nativeN2,1);
    nativePsiXi = solver.evaluatePhysicalDerivative(nativePsi,zNative,1);
    nativePsiXiXi = solver.evaluatePhysicalDerivative(nativePsi,zNative,2);
    residual = -currentKappa^2*nativePsiValues ...
        +(f^2./nativeN2).*nativePsiXiXi ...
        -(f^2.*nativeN2Xi./nativeN2.^2).*nativePsiXi;
    interior = 2:(numel(zNative)-1);
    residualScale = max(currentKappa^2*norm(nativePsiValues(interior)) ...
        +norm((f^2./nativeN2(interior)).*nativePsiXiXi(interior)) ...
        +norm((f^2.*nativeN2Xi(interior)./nativeN2(interior).^2).*nativePsiXi(interior)),realmin);
    inversionResidual(iMode) = norm(residual(interior))/residualScale;

    profile = struct("kappa",currentKappa,"psi",psi(:,iMode),"psiXi",psiXi, ...
        "eta",eta(:,iMode),"etaXi",etaXi,"pressure",wvt.rho0*f*psi(:,iMode));
    profiles(kappa == currentKappa) = {profile};
end

meanEta = -xi/D;
meanEta(1) = 1;
meanEta(end) = 0;
zeroProfile = struct("kappa",0,"psi",zeros(size(xi)),"psiXi",zeros(size(xi)), ...
    "eta",meanEta,"etaXi",-ones(size(xi))/D,"pressure",zeros(size(xi)));
profiles(kappa == 0) = {zeroProfile};

diagnostics = struct();
diagnostics.solverClass = class(modeBasis.solver);
diagnostics.solverOrder = nEVP;
diagnostics.internalModesEVPCommit = "df86687e91faa31bf65941299062d125a96904b1";
diagnostics.uniquePositiveKappa = positiveKappa;
diagnostics.normalizationBeforeRescaling = normalization;
diagnostics.maximumBottomValueError = max(bottomError);
diagnostics.maximumSurfaceValueError = max(surfaceError);
diagnostics.maximumAPVResidual = max(apvResidual);
diagnostics.maximumInversionResidual = max(inversionResidual);
diagnostics.apvResidual = apvResidual;
diagnostics.inversionResidual = inversionResidual;
end
