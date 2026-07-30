function family = buildBoundaryCompleteVerticalReferenceFamily(problem,kappa,numberOfModes,robinLength,solverOrder,z)
% Build the boundary-complete vertical reference family at fixed kappa.
arguments
    problem (1,1) WVTerrainEnergyGalerkin
    kappa (1,1) double {mustBePositive}
    numberOfModes (1,1) double {mustBeInteger,mustBePositive}
    robinLength (1,1) double {mustBeReal}
    solverOrder (1,1) double {mustBeInteger,mustBePositive}
    z (:,1) double {mustBeReal,mustBeFinite}
end

wvt = problem.originatingTransform;
D = wvt.Lz;
f = wvt.f;
g = wvt.g;
N2 = wvt.N2Function;
if robinLength == 0
    error("WVTerrainEnergyGalerkin:ZeroRobinLength", ...
        "The signed Robin length must be nonzero or Inf.")
end

solver = IMSolverSpectral(nEVP=solverOrder);
waveEVP = IMInternalModes.waveModesAtWavenumber( ...
    N2=N2,zDomain=[-D 0],k=kappa,f0=f,g=g);
waveBasis = solver.solveEVP(waveEVP,nModes=numberOfModes);
waveG = evaluateBasis(waveBasis,"G",z);
waveF = evaluateBasis(waveBasis,"F",z);
waveH = reshape(waveBasis.h,1,[]);
waveOmega = sqrt(f^2+g*kappa^2*waveH);

p = @(zValue,~) 1./N2(zValue);
if isinf(robinLength)
    bottomBoundary = IMBoundaryCondition.neumann();
else
    N2Bottom = N2(-D);
    pBottom = 1/N2Bottom;
    bottomBoundary = IMBoundaryCondition.robin(pBottom/robinLength,1);
end
robinEVP = IMInternalModes(name="terrainEnergyRobinGeostrophic", ...
    formulation="F",N2=N2,zDomain=[-D 0], ...
    p=p,q=@(zValue,~) zeros(size(zValue)), ...
    r=@(zValue,~) ones(size(zValue))/g, ...
    surfaceBoundary=IMBoundaryCondition.neumann(), ...
    bottomBoundary=bottomBoundary,f0=f,g=g);
robinBasis = solver.solveEVP(robinEVP,nModes=numberOfModes);
robinF = evaluateBasis(robinBasis,"F",z,normalization="unity");
robinFz = evaluateBasis(robinBasis,"uz",z,normalization="unity");
robinEta = -(f./N2(z)).*robinFz;
robinEigenvalue = reshape(robinBasis.eigenvalues,1,[]);
robinModeNumber = reshape(robinBasis.modeNumber,1,[]);
numberOfNegativeRobinModes = nnz(robinEigenvalue < 0);
expectedNegativeRobinModes = double(isfinite(robinLength) && robinLength < 0);

N2Bottom = N2(-D);
bottomProblem = IMSurfaceGeostrophicModes.atWavenumber( ...
    N2=N2,zDomain=[-D 0],f0=f,g=g,k=kappa, ...
    g0=Inf,gd=N2Bottom,surfaceAnomaly="noFreeSurface");
bottomBasis = solver.solveSurfaceGeostrophicModes(bottomProblem);
bottomPsi = evaluateBasis(bottomBasis,"F",z);
bottomEta = (f/g)*evaluateBasis(bottomBasis,"G",z);
bottomEndpointEta = (f/g)*evaluateBasis( ...
    bottomBasis,"G",[-D;0]);
normalization = bottomEndpointEta(1);
if abs(normalization) <= 100*eps
    error("WVTerrainEnergyGalerkin:BottomInversionNormalizationFailure", ...
        "The zero-APV bottom inversion has a vanishing bottom displacement.")
end
bottomPsi = bottomPsi/normalization;
bottomEta = bottomEta/normalization;
bottomEndpointEta = bottomEndpointEta/normalization;

waveEndpointG = evaluateBasis(waveBasis,"G",[-D;0]);
waveEndpointDefect = max(abs(waveEndpointG),[],"all") ...
    /max(max(abs(waveG),[],"all"),realmin);
robinEndpointF = evaluateBasis( ...
    robinBasis,"F",[-D;0],normalization="unity");
robinEndpointFz = evaluateBasis( ...
    robinBasis,"uz",[-D;0],normalization="unity");
robinScale = max(abs(robinFz),[],"all") ...
    +max(abs(robinF),[],"all")/D;
surfaceRobinDefect = max(abs(robinEndpointFz(2,:))) ...
    /max(robinScale,realmin);
if isinf(robinLength)
    bottomRobinResidual = robinEndpointFz(1,:);
else
    bottomRobinResidual = robinEndpointFz(1,:) ...
        -robinEndpointF(1,:)/robinLength;
end
bottomRobinScale = robinScale;
bottomRobinDefect = max(abs(bottomRobinResidual)) ...
    /max(bottomRobinScale,realmin);
bottomEtaXi = -(kappa^2/f)*bottomPsi;
bottomAPV = -kappa^2*bottomPsi-f*bottomEtaXi;
bottomAPVDefect = norm(bottomAPV) ...
    /max(kappa^2*norm(bottomPsi)+abs(f)*norm(bottomEtaXi),realmin);

family = struct;
family.kappa = kappa;
family.numberOfModes = numberOfModes;
family.robinLength = robinLength;
family.solverOrder = solverOrder;
family.wave = struct("F",waveF,"G",waveG,"h",waveH, ...
    "omega",waveOmega,"modeNumber",waveBasis.modeNumber, ...
    "eigenvalue",waveBasis.eigenvalues);
family.robin = struct("F",robinF,"eta",robinEta, ...
    "eigenvalue",robinEigenvalue,"modeNumber",robinModeNumber, ...
    "numberOfNegativeModes",numberOfNegativeRobinModes, ...
    "expectedNegativeModes",expectedNegativeRobinModes);
family.bottom = struct("psi",bottomPsi,"eta",bottomEta);
family.diagnostics = struct( ...
    "waveEndpointDefect",waveEndpointDefect, ...
    "surfaceRobinDefect",surfaceRobinDefect, ...
    "bottomRobinDefect",bottomRobinDefect, ...
    "negativeRobinModeCountDefect", ...
    abs(numberOfNegativeRobinModes-expectedNegativeRobinModes), ...
    "bottomAPVDefect",bottomAPVDefect, ...
    "bottomValueDefect",abs(bottomEndpointEta(1)-1), ...
    "bottomSurfaceValueDefect",abs(bottomEndpointEta(2)));
end

function values = evaluateBasis(basis,variable,z,options)
arguments
    basis
    variable (1,1) string {mustBeMember(variable,["F","G","uz"])}
    z (:,1) double
    options.normalization = []
end

[zDescending,order] = sort(z,"descend");
inverseOrder = zeros(size(order));
inverseOrder(order) = 1:numel(order);
if isempty(options.normalization)
    switch variable
        case "F"
            values = basis.F(zDescending);
        case "G"
            values = basis.G(zDescending);
        case "uz"
            values = basis.uz(zDescending);
    end
else
    switch variable
        case "F"
            values = basis.F( ...
                zDescending,normalization=options.normalization);
        case "G"
            values = basis.G( ...
                zDescending,normalization=options.normalization);
        case "uz"
            values = basis.uz( ...
                zDescending,normalization=options.normalization);
    end
end
values = values(inverseOrder,:);
end
