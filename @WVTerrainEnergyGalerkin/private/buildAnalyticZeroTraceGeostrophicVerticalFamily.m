function family = buildAnalyticZeroTraceGeostrophicVerticalFamily( ...
    problem,kappa,numberOfModes,z,weight)
% Build exact constant-N mixed Dirichlet--Neumann coordinates.
arguments
    problem (1,1) WVTerrainEnergyGalerkin
    kappa (1,1) double {mustBePositive}
    numberOfModes (1,1) double {mustBeInteger,mustBePositive}
    z (:,1) double {mustBeReal,mustBeFinite}
    weight (:,1) double {mustBeReal,mustBeFinite,mustBePositive}
end

wvt = problem.originatingTransform;
D = wvt.Lz;
f = wvt.f;
g = wvt.g;
N2Values = wvt.N2Function(z);
if isscalar(N2Values)
    N2Values = repmat(N2Values,size(z));
end
N2Value = mean(N2Values);
constantDefect = max(abs(N2Values-N2Value),[],"all") ...
    /max(abs(N2Value),realmin);
if constantDefect > 1e-12
    error("WVTerrainEnergyGalerkin:AnalyticGlobalTangentRequiresConstantN", ...
        "The analytic mixed-boundary family requires constant N2.")
end

modeNumber = 1:numberOfModes;
verticalWavenumber = (modeNumber-0.5)*pi/D;
mu = kappa^2+(f^2/N2Value)*verticalWavenumber.^2;
normalization = sqrt(2./(D*mu));
argument = (z+D)*verticalWavenumber;
F = sin(argument).*normalization;
Fz = cos(argument).*(verticalWavenumber.*normalization);
Fzz = -F.*verticalWavenumber.^2;
eta = -(f/N2Value)*Fz;
eigenvalue = (g/N2Value)*verticalWavenumber.^2;

energy = F'*(weight.*(kappa^2*F)) ...
    +Fz'*(weight.*((f^2/N2Value)*Fz));
q = -F.*mu;
enstrophy = q'*(weight.*q);
identity = eye(numberOfModes);
endpointScale = max(abs(F),[],"all") ...
    +D*max(abs(Fz),[],"all");
bottomValues = sin(zeros(1,numberOfModes)).*normalization;
surfaceDerivatives = cos(D*verticalWavenumber) ...
    .*(verticalWavenumber.*normalization);
odeResidual = Fzz+F.*verticalWavenumber.^2;

family = struct;
family.F = F;
family.Fz = Fz;
family.eta = eta;
family.eigenvalue = eigenvalue;
family.mu = mu;
family.modeNumber = modeNumber;
family.diagnostics = struct( ...
    "bottomEndpointDefect",max(abs(bottomValues),[],"all") ...
        /max(endpointScale,realmin), ...
    "surfaceEndpointDefect",D*max(abs(surfaceDerivatives),[],"all") ...
        /max(endpointScale,realmin), ...
    "differentialEquationDefect",norm(odeResidual,"fro") ...
        /max(norm(Fzz,"fro") ...
        +norm(F.*verticalWavenumber.^2,"fro"),realmin), ...
    "normalizationDefect",norm(energy-identity,"fro") ...
        /max(norm(identity,"fro"),realmin), ...
    "energyDiagonalizationDefect",offDiagonalDefect(energy), ...
    "enstrophyDiagonalizationDefect",offDiagonalDefect(enstrophy), ...
    "constantStratificationDefect",constantDefect, ...
    "constantStratificationAnalyticDefect",0, ...
    "solverOrder",NaN);
end

function value = offDiagonalDefect(matrix)
diagonal = diag(diag(matrix));
value = norm(matrix-diagonal,"fro")/max(norm(matrix,"fro"),realmin);
end
