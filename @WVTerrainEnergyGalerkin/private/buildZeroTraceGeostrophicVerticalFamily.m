function family = buildZeroTraceGeostrophicVerticalFamily( ...
    problem,kappa,numberOfModes,solverOrder,z,weight)
% Build mixed Dirichlet--Neumann geostrophic vertical coordinates.
arguments
    problem (1,1) WVTerrainEnergyGalerkin
    kappa (1,1) double {mustBeNonnegative}
    numberOfModes (1,1) double {mustBeInteger,mustBePositive}
    solverOrder (1,1) double {mustBeInteger,mustBePositive}
    z (:,1) double {mustBeReal,mustBeFinite}
    weight (:,1) double {mustBeReal,mustBeFinite,mustBePositive}
end

wvt = problem.originatingTransform;
D = wvt.Lz;
f = wvt.f;
g = wvt.g;
N2 = wvt.N2Function;
solver = IMSolverSpectral(nEVP=solverOrder);
evp = IMInternalModes(name="terrainEnergyZeroTraceGeostrophic", ...
    formulation="F",N2=N2,zDomain=[-D 0], ...
    p=@(zValue,~)1./N2(zValue), ...
    q=@(zValue,~)zeros(size(zValue)), ...
    r=@(zValue,~)ones(size(zValue))/g, ...
    surfaceBoundary=IMBoundaryCondition.neumann(), ...
    bottomBoundary=IMBoundaryCondition.dirichlet(),f0=f,g=g);
basis = solver.solveEVP(evp,nModes=numberOfModes);
F = evaluateBasis(basis,"F",z);
Fz = evaluateBasis(basis,"uz",z);
eta = -(f./N2(z)).*Fz;
eigenvalue = reshape(basis.eigenvalues,1,[]);
mu = kappa^2+(f^2/g)*eigenvalue;

endpointF = evaluateBasis(basis,"F",[-D;0]);
endpointFz = evaluateBasis(basis,"uz",[-D;0]);
endpointScale = max(abs(F),[],"all") ...
    +D*max(abs(Fz),[],"all");
bottomDefect = max(abs(endpointF(1,:))) ...
    /max(endpointScale,realmin);
surfaceDefect = D*max(abs(endpointFz(2,:))) ...
    /max(endpointScale,realmin);

N2Values = N2(z);
if isscalar(N2Values)
    N2Values = repmat(N2Values,size(z));
end
energy = F'*(weight.*(kappa^2*F)) ...
    +Fz'*(weight.*((f^2./N2Values).*Fz));
q = -F.*mu;
enstrophy = q'*(weight.*q);
energyDiagonalization = offDiagonalDefect(energy);
enstrophyDiagonalization = offDiagonalDefect(enstrophy);

analyticDefect = NaN;
values = N2Values(:);
if max(values)-min(values) <=100*eps*max(values)
    analytic = zeros(size(F));
    for iMode = 1:numberOfModes
        verticalWavenumber = (iMode-0.5)*pi/D;
        analytic(:,iMode) = sin(verticalWavenumber*(z+D));
    end
    analyticDefect = subspaceDefect(F,analytic,weight);
end

family = struct;
family.F = F;
family.Fz = Fz;
family.eta = eta;
family.eigenvalue = eigenvalue;
family.mu = mu;
family.modeNumber = reshape(basis.modeNumber,1,[]);
family.diagnostics = struct("bottomEndpointDefect",bottomDefect, ...
    "surfaceEndpointDefect",surfaceDefect, ...
    "energyDiagonalizationDefect",energyDiagonalization, ...
    "enstrophyDiagonalizationDefect",enstrophyDiagonalization, ...
    "constantStratificationAnalyticDefect",analyticDefect, ...
    "solverOrder",solverOrder);
end

function values = evaluateBasis(basis,variable,z)
[zDescending,order] = sort(z,"descend");
inverseOrder = zeros(size(order));
inverseOrder(order) = 1:numel(order);
if variable == "F"
    values = basis.F(zDescending,normalization="unity");
else
    values = basis.uz(zDescending,normalization="unity");
end
values = values(inverseOrder,:);
end

function value = offDiagonalDefect(matrix)
diagonal = diag(diag(matrix));
value = norm(matrix-diagonal,"fro")/max(norm(matrix,"fro"),realmin);
end

function value = subspaceDefect(first,second,weight)
[first,~,~] = svd(sqrt(weight).*first,"econ");
[second,~,~] = svd(sqrt(weight).*second,"econ");
firstProjector = first*first';
secondProjector = second*second';
value = norm(firstProjector-secondProjector,"fro") ...
    /max(norm(firstProjector,"fro")+norm(secondProjector,"fro"),realmin);
end
