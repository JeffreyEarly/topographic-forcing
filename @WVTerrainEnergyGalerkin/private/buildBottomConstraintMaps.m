function [B,R,diagnostics] = buildBottomConstraintMaps(problem)
% Build resolved bottom-value and bottom-evolution maps.

forms = problem.finiteTerrainForms;
wvt = problem.originatingTransform;
Nx = forms.oversampledSize(1);
Ny = forms.oversampledSize(2);
nXY = Nx*Ny;
nState = height(problem.stateLayout);
nHorizontal = height(problem.horizontalLayout);
[x,y] = ndgrid((0:Nx-1)'*wvt.Lx/Nx,(0:Ny-1)'*wvt.Ly/Ny);
phase = exp(1i*(x(:)*problem.horizontalLayout.k.'+y(:)*problem.horizontalLayout.l.'));

B = zeros(nHorizontal,nState);
uHatBottom = zeros(nXY,nState);
vHatBottom = zeros(nXY,nState);
for iK = 1:nHorizontal
    rows = find(problem.stateLayout.horizontalIndex == iK);
    block = problem.basisBlocks{iK};
    uHatBottom(:,rows) = phase(:,iK)*block.uHat(1,:);
    vHatBottom(:,rows) = phase(:,iK)*block.vHat(1,:);
    bottomRow = rows(problem.stateLayout.component(rows) == "etaB");
    if numel(bottomRow) ~= 1
        error("WVTerrainEnergyGalerkin:InvalidBottomCoordinate", ...
            "Every retained horizontal mode must contain exactly one bottom-displacement coordinate.")
    end
    B(iK,bottomRow) = 1;
end

gradLnGammaX = forms.oversampledGradLnGammaX(:);
gradLnGammaY = forms.oversampledGradLnGammaY(:);
bottomTendencyGrid = -wvt.Lz*(gradLnGammaX.*uHatBottom+gradLnGammaY.*vHatBottom);
R = phase'*bottomTendencyGrid/nXY;
resolvedBottomTendencyGrid = phase*R;
projectionResidual = norm(bottomTendencyGrid-resolvedBottomTendencyGrid,"fro")/max(norm(bottomTendencyGrid,"fro"),realmin);
phaseGramDefect = norm(phase'*phase/nXY-eye(nHorizontal),"fro")/sqrt(nHorizontal);

diagnostics = struct;
diagnostics.oversampledBottomTendencyMatrix = bottomTendencyGrid;
diagnostics.projectionResidual = projectionResidual;
diagnostics.phaseGramDefect = phaseGramDefect;
diagnostics.oversampledSize = [Nx Ny];
end
