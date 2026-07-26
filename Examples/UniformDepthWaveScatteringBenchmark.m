function result = UniformDepthWaveScatteringBenchmark(options)
% Validate autonomous scattering with the uniform-depth frequency shift.
%
% A spatially uniform topographic height changes the physical depth from
% $$D$$ to $$D-h_0$$. The first-order scattering operator is compared with
% the exact constant-stratification dispersion relation and with an
% adaptive WaveVortexModel integration.
%
% - Topic: Examples
% - Declaration: result = UniformDepthWaveScatteringBenchmark(options)
% - Parameter options.topographicHeights: positive uniform heights in meters
% - Parameter options.resolution: transform resolution $$[N_x,N_y,N_z]$$
% - Parameter options.relativeTolerance: adaptive integration tolerance
% - Parameter options.shouldMakeFigures: whether to create the diagnostic figure
% - Returns result: frequency, convergence, integration, and figure diagnostics
arguments
    options.topographicHeights (1,:) double {mustBePositive,mustBeFinite} = [50 25 12.5]
    options.resolution (1,3) double {mustBeInteger,mustBePositive} = [8 4 9]
    options.relativeTolerance (1,1) double {mustBePositive,mustBeFinite} = 1e-9
    options.shouldMakeFigures (1,1) logical = true
end

domainSize = [20e3 20e3 2e3];
N2 = 2e-5;
latitude = 45;
kMode = 2;
lMode = 0;
j = 1;
wvt = WVTransformBoussinesq(domainSize,options.resolution,N2=@(z)N2*ones(size(z)),latitude=latitude,shouldAntialias=false);
wvt.t0 = 0;
wvt.t = 0;
[plusIndices,minusIndices] = horizontalWaveBlock(wvt,kMode,lMode);
omega0 = wvt.Omega(plusIndices(wvt.J(plusIndices) == j));
k = 2*pi*kMode/wvt.Lx;
m0 = j*pi/wvt.Lz;
dOmegaDDepth = m0^2*k^2*(N2-wvt.f^2)/(omega0*wvt.Lz*(k^2+m0^2)^2);

topographicHeights = sort(options.topographicHeights,"descend");
operatorFrequency = zeros(size(topographicHeights));
exactFrequency = zeros(size(topographicHeights));
firstOrderFrequency = omega0-topographicHeights*dOmegaDDepth;
generators = cell(size(topographicHeights));
for iHeight = 1:numel(topographicHeights)
    [generators{iHeight},operatorFrequency(iHeight)] = physicalGenerator(wvt,topographicHeights(iHeight),plusIndices,minusIndices,j);
    physicalDepth = wvt.Lz-topographicHeights(iHeight);
    verticalWavenumber = j*pi/physicalDepth;
    exactFrequency(iHeight) = sqrt((N2*k^2+wvt.f^2*verticalWavenumber^2)/(k^2+verticalWavenumber^2));
end
frequencyError = abs(operatorFrequency-exactFrequency);
convergenceFit = polyfit(log(topographicHeights),log(frequencyError),1);

generator = generators{1};
[eigenvectors,eigenvalues] = eig(generator,"vector");
[~,iEigenmode] = min(abs(imag(eigenvalues)-operatorFrequency(1)));
initialPhysicalCoefficients = eigenvectors(:,iEigenmode);
initialPhysicalCoefficients = 0.01*initialPhysicalCoefficients/max(abs(initialPhysicalCoefficients));
wvt.Ap(:) = 0;
wvt.Am(:) = 0;
wvt.A0(:) = 0;
wvt.Ap(plusIndices) = initialPhysicalCoefficients(1:numel(plusIndices));
wvt.Am(minusIndices) = initialPhysicalCoefficients(numel(plusIndices)+1:end);
initialAp = wvt.Ap;
initialAm = wvt.Am;
forcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=topographicHeights(1)*ones(wvt.Nx,wvt.Ny));
wvt.removeAllForcing();
wvt.addForcing(forcing);
model = WVModel(wvt);
model.setupIntegrator(integratorType="adaptive",absTolerance=options.relativeTolerance,relTolerance=options.relativeTolerance);
finalTime = 4*2*pi/operatorFrequency(1);
warningState = warning;
warningCleanup = onCleanup(@()warning(warningState));
warning("off","all")
model.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
clear warningCleanup

exactPhysicalCoefficients = expm(generator*finalTime)*initialPhysicalCoefficients;
exactAp = complex(zeros(size(wvt.Ap)));
exactAm = complex(zeros(size(wvt.Am)));
exactAp(plusIndices) = exactPhysicalCoefficients(1:numel(plusIndices)).*exp(-wvt.iOmega(plusIndices)*finalTime);
exactAm(minusIndices) = exactPhysicalCoefficients(numel(plusIndices)+1:end).*exp(wvt.iOmega(minusIndices)*finalTime);
coefficientError = sqrt(sum(wvt.Apm_TE_factor(:).*(abs(wvt.Ap(:)-exactAp(:)).^2+abs(wvt.Am(:)-exactAm(:)).^2))/sum(wvt.Apm_TE_factor(:).*(abs(exactAp(:)).^2+abs(exactAm(:)).^2)));
balancedEnergy = sum(wvt.A0_TE_factor(:).*abs(wvt.A0(:)).^2);

wvt.Ap = initialAp;
wvt.Am = initialAm;
wvt.A0(:) = 0;
wvt.t = 0;
[Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
[Fu,Fv,~,Feta] = wvt.transformWaveVortexToUVWEta(Fp,Fm,F0,wvt.t);
qgpvTendency = wvt.diffX(Fv)-wvt.diffY(Fu)-wvt.f*wvt.diffZG(Feta);

figureHandles = gobjects(0);
if options.shouldMakeFigures
    figureHandles(1) = figure(Name="Uniform-depth wave scattering",Color="w",Position=[100 100 1150 480]);
    layout = tiledlayout(figureHandles(1),1,2,TileSpacing="compact",Padding="compact");
    axesHandle = nexttile(layout,1);
    plot(axesHandle,topographicHeights,1e6*(operatorFrequency-omega0),"o-",LineWidth=1.5,MarkerFaceColor=[0.1 0.4 0.8])
    hold(axesHandle,"on")
    plot(axesHandle,topographicHeights,1e6*(exactFrequency-omega0),"--",LineWidth=1.5)
    plot(axesHandle,topographicHeights,1e6*(firstOrderFrequency-omega0),":",LineWidth=1.5)
    grid(axesHandle,"on")
    xlabel(axesHandle,"uniform topographic height h_0 (m)")
    ylabel(axesHandle,"\omega-\omega_D (10^{-6} s^{-1})")
    legend(axesHandle,"scattering operator","exact depth D-h_0","linear derivative",Location="best")
    title(axesHandle,"Frequency correction")

    axesHandle = nexttile(layout,2);
    loglog(axesHandle,topographicHeights,frequencyError,"o-",LineWidth=1.5,MarkerFaceColor=[0.85 0.33 0.10])
    hold(axesHandle,"on")
    referenceError = frequencyError(end)*(topographicHeights/topographicHeights(end)).^2;
    loglog(axesHandle,topographicHeights,referenceError,"--",LineWidth=1.2)
    grid(axesHandle,"on")
    xlabel(axesHandle,"uniform topographic height h_0 (m)")
    ylabel(axesHandle,"frequency error (s^{-1})")
    legend(axesHandle,"measured error","O(h_0^2)",Location="best")
    title(axesHandle,sprintf("Observed exponent %.3f",convergenceFit(1)))
    title(layout,"Uniform-depth first-order scattering benchmark")
end

configuration = struct(domainSize=domainSize,resolution=options.resolution,N2=N2,latitude=latitude,kMode=kMode,lMode=lMode,j=j,finalTime=finalTime);
result = struct(configuration=configuration,topographicHeights=topographicHeights,omega0=omega0,operatorFrequency=operatorFrequency,exactFrequency=exactFrequency,firstOrderFrequency=firstOrderFrequency,frequencyError=frequencyError,convergenceExponent=convergenceFit(1),coefficientError=coefficientError,balancedEnergy=balancedEnergy,qgpvTendencyNorm=norm(qgpvTendency(:)),figureHandles=figureHandles);
disp(table(topographicHeights(:),operatorFrequency(:),exactFrequency(:),frequencyError(:),VariableNames=["Height","OperatorFrequency","ExactFrequency","Error"]))
end

function [generator,frequency] = physicalGenerator(wvt,topographicHeight,plusIndices,minusIndices,j)
forcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=topographicHeight*ones(wvt.Nx,wvt.Ny));
numberOfPlusModes = numel(plusIndices);
numberOfModes = numberOfPlusModes+numel(minusIndices);
coupling = complex(zeros(numberOfModes));
for iInput = 1:numberOfModes
    wvt.Ap(:) = 0;
    wvt.Am(:) = 0;
    wvt.A0(:) = 0;
    wvt.t = 0;
    if iInput <= numberOfPlusModes
        wvt.Ap(plusIndices(iInput)) = 1;
    else
        wvt.Am(minusIndices(iInput-numberOfPlusModes)) = 1;
    end
    [Fp,Fm] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
    coupling(:,iInput) = [Fp(plusIndices); Fm(minusIndices)];
end
signedFrequency = [wvt.Omega(plusIndices); -wvt.Omega(minusIndices)];
generator = diag(1i*signedFrequency)+coupling;
eigenvalues = eig(generator);
positiveFrequency = sort(imag(eigenvalues(imag(eigenvalues) > 0)),"descend");
frequency = positiveFrequency(j);
end

function [plusIndices,minusIndices] = horizontalWaveBlock(wvt,kMode,lMode)
k = 2*pi*kMode/wvt.Lx;
l = 2*pi*lMode/wvt.Ly;
tolerance = 100*eps(max([1 abs(k) abs(l)]));
horizontalMask = abs(wvt.K-k) <= tolerance & abs(wvt.L-l) <= tolerance;
plusIndices = find(horizontalMask & wvt.waveComponent.maskAp);
minusIndices = find(horizontalMask & wvt.waveComponent.maskAm);
if isempty(plusIndices) || isempty(minusIndices)
    error("UniformDepthWaveScatteringBenchmark:MissingWaveBlock", "The requested horizontal wave block is not retained by the transform.")
end
end
