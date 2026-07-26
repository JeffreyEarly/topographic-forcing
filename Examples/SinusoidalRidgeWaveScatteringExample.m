function result = SinusoidalRidgeWaveScatteringExample(options)
% Scatter one internal wave from a sinusoidal ridge.
%
% The example evolves one incident wave with the autonomous first-order
% bottom scattering forcing, writes standard WaveVortexModel output, and
% reconstructs the bottom displacement from the saved accepted states.
%
% - Topic: Examples
% - Declaration: result = SinusoidalRidgeWaveScatteringExample(options)
% - Parameter options.resolution: transform resolution $$[N_x,N_y,N_z]$$
% - Parameter options.topographicHeight: sinusoidal-ridge amplitude in meters
% - Parameter options.numberOfWavePeriods: duration in incident-wave periods
% - Parameter options.numberOfOutputTimes: number of standard output samples
% - Parameter options.relativeTolerance: adaptive integration tolerance
% - Parameter options.shouldAntialias: antialias setting for the transform
% - Parameter options.outputPath: optional retained WaveVortexModel NetCDF path
% - Parameter options.shouldMakeFigures: whether to create diagnostic figures
% - Returns result: configuration, spectral, energetic, bottom, and figure diagnostics
arguments
    options.resolution (1,3) double {mustBeInteger,mustBePositive} = [32 4 17]
    options.topographicHeight (1,1) double {mustBePositive,mustBeFinite} = 50
    options.numberOfWavePeriods (1,1) double {mustBePositive,mustBeFinite} = 4
    options.numberOfOutputTimes (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.numberOfOutputTimes,3)} = 129
    options.relativeTolerance (1,1) double {mustBePositive,mustBeFinite} = 1e-9
    options.shouldAntialias (1,1) logical = true
    options.outputPath (1,1) string = ""
    options.shouldMakeFigures (1,1) logical = true
end

domainSize = [20e3 20e3 2e3];
N2 = 2e-5;
latitude = 45;
kMode = 2;
lMode = 0;
j = 1;
terrainMode = 1;
wvt = WVTransformBoussinesq(domainSize,options.resolution,N2=@(z)N2*ones(size(z)),latitude=latitude,shouldAntialias=options.shouldAntialias);
wvt.t0 = 0;
wvt.t = 0;
[x,~] = ndgrid(wvt.x,wvt.y);
topographicHeight = options.topographicHeight*cos(2*pi*terrainMode*x/wvt.Lx);
wvt.initWithWaveModes(kMode=kMode,lMode=lMode,j=j,phi=0,u=0.01,sign=1);
wvt.A0(:) = 0;
incidentIndex = modeIndex(wvt,kMode,lMode,j,"plus");
incidentFrequency = wvt.Omega(incidentIndex);
wavePeriod = 2*pi/incidentFrequency;
finalTime = options.numberOfWavePeriods*wavePeriod;

forcingName = "sinusoidal-ridge wave scattering";
forcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=topographicHeight,name=forcingName);
wvt.removeAllForcing();
wvt.addForcing(forcing);
model = WVModel(wvt);
model.setupIntegrator(integratorType="adaptive",absTolerance=options.relativeTolerance,relTolerance=options.relativeTolerance);

isTemporaryOutput = strlength(options.outputPath) == 0;
if isTemporaryOutput
    outputPath = string(tempname)+".nc";
else
    outputPath = options.outputPath;
end
outputCleanup = onCleanup(@()deleteTemporaryOutput(outputPath,isTemporaryOutput));
outputInterval = finalTime/(options.numberOfOutputTimes-1);
model.createNetCDFFileForModelOutput(outputPath,outputInterval=outputInterval,shouldOverwriteExisting=true);
warningState = warning;
warningCleanup = onCleanup(@()warning(warningState));
warning("off","all")
model.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
model.closeNetCDFFile();
clear warningCleanup

bottomDiagnostics = WVBottomWaveScatteringForcing.bottomDisplacementFromFile(outputPath,forcingName=forcingName);
[savedDiagnostics,wvtSaved] = savedWaveDiagnostics(outputPath,forcingName,topographicHeight);
time = bottomDiagnostics.time;
[horizontalMode,energyByHorizontalMode] = horizontalEnergySpectrum(wvtSaved,savedDiagnostics.finalAp,savedDiagnostics.finalAm);
sidebandModes = [kMode-terrainMode kMode+terrainMode];
sidebandEnergy = zeros(numel(time),numel(sidebandModes));
for iSideband = 1:numel(sidebandModes)
    mask = abs(wvtSaved.K-2*pi*sidebandModes(iSideband)/wvtSaved.Lx) <= 100*eps(max(1,abs(2*pi*sidebandModes(iSideband)/wvtSaved.Lx))) & abs(wvtSaved.L) <= eps;
    sidebandEnergy(:,iSideband) = sum(savedDiagnostics.modalEnergy(:,mask),2);
end

figureHandles = gobjects(0);
if options.shouldMakeFigures
    figureHandles = createFigures(wvtSaved,time,wavePeriod,topographicHeight,savedDiagnostics,bottomDiagnostics,horizontalMode,energyByHorizontalMode,sidebandModes,sidebandEnergy);
end

configuration = struct(domainSize=domainSize,resolution=options.resolution,N2=N2,latitude=latitude,kMode=kMode,lMode=lMode,j=j,terrainMode=terrainMode,topographicHeight=options.topographicHeight,wavePeriod=wavePeriod,finalTime=finalTime,relativeTolerance=options.relativeTolerance,shouldAntialias=options.shouldAntialias);
retainedOutputPath = "";
if ~isTemporaryOutput
    retainedOutputPath = outputPath;
end
result = struct(configuration=configuration,time=time,topographicHeight=topographicHeight,waveEnergy=savedDiagnostics.waveEnergy,firstOrderEnergy=savedDiagnostics.firstOrderEnergy,balancedEnergy=savedDiagnostics.balancedEnergy,qgpvNorm=savedDiagnostics.qgpvNorm,horizontalMode=horizontalMode,energyByHorizontalMode=energyByHorizontalMode,sidebandModes=sidebandModes,sidebandEnergy=sidebandEnergy,bottomVelocity=bottomDiagnostics.bottomVelocity,bottomDisplacement=bottomDiagnostics.bottomDisplacement,bottomQuadratureMethod=bottomDiagnostics.quadratureMethod,outputPath=retainedOutputPath,figureHandles=figureHandles);
clear outputCleanup
end

function [diagnostics,wvt] = savedWaveDiagnostics(path,forcingName,topographicHeight)
[wvt,ncfile] = WVTransformBoussinesq.waveVortexTransformFromFile(char(path),iTime=1,shouldReadOnly=true);
fileCleanup = onCleanup(@()ncfile.close());
group = ncfile.groupWithName("wave-vortex");
time = reshape(group.readVariables("t"),[],1);
numberOfTimes = numel(time);
modalEnergy = zeros(numberOfTimes,numel(wvt.Ap));
waveEnergy = zeros(numberOfTimes,1);
firstOrderEnergy = zeros(numberOfTimes,1);
balancedEnergy = zeros(numberOfTimes,1);
qgpvNorm = zeros(numberOfTimes,1);
forcing = wvt.forcingWithName(char(forcingName));
for iTime = 1:numberOfTimes
    [wvt.Ap,wvt.Am,wvt.A0] = group.readVariablesAtIndexAlongDimension("t",iTime,"Ap","Am","A0");
    wvt.t = time(iTime);
    modalEnergy(iTime,:) = reshape(wvt.Apm_TE_factor.*(abs(wvt.Ap).^2+abs(wvt.Am).^2),1,[]);
    waveEnergy(iTime) = sum(modalEnergy(iTime,:));
    balancedEnergy(iTime) = sum(wvt.A0_TE_factor(:).*abs(wvt.A0(:)).^2);
    [~,bottomFields] = forcing.bottomVelocityFromWaveState(wvt);
    firstOrderEnergy(iTime) = waveEnergy(iTime)-0.5*mean(topographicHeight.*(bottomFields.u.^2+bottomFields.v.^2),"all");
    qgpvNorm(iTime) = norm(wvt.qgpv(:));
end
finalAp = wvt.Ap;
finalAm = wvt.Am;
diagnostics = struct(modalEnergy=modalEnergy,waveEnergy=waveEnergy,firstOrderEnergy=firstOrderEnergy,balancedEnergy=balancedEnergy,qgpvNorm=qgpvNorm,finalAp=finalAp,finalAm=finalAm);
clear fileCleanup
end

function [horizontalMode,energyByHorizontalMode] = horizontalEnergySpectrum(wvt,Ap,Am)
modalEnergy = wvt.Apm_TE_factor.*(abs(Ap).^2+abs(Am).^2);
retained = wvt.waveComponent.maskAp | wvt.waveComponent.maskAm;
modeNumber = round(wvt.K/wvt.dk);
horizontalMode = unique(modeNumber(retained & abs(wvt.L) <= eps));
energyByHorizontalMode = zeros(size(horizontalMode));
for iMode = 1:numel(horizontalMode)
    energyByHorizontalMode(iMode) = sum(modalEnergy(retained & modeNumber == horizontalMode(iMode) & abs(wvt.L) <= eps));
end
end

function figureHandles = createFigures(wvt,time,wavePeriod,topographicHeight,diagnostics,bottomDiagnostics,horizontalMode,energyByHorizontalMode,sidebandModes,sidebandEnergy)
figureHandles(1) = figure(Name="Sinusoidal-ridge wave scattering",Color="w",Position=[100 100 1200 780]);
layout = tiledlayout(figureHandles(1),2,2,TileSpacing="compact",Padding="compact");

axesHandle = nexttile(layout,1);
plot(axesHandle,wvt.x/1e3,topographicHeight(:,1),LineWidth=1.6)
grid(axesHandle,"on")
xlabel(axesHandle,"x (km)")
ylabel(axesHandle,"h (m)")
title(axesHandle,"Sinusoidal ridge")

axesHandle = nexttile(layout,2);
stem(axesHandle,horizontalMode,energyByHorizontalMode,"filled")
grid(axesHandle,"on")
xlabel(axesHandle,"zonal mode k")
ylabel(axesHandle,"final wave energy")
title(axesHandle,"Incident mode and scattered sidebands")

axesHandle = nexttile(layout,3);
plot(axesHandle,time/wavePeriod,sidebandEnergy,LineWidth=1.4)
grid(axesHandle,"on")
xlabel(axesHandle,"time / incident period")
ylabel(axesHandle,"sideband energy")
legend(axesHandle,compose("k=%d",sidebandModes),Location="best")
title(axesHandle,"First terrain sidebands")

axesHandle = nexttile(layout,4);
plot(axesHandle,time/wavePeriod,diagnostics.waveEnergy-diagnostics.waveEnergy(1),LineWidth=1.4)
hold(axesHandle,"on")
plot(axesHandle,time/wavePeriod,diagnostics.firstOrderEnergy-diagnostics.firstOrderEnergy(1),"--",LineWidth=1.4)
grid(axesHandle,"on")
xlabel(axesHandle,"time / incident period")
ylabel(axesHandle,"energy change")
legend(axesHandle,"flat wave energy","first-order physical energy",Location="best")
title(axesHandle,"Autonomous energy exchange")
title(layout,"First-order scattering from a sinusoidal ridge")

figureHandles(2) = figure(Name="Bottom scattering diagnostics",Color="w",Position=[150 150 1200 480]);
layout = tiledlayout(figureHandles(2),1,2,TileSpacing="compact",Padding="compact");
axesHandle = nexttile(layout,1);
imagesc(axesHandle,wvt.x/1e3,wvt.y/1e3,bottomDiagnostics.bottomVelocity(:,:,end).')
axis(axesHandle,"xy","image")
xlabel(axesHandle,"x (km)")
ylabel(axesHandle,"y (km)")
title(axesHandle,"Final bottom velocity g_b")
colorbar(axesHandle)
axesHandle = nexttile(layout,2);
imagesc(axesHandle,wvt.x/1e3,wvt.y/1e3,bottomDiagnostics.bottomDisplacement(:,:,end).')
axis(axesHandle,"xy","image")
xlabel(axesHandle,"x (km)")
ylabel(axesHandle,"y (km)")
title(axesHandle,"Postprocessed bottom displacement \eta_d")
colorbar(axesHandle)
colormap(figureHandles(2),"turbo")
end

function index = modeIndex(wvt,kMode,lMode,j,branch)
k = 2*pi*kMode/wvt.Lx;
l = 2*pi*lMode/wvt.Ly;
tolerance = 100*eps(max([1 abs(k) abs(l)]));
mask = wvt.J == j & abs(wvt.K-k) <= tolerance & abs(wvt.L-l) <= tolerance;
if branch == "plus"
    mask = mask & wvt.waveComponent.maskAp;
else
    mask = mask & wvt.waveComponent.maskAm;
end
index = find(mask,1);
if isempty(index)
    error("SinusoidalRidgeWaveScatteringExample:MissingWaveMode", "The requested incident mode is not retained by the transform.")
end
end

function deleteTemporaryOutput(path,shouldDelete)
if shouldDelete && isfile(path)
    delete(path)
end
end
