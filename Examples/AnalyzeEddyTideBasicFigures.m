function [figureHandles,summary] = AnalyzeEddyTideBasicFigures(eddyFile,noEddyFile,options)
% Create the basic paired figures for a completed eddy-tide experiment.
%
% The analysis compares wave and geostrophic energy, domain-maximum wave
% speed, final wave-energy spectra, and final surface vertical vorticity.
% Energy is obtained through `WVDiagnostics`; the other quantities are
% reconstructed from the saved wave-vortex coefficients.
%
% ```matlab
% [figures,summary] = AnalyzeEddyTideBasicFigures(eddyFile,noEddyFile);
% ```
%
% - Topic: Examples
% - Declaration: [figureHandles,summary] = AnalyzeEddyTideBasicFigures(eddyFile,noEddyFile,options)
% - Parameter eddyFile: completed initial-eddy simulation NetCDF path
% - Parameter noEddyFile: completed no-initial-eddy simulation NetCDF path
% - Parameter options.figureVisible: figure visibility, `"on"` or `"off"`
% - Parameter options.shouldExport: whether to export PNG figures and a MAT summary
% - Parameter options.exportDirectory: output directory, default beside `eddyFile`
% - Parameter options.exportPrefix: common exported-file prefix
% - Parameter options.exportResolution: exported PNG resolution in dots per inch
% - Parameter options.shouldOverwriteExisting: whether to replace existing analysis files
% - Returns figureHandles: energy, speed, spectrum, and vorticity figures
% - Returns summary: paired diagnostics, extrema, spectra, vorticity, and output paths
arguments (Input)
    eddyFile (1,1) string {mustBeFile}
    noEddyFile (1,1) string {mustBeFile}
    options.figureVisible (1,1) string {mustBeMember(options.figureVisible,["on" "off"])} = "on"
    options.shouldExport (1,1) logical = true
    options.exportDirectory (1,1) string = ""
    options.exportPrefix (1,1) string = "eddy-tide-basic"
    options.exportResolution (1,1) double {mustBeInteger,mustBePositive} = 300
    options.shouldOverwriteExisting (1,1) logical = false
end
arguments (Output)
    figureHandles (4,1) matlab.ui.Figure
    summary (1,1) struct
end

if exist("WVDiagnostics","class") ~= 8
    error("AnalyzeEddyTideBasicFigures:WVDiagnosticsNotFound", "WaveVortexModelDiagnostics is required. Install it with OceanKit or add its authoring repository to the MATLAB path.")
end

exportDirectory = options.exportDirectory;
if strlength(exportDirectory) == 0
    exportDirectory = string(fileparts(eddyFile));
end
figurePaths = [
    fullfile(exportDirectory,options.exportPrefix+"-energy.png")
    fullfile(exportDirectory,options.exportPrefix+"-maximum-wave-speed.png")
    fullfile(exportDirectory,options.exportPrefix+"-final-wave-energy-spectrum.png")
    fullfile(exportDirectory,options.exportPrefix+"-final-surface-vorticity.png")
    ];
summaryPath = fullfile(exportDirectory,options.exportPrefix+"-summary.mat");
if options.shouldExport
    validateExportPaths([figurePaths; summaryPath],options.shouldOverwriteExisting)
end

[figureHandles(1),energy] = AnalyzeEddyTideEnergy(eddyFile,noEddyFile,figureVisible=options.figureVisible,shouldExport=false);
eddy = diagnosticsForCase(eddyFile,"initial eddy");
noEddy = diagnosticsForCase(noEddyFile,"no initial eddy");
validatePairedCases(eddy,noEddy,energy)

figureHandles(2) = createSpeedFigure(eddy,noEddy,options.figureVisible);
figureHandles(3) = createSpectrumFigure(eddy,noEddy,options.figureVisible);
figureHandles(4) = createVorticityFigure(eddy,noEddy,options.figureVisible);

comparison = struct(waveVorticityColorLimit=max(abs([eddy.vorticity.waveZetaOverF(:); noEddy.vorticity.waveZetaOverF(:)])), ...
    geostrophicVorticityColorLimit=max(abs([eddy.vorticity.geostrophicZetaOverF(:); noEddy.vorticity.geostrophicZetaOverF(:)])));
summary = struct(files=[eddyFile; noEddyFile],energy=energy,eddy=eddy,noEddy=noEddy,comparison=comparison,figurePaths=strings(4,1),summaryPath="");
if options.shouldExport
    for iFigure = 1:numel(figureHandles)
        exportgraphics(figureHandles(iFigure),figurePaths(iFigure),Resolution=options.exportResolution)
    end
    summary.figurePaths = figurePaths;
    summary.summaryPath = summaryPath;
    save(summaryPath,"summary")
end
end

function result = diagnosticsForCase(file,label)
[wvt,ncfile] = WVTransform.waveVortexTransformFromFile(char(file),iTime=1,shouldReadOnly=true);
cleanup = onCleanup(@()ncfile.close());
time = ncfile.readVariables("wave-vortex/t");
time = time(:);
numberOfTimes = numel(time);
maximumHorizontalSpeed = zeros(numberOfTimes,1);
maximumThreeDimensionalSpeed = zeros(numberOfTimes,1);
maximumVerticalSpeed = zeros(numberOfTimes,1);
generationRate = zeros(numberOfTimes,1);
generation = wvt.forcingWithName("bottom wave generation");
zeroWave = complex(zeros(size(wvt.Ap)));
zeroVortex = complex(zeros(size(wvt.A0)));

for iTime = 1:numberOfTimes
    wvt.initFromNetCDFFile(ncfile,iTime=iTime);
    [Fp,Fm] = generation.addSpectralForcing(wvt,zeroWave,zeroWave,zeroVortex);
    generationRate(iTime) = 2*sum(wvt.Apm_TE_factor.*real(Fp.*conj(wvt.Ap)+Fm.*conj(wvt.Am)),"all");

    waveU = wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*wvt.Apt+wvt.UAm.*wvt.Amt);
    waveV = wvt.transformToSpatialDomainWithF(Apm=wvt.VAp.*wvt.Apt+wvt.VAm.*wvt.Amt);
    waveW = wvt.transformToSpatialDomainWithG(Apm=wvt.WAp.*wvt.Apt+wvt.WAm.*wvt.Amt);
    horizontalSpeedSquared = waveU.^2+waveV.^2;
    maximumHorizontalSpeed(iTime) = sqrt(max(horizontalSpeedSquared,[],"all"));
    maximumThreeDimensionalSpeed(iTime) = sqrt(max(horizontalSpeedSquared+waveW.^2,[],"all"));
    maximumVerticalSpeed(iTime) = max(abs(waveW),[],"all");
end

spectrum = finalWaveSpectrum(wvt);
vorticity = finalSurfaceVorticity(wvt);
[maximumSpeed,maximumSpeedIndex] = max(maximumThreeDimensionalSpeed);
[peakHorizontalFraction,peakHorizontalIndex] = max(spectrum.horizontalFraction);
[peakModeFraction,peakModeIndex] = max(spectrum.modeFraction);
result = struct(label=label,time=time,timeDays=time/86400,outputInterval=representativeOutputInterval(time), ...
    generationRate=generationRate,integratedGeneration=trapz(time,generationRate), ...
    maximumHorizontalSpeed=maximumHorizontalSpeed,maximumThreeDimensionalSpeed=maximumThreeDimensionalSpeed, ...
    maximumVerticalSpeed=maximumVerticalSpeed,maximumSpeed=maximumSpeed, ...
    maximumSpeedTimeDays=time(maximumSpeedIndex)/86400,spectrum=spectrum,vorticity=vorticity, ...
    peakHorizontalWavelengthKilometers=spectrum.horizontalWavelengthKilometers(peakHorizontalIndex), ...
    peakHorizontalFraction=peakHorizontalFraction,peakMode=spectrum.verticalMode(peakModeIndex), ...
    peakModeFraction=peakModeFraction,Lx=wvt.Lx,Ly=wvt.Ly,Lz=wvt.Lz,Nx=wvt.Nx,Ny=wvt.Ny,Nz=wvt.Nz,Nj=wvt.Nj);
clear cleanup
end

function spectrum = finalWaveSpectrum(wvt)
waveEnergy = wvt.Apm_TE_factor.*(abs(wvt.waveComponent.maskAp.*wvt.Ap).^2+abs(wvt.waveComponent.maskAm.*wvt.Am).^2);
radialEnergy = wvt.transformToRadialWavenumber(waveEnergy);
totalEnergy = sum(waveEnergy,"all");
if ~isfinite(totalEnergy) || totalEnergy <= 0
    error("AnalyzeEddyTideBasicFigures:InvalidFinalWaveEnergy", "The final wave energy must be finite and positive.")
end
spectrum = struct(radialEnergy=radialEnergy,totalEnergy=totalEnergy, ...
    horizontalWavelengthKilometers=2*pi./wvt.kRadial/1e3,verticalMode=wvt.j(:), ...
    horizontalFraction=sum(radialEnergy,1)/totalEnergy,modeFraction=sum(radialEnergy,2)/totalEnergy);
end

function vorticity = finalSurfaceVorticity(wvt)
waveU = wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*wvt.Apt+wvt.UAm.*wvt.Amt);
waveV = wvt.transformToSpatialDomainWithF(Apm=wvt.VAp.*wvt.Apt+wvt.VAm.*wvt.Amt);
geostrophicU = wvt.transformToSpatialDomainWithF(A0=wvt.UA0.*wvt.A0t);
geostrophicV = wvt.transformToSpatialDomainWithF(A0=wvt.VA0.*wvt.A0t);
waveZeta = wvt.diffX(waveV)-wvt.diffY(waveU);
geostrophicZeta = wvt.diffX(geostrophicV)-wvt.diffY(geostrophicU);
[~,iSurface] = max(wvt.z);
waveZetaOverF = squeeze(waveZeta(:,:,iSurface))/abs(wvt.f);
geostrophicZetaOverF = squeeze(geostrophicZeta(:,:,iSurface))/abs(wvt.f);
vorticity = struct(xKilometers=wvt.x(:)/1e3,yKilometers=wvt.y(:)/1e3,z=wvt.z(iSurface), ...
    waveZetaOverF=waveZetaOverF,geostrophicZetaOverF=geostrophicZetaOverF, ...
    waveRmsOverF=rms(waveZetaOverF,"all"),waveMaximumOverF=max(abs(waveZetaOverF),[],"all"), ...
    geostrophicRmsOverF=rms(geostrophicZetaOverF,"all"),geostrophicMaximumOverF=max(abs(geostrophicZetaOverF),[],"all"));
end

function outputInterval = representativeOutputInterval(time)
if numel(time) < 2
    outputInterval = NaN;
else
    timeStep = diff(time);
    if max(abs(timeStep-timeStep(1))) > 100*eps(max(time))
        error("AnalyzeEddyTideBasicFigures:UnevenOutputTimes", "Saved model times must be evenly spaced.")
    end
    outputInterval = timeStep(1);
end
end

function validatePairedCases(eddy,noEddy,energy)
if ~isequal(eddy.time,noEddy.time)
    error("AnalyzeEddyTideBasicFigures:TimeMismatch", "The eddy and no-eddy simulations must have identical saved times.")
end
if ~isequal([eddy.Lx eddy.Ly eddy.Lz eddy.Nx eddy.Ny eddy.Nz],[noEddy.Lx noEddy.Ly noEddy.Lz noEddy.Nx noEddy.Ny noEddy.Nz])
    error("AnalyzeEddyTideBasicFigures:DomainMismatch", "The eddy and no-eddy simulations must use the same domain and resolution.")
end
if ~isequal(eddy.time,energy.eddy.time) || ~isequal(noEddy.time,energy.noEddy.time)
    error("AnalyzeEddyTideBasicFigures:DiagnosticsTimeMismatch", "WVDiagnostics and model-output times must agree.")
end
numericValues = [eddy.maximumThreeDimensionalSpeed; eddy.maximumVerticalSpeed; eddy.generationRate; ...
    noEddy.maximumThreeDimensionalSpeed; noEddy.maximumVerticalSpeed; noEddy.generationRate];
if any(~isfinite(numericValues))
    error("AnalyzeEddyTideBasicFigures:NonfiniteDiagnostics", "All reconstructed diagnostics must be finite.")
end
end

function figureHandle = createSpeedFigure(eddy,noEddy,visibility)
figureHandle = figure(Name="Eddy-tide maximum wave speed",Color="w",Visible=visibility,Position=[100 100 1200 700]);
layout = tiledlayout(figureHandle,2,1,TileSpacing="compact",Padding="compact");
axesHandle = nexttile(layout,1);
plot(axesHandle,eddy.timeDays,100*eddy.maximumThreeDimensionalSpeed,LineWidth=1.8,DisplayName=eddy.label)
hold(axesHandle,"on")
plot(axesHandle,noEddy.timeDays,100*noEddy.maximumThreeDimensionalSpeed,"--",LineWidth=1.8,DisplayName=noEddy.label)
grid(axesHandle,"on")
xlim(axesHandle,[eddy.timeDays(1) eddy.timeDays(end)])
ylabel(axesHandle,"maximum wave speed (cm s^{-1})")
title(axesHandle,"Maximum three-dimensional wave velocity")
legend(axesHandle,Location="best")
axesHandle = nexttile(layout,2);
plot(axesHandle,eddy.timeDays,100*eddy.maximumVerticalSpeed,LineWidth=1.8,DisplayName=eddy.label)
hold(axesHandle,"on")
plot(axesHandle,noEddy.timeDays,100*noEddy.maximumVerticalSpeed,"--",LineWidth=1.8,DisplayName=noEddy.label)
grid(axesHandle,"on")
xlim(axesHandle,[eddy.timeDays(1) eddy.timeDays(end)])
xlabel(axesHandle,"time (days)")
ylabel(axesHandle,"maximum |w| (cm s^{-1})")
title(axesHandle,"Maximum vertical wave velocity")
legend(axesHandle,Location="best")
title(layout,sprintf("Domain-maximum wave velocity sampled every %g hours",eddy.outputInterval/3600))
end

function figureHandle = createSpectrumFigure(eddy,noEddy,visibility)
validRadial = isfinite(eddy.spectrum.horizontalWavelengthKilometers) & eddy.spectrum.horizontalWavelengthKilometers > 0;
if ~isequal(eddy.spectrum.horizontalWavelengthKilometers,noEddy.spectrum.horizontalWavelengthKilometers) || ~isequal(eddy.spectrum.verticalMode,noEddy.spectrum.verticalMode)
    error("AnalyzeEddyTideBasicFigures:SpectrumGridMismatch", "The paired spectra must use identical wavelength and vertical-mode grids.")
end
eddyLogFraction = log10(eddy.spectrum.radialEnergy(:,validRadial)/eddy.spectrum.totalEnergy);
noEddyLogFraction = log10(noEddy.spectrum.radialEnergy(:,validRadial)/noEddy.spectrum.totalEnergy);
finiteValues = [eddyLogFraction(isfinite(eddyLogFraction)); noEddyLogFraction(isfinite(noEddyLogFraction))];
if isempty(finiteValues)
    error("AnalyzeEddyTideBasicFigures:EmptySpectrum", "The final paired spectra contain no finite nonzero bins.")
end
colorMaximum = max(finiteValues);
colorLimits = [colorMaximum-5 colorMaximum];
figureHandle = figure(Name="Eddy-tide final wave-energy spectrum",Color="w",Visible=visibility,Position=[100 100 1250 850]);
layout = tiledlayout(figureHandle,2,2,TileSpacing="compact",Padding="compact");
plotJointSpectrum(nexttile(layout,1),eddy,eddyLogFraction,validRadial,colorLimits)
plotJointSpectrum(nexttile(layout,2),noEddy,noEddyLogFraction,validRadial,colorLimits)
axesHandle = nexttile(layout,3);
loglog(axesHandle,eddy.spectrum.horizontalWavelengthKilometers(validRadial),eddy.spectrum.horizontalFraction(validRadial),LineWidth=2,DisplayName=eddy.label)
hold(axesHandle,"on")
loglog(axesHandle,noEddy.spectrum.horizontalWavelengthKilometers(validRadial),noEddy.spectrum.horizontalFraction(validRadial),"--",LineWidth=2,DisplayName=noEddy.label)
set(axesHandle,XDir="reverse")
grid(axesHandle,"on")
xlabel(axesHandle,"horizontal wavelength (km)")
ylabel(axesHandle,"wave-energy fraction per radial bin")
legend(axesHandle,Location="best")
title(axesHandle,"Horizontal spectrum")
axesHandle = nexttile(layout,4);
semilogy(axesHandle,eddy.spectrum.verticalMode,eddy.spectrum.modeFraction,"o-",LineWidth=2,DisplayName=eddy.label)
hold(axesHandle,"on")
semilogy(axesHandle,noEddy.spectrum.verticalMode,noEddy.spectrum.modeFraction,"s--",LineWidth=2,DisplayName=noEddy.label)
grid(axesHandle,"on")
xlabel(axesHandle,"vertical mode")
ylabel(axesHandle,"wave-energy fraction")
legend(axesHandle,Location="best")
title(axesHandle,"Vertical-mode spectrum")
colormap(figureHandle,"turbo")
title(layout,sprintf("Wave-energy spectrum at day %g",eddy.timeDays(end)))
end

function plotJointSpectrum(axesHandle,result,logFraction,validRadial,colorLimits)
surface(axesHandle,result.spectrum.horizontalWavelengthKilometers(validRadial),result.spectrum.verticalMode,logFraction,EdgeColor="none")
view(axesHandle,2)
set(axesHandle,XScale="log",XDir="reverse")
clim(axesHandle,colorLimits)
xlabel(axesHandle,"horizontal wavelength (km)")
ylabel(axesHandle,"vertical mode")
title(axesHandle,result.label)
colorbar(axesHandle)
end

function figureHandle = createVorticityFigure(eddy,noEddy,visibility)
waveLimit = max(abs([eddy.vorticity.waveZetaOverF(:); noEddy.vorticity.waveZetaOverF(:)]));
geostrophicLimit = max(abs([eddy.vorticity.geostrophicZetaOverF(:); noEddy.vorticity.geostrophicZetaOverF(:)]));
figureHandle = figure(Name="Eddy-tide final surface vorticity",Color="w",Visible=visibility,Position=[100 100 1320 900]);
layout = tiledlayout(figureHandle,2,2,TileSpacing="compact",Padding="compact");
plotVorticityPanel(nexttile(layout),eddy.vorticity,eddy.vorticity.waveZetaOverF,waveLimit,sprintf("%s: wave, max |\\zeta/f| = %.2f",eddy.label,eddy.vorticity.waveMaximumOverF))
plotVorticityPanel(nexttile(layout),eddy.vorticity,eddy.vorticity.geostrophicZetaOverF,geostrophicLimit,sprintf("%s: geostrophic, max |\\zeta/f| = %.2f",eddy.label,eddy.vorticity.geostrophicMaximumOverF))
plotVorticityPanel(nexttile(layout),noEddy.vorticity,noEddy.vorticity.waveZetaOverF,waveLimit,sprintf("%s: wave, max |\\zeta/f| = %.2f",noEddy.label,noEddy.vorticity.waveMaximumOverF))
plotVorticityPanel(nexttile(layout),noEddy.vorticity,noEddy.vorticity.geostrophicZetaOverF,geostrophicLimit,sprintf("%s: geostrophic, max |\\zeta/f| = %.2f",noEddy.label,noEddy.vorticity.geostrophicMaximumOverF))
colormap(figureHandle,blueWhiteRed(256))
title(layout,sprintf("Surface vertical vorticity at day %g",eddy.timeDays(end)))
xlabel(layout,"x (km)")
ylabel(layout,"y (km)")
end

function plotVorticityPanel(axesHandle,vorticity,field,colorLimit,panelTitle)
imagesc(axesHandle,vorticity.xKilometers,vorticity.yKilometers,field.')
axis(axesHandle,"xy","image")
if colorLimit > 0
    clim(axesHandle,[-colorLimit colorLimit])
end
title(axesHandle,panelTitle)
colorbar(axesHandle)
end

function colormapValues = blueWhiteRed(numberOfColors)
lowerCount = floor(numberOfColors/2);
upperCount = numberOfColors-lowerCount;
blue = [0.085 0.286 0.621];
white = [0.98 0.98 0.98];
red = [0.706 0.016 0.150];
lower = [linspace(blue(1),white(1),lowerCount).',linspace(blue(2),white(2),lowerCount).',linspace(blue(3),white(3),lowerCount).'];
upper = [linspace(white(1),red(1),upperCount).',linspace(white(2),red(2),upperCount).',linspace(white(3),red(3),upperCount).'];
colormapValues = [lower; upper];
end

function validateExportPaths(paths,shouldOverwriteExisting)
for path = paths(:).'
    directory = fileparts(path);
    if strlength(directory) > 0 && ~isfolder(directory)
        mkdir(directory)
    end
    if isfile(path) && ~shouldOverwriteExisting
        error("AnalyzeEddyTideBasicFigures:ExportFileExists", "The analysis output '%s' already exists. Set shouldOverwriteExisting=true to replace it.",path)
    end
end
end
