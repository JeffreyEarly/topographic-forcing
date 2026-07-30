function [figureHandle,analysis,diagnosticsFile] = AnalyzeEddyTideAdaptiveDamping(eddyFile,options)
% Plot azimuthal-mean eddy velocity and adaptive-damping energy removal.
%
% The adaptive-damping panel shows the exact signed spatial work density
% from WaveVortexModel's spectral damping operator. Its volume integral is
% an energy sink, but the local field is not a non-negative viscous
% dissipation density. See `AnalyzeEddyTideAdaptiveDamping.md`.
%
% ```matlab
% [fig,analysis,diagnosticsFile] = AnalyzeEddyTideAdaptiveDamping(eddyFile);
% ```
%
% - Topic: Examples
% - Declaration: [figureHandle,analysis,diagnosticsFile] = AnalyzeEddyTideAdaptiveDamping(eddyFile,options)
% - Parameter eddyFile: initial-eddy simulation NetCDF path
% - Parameter options.timeRangeDays: inclusive averaging interval in days
% - Parameter options.radialBinCount: number of equal-width radial bins
% - Parameter options.eddyCenter: eddy center `[x y]` in meters, or `[NaN NaN]` for the domain center
% - Parameter options.eddyVelocityThreshold: eddy-aligned azimuthal-velocity contour in meters per second
% - Parameter options.diagnosticsStride: model-output stride used when creating or appending diagnostics
% - Parameter options.figureVisible: figure visibility, `"on"` or `"off"`
% - Parameter options.shouldExport: whether to export a PNG
% - Parameter options.exportPath: optional PNG path, default beside `eddyFile`
% - Parameter options.exportResolution: exported PNG resolution in dots per inch
% - Parameter options.shouldOverwriteExisting: whether to replace an existing PNG
% - Returns figureHandle: two-panel velocity and adaptive-damping figure
% - Returns analysis: averaged fields, coordinates, sampling metadata, and energy-flux validation
% - Returns diagnosticsFile: standard WVDiagnostics path used for validation
arguments (Input)
    eddyFile (1,1) string {mustBeFile}
    options.timeRangeDays (1,2) double {mustBeFinite,mustBeNonnegative} = [25 30]
    options.radialBinCount (1,1) double {mustBeInteger,mustBePositive} = 64
    options.eddyCenter (1,2) double = [NaN NaN]
    options.eddyVelocityThreshold (1,1) double {mustBePositive} = 0.02
    options.diagnosticsStride (1,1) double {mustBeInteger,mustBePositive} = 1
    options.figureVisible (1,1) string {mustBeMember(options.figureVisible,["on" "off"])} = "on"
    options.shouldExport (1,1) logical = true
    options.exportPath (1,1) string = ""
    options.exportResolution (1,1) double {mustBeInteger,mustBePositive} = 300
    options.shouldOverwriteExisting (1,1) logical = false
end
arguments (Output)
    figureHandle matlab.ui.Figure
    analysis (1,1) struct
    diagnosticsFile (1,1) string
end

if exist("WVDiagnostics","class") ~= 8
    error("AnalyzeEddyTideAdaptiveDamping:WVDiagnosticsNotFound", "WaveVortexModelDiagnostics is required. Install it with OceanKit or add its authoring repository to the MATLAB path.")
end
if options.timeRangeDays(2) <= options.timeRangeDays(1)
    error("AnalyzeEddyTideAdaptiveDamping:InvalidTimeRange", "timeRangeDays must contain a strictly increasing start and end time.")
end
if ~(all(isnan(options.eddyCenter)) || all(isfinite(options.eddyCenter)))
    error("AnalyzeEddyTideAdaptiveDamping:InvalidEddyCenter", "eddyCenter must contain two finite coordinates or `[NaN NaN]`.")
end

exportPath = prepareExportPath(eddyFile,options);
diagnostics = WVDiagnostics(char(eddyFile));
diagnosticsCleanup = onCleanup(@()diagnostics.close());
ensureDiagnosticsAreCurrent(diagnostics,options.diagnosticsStride)
diagnosticsFile = string(diagnostics.diagpath);

wvt = diagnostics.wvt;
if ~wvt.hasForcingWithName("adaptive damping")
    error("AnalyzeEddyTideAdaptiveDamping:AdaptiveDampingNotFound", "The simulation must contain a forcing named 'adaptive damping'.")
end

eddyCenter = options.eddyCenter;
if all(isnan(eddyCenter))
    eddyCenter = [wvt.Lx/2 wvt.Ly/2];
end
if eddyCenter(1) < 0 || eddyCenter(1) >= wvt.Lx || eddyCenter(2) < 0 || eddyCenter(2) >= wvt.Ly
    error("AnalyzeEddyTideAdaptiveDamping:EddyCenterOutsideDomain", "eddyCenter must lie inside the periodic horizontal domain.")
end

modelTime = diagnostics.t_wv(:);
timeRange = 86400*options.timeRangeDays;
timeTolerance = max(1e-9,100*eps(max(abs(modelTime))));
selectedIndices = find(modelTime >= timeRange(1)-timeTolerance & modelTime <= timeRange(2)+timeTolerance);
if numel(selectedIndices) < 2
    error("AnalyzeEddyTideAdaptiveDamping:InsufficientTimeRecords", "The requested averaging interval must contain at least two saved model records.")
end
selectedTime = modelTime(selectedIndices);
if abs(selectedTime(1)-timeRange(1)) > timeTolerance || abs(selectedTime(end)-timeRange(2)) > timeTolerance
    error("AnalyzeEddyTideAdaptiveDamping:MissingTimeRangeEndpoint", "Both endpoints of timeRangeDays must be present in the saved model times.")
end
if any(diff(selectedTime) <= 0)
    error("AnalyzeEddyTideAdaptiveDamping:InvalidTimeAxis", "Saved model times in the averaging interval must be strictly increasing.")
end
timeWeights = trapezoidalMeanWeights(selectedTime);

diagnosticTime = diagnostics.t_diag(:);
[hasDiagnosticTime,diagnosticIndices] = ismember(selectedTime,diagnosticTime);
if any(~hasDiagnosticTime)
    error("AnalyzeEddyTideAdaptiveDamping:DiagnosticsTimeMismatch", "The diagnostics file must contain every selected model-output time.")
end
[exactFluxes,exactFluxTime] = diagnostics.exactEnergyFluxesOverTime(timeIndices=diagnosticIndices);
if ~isequal(exactFluxTime(:),selectedTime)
    error("AnalyzeEddyTideAdaptiveDamping:DiagnosticsTimeMismatch", "The selected diagnostics and model-output times must agree.")
end
forcingNames = string({exactFluxes.fancyName});
dampingIndex = forcingNames == "adaptive damping";
if nnz(dampingIndex) ~= 1
    error("AnalyzeEddyTideAdaptiveDamping:InvalidAdaptiveDampingDiagnostics", "The diagnostics file must contain exactly one adaptive-damping energy flux.")
end
diagnosticWork = exactFluxes(dampingIndex).te(:);

[radius,theta,radiusMaximum] = polarCoordinates(wvt,eddyCenter);
orientationFactor = eddyOrientationFactor(diagnostics,radius,theta,radiusMaximum);
[meanAzimuthalVelocity,meanWorkDensity,spatialWork] = averageSpatialFields(diagnostics,selectedIndices,timeWeights,theta);
[radialCenters,radialEdges,radialCounts,azimuthalVelocity] = azimuthalAverage(meanAzimuthalVelocity,radius,radiusMaximum,options.radialBinCount);
[~,~,~,energyRemovalRate] = azimuthalAverage(-wvt.rho0*meanWorkDensity,radius,radiusMaximum,options.radialBinCount);
azimuthalVelocity = orientationFactor*azimuthalVelocity;

[relativeFluxError,absoluteFluxError] = validateEnergyFlux(spatialWork,diagnosticWork);
tidalPeriod = NaN;
if wvt.hasForcingWithName("bottom wave generation")
    generation = wvt.forcingWithName("bottom wave generation");
    tidalPeriod = 2*pi/generation.frequency;
    if max(diff(selectedTime)) > tidalPeriod/4
        warning("AnalyzeEddyTideAdaptiveDamping:UnderresolvedTidalCadence", "The saved output interval exceeds one quarter of the forcing period. Interpret the result as a sampled multi-cycle mean, not a resolved tidal-phase average.")
    end
end

hasThresholdContour = min(azimuthalVelocity,[],"all") <= options.eddyVelocityThreshold && max(azimuthalVelocity,[],"all") >= options.eddyVelocityThreshold;
if ~hasThresholdContour
    warning("AnalyzeEddyTideAdaptiveDamping:ThresholdContourMissing", "The eddy-aligned azimuthal-velocity field does not cross the requested threshold, so no threshold contour will be visible.")
end

analysis = struct( ...
    modelFile=eddyFile, ...
    diagnosticsFile=diagnosticsFile, ...
    figurePath="", ...
    time=selectedTime, ...
    timeDays=selectedTime/86400, ...
    timeRangeDays=options.timeRangeDays, ...
    timeWeights=timeWeights, ...
    eddyCenter=eddyCenter, ...
    eddyOrientationFactor=orientationFactor, ...
    eddyVelocityThreshold=options.eddyVelocityThreshold, ...
    radius=radialCenters, ...
    radialEdges=radialEdges, ...
    radialCounts=radialCounts, ...
    depth=wvt.z(:), ...
    azimuthalVelocity=azimuthalVelocity, ...
    adaptiveDampingEnergyRemovalRate=energyRemovalRate, ...
    sampling=struct(outputInterval=max(diff(selectedTime)),tidalPeriod=tidalPeriod), ...
    validation=struct( ...
        spatialWork=spatialWork, ...
        diagnosticWork=diagnosticWork, ...
        absoluteFluxError=absoluteFluxError, ...
        relativeFluxError=relativeFluxError, ...
        maximumRelativeFluxError=max(relativeFluxError), ...
        timeMeanSpatialWork=sum(timeWeights.*spatialWork), ...
        timeMeanDiagnosticWork=sum(timeWeights.*diagnosticWork), ...
        timeMeanEnergyRemovalPerArea=-wvt.rho0*sum(timeWeights.*spatialWork), ...
        hasThresholdContour=hasThresholdContour));

figureHandle = makeFigure(analysis,options.figureVisible);
if options.shouldExport
    exportgraphics(figureHandle,exportPath,Resolution=options.exportResolution)
    analysis.figurePath = exportPath;
end
clear diagnosticsCleanup
end

function exportPath = prepareExportPath(eddyFile,options)
exportPath = "";
if ~options.shouldExport
    return
end

exportPath = options.exportPath;
if strlength(exportPath) == 0
    exportPath = fullfile(fileparts(eddyFile),"eddy-tide-adaptive-damping.png");
elseif ~endsWith(exportPath,".png",IgnoreCase=true)
    exportPath = exportPath+".png";
end
if isfile(exportPath) && ~options.shouldOverwriteExisting
    error("AnalyzeEddyTideAdaptiveDamping:ExportFileExists", "The figure '%s' already exists. Set shouldOverwriteExisting=true to replace it.",exportPath)
end
exportDirectory = fileparts(exportPath);
if strlength(exportDirectory) > 0 && ~isfolder(exportDirectory)
    mkdir(exportDirectory)
end
end

function ensureDiagnosticsAreCurrent(diagnostics,diagnosticsStride)
if ~isfile(diagnostics.diagpath)
    diagnostics.createDiagnosticsFile(stride=diagnosticsStride);
elseif diagnostics.t_diag(end) < diagnostics.t_wv(end)
    if numel(diagnostics.t_diag) < 2
        error("AnalyzeEddyTideAdaptiveDamping:IncompleteDiagnostics", "The existing diagnostics file '%s' has only one incomplete record and cannot be appended safely.",diagnostics.diagpath)
    end
    diagnostics.createDiagnosticsFile(stride=diagnosticsStride);
end
end

function weights = trapezoidalMeanWeights(time)
intervals = diff(time);
weights = zeros(size(time));
weights(1) = intervals(1)/2;
weights(end) = intervals(end)/2;
if numel(time) > 2
    weights(2:end-1) = (intervals(1:end-1)+intervals(2:end))/2;
end
weights = weights/(time(end)-time(1));
end

function [radius,theta,radiusMaximum] = polarCoordinates(wvt,eddyCenter)
deltaX = mod(wvt.X(:,:,1)-eddyCenter(1)+wvt.Lx/2,wvt.Lx)-wvt.Lx/2;
deltaY = mod(wvt.Y(:,:,1)-eddyCenter(2)+wvt.Ly/2,wvt.Ly)-wvt.Ly/2;
radius = hypot(deltaX,deltaY);
theta = atan2(deltaY,deltaX);
radiusMaximum = min(wvt.Lx,wvt.Ly)/2;
end

function orientationFactor = eddyOrientationFactor(diagnostics,radius,theta,radiusMaximum)
diagnostics.iTime = 1;
wvt = diagnostics.wvt;
[~,surfaceIndex] = max(wvt.z);
standardAzimuthalVelocity = -wvt.u_g.*sin(theta)+wvt.v_g.*cos(theta);
surfaceVelocity = standardAzimuthalVelocity(:,:,surfaceIndex);
eddyMask = radius > 0 & radius <= radiusMaximum;
surfaceVelocity = surfaceVelocity(eddyMask);
[dominantSpeed,dominantIndex] = max(abs(surfaceVelocity));
if ~isfinite(dominantSpeed) || dominantSpeed == 0
    error("AnalyzeEddyTideAdaptiveDamping:UndefinedEddyOrientation", "The initial geostrophic flow does not define a finite, nonzero azimuthal-velocity orientation.")
end
orientationFactor = sign(surfaceVelocity(dominantIndex));
end

function [meanAzimuthalVelocity,meanWorkDensity,spatialWork] = averageSpatialFields(diagnostics,selectedIndices,timeWeights,theta)
meanAzimuthalVelocity = [];
meanWorkDensity = [];
spatialWork = zeros(numel(selectedIndices),1);
for iTime = 1:numel(selectedIndices)
    diagnostics.iTime = selectedIndices(iTime);
    wvt = diagnostics.wvt;
    if wvt.isHydrostatic
        [Fu,Fv,Feta] = wvt.spatialFluxForForcingWithName("adaptive damping");
        workDensity = wvt.u.*Fu+wvt.v.*Fv+wvt.eta_true.*reshape(wvt.N2,1,1,[]).*Feta;
    else
        [Fu,Fv,Fw,Feta] = wvt.spatialFluxForForcingWithName("adaptive damping");
        workDensity = wvt.u.*Fu+wvt.v.*Fv+wvt.w.*Fw+wvt.eta_true.*reshape(wvt.N2,1,1,[]).*Feta;
    end
    standardAzimuthalVelocity = -wvt.u_g.*sin(theta)+wvt.v_g.*cos(theta);
    if isempty(meanAzimuthalVelocity)
        meanAzimuthalVelocity = zeros(size(standardAzimuthalVelocity));
        meanWorkDensity = zeros(size(workDensity));
    end
    meanAzimuthalVelocity = meanAzimuthalVelocity+timeWeights(iTime)*standardAzimuthalVelocity;
    meanWorkDensity = meanWorkDensity+timeWeights(iTime)*workDensity;
    spatialWork(iTime) = sum(mean(mean(reshape(wvt.z_int,1,1,[]).*workDensity,1),2),3);
end
end

function [radialCenters,radialEdges,radialCounts,radialMean] = azimuthalAverage(field,radius,radiusMaximum,radialBinCount)
radialEdges = linspace(0,radiusMaximum,radialBinCount+1).';
radialCenters = (radialEdges(1:end-1)+radialEdges(2:end))/2;
radialCounts = zeros(radialBinCount,1);
radialMean = NaN(radialBinCount,size(field,3));
radius = radius(:);
field = reshape(field,[],size(field,3));
for iRadius = 1:radialBinCount
    if iRadius == radialBinCount
        radialMask = radius > 0 & radius >= radialEdges(iRadius) & radius <= radialEdges(iRadius+1);
    else
        radialMask = radius > 0 & radius >= radialEdges(iRadius) & radius < radialEdges(iRadius+1);
    end
    radialCounts(iRadius) = nnz(radialMask);
    if radialCounts(iRadius) > 0
        radialMean(iRadius,:) = mean(field(radialMask,:),1);
    end
end
if any(radialCounts == 0) || any(~isfinite(radialMean),"all")
    error("AnalyzeEddyTideAdaptiveDamping:EmptyRadialBin", "Every radial bin must contain finite model-grid values. Reduce radialBinCount for this resolution.")
end
end

function [relativeError,absoluteError] = validateEnergyFlux(spatialWork,diagnosticWork)
absoluteError = abs(spatialWork-diagnosticWork);
scale = max(abs(spatialWork),abs(diagnosticWork));
relativeError = absoluteError./max(scale,realmin);
if any(absoluteError > 1e-14 & relativeError > 1e-10)
    error("AnalyzeEddyTideAdaptiveDamping:EnergyFluxMismatch", "The spatially integrated adaptive-damping work does not match the exact diagnostics energy flux.")
end
if any(diagnosticWork > 1e-14)
    error("AnalyzeEddyTideAdaptiveDamping:PositiveGlobalDampingWork", "Adaptive damping must be a non-positive global energy tendency within numerical tolerance.")
end
end

function figureHandle = makeFigure(analysis,figureVisible)
figureHandle = figure(Name="Eddy-tide adaptive damping",Color="w",Visible=figureVisible,Units="inches",Position=[1 1 10.5 4.6]);
layout = tiledlayout(figureHandle,1,2,TileSpacing="compact",Padding="compact");
radiusKilometers = analysis.radius/1e3;
depthKilometers = analysis.depth/1e3;
velocityCentimetersPerSecond = 100*analysis.azimuthalVelocity;
energyRemovalMicrowatts = 1e6*analysis.adaptiveDampingEnergyRemovalRate;
divergingMap = blueWhiteRed(257);

velocityAxes = nexttile(layout,1);
surface(velocityAxes,radiusKilometers,depthKilometers,zeros(size(velocityCentimetersPerSecond.')),velocityCentimetersPerSecond.',EdgeColor="none")
view(velocityAxes,2)
colormap(velocityAxes,divergingMap)
velocityLimit = max(abs(velocityCentimetersPerSecond),[],"all");
clim(velocityAxes,symmetricLimits(velocityLimit))
hold(velocityAxes,"on")
contour(velocityAxes,radiusKilometers,depthKilometers,velocityCentimetersPerSecond.',100*[analysis.eddyVelocityThreshold analysis.eddyVelocityThreshold],"k",LineWidth=2)
formatAxes(velocityAxes,radiusKilometers,depthKilometers)
xlabel(velocityAxes,"Radius (km)")
ylabel(velocityAxes,"Depth (km)")
title(velocityAxes,"(a) Eddy-aligned azimuthal velocity")
velocityColorbar = colorbar(velocityAxes);
velocityColorbar.Label.String = "cm s^{-1}";

dampingAxes = nexttile(layout,2);
surface(dampingAxes,radiusKilometers,depthKilometers,zeros(size(energyRemovalMicrowatts.')),energyRemovalMicrowatts.',EdgeColor="none")
view(dampingAxes,2)
colormap(dampingAxes,divergingMap)
dampingLimit = max(abs(energyRemovalMicrowatts),[],"all");
clim(dampingAxes,symmetricLimits(dampingLimit))
hold(dampingAxes,"on")
contour(dampingAxes,radiusKilometers,depthKilometers,velocityCentimetersPerSecond.',100*[analysis.eddyVelocityThreshold analysis.eddyVelocityThreshold],"k",LineWidth=2)
formatAxes(dampingAxes,radiusKilometers,depthKilometers)
xlabel(dampingAxes,"Radius (km)")
ylabel(dampingAxes,"Depth (km)")
title(dampingAxes,"(b) Signed adaptive-damping energy removal")
dampingColorbar = colorbar(dampingAxes);
dampingColorbar.Label.String = "\muW m^{-3}";

title(layout,sprintf("Ein10 eddy, days %g-%g mean; thick contour: u_{\\theta,eddy} = %.2f m s^{-1}",analysis.timeRangeDays(1),analysis.timeRangeDays(2),analysis.eddyVelocityThreshold))
end

function formatAxes(axesHandle,radiusKilometers,depthKilometers)
box(axesHandle,"on")
axesHandle.Layer = "top";
axesHandle.YDir = "normal";
xlim(axesHandle,[0 max(radiusKilometers)])
ylim(axesHandle,[min(depthKilometers) max(depthKilometers)])
end

function limits = symmetricLimits(maximumMagnitude)
if ~isfinite(maximumMagnitude) || maximumMagnitude <= 0
    maximumMagnitude = 1;
end
limits = [-maximumMagnitude maximumMagnitude];
end

function map = blueWhiteRed(colorCount)
positions = [0 0.5 1];
anchors = [0.230 0.299 0.754; 0.94 0.94 0.94; 0.706 0.016 0.150];
map = interp1(positions,anchors,linspace(0,1,colorCount));
end
