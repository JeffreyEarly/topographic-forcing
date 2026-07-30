function [figureHandle,series,diagnosticsFiles] = AnalyzeEddyTideEnergyEnstrophy(eddyFile,noEddyFile,options)
% Compare normalized energy reservoirs and APV potential enstrophy.
%
% Standard `WVDiagnostics` files are created or appended from the model
% outputs when needed. Both simulations are normalized by the initial-eddy
% simulation's initial total quadratic energy and exact APV potential
% enstrophy.
%
% ```matlab
% [fig,series,diagnosticsFiles] = AnalyzeEddyTideEnergyEnstrophy(eddyFile,noEddyFile);
% ```
%
% - Topic: Examples
% - Declaration: [figureHandle,series,diagnosticsFiles] = AnalyzeEddyTideEnergyEnstrophy(eddyFile,noEddyFile,options)
% - Parameter eddyFile: completed initial-eddy simulation NetCDF path
% - Parameter noEddyFile: completed no-initial-eddy simulation NetCDF path
% - Parameter options.diagnosticsStride: model-output stride used when creating or appending diagnostics
% - Parameter options.figureVisible: figure visibility, `"on"` or `"off"`
% - Parameter options.shouldExport: whether to export a PNG
% - Parameter options.exportPath: optional PNG path, default beside `eddyFile`
% - Parameter options.exportResolution: exported PNG resolution in dots per inch
% - Parameter options.shouldOverwriteExisting: whether to replace an existing PNG
% - Returns figureHandle: two-panel energy and enstrophy figure
% - Returns series: raw and normalized paired time series
% - Returns diagnosticsFiles: initial-eddy and no-initial-eddy diagnostics paths
arguments (Input)
    eddyFile (1,1) string {mustBeFile}
    noEddyFile (1,1) string {mustBeFile}
    options.diagnosticsStride (1,1) double {mustBeInteger,mustBePositive} = 1
    options.figureVisible (1,1) string {mustBeMember(options.figureVisible,["on" "off"])} = "on"
    options.shouldExport (1,1) logical = true
    options.exportPath (1,1) string = ""
    options.exportResolution (1,1) double {mustBeInteger,mustBePositive} = 300
    options.shouldOverwriteExisting (1,1) logical = false
end
arguments (Output)
    figureHandle matlab.ui.Figure
    series (1,1) struct
    diagnosticsFiles (2,1) string
end

if exist("WVDiagnostics","class") ~= 8
    error("AnalyzeEddyTideEnergyEnstrophy:WVDiagnosticsNotFound", "WaveVortexModelDiagnostics is required. Install it with OceanKit or add its authoring repository to the MATLAB path.")
end

exportPath = prepareExportPath(eddyFile,options);
[eddySeries,eddyDiagnosticsFile] = seriesFromSimulation(eddyFile,options.diagnosticsStride);
[controlSeries,controlDiagnosticsFile] = seriesFromSimulation(noEddyFile,options.diagnosticsStride);
diagnosticsFiles = [eddyDiagnosticsFile; controlDiagnosticsFile];

if ~isequal(eddySeries.time,controlSeries.time)
    error("AnalyzeEddyTideEnergyEnstrophy:TimeAxesDiffer", "The initial-eddy and no-initial-eddy diagnostics must have identical time axes.")
end

energyNormalization = eddySeries.raw.energy.total(1);
enstrophyNormalization = eddySeries.raw.apvEnstrophy(1);
if ~isfinite(energyNormalization) || energyNormalization <= 0
    error("AnalyzeEddyTideEnergyEnstrophy:InvalidEnergyNormalization", "The initial-eddy simulation's initial total quadratic energy must be finite and positive.")
end
if ~isfinite(enstrophyNormalization) || enstrophyNormalization <= 0
    error("AnalyzeEddyTideEnergyEnstrophy:InvalidEnstrophyNormalization", "The initial-eddy simulation's initial APV potential enstrophy must be finite and positive.")
end

eddySeries.normalized = normalizeSeries(eddySeries.raw,energyNormalization,enstrophyNormalization);
controlSeries.normalized = normalizeSeries(controlSeries.raw,energyNormalization,enstrophyNormalization);
validateNormalizedSeries(eddySeries.normalized,"initial-eddy");
validateNormalizedSeries(controlSeries.normalized,"no-initial-eddy");

series = struct( ...
    normalization=struct(energy=energyNormalization,enstrophy=enstrophyNormalization), ...
    eddy=eddySeries,control=controlSeries,figurePath="");

figureHandle = makeFigure(series,options.figureVisible);
if options.shouldExport
    exportgraphics(figureHandle,exportPath,Resolution=options.exportResolution)
    series.figurePath = exportPath;
end
end

function exportPath = prepareExportPath(eddyFile,options)
exportPath = "";
if ~options.shouldExport
    return
end

exportPath = options.exportPath;
if strlength(exportPath) == 0
    exportPath = fullfile(fileparts(eddyFile),"eddy-tide-energy-enstrophy.png");
elseif ~endsWith(exportPath,".png",IgnoreCase=true)
    exportPath = exportPath+".png";
end
if isfile(exportPath) && ~options.shouldOverwriteExisting
    error("AnalyzeEddyTideEnergyEnstrophy:ExportFileExists", "The figure '%s' already exists. Set shouldOverwriteExisting=true to replace it.",exportPath)
end
exportDirectory = fileparts(exportPath);
if strlength(exportDirectory) > 0 && ~isfolder(exportDirectory)
    mkdir(exportDirectory)
end
end

function [simulationSeries,diagnosticsFile] = seriesFromSimulation(simulationFile,diagnosticsStride)
diagnostics = WVDiagnostics(char(simulationFile));
diagnosticsCleanup = onCleanup(@()diagnostics.close());
if ~isfile(diagnostics.diagpath)
    diagnostics.createDiagnosticsFile(stride=diagnosticsStride);
elseif diagnostics.t_diag(end) < diagnostics.t_wv(end)
    if numel(diagnostics.t_diag) < 2
        error("AnalyzeEddyTideEnergyEnstrophy:IncompleteDiagnostics", "The existing diagnostics file '%s' has only one incomplete record and cannot be appended safely.",diagnostics.diagpath)
    end
    diagnostics.createDiagnosticsFile(stride=diagnosticsStride);
end

[geostrophic,geostrophicKinetic,geostrophicPotential,wave,kinetic,potentialQuadratic,apvEnstrophy] = ...
    diagnostics.diagfile.readVariables("E_g","KE_g","PE_g","E_w","ke","pe_quadratic","enstrophy_apv");
time = diagnostics.t_diag(:);
raw = struct( ...
    energy=struct( ...
        total=kinetic(:)+potentialQuadratic(:), ...
        wave=wave(:), ...
        geostrophic=geostrophic(:), ...
        geostrophicKinetic=geostrophicKinetic(:), ...
        geostrophicPotential=geostrophicPotential(:)), ...
    apvEnstrophy=apvEnstrophy(:));
validateRawSeries(time,raw,diagnostics.diagpath)

diagnosticsFile = string(diagnostics.diagpath);
simulationSeries = struct(time=time,timeDays=time/86400,raw=raw);
clear diagnosticsCleanup
end

function validateRawSeries(time,raw,diagnosticsPath)
energyNames = string(fieldnames(raw.energy));
values = time;
for iEnergy = 1:numel(energyNames)
    energy = raw.energy.(energyNames(iEnergy));
    if numel(energy) ~= numel(time)
        error("AnalyzeEddyTideEnergyEnstrophy:InvalidDiagnosticsSize", "The diagnostics variable '%s' in '%s' must have the same number of records as time.",energyNames(iEnergy),diagnosticsPath)
    end
    values = [values; energy]; %#ok<AGROW>
end
if numel(raw.apvEnstrophy) ~= numel(time)
    error("AnalyzeEddyTideEnergyEnstrophy:InvalidDiagnosticsSize", "APV potential enstrophy in '%s' must have the same number of records as time.",diagnosticsPath)
end
values = [values; raw.apvEnstrophy];
if any(~isfinite(values))
    error("AnalyzeEddyTideEnergyEnstrophy:NonfiniteDiagnostics", "Diagnostics time, energy, and APV potential enstrophy variables in '%s' must be finite.",diagnosticsPath)
end
end

function normalized = normalizeSeries(raw,energyNormalization,enstrophyNormalization)
normalizedEnergy = struct();
energyNames = string(fieldnames(raw.energy));
for iEnergy = 1:numel(energyNames)
    normalizedEnergy.(energyNames(iEnergy)) = raw.energy.(energyNames(iEnergy))/energyNormalization;
end
normalized = struct(energy=normalizedEnergy,apvEnstrophy=raw.apvEnstrophy/enstrophyNormalization);
end

function validateNormalizedSeries(normalized,caseName)
energyNames = string(fieldnames(normalized.energy));
values = normalized.apvEnstrophy;
for iEnergy = 1:numel(energyNames)
    values = [values; normalized.energy.(energyNames(iEnergy))]; %#ok<AGROW>
end
if any(~isfinite(values))
    error("AnalyzeEddyTideEnergyEnstrophy:NonfiniteNormalizedSeries", "All normalized %s energy and enstrophy series must be finite.",caseName)
end
end

function figureHandle = makeFigure(series,figureVisible)
figureHandle = figure( ...
    Name="Eddy-tide energy and enstrophy", ...
    Color="w", ...
    Visible=figureVisible, ...
    Units="inches", ...
    Position=[1 1 7.5 7.5]);
layout = tiledlayout(figureHandle,2,1,TileSpacing="compact",Padding="compact");

constituentNames = ["total" "wave" "geostrophic" "geostrophicKinetic" "geostrophicPotential"];
constituentLabels = [ ...
    "Total Energy $\mathcal{E}$" ...
    "Wave Energy $\mathcal{E}_w$" ...
    "Geostrophic Energy $\mathcal{E}_g$" ...
    "Geostrophic Kinetic $\mathcal{K}_g$" ...
    "Geostrophic Potential $\mathcal{P}_g$"];
constituentColors = [ ...
    0.0000 0.0000 0.0000
    0.0000 0.4470 0.7410
    0.8500 0.3250 0.0980
    0.9290 0.6940 0.1250
    0.4940 0.1840 0.5560];

energyAxes = nexttile(layout,1);
hold(energyAxes,"on")
for iConstituent = 1:numel(constituentNames)
    plot(energyAxes,series.eddy.timeDays,series.eddy.normalized.energy.(constituentNames(iConstituent)), ...
        Color=constituentColors(iConstituent,:),LineStyle="-",LineWidth=1.8,HandleVisibility="off")
    plot(energyAxes,series.control.timeDays,series.control.normalized.energy.(constituentNames(iConstituent)), ...
        Color=constituentColors(iConstituent,:),LineStyle="--",LineWidth=1.8,HandleVisibility="off")
end
legendHandles = gobjects(7,1);
for iConstituent = 1:numel(constituentNames)
    legendHandles(iConstituent) = plot(energyAxes,NaN,NaN, ...
        Color=constituentColors(iConstituent,:),LineWidth=1.8,DisplayName=constituentLabels(iConstituent));
end
legendHandles(6) = plot(energyAxes,NaN,NaN,Color=[0.35 0.35 0.35],LineStyle="-",LineWidth=1.8,DisplayName="Initial eddy");
legendHandles(7) = plot(energyAxes,NaN,NaN,Color=[0.35 0.35 0.35],LineStyle="--",LineWidth=1.8,DisplayName="No initial eddy");
grid(energyAxes,"on")
xlim(energyAxes,timeLimits(series.eddy.timeDays))
ylim(energyAxes,dataLimits(energyValues(series)))
ylabel(energyAxes,"Normalized Energy")
title(energyAxes,"Energy reservoirs")
legend(energyAxes,legendHandles,Interpreter="latex",Location="northwest",NumColumns=2)

enstrophyAxes = nexttile(layout,2);
hold(enstrophyAxes,"on")
enstrophyColor = [0.3010 0.7450 0.9330];
plot(enstrophyAxes,series.eddy.timeDays,series.eddy.normalized.apvEnstrophy, ...
    Color=enstrophyColor,LineStyle="-",LineWidth=1.8,DisplayName="Initial eddy")
plot(enstrophyAxes,series.control.timeDays,series.control.normalized.apvEnstrophy, ...
    Color=enstrophyColor,LineStyle="--",LineWidth=1.8,DisplayName="No initial eddy")
grid(enstrophyAxes,"on")
xlim(enstrophyAxes,timeLimits(series.eddy.timeDays))
ylim(enstrophyAxes,dataLimits([series.eddy.normalized.apvEnstrophy; series.control.normalized.apvEnstrophy]))
xlabel(enstrophyAxes,"Time (days)")
ylabel(enstrophyAxes,"Normalized APV Enstrophy")
title(enstrophyAxes,"Potential enstrophy $\mathcal{Y}$",Interpreter="latex")
legend(enstrophyAxes,Interpreter="latex",Location="northwest")
title(layout,"Topographically forced eddy-tide experiment")
end

function values = energyValues(series)
values = [];
energyNames = string(fieldnames(series.eddy.normalized.energy));
for iEnergy = 1:numel(energyNames)
    values = [values; series.eddy.normalized.energy.(energyNames(iEnergy)); ...
        series.control.normalized.energy.(energyNames(iEnergy))]; %#ok<AGROW>
end
end

function limits = timeLimits(timeDays)
limits = [min(timeDays) max(timeDays)];
if limits(1) == limits(2)
    limits = limits+[-0.5 0.5];
end
end

function limits = dataLimits(values)
minimumValue = min(values);
maximumValue = max(values);
span = maximumValue-minimumValue;
if span == 0
    span = max(abs(maximumValue),1);
end
padding = 0.04*span;
limits = [min(0,minimumValue-padding) maximumValue+padding];
end
