function [figureHandle,energy,diagnosticsFiles] = AnalyzeEddyTideEnergy(eddyFile,noEddyFile,options)
% Compare eddy and no-eddy wave and geostrophic energy.
%
% Standard `WVDiagnostics` files are created from the completed model
% outputs when needed. The plotted wave and geostrophic reservoirs are all
% normalized by the eddy simulation's initial geostrophic energy.
%
% ```matlab
% [fig,energy,diagnosticsFiles] = AnalyzeEddyTideEnergy(eddyFile,noEddyFile);
% ```
%
% - Topic: Examples
% - Declaration: [figureHandle,energy,diagnosticsFiles] = AnalyzeEddyTideEnergy(eddyFile,noEddyFile,options)
% - Parameter eddyFile: completed eddy simulation NetCDF path
% - Parameter noEddyFile: completed no-eddy simulation NetCDF path
% - Parameter options.diagnosticsStride: model-output stride used when creating new diagnostics
% - Parameter options.figureVisible: figure visibility, `"on"` or `"off"`
% - Parameter options.shouldExport: whether to export a PNG
% - Parameter options.exportPath: optional PNG path, default beside `eddyFile`
% - Parameter options.exportResolution: exported PNG resolution in dots per inch
% - Parameter options.shouldOverwriteExisting: whether to replace an existing PNG
% - Returns figureHandle: energy-comparison figure
% - Returns energy: raw and normalized energy time series
% - Returns diagnosticsFiles: eddy and no-eddy diagnostics paths
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
    energy (1,1) struct
    diagnosticsFiles (2,1) string
end

if exist("WVDiagnostics","class") ~= 8
    error("AnalyzeEddyTideEnergy:WVDiagnosticsNotFound", "WaveVortexModelDiagnostics is required. Install it with OceanKit or add its authoring repository to the MATLAB path.")
end

exportPath = "";
if options.shouldExport
    exportPath = options.exportPath;
    if strlength(exportPath) == 0
        exportPath = fullfile(fileparts(eddyFile),"eddy-tide-energy-comparison.png");
    elseif ~endsWith(exportPath,".png",IgnoreCase=true)
        exportPath = exportPath+".png";
    end
    if isfile(exportPath) && ~options.shouldOverwriteExisting
        error("AnalyzeEddyTideEnergy:ExportFileExists", "The figure '%s' already exists. Set shouldOverwriteExisting=true to replace it.",exportPath)
    end
    exportDirectory = fileparts(exportPath);
    if strlength(exportDirectory) > 0 && ~isfolder(exportDirectory)
        mkdir(exportDirectory)
    end
end

[eddyEnergy,eddyDiagnosticsFile] = energyFromSimulation(eddyFile,options.diagnosticsStride);
[noEddyEnergy,noEddyDiagnosticsFile] = energyFromSimulation(noEddyFile,options.diagnosticsStride);
diagnosticsFiles = [eddyDiagnosticsFile; noEddyDiagnosticsFile];

normalization = eddyEnergy.geostrophic(1);
if ~isfinite(normalization) || normalization <= 0
    error("AnalyzeEddyTideEnergy:InvalidNormalization", "The eddy simulation's initial geostrophic energy must be finite and positive.")
end

eddyEnergy.waveNormalized = eddyEnergy.wave/normalization;
eddyEnergy.geostrophicNormalized = eddyEnergy.geostrophic/normalization;
noEddyEnergy.waveNormalized = noEddyEnergy.wave/normalization;
noEddyEnergy.geostrophicNormalized = noEddyEnergy.geostrophic/normalization;
energy = struct(normalization=normalization,eddy=eddyEnergy,noEddy=noEddyEnergy,figurePath="");

figureHandle = figure(Name="Eddy-tide energy comparison",Color="w",Visible=options.figureVisible);
layout = tiledlayout(figureHandle,2,1,TileSpacing="compact",Padding="compact");

waveAxes = nexttile(layout,1);
plot(waveAxes,eddyEnergy.timeDays,eddyEnergy.waveNormalized,LineWidth=2,DisplayName="eddy")
hold(waveAxes,"on")
plot(waveAxes,noEddyEnergy.timeDays,noEddyEnergy.waveNormalized,"--",LineWidth=2,DisplayName="no eddy")
grid(waveAxes,"on")
ylabel(waveAxes,"E_w / E_{g,eddy}(0)")
title(waveAxes,"Wave energy")
legend(waveAxes,Location="best")

geostrophicAxes = nexttile(layout,2);
plot(geostrophicAxes,eddyEnergy.timeDays,eddyEnergy.geostrophicNormalized,LineWidth=2,DisplayName="eddy")
hold(geostrophicAxes,"on")
plot(geostrophicAxes,noEddyEnergy.timeDays,noEddyEnergy.geostrophicNormalized,"--",LineWidth=2,DisplayName="no eddy")
grid(geostrophicAxes,"on")
xlabel(geostrophicAxes,"time (days)")
ylabel(geostrophicAxes,"E_g / E_{g,eddy}(0)")
title(geostrophicAxes,"Geostrophic energy")
legend(geostrophicAxes,Location="best")
title(layout,"Topographically forced eddy-tide experiment")

if options.shouldExport
    exportgraphics(figureHandle,exportPath,Resolution=options.exportResolution)
    energy.figurePath = exportPath;
end
end

function [energy,diagnosticsFile] = energyFromSimulation(simulationFile,diagnosticsStride)
diagnostics = WVDiagnostics(char(simulationFile));
diagnosticsCleanup = onCleanup(@()diagnostics.close());
if ~isfile(diagnostics.diagpath)
    diagnostics.createDiagnosticsFile(stride=diagnosticsStride);
elseif diagnostics.t_diag(end) < diagnostics.t_wv(end)
    if numel(diagnostics.t_diag) < 2
        error("AnalyzeEddyTideEnergy:IncompleteDiagnostics", "The existing diagnostics file '%s' has only one incomplete record and cannot be appended safely.",diagnostics.diagpath)
    end
    diagnostics.createDiagnosticsFile(stride=diagnosticsStride);
end
[reservoirs,time] = diagnostics.quadraticEnergyOverTime(energyReservoirs=[EnergyReservoir.wave EnergyReservoir.geostrophic]);
waveEnergy = reservoirs(1).energy(:);
geostrophicEnergy = reservoirs(2).energy(:);
time = time(:);
if numel(waveEnergy) ~= numel(time) || numel(geostrophicEnergy) ~= numel(time)
    error("AnalyzeEddyTideEnergy:InvalidDiagnosticsSize", "Diagnostics energy and time variables must have the same number of records.")
end
if any(~isfinite([time; waveEnergy; geostrophicEnergy]))
    error("AnalyzeEddyTideEnergy:NonfiniteDiagnostics", "Diagnostics time and energy variables must be finite.")
end
diagnosticsFile = string(diagnostics.diagpath);
energy = struct(time=time,timeDays=time/86400,wave=waveEnergy,geostrophic=geostrophicEnergy);
clear diagnosticsCleanup
end
