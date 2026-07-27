function [figureHandles,budget,diagnosticsFiles] = AnalyzeEddyTideBudgets(eddyFile,noEddyFile,options)
% Diagnose forcing-resolved eddy-tide energy and enstrophy budgets.
%
% Standard `WVDiagnostics` files are created or extended as needed. The
% returned budget contains state histories, forcing and nonlinear-triad
% rates, cumulative integrals, closure metrics, and eddy-minus-control
% differences.
%
% ```matlab
% [figures,budget,diagnosticsFiles] = AnalyzeEddyTideBudgets(eddyFile,noEddyFile);
% ```
%
% - Topic: Examples
% - Declaration: [figureHandles,budget,diagnosticsFiles] = AnalyzeEddyTideBudgets(eddyFile,noEddyFile,options)
% - Parameter eddyFile: completed eddy simulation NetCDF path
% - Parameter noEddyFile: completed no-eddy simulation NetCDF path
% - Parameter options.diagnosticsStride: model-output stride used when creating diagnostics
% - Parameter options.figureVisible: figure visibility, `"on"` or `"off"`
% - Parameter options.shouldExport: whether to export PNG figures
% - Parameter options.exportDirectory: output directory, default beside `eddyFile`
% - Parameter options.exportPrefix: prefix for exported figure filenames
% - Parameter options.shouldRequireClosure: whether failed closure tolerances raise an error
% - Returns figureHandles: reservoir, forcing-budget, and triad-budget figures
% - Returns budget: forcing-resolved diagnostics and closure metrics
% - Returns diagnosticsFiles: eddy and no-eddy WVDiagnostics paths
arguments (Input)
    eddyFile (1,1) string {mustBeFile}
    noEddyFile (1,1) string {mustBeFile}
    options.diagnosticsStride (1,1) double {mustBeInteger,mustBePositive} = 1
    options.figureVisible (1,1) string {mustBeMember(options.figureVisible,["on" "off"])} = "on"
    options.shouldExport (1,1) logical = true
    options.exportDirectory (1,1) string = ""
    options.exportPrefix (1,1) string = "eddy-tide"
    options.exportResolution (1,1) double {mustBeInteger,mustBePositive} = 300
    options.shouldOverwriteExisting (1,1) logical = false
    options.shouldRequireClosure (1,1) logical = false
    options.quadraticClosureTolerance (1,1) double {mustBePositive,mustBeFinite} = 0.05
    options.exactClosureTolerance (1,1) double {mustBePositive,mustBeFinite} = 0.10
end

if exist("WVDiagnostics","class") ~= 8
    error("AnalyzeEddyTideBudgets:WVDiagnosticsNotFound", ...
        "WaveVortexModelDiagnostics is required. Install it with OceanKit or add its authoring repository to the MATLAB path.")
end

exportDirectory = options.exportDirectory;
if strlength(exportDirectory) == 0
    exportDirectory = string(fileparts(eddyFile));
end
figurePaths = [
    fullfile(exportDirectory,options.exportPrefix+"-reservoirs.png")
    fullfile(exportDirectory,options.exportPrefix+"-forcing-budgets.png")
    fullfile(exportDirectory,options.exportPrefix+"-triad-budgets.png")
    ];
if options.shouldExport
    validateExportPaths(figurePaths,options.shouldOverwriteExisting)
end

[eddy,eddyDiagnosticsFile] = diagnosticsForSimulation(eddyFile,options.diagnosticsStride);
[control,controlDiagnosticsFile] = diagnosticsForSimulation(noEddyFile,options.diagnosticsStride);
diagnosticsFiles = [eddyDiagnosticsFile; controlDiagnosticsFile];
if ~isequal(eddy.time,control.time)
    error("AnalyzeEddyTideBudgets:TimeMismatch", ...
        "The eddy and no-eddy diagnostics must use identical output times.")
end

budget = struct(eddy=eddy,control=control,excess=subtractCases(eddy,control));
budget.normalization = struct(energy=eddy.state.geostrophic(1), ...
    enstrophy=eddy.state.quadraticEnstrophy(1));
if ~isfinite(budget.normalization.energy) || budget.normalization.energy <= 0
    error("AnalyzeEddyTideBudgets:InvalidEnergyNormalization", ...
        "The eddy simulation's initial geostrophic energy must be finite and positive.")
end
if ~isfinite(budget.normalization.enstrophy) || budget.normalization.enstrophy <= 0
    error("AnalyzeEddyTideBudgets:InvalidEnstrophyNormalization", ...
        "The eddy simulation's initial quadratic potential enstrophy must be finite and positive.")
end
budget.acceptance = evaluateAcceptance(eddy,control,options);
budget.figurePaths = strings(3,1);
if options.shouldRequireClosure && ~budget.acceptance.passed
    error("AnalyzeEddyTideBudgets:ClosureFailure", ...
        "One or more budgets exceed the requested closure tolerance.")
end

figureHandles = gobjects(3,1);
figureHandles(1) = createReservoirFigure(budget,options.figureVisible);
figureHandles(2) = createForcingFigure(budget,options.figureVisible);
figureHandles(3) = createTriadFigure(budget,options.figureVisible);
if options.shouldExport
    for iFigure = 1:numel(figureHandles)
        exportgraphics(figureHandles(iFigure),figurePaths(iFigure),Resolution=options.exportResolution)
    end
    budget.figurePaths = figurePaths;
end
end

function [result,diagnosticsFile] = diagnosticsForSimulation(simulationFile,diagnosticsStride)
diagnostics = WVDiagnostics(char(simulationFile));
cleanup = onCleanup(@()closeDiagnosticsFiles(diagnostics));
if ~isfile(diagnostics.diagpath)
    diagnostics.createDiagnosticsFile(stride=diagnosticsStride);
elseif diagnostics.t_diag(end) < diagnostics.t_wv(end)
    if numel(diagnostics.t_diag) < 2
        error("AnalyzeEddyTideBudgets:IncompleteDiagnostics", ...
            "The diagnostics file '%s' has only one incomplete record and cannot be appended safely.",diagnostics.diagpath)
    end
    diagnostics.createDiagnosticsFile(stride=diagnosticsStride);
end

energyReservoirs = [EnergyReservoir.geostrophic EnergyReservoir.mda ...
    EnergyReservoir.geostrophic_mda EnergyReservoir.wave EnergyReservoir.total];
[reservoirs,time] = diagnostics.quadraticEnergyOverTime(energyReservoirs=energyReservoirs);
[exactEnergy,exactEnergyTime] = diagnostics.exactEnergyOverTime;
[quadraticEnstrophy,quadraticEnstrophyTime] = diagnostics.quadraticEnstrophyOverTime;
[exactEnstrophy,exactEnstrophyTime] = diagnostics.exactEnstrophyOverTime;
[energyForcing,energyForcingTime] = diagnostics.quadraticEnergyFluxesOverTime(energyReservoirs=energyReservoirs);
[exactEnergyForcing,exactEnergyForcingTime] = diagnostics.exactEnergyFluxesOverTime;
[enstrophyForcing,enstrophyForcingTime] = diagnostics.quadraticEnstrophyFluxesOverTime;
[exactEnstrophyForcing,exactEnstrophyForcingTime] = diagnostics.exactEnstrophyFluxesOverTime;
[energyTriads,energyTriadTime] = diagnostics.quadraticEnergyTriadFluxesOverTime(energyReservoirs=energyReservoirs);
[enstrophyTriads,enstrophyTriadTime] = diagnostics.quadraticEnstrophyTriadFluxesOverTime;
time = time(:);
allTimes = {exactEnergyTime energyForcingTime exactEnergyForcingTime ...
    quadraticEnstrophyTime exactEnstrophyTime enstrophyForcingTime ...
    exactEnstrophyForcingTime energyTriadTime enstrophyTriadTime};
if ~all(cellfun(@(value)isequal(time,value(:)),allTimes))
    error("AnalyzeEddyTideBudgets:DiagnosticsTimeMismatch", ...
        "WVDiagnostics accessors returned inconsistent time axes.")
end

generation = diagnostics.wvt.forcingWithName("bottom wave generation");
tidalPeriod = 2*pi/generation.frequency;
if numel(time) > 1 && max(diff(time)) > tidalPeriod/4
    warning("AnalyzeEddyTideBudgets:UnderresolvedFluxCadence", ...
        "The output interval exceeds one quarter tidal period; fluxes may be temporally aliased.")
end

state = struct(geostrophic=reservoirs(1).energy(:),mda=reservoirs(2).energy(:), ...
    balanced=reservoirs(3).energy(:),wave=reservoirs(4).energy(:), ...
    quadraticTotal=reservoirs(5).energy(:),exactTotal=exactEnergy(:), ...
    quadraticEnstrophy=quadraticEnstrophy(:),exactEnstrophy=exactEnstrophy(:));
energyFields = ["te_g" "te_mda" "te_gmda" "te_wave" "te_quadratic"];
forcing = struct(energy=integrateFluxes(energyForcing,time,energyFields), ...
    exactEnergy=integrateFluxes(exactEnergyForcing,time,"te"), ...
    enstrophy=integrateFluxes(enstrophyForcing,time,"Z0"), ...
    exactEnstrophy=integrateFluxes(exactEnstrophyForcing,time,"Z0"));
triads = struct(energy=integrateFluxes(energyTriads,time,energyFields), ...
    enstrophy=integrateFluxes(enstrophyTriads,time,"Z0"));
closure = struct( ...
    geostrophic=closureForField(state.geostrophic,forcing.energy,"te_g",time), ...
    mda=closureForField(state.mda,forcing.energy,"te_mda",time), ...
    balanced=closureForField(state.balanced,forcing.energy,"te_gmda",time), ...
    wave=closureForField(state.wave,forcing.energy,"te_wave",time), ...
    quadraticTotal=closureForField(state.quadraticTotal,forcing.energy,"te_quadratic",time), ...
    exactTotal=closureForField(state.exactTotal,forcing.exactEnergy,"te",time), ...
    quadraticEnstrophy=closureForField(state.quadraticEnstrophy,forcing.enstrophy,"Z0",time), ...
    exactEnstrophy=closureForField(state.exactEnstrophy,forcing.exactEnstrophy,"Z0",time));
result = struct(time=time,timeDays=time/86400,tidalPeriod=tidalPeriod, ...
    state=state,forcing=forcing,triads=triads,closure=closure);
diagnosticsFile = string(diagnostics.diagpath);
clear cleanup
end

function fluxes = integrateFluxes(fluxes,time,fields)
for iFlux = 1:numel(fluxes)
    for field = fields
        rate = real(fluxes(iFlux).(field));
        rate = rate(:);
        cumulative = cumtrapz(time,rate);
        fluxes(iFlux).(field) = rate;
        fluxes(iFlux).(field+"Cumulative") = cumulative;
        fluxes(iFlux).(field+"Integrated") = cumulative(end);
    end
end
end

function result = closureForField(state,fluxes,field,time)
observed = state(end)-state(1);
integrated = 0;
gross = 0;
for iFlux = 1:numel(fluxes)
    integrated = integrated+fluxes(iFlux).(field+"Integrated");
    gross = gross+trapz(time,abs(fluxes(iFlux).(field)));
end
stateScale = max([abs(observed); max(state)-min(state); eps]);
result = struct(observed=observed,integrated=integrated,residual=integrated-observed, ...
    relativeStateResidual=abs(integrated-observed)/stateScale, ...
    relativeGrossResidual=abs(integrated-observed)/max(gross,eps),grossThroughput=gross);
end

function acceptance = evaluateAcceptance(eddy,control,options)
quadraticNames = ["geostrophic" "mda" "balanced" "wave" "quadraticTotal" "quadraticEnstrophy"];
exactNames = ["exactTotal" "exactEnstrophy"];
quadraticPassed = true;
exactPassed = true;
for caseValue = {eddy control}
    value = caseValue{1};
    for name = quadraticNames
        quadraticPassed = quadraticPassed && ...
            value.closure.(name).relativeStateResidual <= options.quadraticClosureTolerance;
    end
    for name = exactNames
        exactPassed = exactPassed && ...
            value.closure.(name).relativeGrossResidual <= options.exactClosureTolerance;
    end
end
acceptance = struct(passed=quadraticPassed && exactPassed, ...
    quadraticPassed=quadraticPassed,exactPassed=exactPassed, ...
    quadraticTolerance=options.quadraticClosureTolerance, ...
    exactTolerance=options.exactClosureTolerance);
end

function excess = subtractCases(eddy,control)
excess = struct(time=eddy.time,timeDays=eddy.timeDays,tidalPeriod=eddy.tidalPeriod);
stateNames = string(fieldnames(eddy.state));
for name = stateNames(:).'
    excess.state.(name) = eddy.state.(name)-control.state.(name);
end
excess.forcing = struct( ...
    energy=subtractFluxes(eddy.forcing.energy,control.forcing.energy), ...
    exactEnergy=subtractFluxes(eddy.forcing.exactEnergy,control.forcing.exactEnergy), ...
    enstrophy=subtractFluxes(eddy.forcing.enstrophy,control.forcing.enstrophy), ...
    exactEnstrophy=subtractFluxes(eddy.forcing.exactEnstrophy,control.forcing.exactEnstrophy));
excess.triads = struct( ...
    energy=subtractFluxes(eddy.triads.energy,control.triads.energy), ...
    enstrophy=subtractFluxes(eddy.triads.enstrophy,control.triads.enstrophy));
end

function difference = subtractFluxes(left,right)
if ~isequal(string({left.name}),string({right.name}))
    error("AnalyzeEddyTideBudgets:FluxMismatch", ...
        "The eddy and no-eddy diagnostics contain different forcing or triad names.")
end
difference = left;
for iFlux = 1:numel(left)
    fields = string(fieldnames(left(iFlux)));
    for field = fields(:).'
        if isnumeric(left(iFlux).(field))
            difference(iFlux).(field) = left(iFlux).(field)-right(iFlux).(field);
        end
    end
end
end

function figureHandle = createReservoirFigure(budget,visibility)
figureHandle = figure(Name="Eddy-tide reservoirs",Color="w",Visible=visibility);
layout = tiledlayout(figureHandle,2,2,TileSpacing="compact",Padding="compact");
plotReservoir(nexttile(layout,1),budget,"wave",budget.normalization.energy,"Wave energy")
plotReservoir(nexttile(layout,2),budget,"geostrophic",budget.normalization.energy,"Geostrophic energy")
plotReservoir(nexttile(layout,3),budget,"mda",budget.normalization.energy,"Mean-density-anomaly energy")
plotReservoir(nexttile(layout,4),budget,"quadraticEnstrophy",budget.normalization.enstrophy,"Quadratic potential enstrophy")
title(layout,"Generation-only eddy-tide reservoirs")
end

function plotReservoir(axesHandle,budget,field,normalization,titleText)
plot(axesHandle,budget.eddy.timeDays,budget.eddy.state.(field)/normalization,LineWidth=1.8,DisplayName="eddy")
hold(axesHandle,"on")
plot(axesHandle,budget.control.timeDays,budget.control.state.(field)/normalization,"--",LineWidth=1.8,DisplayName="no initial eddy")
plot(axesHandle,budget.excess.timeDays,budget.excess.state.(field)/normalization,":",LineWidth=1.8,DisplayName="eddy - control")
grid(axesHandle,"on")
xlabel(axesHandle,"time (days)")
ylabel(axesHandle,"normalized reservoir")
title(axesHandle,titleText)
legend(axesHandle,Location="best")
end

function figureHandle = createForcingFigure(budget,visibility)
figureHandle = figure(Name="Eddy-tide forcing budgets",Color="w",Visible=visibility);
layout = tiledlayout(figureHandle,2,3,TileSpacing="compact",Padding="compact");
cases = {budget.eddy budget.control};
caseNames = ["eddy" "no initial eddy"];
for iCase = 1:2
    plotCumulativeFlux(nexttile(layout,(iCase-1)*3+1),cases{iCase}, ...
        cases{iCase}.forcing.energy,"te_wave","wave",caseNames(iCase)+" wave")
    plotCumulativeFlux(nexttile(layout,(iCase-1)*3+2),cases{iCase}, ...
        cases{iCase}.forcing.energy,"te_gmda","balanced",caseNames(iCase)+" balanced")
    plotCumulativeFlux(nexttile(layout,(iCase-1)*3+3),cases{iCase}, ...
        cases{iCase}.forcing.enstrophy,"Z0","quadraticEnstrophy",caseNames(iCase)+" enstrophy")
end
title(layout,"Cumulative forcing-resolved budgets")
end

function figureHandle = createTriadFigure(budget,visibility)
figureHandle = figure(Name="Eddy-tide nonlinear triads",Color="w",Visible=visibility);
layout = tiledlayout(figureHandle,2,2,TileSpacing="compact",Padding="compact");
cases = {budget.eddy budget.control};
caseNames = ["eddy" "no initial eddy"];
for iCase = 1:2
    plotCumulativeFlux(nexttile(layout,(iCase-1)*2+1),cases{iCase}, ...
        cases{iCase}.triads.energy,"te_gmda","balanced",caseNames(iCase)+" balanced energy")
    plotCumulativeFlux(nexttile(layout,(iCase-1)*2+2),cases{iCase}, ...
        cases{iCase}.triads.enstrophy,"Z0","quadraticEnstrophy",caseNames(iCase)+" enstrophy")
end
title(layout,"Cumulative nonlinear-triad pathways")
end

function plotCumulativeFlux(axesHandle,caseValue,fluxes,field,stateField,titleText)
colors = lines(numel(fluxes));
hold(axesHandle,"on")
for iFlux = 1:numel(fluxes)
    plot(axesHandle,caseValue.timeDays,fluxes(iFlux).(field+"Cumulative"), ...
        Color=colors(iFlux,:),LineWidth=1.4,DisplayName=string(fluxes(iFlux).fancyName))
end
observed = caseValue.state.(stateField)-caseValue.state.(stateField)(1);
plot(axesHandle,caseValue.timeDays,observed,"k--",LineWidth=2,DisplayName="observed change")
grid(axesHandle,"on")
xlabel(axesHandle,"time (days)")
ylabel(axesHandle,"cumulative flux")
title(axesHandle,titleText)
legend(axesHandle,Location="best")
end

function validateExportPaths(paths,shouldOverwrite)
for path = paths(:).'
    directory = string(fileparts(path));
    if strlength(directory) > 0 && ~isfolder(directory)
        mkdir(directory)
    end
    if isfile(path) && ~shouldOverwrite
        error("AnalyzeEddyTideBudgets:ExportFileExists", ...
            "The figure '%s' already exists. Set shouldOverwriteExisting=true to replace it.",path)
    end
end
end

function closeDiagnosticsFiles(diagnostics)
if ~isempty(diagnostics.diagfile)
    diagnostics.diagfile.close();
end
if ~isempty(diagnostics.wvfile)
    diagnostics.wvfile.close();
end
end
