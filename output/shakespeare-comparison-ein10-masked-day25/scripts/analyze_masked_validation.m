function analyze_masked_validation
scriptDirectory = fileparts(mfilename("fullpath"));
runDirectory = fileparts(scriptDirectory);
checkpointDirectory = fullfile(runDirectory,"checkpoints","day-0025");
repositoryRoot = fileparts(fileparts(fileparts(scriptDirectory)));
workspaceRoot = fileparts(repositoryRoot);
addpath(repositoryRoot)
addpath(fullfile(repositoryRoot,"Examples"))
addpath(fullfile(workspaceRoot,"wave-vortex-model"))
addpath(fullfile(workspaceRoot,"wave-vortex-model-diagnostics"))

newEddyFile = fullfile(runDirectory,"eddy-tide-topographic-eddy-generation-only-Lxy500km-Nxy256-Ein10-masked-day25-hrms100m-lmin6km-seed2023.nc");
newControlFile = fullfile(runDirectory,"eddy-tide-topographic-no-eddy-generation-only-Lxy500km-Nxy256-Ein10-masked-day25-hrms100m-lmin6km-seed2023.nc");
oldEddyFile = fullfile(repositoryRoot,"output","eddy-tide-topographic-eddy-generation-only-Lxy500km-Nxy256-Ein10-hrms100m-lmin6km-seed2023.nc");
oldControlFile = fullfile(repositoryRoot,"output","eddy-tide-topographic-no-eddy-generation-only-Lxy500km-Nxy256-Ein10-hrms100m-lmin6km-seed2023.nc");
oldSummaryFile = fullfile(repositoryRoot,"output","checkpoints","day-0025","eddy-tide-Lxy500km-Nxy256-Ein10-lmin6km-day25-summary.mat");
prefix = "eddy-tide-Lxy500km-Nxy256-Ein10-masked-day25";
expectedTime = (0:6*3600:25*86400).';

validation = validateProductionPair(newEddyFile,newControlFile,oldEddyFile,oldControlFile,expectedTime);

[basicFigures,basicSummary] = AnalyzeEddyTideBasicFigures(newEddyFile,newControlFile, ...
    figureVisible="off",shouldExport=true,exportDirectory=checkpointDirectory, ...
    exportPrefix=prefix,exportResolution=500,shouldOverwriteExisting=true);
close(basicFigures)

[newBudgetFigures,newBudget,newDiagnosticsFiles] = AnalyzeEddyTideBudgets(newEddyFile,newControlFile, ...
    diagnosticsStride=1,figureVisible="off",shouldExport=false);
close(newBudgetFigures)
[oldBudgetFigures,oldBudget,oldDiagnosticsFiles] = AnalyzeEddyTideBudgets(oldEddyFile,oldControlFile, ...
    diagnosticsStride=1,figureVisible="off",shouldExport=false);
close(oldBudgetFigures)
assert(isequal(newBudget.eddy.time,expectedTime))
assert(isequal(newBudget.control.time,expectedTime))
assert(isequal(oldBudget.eddy.time(1:numel(expectedTime)),expectedTime))
assert(isequal(oldBudget.control.time(1:numel(expectedTime)),expectedTime))

oldSpatialPath = fullfile(checkpointDirectory,prefix+"-unmasked-eddy-days20-25-adaptive-damping.png");
newSpatialPath = fullfile(checkpointDirectory,prefix+"-masked-eddy-days20-25-adaptive-damping.png");
[oldSpatialFigure,oldSpatial,oldSpatialDiagnostics] = AnalyzeEddyTideAdaptiveDamping(oldEddyFile, ...
    timeRangeDays=[20 25],figureVisible="off",shouldExport=true,exportPath=oldSpatialPath, ...
    exportResolution=500,shouldOverwriteExisting=true);
close(oldSpatialFigure)
[newSpatialFigure,newSpatial,newSpatialDiagnostics] = AnalyzeEddyTideAdaptiveDamping(newEddyFile, ...
    timeRangeDays=[20 25],figureVisible="off",shouldExport=true,exportPath=newSpatialPath, ...
    exportResolution=500,shouldOverwriteExisting=true);
close(newSpatialFigure)

oldSeries = pairedSeries(oldBudget,numel(expectedTime));
newSeries = pairedSeries(newBudget,numel(expectedTime));
assert(isequal(oldSeries.time,newSeries.time))
energyComparisonPath = fullfile(checkpointDirectory,prefix+"-old-versus-masked-absolute-energy.png");
energyFigure = makeEnergyComparison(oldSeries,newSeries);
exportgraphics(energyFigure,energyComparisonPath,Resolution=500)
close(energyFigure)

budgetComparisonPath = fullfile(checkpointDirectory,prefix+"-old-versus-masked-generation-damping-budget.png");
budgetFigure = makeBudgetComparison(oldSeries,newSeries);
exportgraphics(budgetFigure,budgetComparisonPath,Resolution=500)
close(budgetFigure)

spatialComparisonPath = fullfile(checkpointDirectory,prefix+"-old-versus-masked-spatial-damping.png");
spatialFigure = makeSpatialComparison(oldSpatial,newSpatial);
exportgraphics(spatialFigure,spatialComparisonPath,Resolution=500)
close(spatialFigure)

loadedOldSummary = load(oldSummaryFile,"summary");
oldBasicSummary = loadedOldSummary.summary;
metrics = controlMetrics(oldSeries.control,newSeries.control,oldBasicSummary.noEddy,basicSummary.noEddy);
metrics.eddy = caseChangeMetrics(oldSeries.eddy,newSeries.eddy,oldBasicSummary.eddy,basicSummary.eddy);
metrics.spatialValidation = struct( ...
    oldMaximumRelativeFluxError=oldSpatial.validation.maximumRelativeFluxError, ...
    maskedMaximumRelativeFluxError=newSpatial.validation.maximumRelativeFluxError, ...
    oldTimeMeanDampingWork=oldSpatial.validation.timeMeanDiagnosticWork, ...
    maskedTimeMeanDampingWork=newSpatial.validation.timeMeanDiagnosticWork);

allDampingRates = [oldSeries.eddy.dampingRate; oldSeries.control.dampingRate; ...
    newSeries.eddy.dampingRate; newSeries.control.dampingRate];
dampingTolerance = 1e-14;
assert(all(allDampingRates <= dampingTolerance))
assert(all(isfinite(structNumericValues(oldSeries))))
assert(all(isfinite(structNumericValues(newSeries))))

diagnosticsFiles = struct(new=newDiagnosticsFiles,old=oldDiagnosticsFiles, ...
    oldSpatial=oldSpatialDiagnostics,newSpatial=newSpatialDiagnostics);
figurePaths = [basicSummary.figurePaths; energyComparisonPath; budgetComparisonPath; ...
    oldSpatialPath; newSpatialPath; spatialComparisonPath];
summaryPath = fullfile(checkpointDirectory,prefix+"-validation-summary.mat");
save(summaryPath,"validation","oldSeries","newSeries","metrics","basicSummary", ...
    "oldSpatial","newSpatial","diagnosticsFiles","figurePaths","-v7.3")

filesToReopen = [newEddyFile; newControlFile; newDiagnosticsFiles(:)];
for file = filesToReopen(:).'
    writableFile = NetCDFFile(char(file),shouldReadOnly=false);
    writableFile.close();
end

completion = struct(completedAt=string(datetime("now","TimeZone","local")), ...
    summaryPath=summaryPath,figurePaths=figurePaths,metrics=metrics);
save(fullfile(runDirectory,"analysis-success.mat"),"completion","-v7.3")
end

function validation = validateProductionPair(newEddyFile,newControlFile,oldEddyFile,oldControlFile,expectedTime)
newEddy = inspectFile(newEddyFile,expectedTime);
newControl = inspectFile(newControlFile,expectedTime);
oldEddy = inspectFile(oldEddyFile,[]);
oldControl = inspectFile(oldControlFile,[]);
assert(isequal(newEddy.time,newControl.time))
assert(isequal(newEddy.terrain,newControl.terrain))
assert(isequal(newEddy.terrain,oldEddy.terrain))
assert(isequal(newControl.terrain,oldControl.terrain))
assert(isequal(newEddy.domain,newControl.domain))
assert(isequal(newEddy.domain,oldEddy.domain))
assert(isequal(newControl.domain,oldControl.domain))
assert(isequal(newEddy.velocityAmplitude,newControl.velocityAmplitude))
assert(isequal(newEddy.velocityAmplitude,oldEddy.velocityAmplitude))
assert(isequal(newControl.velocityAmplitude,oldControl.velocityAmplitude))
validation = struct(newEddy=newEddy,newControl=newControl,oldEddy=oldEddy,oldControl=oldControl);
end

function result = inspectFile(path,expectedTime)
[wvt,ncfile] = WVTransform.waveVortexTransformFromFile(char(path),iTime=Inf,shouldReadOnly=true);
cleanup = onCleanup(@()ncfile.close());
time = ncfile.readVariables("wave-vortex/t");
if ~isempty(expectedTime)
    assert(isequal(time(:),expectedTime))
end
assert(all(isfinite([wvt.Ap(:); wvt.Am(:); wvt.A0(:)])))
generation = wvt.forcingWithName("bottom wave generation");
adaptiveDamping = wvt.forcingWithName("adaptive damping");
[mask,components] = generation.spectralGenerationMask();
expectedMask = (logical(wvt.waveComponent.maskAp) | logical(wvt.waveComponent.maskAm)) & adaptiveDamping.damp == 0;
if generation.shouldAvoidAdaptiveDamping
    assert(isequal(mask,expectedMask))
    assert(isequal(components.adaptiveDamping,adaptiveDamping.damp == 0))
    savedTime = wvt.t;
    wvt.t = max(savedTime,generation.rampDuration);
    [Fp,Fm] = generation.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
    dampingRegion = adaptiveDamping.damp ~= 0;
    assert(~any(Fp(dampingRegion),"all") && ~any(Fm(dampingRegion),"all"))
    assert(norm([Fp(~dampingRegion); Fm(~dampingRegion)]) > 0)
    wvt.t = savedTime;
end
result = struct(path=string(path),time=time(:),recordCount=numel(time), ...
    domain=[wvt.Lx wvt.Ly wvt.Lz wvt.Nx wvt.Ny wvt.Nz wvt.Nj], ...
    terrain=generation.topographicHeight,velocityAmplitude=generation.barotropicVelocityAmplitude, ...
    shouldAvoidAdaptiveDamping=generation.shouldAvoidAdaptiveDamping,maskCount=nnz(mask));
clear cleanup
end

function result = pairedSeries(budget,numberOfRecords)
result = struct(time=budget.eddy.time(1:numberOfRecords), ...
    timeDays=budget.eddy.timeDays(1:numberOfRecords), ...
    eddy=caseSeries(budget.eddy,numberOfRecords), ...
    control=caseSeries(budget.control,numberOfRecords));
end

function result = caseSeries(caseValue,numberOfRecords)
generation = forcingSeries(caseValue.forcing.exactEnergy,"bottom wave generation",numberOfRecords);
damping = forcingSeries(caseValue.forcing.exactEnergy,"adaptive damping",numberOfRecords);
result = struct( ...
    wave=caseValue.state.wave(1:numberOfRecords), ...
    geostrophic=caseValue.state.geostrophic(1:numberOfRecords), ...
    total=caseValue.state.quadraticTotal(1:numberOfRecords), ...
    generationRate=generation.rate,generationCumulative=generation.cumulative, ...
    dampingRate=damping.rate,dampingCumulative=damping.cumulative);
end

function result = forcingSeries(fluxes,fancyName,numberOfRecords)
names = string({fluxes.fancyName});
index = find(names == fancyName);
assert(isscalar(index))
result = struct(rate=fluxes(index).te(1:numberOfRecords), ...
    cumulative=fluxes(index).teCumulative(1:numberOfRecords));
end

function figureHandle = makeEnergyComparison(oldSeries,newSeries)
figureHandle = figure(Name="Old versus masked absolute energy",Color="w",Visible="off", ...
    Units="inches",Position=[1 1 12 7]);
layout = tiledlayout(figureHandle,2,3,TileSpacing="compact",Padding="compact");
caseNames = ["Eddy" "No initial eddy"];
oldCases = {oldSeries.eddy oldSeries.control};
newCases = {newSeries.eddy newSeries.control};
fields = ["wave" "geostrophic" "total"];
titles = ["Wave energy" "Geostrophic energy" "Total quadratic energy"];
for iCase = 1:2
    for iField = 1:3
        axesHandle = nexttile(layout,(iCase-1)*3+iField);
        plot(axesHandle,oldSeries.timeDays,oldCases{iCase}.(fields(iField)), ...
            Color=[0.3 0.3 0.3],LineWidth=1.8,DisplayName="unmasked")
        hold(axesHandle,"on")
        plot(axesHandle,newSeries.timeDays,newCases{iCase}.(fields(iField)), ...
            Color=[0 0.45 0.74],LineWidth=1.8,DisplayName="masked")
        grid(axesHandle,"on")
        xlim(axesHandle,[0 25])
        title(axesHandle,caseNames(iCase)+": "+titles(iField))
        if iCase == 2
            xlabel(axesHandle,"Time (days)")
        end
        if iField == 1
            ylabel(axesHandle,"Energy (m^3 s^{-2})")
        end
        legend(axesHandle,Location="best")
    end
end
title(layout,"Ein10 day-25 absolute energy: original versus damping-excluded generation")
end

function figureHandle = makeBudgetComparison(oldSeries,newSeries)
figureHandle = figure(Name="Old versus masked energy budgets",Color="w",Visible="off", ...
    Units="inches",Position=[1 1 10.5 7]);
layout = tiledlayout(figureHandle,2,2,TileSpacing="compact",Padding="compact");
caseNames = ["Eddy" "No initial eddy"];
oldCases = {oldSeries.eddy oldSeries.control};
newCases = {newSeries.eddy newSeries.control};
for iCase = 1:2
    axesHandle = nexttile(layout,(iCase-1)*2+1);
    plotFluxComparison(axesHandle,oldSeries.timeDays,oldCases{iCase}.generationCumulative, ...
        newCases{iCase}.generationCumulative,"Cumulative generation",caseNames(iCase))
    axesHandle = nexttile(layout,(iCase-1)*2+2);
    plotFluxComparison(axesHandle,oldSeries.timeDays,oldCases{iCase}.dampingCumulative, ...
        newCases{iCase}.dampingCumulative,"Cumulative adaptive damping",caseNames(iCase))
end
title(layout,"Ein10 day-25 exact energy work: original versus damping-excluded generation")
end

function plotFluxComparison(axesHandle,timeDays,oldValue,newValue,titleText,caseName)
plot(axesHandle,timeDays,oldValue,Color=[0.3 0.3 0.3],LineWidth=1.8,DisplayName="unmasked")
hold(axesHandle,"on")
plot(axesHandle,timeDays,newValue,Color=[0 0.45 0.74],LineWidth=1.8,DisplayName="masked")
grid(axesHandle,"on")
xlim(axesHandle,[0 25])
xlabel(axesHandle,"Time (days)")
ylabel(axesHandle,"Cumulative work (m^3 s^{-2})")
title(axesHandle,caseName+": "+titleText)
legend(axesHandle,Location="best")
end

function figureHandle = makeSpatialComparison(oldAnalysis,newAnalysis)
assert(isequal(oldAnalysis.radius,newAnalysis.radius))
assert(isequal(oldAnalysis.depth,newAnalysis.depth))
radius = oldAnalysis.radius/1e3;
depth = oldAnalysis.depth/1e3;
oldVelocity = 100*oldAnalysis.azimuthalVelocity;
newVelocity = 100*newAnalysis.azimuthalVelocity;
oldDamping = 1e6*oldAnalysis.adaptiveDampingEnergyRemovalRate;
newDamping = 1e6*newAnalysis.adaptiveDampingEnergyRemovalRate;
velocityLimit = max(abs([oldVelocity(:); newVelocity(:)]));
dampingLimit = max(abs([oldDamping(:); newDamping(:)]));
figureHandle = figure(Name="Old versus masked spatial damping",Color="w",Visible="off", ...
    Units="inches",Position=[1 1 12 8]);
layout = tiledlayout(figureHandle,2,2,TileSpacing="compact",Padding="compact");
plotSpatialPanel(nexttile(layout,1),radius,depth,oldVelocity,velocityLimit, ...
    "(a) Unmasked eddy-aligned velocity","cm s^{-1}")
plotSpatialPanel(nexttile(layout,2),radius,depth,newVelocity,velocityLimit, ...
    "(b) Masked eddy-aligned velocity","cm s^{-1}")
plotSpatialPanel(nexttile(layout,3),radius,depth,oldDamping,dampingLimit, ...
    "(c) Unmasked signed damping removal","\muW m^{-3}")
plotSpatialPanel(nexttile(layout,4),radius,depth,newDamping,dampingLimit, ...
    "(d) Masked signed damping removal","\muW m^{-3}")
title(layout,"Ein10 eddy, days 20-25 mean; common color limits within each row")
end

function plotSpatialPanel(axesHandle,radius,depth,value,colorLimit,titleText,colorbarText)
surface(axesHandle,radius,depth,zeros(size(value.')),value.',EdgeColor="none")
view(axesHandle,2)
colormap(axesHandle,blueWhiteRed(257))
clim(axesHandle,[-colorLimit colorLimit])
box(axesHandle,"on")
axesHandle.Layer = "top";
axesHandle.YDir = "normal";
xlim(axesHandle,[0 max(radius)])
ylim(axesHandle,[min(depth) max(depth)])
xlabel(axesHandle,"Radius (km)")
ylabel(axesHandle,"Depth (km)")
title(axesHandle,titleText)
colorbarHandle = colorbar(axesHandle);
colorbarHandle.Label.String = colorbarText;
end

function metrics = controlMetrics(oldCase,newCase,oldBasic,newBasic)
metrics = caseChangeMetrics(oldCase,newCase,oldBasic,newBasic);
metrics.finalGeostrophicEnergy = newCase.geostrophic(end);
metrics.finalWaveEnergy = newCase.wave(end);
metrics.geostrophicToWaveRatio = newCase.geostrophic(end)/newCase.wave(end);
metrics.geostrophicFractionOfWaveAndGeostrophic = newCase.geostrophic(end)/(newCase.geostrophic(end)+newCase.wave(end));
end

function metrics = caseChangeMetrics(oldCase,newCase,oldBasic,newBasic)
oldMaximumHorizontalSpeed = max(oldBasic.maximumHorizontalSpeed);
newMaximumHorizontalSpeed = max(newBasic.maximumHorizontalSpeed);
oldMaximumSpeed = oldBasic.maximumSpeed;
newMaximumSpeed = newBasic.maximumSpeed;
metrics = struct( ...
    finalGeostrophicEnergyChange=change(oldCase.geostrophic(end),newCase.geostrophic(end)), ...
    finalWaveEnergyChange=change(oldCase.wave(end),newCase.wave(end)), ...
    maximumHorizontalWaveSpeedChange=change(oldMaximumHorizontalSpeed,newMaximumHorizontalSpeed), ...
    maximumThreeDimensionalWaveSpeedChange=change(oldMaximumSpeed,newMaximumSpeed), ...
    cumulativeGenerationChange=change(oldCase.generationCumulative(end),newCase.generationCumulative(end)), ...
    cumulativeDampingChange=change(oldCase.dampingCumulative(end),newCase.dampingCumulative(end)));
end

function result = change(oldValue,newValue)
result = struct(old=oldValue,masked=newValue,absolute=newValue-oldValue, ...
    fractional=(newValue-oldValue)/max(abs(oldValue),realmin));
end

function values = structNumericValues(value)
values = [];
fields = fieldnames(value);
for iField = 1:numel(fields)
    fieldValue = value.(fields{iField});
    if isnumeric(fieldValue)
        values = [values; fieldValue(:)]; %#ok<AGROW>
    elseif isstruct(fieldValue)
        values = [values; structNumericValues(fieldValue)]; %#ok<AGROW>
    end
end
end

function map = blueWhiteRed(numberOfColors)
anchors = [0.230 0.299 0.754; 1 1 1; 0.706 0.016 0.150];
map = interp1([0 0.5 1],anchors,linspace(0,1,numberOfColors),"linear");
end
