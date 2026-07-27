function [benchmark,figureHandle] = GoffScatteringPerturbationBenchmark(inputFile,options)
% Test autonomous scattering across scaled Goff terrain amplitudes.
%
% The benchmark reads wave states from a generation-only simulation,
% evaluates the first-order bottom-energy identity at every selected time,
% and performs short scattering-only integrations from one representative
% state. The terrain realization is held fixed while its RMS height is
% scaled.
%
% - Topic: Examples
% - Declaration: [benchmark,figureHandle] = GoffScatteringPerturbationBenchmark(inputFile,options)
% - Parameter inputFile: restartable generation-only model output
% - Parameter options.rmsHeights: terrain RMS heights in meters
% - Parameter options.timeIndices: saved records used for instantaneous diagnostics
% - Parameter options.referenceTime: initial time for scattering-only integrations
% - Parameter options.numberOfIntegrationPeriods: number of M2 periods to integrate
% - Parameter options.relativeTolerance: adaptive integration tolerance
% - Parameter options.residualRatioTolerance: admissible RMS residual/work ratio
% - Parameter options.energyDriftTolerance: admissible fractional first-order energy drift
% - Parameter options.shouldMakeFigure: whether to create a summary figure
% - Parameter options.figureVisible: figure visibility, `"on"` or `"off"`
% - Returns benchmark: scaling, residual, drift, and acceptance diagnostics
% - Returns figureHandle: summary figure, or an empty graphics array
arguments (Input)
    inputFile (1,1) string {mustBeFile}
    options.rmsHeights (1,:) double {mustBePositive,mustBeFinite} = [100 50 25 12.5 6.25 3.125]
    options.timeIndices (1,:) double = Inf
    options.referenceTime (1,1) double {mustBeNonnegative,mustBeFinite} = 5*86400
    options.numberOfIntegrationPeriods (1,1) double {mustBeNonnegative,mustBeFinite} = 5
    options.relativeTolerance (1,1) double {mustBePositive,mustBeFinite} = 1e-6
    options.residualRatioTolerance (1,1) double {mustBePositive,mustBeFinite} = 0.10
    options.energyDriftTolerance (1,1) double {mustBePositive,mustBeFinite} = 0.05
    options.shouldMakeFigure (1,1) logical = true
    options.figureVisible (1,1) string {mustBeMember(options.figureVisible,["on" "off"])} = "on"
end

if numel(options.rmsHeights) < 3 || numel(unique(options.rmsHeights)) < 3
    error("GoffScatteringPerturbationBenchmark:InsufficientTerrainScales", ...
        "At least three distinct RMS heights are required to estimate scaling exponents.")
end

[wvt,ncfile] = WVTransform.waveVortexTransformFromFile(char(inputFile),iTime=1,shouldReadOnly=true);
fileCleanup = onCleanup(@()ncfile.close());
time = ncfile.readVariables("wave-vortex/t");
generation = wvt.forcingWithName("bottom wave generation");
if isempty(generation)
    error("GoffScatteringPerturbationBenchmark:MissingGenerationForcing", ...
        "The input file must contain bottom wave generation.")
end
tidalPeriod = 2*pi/generation.frequency;
terrainShape = generation.topographicHeight/rms(generation.topographicHeight,"all");
if isscalar(options.timeIndices) && isinf(options.timeIndices)
    timeIndices = 1:numel(time);
else
    timeIndices = options.timeIndices;
end
if any(timeIndices < 1 | timeIndices > numel(time) | timeIndices ~= round(timeIndices))
    error("GoffScatteringPerturbationBenchmark:InvalidTimeIndices", ...
        "timeIndices must identify saved model records.")
end

numberOfHeights = numel(options.rmsHeights);
numberOfTimes = numel(timeIndices);
flatWork = zeros(numberOfTimes,numberOfHeights);
freeBottomCorrection = zeros(numberOfTimes,numberOfHeights);
autonomousBottomCorrection = zeros(numberOfTimes,numberOfHeights);
autonomousResidual = zeros(numberOfTimes,numberOfHeights);
scatteringForcings = cell(numberOfHeights,1);
for iHeight = 1:numberOfHeights
    scatteringForcings{iHeight} = WVBottomWaveScatteringForcing(wvt, ...
        topographicHeight=options.rmsHeights(iHeight)*terrainShape, ...
        name="benchmark scattering");
end

for iTime = 1:numberOfTimes
    wvt.initFromNetCDFFile(ncfile,iTime=timeIndices(iTime));
    Ap = wvt.Ap;
    Am = wvt.Am;
    for iHeight = 1:numberOfHeights
        scattering = scatteringForcings{iHeight};
        [~,bottom] = scattering.bottomVelocityFromWaveState(wvt);
        [Fp,Fm] = scattering.addSpectralForcing(wvt, ...
            complex(zeros(size(Ap))),complex(zeros(size(Am))), ...
            complex(zeros(size(wvt.A0))));
        flatWork(iTime,iHeight) = 2*sum(wvt.Apm_TE_factor.* ...
            real(Fp.*conj(Ap)+Fm.*conj(Am)),"all");

        wvt.Ap = wvt.iOmega.*Ap;
        wvt.Am = -wvt.iOmega.*Am;
        [~,freeTendency] = scattering.bottomVelocityFromWaveState(wvt);
        freeBottomCorrection(iTime,iHeight) = mean(scattering.topographicHeight.* ...
            (bottom.u.*freeTendency.u+bottom.v.*freeTendency.v),"all");

        wvt.Ap = Fp;
        wvt.Am = Fm;
        [~,scatteringTendency] = scattering.bottomVelocityFromWaveState(wvt);
        autonomousBottomCorrection(iTime,iHeight) = mean(scattering.topographicHeight.* ...
            (bottom.u.*scatteringTendency.u+bottom.v.*scatteringTendency.v),"all");
        autonomousResidual(iTime,iHeight) = flatWork(iTime,iHeight) ...
            -freeBottomCorrection(iTime,iHeight)-autonomousBottomCorrection(iTime,iHeight);
        wvt.Ap = Ap;
        wvt.Am = Am;
    end
end

workMagnitude = rms(flatWork,1);
residualMagnitude = rms(autonomousResidual,1);
residualRatio = residualMagnitude./max(workMagnitude,eps);
workExponent = scalingExponent(options.rmsHeights,workMagnitude);
residualExponent = scalingExponent(options.rmsHeights,residualMagnitude);
strictIdentityError = max(abs(flatWork-freeBottomCorrection),[],"all") ...
    /max(max(abs(flatWork),[],"all"),eps);

[~,referenceIndex] = min(abs(time-options.referenceTime));
energyDrift = zeros(1,numberOfHeights);
fractionalEnergyDrift = zeros(1,numberOfHeights);
if options.numberOfIntegrationPeriods > 0
    for iHeight = 1:numberOfHeights
        [integrationTransform,integrationFile] = WVTransform.waveVortexTransformFromFile( ...
            char(inputFile),iTime=referenceIndex,shouldReadOnly=true);
        integrationCleanup = onCleanup(@()integrationFile.close());
        integrationGeneration = integrationTransform.forcingWithName("bottom wave generation");
        integrationTerrain = options.rmsHeights(iHeight) ...
            *integrationGeneration.topographicHeight/rms(integrationGeneration.topographicHeight,"all");
        scattering = WVBottomWaveScatteringForcing(integrationTransform, ...
            topographicHeight=integrationTerrain,name="benchmark scattering");
        integrationTransform.removeAllForcing();
        integrationTransform.addForcing(scattering);
        initialEnergy = firstOrderWaveEnergy(integrationTransform,scattering);
        model = WVModel(integrationTransform);
        model.setupIntegrator(integratorType="adaptive", ...
            absTolerance=options.relativeTolerance,relTolerance=options.relativeTolerance);
        model.integrateToTime(integrationTransform.t+options.numberOfIntegrationPeriods*tidalPeriod, ...
            shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
        finalEnergy = firstOrderWaveEnergy(integrationTransform,scattering);
        energyDrift(iHeight) = finalEnergy-initialEnergy;
        fractionalEnergyDrift(iHeight) = abs(energyDrift(iHeight))/max(abs(initialEnergy),eps);
        clear integrationCleanup
    end
end
if all(energyDrift ~= 0)
    driftExponent = scalingExponent(options.rmsHeights,abs(energyDrift));
else
    driftExponent = NaN;
end

isAdmissible = residualRatio <= options.residualRatioTolerance ...
    & fractionalEnergyDrift <= options.energyDriftTolerance;
if any(isAdmissible)
    maximumAdmissibleHeight = max(options.rmsHeights(isAdmissible));
else
    maximumAdmissibleHeight = NaN;
end
benchmark = struct(inputFile=inputFile,time=time(timeIndices),timeIndices=timeIndices(:), ...
    referenceTime=time(referenceIndex),tidalPeriod=tidalPeriod,rmsHeights=options.rmsHeights, ...
    flatWork=flatWork,freeBottomCorrection=freeBottomCorrection, ...
    autonomousBottomCorrection=autonomousBottomCorrection, ...
    autonomousResidual=autonomousResidual,workMagnitude=workMagnitude, ...
    residualMagnitude=residualMagnitude,residualRatio=residualRatio, ...
    energyDrift=energyDrift,fractionalEnergyDrift=fractionalEnergyDrift, ...
    workExponent=workExponent,residualExponent=residualExponent, ...
    driftExponent=driftExponent,strictIdentityError=strictIdentityError, ...
    isAdmissible=isAdmissible,maximumAdmissibleHeight=maximumAdmissibleHeight, ...
    residualRatioTolerance=options.residualRatioTolerance, ...
    energyDriftTolerance=options.energyDriftTolerance);

figureHandle = gobjects(0);
if options.shouldMakeFigure
    figureHandle = createFigure(benchmark,options.figureVisible);
end
clear fileCleanup
end

function exponent = scalingExponent(heights,values)
fit = polyfit(log(heights(:)),log(values(:)),1);
exponent = fit(1);
end

function energy = firstOrderWaveEnergy(wvt,scattering)
[~,bottom] = scattering.bottomVelocityFromWaveState(wvt);
flatEnergy = sum(wvt.Apm_TE_factor.*(abs(wvt.Ap).^2+abs(wvt.Am).^2),"all");
bottomCorrection = 0.5*mean(scattering.topographicHeight.* ...
    (bottom.u.^2+bottom.v.^2),"all");
energy = flatEnergy-bottomCorrection;
end

function figureHandle = createFigure(benchmark,visibility)
figureHandle = figure(Name="Goff scattering perturbation benchmark", ...
    Color="w",Visible=visibility);
layout = tiledlayout(figureHandle,1,2,TileSpacing="compact",Padding="compact");
axesHandle = nexttile(layout,1);
loglog(axesHandle,benchmark.rmsHeights,benchmark.workMagnitude,"o-", ...
    LineWidth=1.8,DisplayName=sprintf("first-order work, p=%.2f",benchmark.workExponent))
hold(axesHandle,"on")
loglog(axesHandle,benchmark.rmsHeights,benchmark.residualMagnitude,"s-", ...
    LineWidth=1.8,DisplayName=sprintf("autonomous residual, p=%.2f",benchmark.residualExponent))
grid(axesHandle,"on")
xlabel(axesHandle,"terrain RMS height (m)")
ylabel(axesHandle,"RMS energy tendency")
legend(axesHandle,Location="best")

axesHandle = nexttile(layout,2);
loglog(axesHandle,benchmark.rmsHeights,benchmark.residualRatio,"o-", ...
    LineWidth=1.8,DisplayName="residual / work")
hold(axesHandle,"on")
loglog(axesHandle,benchmark.rmsHeights,benchmark.fractionalEnergyDrift,"s-", ...
    LineWidth=1.8,DisplayName="five-period energy drift")
yline(axesHandle,benchmark.residualRatioTolerance,":",DisplayName="residual threshold")
yline(axesHandle,benchmark.energyDriftTolerance,"--",DisplayName="drift threshold")
grid(axesHandle,"on")
xlabel(axesHandle,"terrain RMS height (m)")
ylabel(axesHandle,"fraction")
legend(axesHandle,Location="best")
title(layout,"Autonomous first-order scattering validity")
end
