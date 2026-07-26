function benchmark = SinusoidalRidgeWaveGenerationBenchmark(options)
% Validate prescribed bottom wave generation against its exact modal response.
%
% A uniform M2 current flows over a sinusoidal ridge for one tidal period.
% Adaptive WaveVortexModel integrations are compared with the exact
% interaction-representation coefficients. By default, the example also
% shows the convergence, energy budget, modal distribution, and generated
% wave fields.
%
% - Topic: Examples
% - Declaration: benchmark = SinusoidalRidgeWaveGenerationBenchmark(options)
% - Parameter options.relativeTolerances: adaptive relative and absolute tolerances
% - Parameter options.shouldAntialias: antialias setting for the transform
% - Parameter options.resolution: transform resolution $$[N_x,N_y,N_z]$$
% - Parameter options.shouldMakeFigures: whether to create the two diagnostic figures
% - Returns benchmark: exact-reference, adaptive-integration, energy, modal, and figure diagnostics
arguments
    options.relativeTolerances (1,:) double {mustBePositive,mustBeFinite} = [1e-6 1e-8 1e-10]
    options.shouldAntialias (1,1) logical = true
    options.resolution (1,3) double {mustBeInteger,mustBePositive} = [32 4 17]
    options.shouldMakeFigures (1,1) logical = true
end

domainSize = [20e3 20e3 2e3];
N2 = 2e-5;
latitude = 45;
terrainAmplitude = 50;
barotropicVelocityAmplitude = [0.05; 0];
frequency = 2*pi/(12.4206012*3600);
tidalPeriod = 2*pi/frequency;

wvt = WVTransformBoussinesq(domainSize,options.resolution,N2=@(z) N2*ones(size(z)),latitude=latitude,shouldAntialias=options.shouldAntialias);
wvt.t0 = 0;
wvt.t = 0;
x = reshape(wvt.x,[],1);
topographicHeight = terrainAmplitude*cos(2*pi*x/wvt.Lx).*ones(1,wvt.Ny);
forcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=topographicHeight,barotropicVelocityAmplitude=barotropicVelocityAmplitude,frequency=frequency,rampDuration=0,startTime=0);
[responsePlusX,responseMinusX] = unitCurrentResponse(wvt,topographicHeight,frequency,[1; 0]);
[responsePlusY,responseMinusY] = unitCurrentResponse(wvt,topographicHeight,frequency,[0; 1]);

numberOfDiagnosticTimes = 401;
time = linspace(0,tidalPeriod,numberOfDiagnosticTimes);
plusEnergy = zeros(size(time));
minusEnergy = zeros(size(time));
sourcePower = zeros(size(time));
for iTime = 1:numberOfDiagnosticTimes
    [ApReference,AmReference] = exactCoefficients(time(iTime),wvt.Omega,frequency,barotropicVelocityAmplitude,responsePlusX,responsePlusY,responseMinusX,responseMinusY);
    plusEnergy(iTime) = sum(wvt.Apm_TE_factor(:).*abs(ApReference(:)).^2);
    minusEnergy(iTime) = sum(wvt.Apm_TE_factor(:).*abs(AmReference(:)).^2);
    sourcePower(iTime) = exactBottomPower(time(iTime),wvt,forcing,frequency,barotropicVelocityAmplitude,responsePlusX,responsePlusY,responseMinusX,responseMinusY);
end
referenceEnergy = plusEnergy+minusEnergy;
cumulativeBottomWork = cumtrapz(time,sourcePower);
[ApReference,AmReference] = exactCoefficients(tidalPeriod,wvt.Omega,frequency,barotropicVelocityAmplitude,responsePlusX,responsePlusY,responseMinusX,responseMinusY);
referenceFinalEnergy = referenceEnergy(end);
quadratureAbsoluteTolerance = max(1e-15,1e-13*referenceFinalEnergy);
integratedBottomWork = integral(@(t) exactBottomPower(t,wvt,forcing,frequency,barotropicVelocityAmplitude,responsePlusX,responsePlusY,responseMinusX,responseMinusY),0,tidalPeriod,ArrayValued=true,RelTol=1e-12,AbsTol=quadratureAbsoluteTolerance);
referenceEnergyWorkError = abs(integratedBottomWork-referenceFinalEnergy)/max(referenceFinalEnergy,eps);

relativeTolerances = options.relativeTolerances;
coefficientError = zeros(size(relativeTolerances));
numericalEnergy = zeros(size(relativeTolerances));
energyWorkError = zeros(size(relativeTolerances));
balancedEnergyFraction = zeros(size(relativeTolerances));
numberOfFluxComputations = zeros(size(relativeTolerances));
finalPlusEnergy = zeros(size(relativeTolerances));
finalMinusEnergy = zeros(size(relativeTolerances));
finalAp = cell(size(relativeTolerances));
finalAm = cell(size(relativeTolerances));
finalA0 = cell(size(relativeTolerances));
referenceNorm = sqrt(referenceFinalEnergy);

warningState = warning;
warningCleanup = onCleanup(@()warning(warningState));
warning("off","all")
for iRun = 1:numel(relativeTolerances)
    wvt.Ap(:) = 0;
    wvt.Am(:) = 0;
    wvt.A0(:) = 0;
    wvt.t0 = 0;
    wvt.t = 0;
    wvt.removeAllForcing();
    wvt.addForcing(forcing);
    model = WVModel(wvt);
    model.setupIntegrator(integratorType="adaptive",absTolerance=relativeTolerances(iRun),relTolerance=relativeTolerances(iRun));
    model.integrateToTime(tidalPeriod,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);

    differenceEnergy = sum(wvt.Apm_TE_factor(:).*(abs(wvt.Ap(:)-ApReference(:)).^2+abs(wvt.Am(:)-AmReference(:)).^2));
    coefficientError(iRun) = sqrt(differenceEnergy)/max(referenceNorm,eps);
    finalPlusEnergy(iRun) = sum(wvt.Apm_TE_factor(:).*abs(wvt.Ap(:)).^2);
    finalMinusEnergy(iRun) = sum(wvt.Apm_TE_factor(:).*abs(wvt.Am(:)).^2);
    numericalEnergy(iRun) = finalPlusEnergy(iRun)+finalMinusEnergy(iRun);
    balancedEnergy = sum(wvt.A0_TE_factor(:).*abs(wvt.A0(:)).^2);
    balancedEnergyFraction(iRun) = balancedEnergy/max(numericalEnergy(iRun),eps);
    energyWorkError(iRun) = abs(numericalEnergy(iRun)-integratedBottomWork)/max(abs(integratedBottomWork),eps);
    numberOfFluxComputations(iRun) = model.nFluxComputations;
    finalAp{iRun} = wvt.Ap;
    finalAm{iRun} = wvt.Am;
    finalA0{iRun} = wvt.A0;
end
clear warningCleanup

wvt.Ap = finalAp{end};
wvt.Am = finalAm{end};
wvt.A0 = finalA0{end};
wvt.t = tidalPeriod;
qgpvNorm = norm(wvt.qgpv(:));
waveVelocityScale = sqrt(2*numericalEnergy(end)/wvt.Lz);
normalizedQGPV = qgpvNorm/max(abs(wvt.f)*waveVelocityScale/wvt.Lz,eps);

[sourceFp,sourceFm,~] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
terrainWavenumber = 2*pi/wvt.Lx;
expectedHorizontalSupport = abs(abs(wvt.K)-terrainWavenumber) <= 100*eps(terrainWavenumber) & abs(wvt.L) <= 100*eps(terrainWavenumber);
sourceNorm = norm([sourceFp(:); sourceFm(:)]);
sourceSupportError = norm([sourceFp(~expectedHorizontalSupport); sourceFm(~expectedHorizontalSupport)])/max(sourceNorm,eps);

waveModeNumbers = unique(wvt.J(wvt.waveComponent.maskAp | wvt.waveComponent.maskAm));
plusEnergyByMode = zeros(size(waveModeNumbers));
minusEnergyByMode = zeros(size(waveModeNumbers));
for iMode = 1:numel(waveModeNumbers)
    modeMask = wvt.J == waveModeNumbers(iMode);
    plusEnergyByMode(iMode) = sum(wvt.Apm_TE_factor(modeMask).*abs(ApReference(modeMask)).^2);
    minusEnergyByMode(iMode) = sum(wvt.Apm_TE_factor(modeMask).*abs(AmReference(modeMask)).^2);
end

[~,iSnapshot] = max(referenceEnergy);
snapshotTime = time(iSnapshot);
[snapshotAp,snapshotAm] = exactCoefficients(snapshotTime,wvt.Omega,frequency,barotropicVelocityAmplitude,responsePlusX,responsePlusY,responseMinusX,responseMinusY);
wvt.Ap = snapshotAp;
wvt.Am = snapshotAm;
wvt.A0(:) = 0;
wvt.t = snapshotTime;
[~,iBottom] = min(wvt.z);
snapshotBottomVelocity = forcing.bottomVelocityAtTime(snapshotTime);
snapshotKinematicPressure = wvt.p(:,:,iBottom)/wvt.rho0;
snapshotWorkDensity = snapshotKinematicPressure.*snapshotBottomVelocity;
snapshotW = wvt.w;
snapshotEta = wvt.eta;

figureHandles = gobjects(0);
if options.shouldMakeFigures
    figureHandles = createFigures(wvt,relativeTolerances,coefficientError,time,referenceEnergy,cumulativeBottomWork,numericalEnergy,plusEnergy,minusEnergy,finalPlusEnergy,finalMinusEnergy,waveModeNumbers,plusEnergyByMode,minusEnergyByMode,topographicHeight,snapshotTime,snapshotWorkDensity,snapshotW,snapshotEta);
end

configuration = struct(domainSize=domainSize,resolution=options.resolution,N2=N2,latitude=latitude,terrainAmplitude=terrainAmplitude,barotropicVelocityAmplitude=barotropicVelocityAmplitude,frequency=frequency,tidalPeriod=tidalPeriod,shouldAntialias=options.shouldAntialias);
benchmark = struct(configuration=configuration,time=time,relativeTolerances=relativeTolerances,coefficientError=coefficientError,numberOfFluxComputations=numberOfFluxComputations,referenceEnergy=referenceEnergy,numericalEnergy=numericalEnergy,sourcePower=sourcePower,cumulativeBottomWork=cumulativeBottomWork,integratedBottomWork=integratedBottomWork,referenceEnergyWorkError=referenceEnergyWorkError,energyWorkError=energyWorkError,plusEnergy=plusEnergy,minusEnergy=minusEnergy,finalPlusEnergy=finalPlusEnergy,finalMinusEnergy=finalMinusEnergy,waveModeNumbers=waveModeNumbers,plusEnergyByMode=plusEnergyByMode,minusEnergyByMode=minusEnergyByMode,balancedEnergyFraction=balancedEnergyFraction,normalizedQGPV=normalizedQGPV,sourceSupportError=sourceSupportError,snapshotTime=snapshotTime,figureHandles=figureHandles);

summary = table(relativeTolerances(:),coefficientError(:),energyWorkError(:),numberOfFluxComputations(:),VariableNames=["Tolerance","CoefficientError","EnergyWorkError","FluxEvaluations"]);
disp(summary)
end

function [responsePlus,responseMinus] = unitCurrentResponse(wvt,topographicHeight,frequency,velocity)
unitForcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=topographicHeight,barotropicVelocityAmplitude=velocity,frequency=frequency,rampDuration=0,startTime=0);
wvt.t = 0;
[responsePlus,responseMinus] = unitForcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
end

function [Ap,Am] = exactCoefficients(t,Omega,frequency,velocityAmplitude,responsePlusX,responsePlusY,responseMinusX,responseMinusY)
plusResponse = cat(3,responsePlusX,responsePlusY);
minusResponse = cat(3,responseMinusX,responseMinusY);
Ap = complex(zeros(size(Omega)));
Am = complex(zeros(size(Omega)));
for iDirection = 1:2
    amplitude = velocityAmplitude(iDirection);
    Ap = Ap+0.5*plusResponse(:,:,iDirection).*(amplitude*oscillatoryIntegral(Omega+frequency,t)+conj(amplitude)*oscillatoryIntegral(Omega-frequency,t));
    Am = Am+0.5*minusResponse(:,:,iDirection).*(amplitude*oscillatoryIntegral(frequency-Omega,t)+conj(amplitude)*oscillatoryIntegral(-frequency-Omega,t));
end
end

function value = oscillatoryIntegral(nu,t)
argument = nu*t;
value = t*(1-0.5i*argument-argument.^2/6);
regular = abs(argument) > sqrt(eps);
value(regular) = -expm1(-1i*argument(regular))./(1i*nu(regular));
end

function power = exactBottomPower(t,wvt,forcing,frequency,velocityAmplitude,responsePlusX,responsePlusY,responseMinusX,responseMinusY)
power = zeros(size(t));
[~,iBottom] = min(wvt.z);
for iTime = 1:numel(t)
    [Ap,Am] = exactCoefficients(t(iTime),wvt.Omega,frequency,velocityAmplitude,responsePlusX,responsePlusY,responseMinusX,responseMinusY);
    wvt.Ap = Ap;
    wvt.Am = Am;
    wvt.A0(:) = 0;
    wvt.t = t(iTime);
    kinematicPressure = wvt.p(:,:,iBottom)/wvt.rho0;
    power(iTime) = mean(kinematicPressure.*forcing.bottomVelocityAtTime(t(iTime)),"all");
end
end

function figureHandles = createFigures(wvt,tolerances,coefficientError,time,energy,cumulativeWork,numericalEnergy,plusEnergy,minusEnergy,numericalPlusEnergy,numericalMinusEnergy,modeNumbers,plusEnergyByMode,minusEnergyByMode,topographicHeight,snapshotTime,workDensity,w,eta)
tidalPeriod = time(end);
figureHandles(1) = figure(Name="Exact sinusoidal-ridge validation",Color="w",Position=[100 100 1250 850]);
layout = tiledlayout(figureHandles(1),2,2,TileSpacing="compact",Padding="compact");

axesHandle = nexttile(layout,1);
loglog(axesHandle,tolerances,coefficientError,"o-",LineWidth=1.5,MarkerFaceColor=[0.1 0.4 0.8])
grid(axesHandle,"on")
xlabel(axesHandle,"adaptive tolerance")
ylabel(axesHandle,"energy-norm coefficient error")
title(axesHandle,"Exact coefficient convergence")

axesHandle = nexttile(layout,2);
plot(axesHandle,time/tidalPeriod,energy,LineWidth=1.6)
hold(axesHandle,"on")
plot(axesHandle,time/tidalPeriod,cumulativeWork,"--",LineWidth=1.4)
scatter(axesHandle,ones(size(numericalEnergy)),numericalEnergy,32,"k","filled",HandleVisibility="off")
grid(axesHandle,"on")
xlabel(axesHandle,"time / M2 period")
ylabel(axesHandle,"energy or work (m^3 s^{-2})")
legend(axesHandle,"exact wave energy","integrated bottom work",Location="best")
title(axesHandle,"Energy supplied by bottom pressure work")

axesHandle = nexttile(layout,3);
plot(axesHandle,time/tidalPeriod,plusEnergy,LineWidth=1.5)
hold(axesHandle,"on")
plot(axesHandle,time/tidalPeriod,minusEnergy,"--",LineWidth=1.5)
scatter(axesHandle,ones(size(numericalPlusEnergy)),numericalPlusEnergy,28,[0 0.4470 0.7410],"filled",HandleVisibility="off")
scatter(axesHandle,ones(size(numericalMinusEnergy)),numericalMinusEnergy,28,[0.8500 0.3250 0.0980],"filled",HandleVisibility="off")
grid(axesHandle,"on")
xlabel(axesHandle,"time / M2 period")
ylabel(axesHandle,"branch energy (m^3 s^{-2})")
legend(axesHandle,"A_+","A_-",Location="best")
title(axesHandle,"Wave-branch response")

axesHandle = nexttile(layout,4);
bar(axesHandle,modeNumbers,[plusEnergyByMode(:) minusEnergyByMode(:)],"stacked")
grid(axesHandle,"on")
xlabel(axesHandle,"vertical mode j")
ylabel(axesHandle,"final energy (m^3 s^{-2})")
legend(axesHandle,"A_+","A_-",Location="best")
title(axesHandle,"Vertical-mode distribution")
title(layout,"Exact sinusoidal-ridge wave generation")

figureHandles(2) = figure(Name="Exact sinusoidal-ridge wave fields",Color="w",Position=[150 100 1250 850]);
layout = tiledlayout(figureHandles(2),2,2,TileSpacing="compact",Padding="compact");

axesHandle = nexttile(layout,1);
plot(axesHandle,wvt.x/1e3,topographicHeight(:,1),LineWidth=1.6)
grid(axesHandle,"on")
xlabel(axesHandle,"x (km)")
ylabel(axesHandle,"h (m)")
title(axesHandle,"Bottom topography")

axesHandle = nexttile(layout,2);
plot(axesHandle,wvt.x/1e3,workDensity(:,1),LineWidth=1.6)
grid(axesHandle,"on")
xlabel(axesHandle,"x (km)")
ylabel(axesHandle,"(p_d/\rho_0)g_b (m^3 s^{-3})")
title(axesHandle,"Local bottom work density")

axesHandle = nexttile(layout,3);
imagesc(axesHandle,wvt.x/1e3,wvt.z/1e3,squeeze(w(:,1,:)).')
axis(axesHandle,"xy")
xlabel(axesHandle,"x (km)")
ylabel(axesHandle,"z (km)")
title(axesHandle,"Vertical velocity w (m s^{-1})")
colorbar(axesHandle)
applySymmetricColorLimits(axesHandle,w)

axesHandle = nexttile(layout,4);
imagesc(axesHandle,wvt.x/1e3,wvt.z/1e3,squeeze(eta(:,1,:)).')
axis(axesHandle,"xy")
xlabel(axesHandle,"x (km)")
ylabel(axesHandle,"z (km)")
title(axesHandle,"Displacement \eta (m)")
colorbar(axesHandle)
applySymmetricColorLimits(axesHandle,eta)
colormap(figureHandles(2),"turbo")
title(layout,sprintf("Generated fields at t/T_{M2} = %.3f",snapshotTime/tidalPeriod))
end

function applySymmetricColorLimits(axesHandle,field)
limit = max(abs(field),[],"all");
if limit == 0
    limit = 1;
end
clim(axesHandle,[-limit limit])
end
