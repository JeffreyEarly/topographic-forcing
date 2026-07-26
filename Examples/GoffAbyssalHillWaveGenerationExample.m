function result = GoffAbyssalHillWaveGenerationExample(options)
% Generate waves from prescribed M2 flow over Goff abyssal hills.
%
% The example uses an exponential stationary stratification profile and a
% deterministic two-dimensional Goff terrain realization. It integrates
% the linear wave-only response with the adaptive WaveVortexModel
% integrator and compares wave-energy growth with bottom pressure work.
%
% - Topic: Examples
% - Declaration: result = GoffAbyssalHillWaveGenerationExample(options)
% - Parameter options.resolution: transform resolution $$[N_x,N_y,N_z]$$
% - Parameter options.minimumWavelength: shortest retained terrain wavelength in meters
% - Parameter options.numberOfOutputTimes: number of diagnostic output times
% - Parameter options.relativeTolerance: adaptive relative and absolute tolerance
% - Parameter options.shouldAntialias: antialias setting for the transform
% - Parameter options.shouldMakeFigures: whether to create the diagnostic figure
% - Returns result: configuration, terrain, spectral, wave, energy, and figure diagnostics
arguments
    options.resolution (1,3) double {mustBeInteger,mustBePositive} = [64 64 17]
    options.minimumWavelength (1,1) double {mustBePositive,mustBeFinite} = 40e3
    options.numberOfOutputTimes (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.numberOfOutputTimes,3)} = 41
    options.relativeTolerance (1,1) double {mustBePositive,mustBeFinite} = 1e-8
    options.shouldAntialias (1,1) logical = true
    options.shouldMakeFigures (1,1) logical = true
end

domainSize = [750e3 750e3 2e3];
latitude = 45;
frequency = 2*pi/(12.4206012*3600);
tidalPeriod = 2*pi/frequency;
barotropicVelocityAmplitude = [0.05; 0];
N2Function = @(z) 2e-5*exp(z/4000);

wvt = WVTransformBoussinesq(domainSize,options.resolution,N2=N2Function,latitude=latitude,shouldAntialias=options.shouldAntialias);
wvt.t0 = 0;
wvt.t = 0;
[virtualDepth,topographicHeight,topographyDiagnostics] = WVBottomWaveGenerationForcing.goffAbyssalHillTopography(wvt,rmsHeight=100,cornerWavenumber=1e-4,minimumWavelength=options.minimumWavelength,randomSeed=2023);
forcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=topographicHeight,barotropicVelocityAmplitude=barotropicVelocityAmplitude,frequency=frequency,rampDuration=0,startTime=0);
wvt.removeAllForcing();
wvt.addForcing(forcing);

model = WVModel(wvt);
model.setupIntegrator(integratorType="adaptive",absTolerance=options.relativeTolerance,relTolerance=options.relativeTolerance);
time = linspace(0,tidalPeriod,options.numberOfOutputTimes);
waveEnergy = zeros(size(time));
balancedEnergy = zeros(size(time));
modalSourcePower = zeros(size(time));
bottomSourcePower = zeros(size(time));
warningState = warning;
warningCleanup = onCleanup(@()warning(warningState));
warning("off","all")
for iTime = 1:numel(time)
    if iTime > 1
        model.integrateToTime(time(iTime),shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
    end
    diagnostics = sourceDiagnostics(wvt,forcing);
    waveEnergy(iTime) = sum(wvt.Apm_TE_factor(:).*(abs(wvt.Ap(:)).^2+abs(wvt.Am(:)).^2));
    balancedEnergy(iTime) = sum(wvt.A0_TE_factor(:).*abs(wvt.A0(:)).^2);
    modalSourcePower(iTime) = diagnostics.modalPower;
    bottomSourcePower(iTime) = diagnostics.bottomPower;
end
clear warningCleanup

cumulativeBottomWork = cumtrapz(time,bottomSourcePower);
energyWorkError = abs(waveEnergy(end)-cumulativeBottomWork(end))/max([waveEnergy(end) abs(cumulativeBottomWork(end)) eps]);
powerAgreement = norm(modalSourcePower-bottomSourcePower)/max(norm(bottomSourcePower),eps);
qgpvNorm = norm(wvt.qgpv(:));
waveVelocityScale = sqrt(2*waveEnergy(end)/wvt.Lz);
normalizedQGPV = qgpvNorm/max(abs(wvt.f)*waveVelocityScale/wvt.Lz,eps);
[modeNumbers,radialWavenumber,energyByModeAndWavenumber] = modalEnergyDistribution(wvt);

figureHandles = gobjects(0);
if options.shouldMakeFigures
    figureHandles = createFigure(wvt,time,tidalPeriod,topographicHeight,topographyDiagnostics,waveEnergy,cumulativeBottomWork,modeNumbers,radialWavenumber,energyByModeAndWavenumber);
end

configuration = struct(domainSize=domainSize,resolution=options.resolution,latitude=latitude,frequency=frequency,tidalPeriod=tidalPeriod,barotropicVelocityAmplitude=barotropicVelocityAmplitude,minimumWavelength=options.minimumWavelength,relativeTolerance=options.relativeTolerance,shouldAntialias=options.shouldAntialias);
result = struct(configuration=configuration,virtualDepth=virtualDepth,topographicHeight=topographicHeight,topographyDiagnostics=topographyDiagnostics,time=time,waveEnergy=waveEnergy,balancedEnergy=balancedEnergy,modalSourcePower=modalSourcePower,bottomSourcePower=bottomSourcePower,cumulativeBottomWork=cumulativeBottomWork,energyWorkError=energyWorkError,powerAgreement=powerAgreement,normalizedQGPV=normalizedQGPV,modeNumbers=modeNumbers,radialWavenumber=radialWavenumber,energyByModeAndWavenumber=energyByModeAndWavenumber,figureHandles=figureHandles);
end

function diagnostics = sourceDiagnostics(wvt,forcing)
[Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
modalPower = 2*sum(wvt.Apm_TE_factor(:).*real(Fp(:).*conj(wvt.Ap(:))+Fm(:).*conj(wvt.Am(:))));
[~,iBottom] = min(wvt.z);
bottomPower = mean((wvt.p(:,:,iBottom)/wvt.rho0).*forcing.bottomVelocityAtTime(wvt.t),"all");
diagnostics = struct(F0=F0,modalPower=modalPower,bottomPower=bottomPower);
end

function [modeNumbers,radialWavenumber,energyDistribution] = modalEnergyDistribution(wvt)
waveMask = wvt.waveComponent.maskAp | wvt.waveComponent.maskAm;
modeNumbers = unique(wvt.J(waveMask));
binWidth = min(wvt.dk,wvt.dl);
radialEdges = (-0.5:ceil(max(wvt.Kh,[],"all")/binWidth)+0.5)*binWidth;
radialWavenumber = (radialEdges(1:end-1)+radialEdges(2:end))/2;
radialBin = discretize(wvt.Kh,radialEdges);
modalEnergy = wvt.Apm_TE_factor.*(abs(wvt.Ap).^2+abs(wvt.Am).^2);
energyDistribution = zeros(numel(modeNumbers),numel(radialWavenumber));
for iMode = 1:numel(modeNumbers)
    for iBin = 1:numel(radialWavenumber)
        indices = waveMask & wvt.J == modeNumbers(iMode) & radialBin == iBin;
        energyDistribution(iMode,iBin) = sum(modalEnergy(indices));
    end
end
end

function figureHandles = createFigure(wvt,time,tidalPeriod,topographicHeight,diagnostics,waveEnergy,cumulativeBottomWork,modeNumbers,radialWavenumber,energyDistribution)
figureHandles = figure(Name="Goff abyssal-hill wave generation",Color="w",Position=[100 100 1250 850]);
layout = tiledlayout(figureHandles,2,2,TileSpacing="compact",Padding="compact");

axesHandle = nexttile(layout,1);
imagesc(axesHandle,wvt.x/1e3,wvt.y/1e3,topographicHeight.')
axis(axesHandle,"xy","image")
xlabel(axesHandle,"x (km)")
ylabel(axesHandle,"y (km)")
title(axesHandle,sprintf("Goff topography, h_{rms} = %.0f m",diagnostics.rmsHeight))
colorbar(axesHandle)

axesHandle = nexttile(layout,2);
validSpectrum = diagnostics.radialWavenumber > 0 & diagnostics.radialModeCount > 0 & isfinite(diagnostics.radialPowerSpectrum) & isfinite(diagnostics.radialTargetPowerSpectrum);
loglog(axesHandle,diagnostics.radialWavenumber(validSpectrum),diagnostics.radialPowerSpectrum(validSpectrum),"o",MarkerSize=4)
hold(axesHandle,"on")
loglog(axesHandle,diagnostics.radialWavenumber(validSpectrum),diagnostics.radialTargetPowerSpectrum(validSpectrum),LineWidth=1.5)
xline(axesHandle,diagnostics.cutoffWavenumber,"--",HandleVisibility="off")
grid(axesHandle,"on")
xlabel(axesHandle,"K (rad m^{-1})")
ylabel(axesHandle,"P_h(K)")
legend(axesHandle,"realized","target",Location="best")
title(axesHandle,"Terrain spectrum")

axesHandle = nexttile(layout,3);
plot(axesHandle,time/tidalPeriod,waveEnergy,LineWidth=1.6)
hold(axesHandle,"on")
plot(axesHandle,time/tidalPeriod,cumulativeBottomWork,"--",LineWidth=1.5)
grid(axesHandle,"on")
xlabel(axesHandle,"time / M2 period")
ylabel(axesHandle,"energy or work (m^3 s^{-2})")
legend(axesHandle,"wave energy","bottom pressure work",Location="best")
title(axesHandle,"Prescribed energy input")

axesHandle = nexttile(layout,4);
imagesc(axesHandle,radialWavenumber/1e-3,modeNumbers,energyDistribution)
axis(axesHandle,"xy")
xlabel(axesHandle,"K (10^{-3} rad m^{-1})")
ylabel(axesHandle,"vertical mode j")
title(axesHandle,"Final wave-energy distribution")
colorbar(axesHandle)
colormap(figureHandles,"turbo")
title(layout,"Variable-stratification Goff abyssal-hill generation")
end
