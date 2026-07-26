function result = GaussianRidgeWaveScatteringExample(options)
% Scatter a mode-one M2 wave packet from a Gaussian ridge.
%
% The example uses the first-order mean-depth bottom condition and an
% adaptive WaveVortexModel integration. Standard NetCDF output is optional;
% the complete sampled coefficient trajectory is always returned.
%
% - Topic: Examples
% - Declaration: result = GaussianRidgeWaveScatteringExample(options)
% - Parameter options.resolution: transform resolution $$[N_x,N_y,N_z]$$
% - Parameter options.heightRatio: ridge height divided by reference depth
% - Parameter options.criticality: maximum ridge slope divided by the M2 ray slope
% - Parameter options.numberOfOutputTimes: number of sampled output times
% - Parameter options.relativeTolerance: adaptive relative and absolute tolerance
% - Parameter options.shouldAntialias: antialias setting for the transform
% - Parameter options.outputPath: optional standard WaveVortexModel NetCDF path
% - Parameter options.shouldOverwriteExisting: whether an existing output file may be replaced
% - Parameter options.shouldMakeFigures: whether to create the two summary figures
% - Returns result: configuration, coefficients, diagnostics, snapshots, output path, and figures
arguments
    options.resolution (1,3) double {mustBeInteger,mustBePositive} = [64 4 9]
    options.heightRatio (1,1) double {mustBePositive,mustBeLessThan(options.heightRatio,1)} = 0.05
    options.criticality (1,1) double {mustBePositive,mustBeFinite} = 0.05
    options.numberOfOutputTimes (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.numberOfOutputTimes,3)} = 121
    options.relativeTolerance (1,1) double {mustBePositive,mustBeFinite} = 1e-8
    options.shouldAntialias (1,1) logical = true
    options.outputPath (1,1) string = ""
    options.shouldOverwriteExisting (1,1) logical = false
    options.shouldMakeFigures (1,1) logical = true
end

if options.heightRatio > 0.1
    warning("GaussianRidgeWaveScatteringExample:LargeTerrainAmplitude", "heightRatio=%.3g is outside the preferred first-order range h_0/D <= 0.1.",options.heightRatio)
end
if strlength(options.outputPath) > 0 && isfile(options.outputPath) && ~options.shouldOverwriteExisting
    error("GaussianRidgeWaveScatteringExample:OutputExists", "The output file '%s' already exists. Set shouldOverwriteExisting=true to replace it.",options.outputPath)
end

[wvt,configuration] = gaussianRidgeScatteringConfiguration(options.resolution,6,options.shouldAntialias);
x = reshape(wvt.x,[],1);
ridgeCenter = 3*configuration.firstModeWavelength;
ridgeHeight = options.heightRatio*configuration.depth;
ridgeWidth = ridgeHeight/(options.criticality*configuration.mu*sqrt(exp(1)));
ridge = ridgeHeight*exp(-0.5*((x-ridgeCenter)/ridgeWidth).^2);
topographicHeight = repmat(ridge,1,wvt.Ny);
maximumSlope = maximumPeriodicSlope(ridge,wvt.Lx);
pointsAcrossFWHM = 2*sqrt(2*log(2))*ridgeWidth/(wvt.Lx/wvt.Nx);

packetCenter = 1.5*configuration.firstModeWavelength;
envelopeWidth = 0.5*configuration.firstModeWavelength;
initializeRightwardWavePacket(wvt,configuration.carrierMode,packetCenter,envelopeWidth);
maximumVelocity = max(hypot(wvt.u,wvt.v),[],"all");
if maximumVelocity == 0
    error("GaussianRidgeWaveScatteringExample:EmptyWavePacket", "The requested transform does not retain any wave-packet modes.")
end
wvt.Ap = wvt.Ap*(1e-2/maximumVelocity);
wvt.Am = wvt.Am*(1e-2/maximumVelocity);
initialAp = wvt.Ap;
initialAm = wvt.Am;
initialA0 = wvt.A0;

forcingName = "Gaussian-ridge wave scattering";
forcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=topographicHeight,name=forcingName);
wvt.removeAllForcing();
wvt.addForcing(forcing);
warningState = warning;
warningCleanup = onCleanup(@()warning(warningState));
warning("off","all")
model = WVModel(wvt);
warning(warningState)
model.setupIntegrator(integratorType="adaptive",absTolerance=options.relativeTolerance,relTolerance=options.relativeTolerance);
finalTime = 3*configuration.firstModeWavelength/configuration.groupVelocity;
outputTimes = linspace(0,finalTime,options.numberOfOutputTimes);
if strlength(options.outputPath) > 0
    model.createNetCDFFileForModelOutput(options.outputPath,outputInterval=outputTimes(2)-outputTimes(1),shouldOverwriteExisting=options.shouldOverwriteExisting);
end

Ap = complex(zeros([size(wvt.Ap) options.numberOfOutputTimes]));
Am = complex(zeros([size(wvt.Am) options.numberOfOutputTimes]));
A0 = complex(zeros([size(wvt.A0) options.numberOfOutputTimes]));
Ap(:,:,1) = wvt.Ap;
Am(:,:,1) = wvt.Am;
A0(:,:,1) = wvt.A0;
try
    for iTime = 2:options.numberOfOutputTimes
        model.integrateToTime(outputTimes(iTime),shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
        Ap(:,:,iTime) = wvt.Ap;
        Am(:,:,iTime) = wvt.Am;
        A0(:,:,iTime) = wvt.A0;
    end
    model.closeNetCDFFile();
catch exception
    model.closeNetCDFFile();
    rethrow(exception)
end
clear warningCleanup

diagnostics = gaussianRidgeScatteringDiagnostics(wvt,forcing,Ap,Am,A0,outputTimes);
flatAp = repmat(initialAp,1,1,options.numberOfOutputTimes);
flatAm = repmat(initialAm,1,1,options.numberOfOutputTimes);
flatA0 = repmat(initialA0,1,1,options.numberOfOutputTimes);
flatForcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny),name="flat reference");
flatDiagnostics = gaussianRidgeScatteringDiagnostics(wvt,flatForcing,flatAp,flatAm,flatA0,outputTimes);
snapshotIndices = unique([1 round((options.numberOfOutputTimes+1)/2) options.numberOfOutputTimes]);
snapshots = scatteringSnapshots(wvt,Ap,Am,A0,outputTimes,snapshotIndices);

configuration.ridgeCenter = ridgeCenter;
configuration.ridgeHeight = ridgeHeight;
configuration.ridgeWidth = ridgeWidth;
configuration.heightRatio = ridgeHeight/configuration.depth;
configuration.requestedCriticality = options.criticality;
configuration.discreteCriticality = maximumSlope/configuration.mu;
configuration.pointsAcrossFWHM = pointsAcrossFWHM;
configuration.packetCenter = packetCenter;
configuration.packetEnvelopeWidth = envelopeWidth;
configuration.maximumInitialHorizontalVelocity = 1e-2;
configuration.finalTime = finalTime;
configuration.relativeTolerance = options.relativeTolerance;
configuration.numberOfOutputTimes = options.numberOfOutputTimes;

figureHandles = gobjects(0);
if options.shouldMakeFigures
    figureHandles = plotGaussianRidgeResult(wvt,configuration,topographicHeight,outputTimes,diagnostics,flatDiagnostics,snapshots);
end
result = struct(configuration=configuration,topographicHeight=topographicHeight,time=reshape(outputTimes,[],1),distance=configuration.groupVelocity*reshape(outputTimes,[],1),coefficients=struct(Ap=Ap,Am=Am,A0=A0),flatReferenceCoefficients=struct(Ap=initialAp,Am=initialAm,A0=initialA0),diagnostics=diagnostics,flatReferenceDiagnostics=flatDiagnostics,snapshots=snapshots,integration=struct(integrator="adaptive ode78",numberOfFluxComputations=model.nFluxComputations,relativeTolerance=options.relativeTolerance),outputPath=options.outputPath,figureHandles=figureHandles);
end

function initializeRightwardWavePacket(wvt,carrierMode,packetCenter,envelopeWidth)
spectralWidth = wvt.Lx/(2*pi*envelopeWidth);
maximumMode = max(round(abs(wvt.K(:))*wvt.Lx/(2*pi)));
minimumMode = max(1,ceil(carrierMode-4*spectralWidth));
maximumPacketMode = min(maximumMode,floor(carrierMode+4*spectralWidth));
kMode = reshape(minimumMode:maximumPacketMode,[],1);
weights = exp(-0.5*((kMode-carrierMode)/spectralWidth).^2);
k = 2*pi*kMode/wvt.Lx;
phase = -k*packetCenter;
wvt.removeAll();
for iMode = 1:numel(kMode)
    wvt.setWaveModes(kMode=kMode(iMode),lMode=0,j=1,phi=phase(iMode),u=1e-4*weights(iMode),sign=-1);
end
end

function maximumSlope = maximumPeriodicSlope(height,Lx)
numberOfPoints = numel(height);
mode = [0:floor(numberOfPoints/2) -ceil(numberOfPoints/2)+1:-1].';
wavenumber = 2*pi*mode/Lx;
slope = real(ifft(1i*wavenumber.*fft(height)));
maximumSlope = max(abs(slope));
end

function snapshots = scatteringSnapshots(wvt,Ap,Am,A0,time,indices)
snapshots = repmat(struct(time=0,u=[],w=[],eta=[]),numel(indices),1);
for iSnapshot = 1:numel(indices)
    index = indices(iSnapshot);
    [u,~,w,eta] = wvt.transformWaveVortexToUVWEta(Ap(:,:,index),Am(:,:,index),A0(:,:,index),time(index));
    snapshots(iSnapshot).time = time(index);
    snapshots(iSnapshot).u = squeeze(u(:,1,:));
    snapshots(iSnapshot).w = squeeze(w(:,1,:));
    snapshots(iSnapshot).eta = squeeze(eta(:,1,:));
end
end

function figureHandles = plotGaussianRidgeResult(wvt,configuration,topographicHeight,time,diagnostics,flatDiagnostics,snapshots)
figureHandles = gobjects(2,1);
figureHandles(1) = figure(Name="Gaussian-ridge wave-scattering fields",Color="w",Position=[100 100 1200 760]);
layout = tiledlayout(figureHandles(1),2,2,TileSpacing="compact",Padding="compact");
axesHandle = nexttile(layout,1);
plot(axesHandle,wvt.x/1e3,topographicHeight(:,1),LineWidth=1.6)
grid(axesHandle,"on")
xlabel(axesHandle,"x (km)")
ylabel(axesHandle,"h (m)")
title(axesHandle,sprintf("Gaussian ridge: h_0/D=%.3f, criticality=%.3f",configuration.heightRatio,configuration.discreteCriticality))
sectionLimit = max(abs(cat(3,snapshots.u)),[],"all");
for iSnapshot = 1:numel(snapshots)
    axesHandle = nexttile(layout,iSnapshot+1);
    imagesc(axesHandle,wvt.x/1e3,wvt.z/1e3,snapshots(iSnapshot).u.'*1e2)
    axis(axesHandle,"xy")
    clim(axesHandle,1e2*sectionLimit*[-1 1])
    xlabel(axesHandle,"x (km)")
    ylabel(axesHandle,"z (km)")
    title(axesHandle,sprintf("u at t/T_{M2}=%.2f",snapshots(iSnapshot).time*configuration.targetFrequency/(2*pi)))
    colorbar(axesHandle)
end
colormap(figureHandles(1),divergingColormap())
title(layout,"First-order Gaussian-ridge wave scattering")

figureHandles(2) = figure(Name="Gaussian-ridge modal scattering",Color="w",Position=[150 150 1150 480]);
layout = tiledlayout(figureHandles(2),1,2,TileSpacing="compact",Padding="compact");
axesHandle = nexttile(layout,1);
normalization = diagnostics.initialEnergy;
yyaxis(axesHandle,"left")
rightwardLine = plot(axesHandle,configuration.groupVelocity*time/1e3,diagnostics.rightwardFirstModeEnergy/normalization,LineWidth=1.5);
hold(axesHandle,"on")
ylabel(axesHandle,"rightward mode 1 / initial energy")
yyaxis(axesHandle,"right")
leftwardLine = plot(axesHandle,configuration.groupVelocity*time/1e3,diagnostics.leftwardFirstModeEnergy/normalization,LineWidth=1.5);
higherModeLine = plot(axesHandle,configuration.groupVelocity*time/1e3,diagnostics.higherModeEnergy/normalization,LineWidth=1.5);
ylabel(axesHandle,"scattered energy / initial energy")
grid(axesHandle,"on")
xlabel(axesHandle,"packet travel distance (km)")
legend(axesHandle,[rightwardLine leftwardLine higherModeLine],["rightward mode 1" "leftward mode 1" "higher modes"],Location="best")
title(axesHandle,"Modal scattering")

axesHandle = nexttile(layout,2);
plot(axesHandle,time*configuration.targetFrequency/(2*pi),diagnostics.relativeFlatEnergyChange,LineWidth=1.5)
hold(axesHandle,"on")
plot(axesHandle,time*configuration.targetFrequency/(2*pi),diagnostics.relativeFirstOrderEnergyChange,LineWidth=1.5)
plot(axesHandle,time*configuration.targetFrequency/(2*pi),flatDiagnostics.relativeFlatEnergyChange,"--",LineWidth=1.2)
grid(axesHandle,"on")
xlabel(axesHandle,"time / M2 period")
ylabel(axesHandle,"relative energy change")
legend(axesHandle,"flat wave energy","first-order physical energy","flat control",Location="best")
title(axesHandle,"Energy diagnostics")
title(layout,"Wave-only modal and energetic response")
end

function colors = divergingColormap()
colors = interp1([0 0.5 1],[0.15 0.30 0.72; 0.98 0.98 0.98; 0.70 0.09 0.17],linspace(0,1,256));
end
