function [model,wvt,filename] = EddyTideTopographicForcingSimulation(options)
% Run the eddy-tide experiment with topographic wave forcing.
%
% The experiment preserves the constant-stratification domain and shallow
% Gaussian eddy used by the JPO2026 minimal simulation, but starts with no
% wave field. A prescribed barotropic M2 current generates waves over a
% deterministic Goff abyssal-hill realization. Autonomous first-order
% topographic scattering is available as an opt-in validation experiment,
% but is not part of the default generation-only calculation.
%
% The default integration is a 30-day validation run. Set `includeEddy` to
% `false` for the matched no-eddy control. Both cases use the same terrain
% when their terrain options and random seed agree.
%
% ```matlab
% [~,~,eddyFile] = EddyTideTopographicForcingSimulation;
% [~,~,controlFile] = EddyTideTopographicForcingSimulation(includeEddy=false);
% ```
%
% - Topic: Examples
% - Declaration: [model,wvt,filename] = EddyTideTopographicForcingSimulation(options)
% - Parameter options.Nxy: horizontal grid resolution
% - Parameter options.horizontalDomainSize: square horizontal-domain size in meters, or `NaN` for four mode-one wavelengths
% - Parameter options.latitude: latitude in degrees north
% - Parameter options.includeEddy: whether to initialize the shallow Gaussian eddy
% - Parameter options.maxT: final integration time in seconds
% - Parameter options.outputInterval: model-output interval in seconds
% - Parameter options.outputDirectory: directory for the restartable NetCDF output
% - Parameter options.outputFilename: optional explicit output filename
% - Parameter options.shouldOverwriteExisting: whether to replace an existing output file
% - Parameter options.barotropicVelocityAmplitude: complex zonal and meridional M2 velocity amplitudes in meters per second
% - Parameter options.tidalPeriod: barotropic tidal period in seconds
% - Parameter options.rampDuration: half-cosine startup-ramp duration, or `NaN` for one tidal period
% - Parameter options.rmsHeight: post-filter RMS terrain height in meters
% - Parameter options.cornerWavenumber: Goff terrain corner wavenumber in radians per meter
% - Parameter options.minimumWavelength: shortest retained terrain wavelength in meters
% - Parameter options.randomSeed: nonnegative terrain random seed
% - Parameter options.shouldUseScattering: whether to add autonomous first-order wave scattering
% - Parameter options.shouldAntialias: transform antialias setting
% - Parameter options.shouldShowIntegrationDiagnostics: whether to print integration progress
% - Returns model: integrated `WVModel` with its output file closed
% - Returns wvt: integrated `WVTransformConstantStratification`
% - Returns filename: restartable NetCDF output path
arguments (Input)
    options.Nxy (1,1) double {mustBeInteger,mustBePositive} = 128
    options.horizontalDomainSize (1,1) double = NaN
    options.latitude (1,1) double {mustBeFinite} = 45
    options.includeEddy (1,1) logical = true
    options.maxT (1,1) double {mustBePositive,mustBeFinite} = 30*86400
    options.outputInterval (1,1) double {mustBePositive,mustBeFinite} = 3600
    options.outputDirectory (1,1) string = defaultOutputDirectory()
    options.outputFilename (1,1) string = ""
    options.shouldOverwriteExisting (1,1) logical = false
    options.barotropicVelocityAmplitude (2,1) double = [0.05; 0]
    options.tidalPeriod (1,1) double {mustBePositive,mustBeFinite} = 12.420602*3600
    options.rampDuration (1,1) double = NaN
    options.rmsHeight (1,1) double {mustBePositive,mustBeFinite} = 100
    options.cornerWavenumber (1,1) double {mustBePositive,mustBeFinite} = 1e-4
    options.minimumWavelength (1,1) double {mustBePositive,mustBeFinite} = 20e3
    options.randomSeed (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 2023
    options.shouldUseScattering (1,1) logical = false
    options.shouldAntialias (1,1) logical = true
    options.shouldShowIntegrationDiagnostics (1,1) logical = true
end
arguments (Output)
    model WVModel
    wvt WVTransformConstantStratification
    filename (1,1) string
end

if any(~isfinite(options.barotropicVelocityAmplitude),"all")
    error("EddyTideTopographicForcingSimulation:InvalidBarotropicVelocityAmplitude", "barotropicVelocityAmplitude must contain finite values.")
end
if ~(isnan(options.horizontalDomainSize) || isfinite(options.horizontalDomainSize) && options.horizontalDomainSize > 0)
    error("EddyTideTopographicForcingSimulation:InvalidHorizontalDomainSize", "horizontalDomainSize must be NaN or a finite positive scalar in meters.")
end
if ~(isnan(options.rampDuration) || isfinite(options.rampDuration) && options.rampDuration >= 0)
    error("EddyTideTopographicForcingSimulation:InvalidRampDuration", "rampDuration must be NaN or a finite nonnegative scalar.")
end
if isnan(options.rampDuration)
    rampDuration = options.tidalPeriod;
else
    rampDuration = options.rampDuration;
end

N0 = sqrt(2e-5);
Lz = 2000;
N2 = @(z) N0*N0*ones(size(z));
frequency = 2*pi/options.tidalPeriod;

internalModes = InternalModesWKBSpectral(N2=N2,zIn=[-Lz 0],latitude=options.latitude);
[~,~,~,kM2] = internalModes.ModesAtFrequency(frequency);
modeOneWavelength = 2*pi/kM2(1);
if isnan(options.horizontalDomainSize)
    Lxy = 4*modeOneWavelength;
else
    Lxy = options.horizontalDomainSize;
end
Nz = WVStratification.verticalResolutionForHorizontalResolution(Lxy,Lz,options.Nxy,N2=N2,latitude=options.latitude);
wvt = WVTransformConstantStratification([Lxy Lxy Lz],[options.Nxy options.Nxy Nz],N0=N0,latitude=options.latitude,isHydrostatic=false,shouldAntialias=options.shouldAntialias);

wvt.addForcing(WVAdaptiveDamping(wvt));
adaptiveDamping = wvt.forcingWithName("adaptive damping");
wvt.removeAll;

if options.includeEddy
    addShallowGaussianEddy(wvt,adaptiveDamping);
end

[~,topographicHeight] = WVBottomWaveGenerationForcing.goffAbyssalHillTopography(wvt,rmsHeight=options.rmsHeight,cornerWavenumber=options.cornerWavenumber,minimumWavelength=options.minimumWavelength,randomSeed=options.randomSeed);
generation = WVBottomWaveGenerationForcing(wvt,topographicHeight=topographicHeight,barotropicVelocityAmplitude=options.barotropicVelocityAmplitude,frequency=frequency,rampDuration=rampDuration,startTime=0,name="bottom wave generation");
wvt.addForcing(generation);
if options.shouldUseScattering
    scattering = WVBottomWaveScatteringForcing(wvt,topographicHeight=topographicHeight,name="bottom wave scattering");
    wvt.addForcing(scattering);
end

if ~isfolder(options.outputDirectory)
    mkdir(options.outputDirectory)
end
outputFilename = options.outputFilename;
if strlength(outputFilename) == 0
    outputFilename = defaultOutputFilename(options);
elseif ~endsWith(outputFilename,".nc",IgnoreCase=true)
    outputFilename = outputFilename+".nc";
end
filename = fullfile(options.outputDirectory,outputFilename);
if isfile(filename) && ~options.shouldOverwriteExisting
    error("EddyTideTopographicForcingSimulation:OutputFileExists", "The model output '%s' already exists. Set shouldOverwriteExisting=true to replace it.",filename)
end

model = WVModel(wvt);
model.createNetCDFFileForModelOutput(char(filename),outputInterval=options.outputInterval,shouldOverwriteExisting=options.shouldOverwriteExisting);
outputCleanup = onCleanup(@()model.closeNetCDFFile());
model.integrateToTime(options.maxT,shouldShowIntegrationDiagnostics=options.shouldShowIntegrationDiagnostics,callback=@(~)[]);
clear outputCleanup
end

function addShallowGaussianEddy(wvt,adaptiveDamping)
x0 = wvt.Lx/2;
y0 = wvt.Ly/2;
Le = 80e3;
He = 300;
U = 0.10;
verticalStructure = @(z) exp(-(z/He/sqrt(2)).^2);
horizontalStructure = @(x,y) exp(-((x-x0)/Le).^2-((y-y0)/Le).^2);
psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*verticalStructure(z).*(horizontalStructure(x,y)-pi*Le*Le/(wvt.Lx*wvt.Ly));
wvt.addGeostrophicStreamfunction(psi);
wvt.A0(wvt.Kh > adaptiveDamping.k_damp) = 0;
end

function filename = defaultOutputFilename(options)
if options.includeEddy
    experimentName = "eddy";
else
    experimentName = "no-eddy";
end
if options.shouldUseScattering
    forcingName = "generation-scattering";
else
    forcingName = "generation-only";
end
currentSpeedCentimetersPerSecond = 100*norm(options.barotropicVelocityAmplitude);
if isnan(options.horizontalDomainSize)
    domainName = "";
else
    domainName = sprintf("-Lxy%gkm",options.horizontalDomainSize/1e3);
end
filename = string(sprintf("eddy-tide-topographic-%s-%s%s-Nxy%d-Ubt%gcms-hrms%gm-lmin%gkm-seed%d.nc",experimentName,forcingName,domainName,options.Nxy,currentSpeedCentimetersPerSecond,options.rmsHeight,options.minimumWavelength/1e3,options.randomSeed));
end

function outputDirectory = defaultOutputDirectory()
[scriptFolder,~,~] = fileparts(mfilename("fullpath"));
outputDirectory = string(fullfile(fileparts(scriptFolder),"output"));
end
