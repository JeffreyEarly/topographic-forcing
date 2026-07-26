function movie = GaussianRidgeWaveScatteringMovie(inputPath,options)
% Render an x-z movie from Gaussian-ridge WaveVortexModel output.
%
% The renderer reads standard `wave-vortex` coefficients, reconstructs a
% fixed-scale first-order mean-depth section, and combines it with spatial
% and modal wave-energy diagnostics.
%
% - Topic: Examples
% - Declaration: movie = GaussianRidgeWaveScatteringMovie(inputPath,options)
% - Parameter inputPath: standard WaveVortexModel NetCDF output
% - Parameter options.videoPath: MPEG-4 path; empty derives it from `inputPath`
% - Parameter options.posterPath: PNG path; empty derives it from `inputPath`
% - Parameter options.forcingName: scattering forcing name; empty requires exactly one
% - Parameter options.field: main section field, one of `u`, `w`, or `eta`
% - Parameter options.iTime: unique increasing output indices; empty selects all
% - Parameter options.frameRate: MPEG-4 frames per second
% - Parameter options.quality: MPEG-4 quality from 0 through 100
% - Parameter options.verticalExaggeration: vertical exaggeration of the x-z section
% - Parameter options.shouldOverwriteExisting: whether existing movie products may be replaced
% - Returns movie: paths, frame metadata, scales, selected indices, and poster frame
arguments (Input)
    inputPath (1,1) string {mustBeFile}
    options.videoPath (1,1) string = ""
    options.posterPath (1,1) string = ""
    options.forcingName (1,1) string = ""
    options.field (1,1) string {mustBeMember(options.field,["u" "w" "eta"])} = "u"
    options.iTime double {mustBeInteger,mustBePositive} = double.empty(0,1)
    options.frameRate (1,1) double {mustBeInteger,mustBePositive} = 15
    options.quality (1,1) double {mustBeInteger,mustBeBetween(options.quality,0,100)} = 95
    options.verticalExaggeration (1,1) double {mustBePositive,mustBeFinite} = 20
    options.shouldOverwriteExisting (1,1) logical = false
end
arguments (Output)
    movie (1,1) struct
end

[inputFolder,inputName] = fileparts(inputPath);
if strlength(options.videoPath) == 0
    videoPath = fullfile(inputFolder,inputName+"-scattering.mp4");
else
    videoPath = options.videoPath;
end
if strlength(options.posterPath) == 0
    posterPath = fullfile(inputFolder,inputName+"-scattering.png");
else
    posterPath = options.posterPath;
end
prepareOutputPath(videoPath,options.shouldOverwriteExisting);
prepareOutputPath(posterPath,options.shouldOverwriteExisting);

[wvt,ncfile] = WVTransformBoussinesq.waveVortexTransformFromFile(char(inputPath),iTime=1,shouldReadOnly=true);
fileCleanup = onCleanup(@()ncfile.close());
if ~ncfile.hasGroupWithName("wave-vortex")
    error("GaussianRidgeWaveScatteringMovie:MissingOutputGroup", "The file does not contain the standard 'wave-vortex' output group.")
end
group = ncfile.groupWithName("wave-vortex");
for variableName = ["t" "Ap" "Am" "A0"]
    if ~group.hasVariableWithName(variableName)
        error("GaussianRidgeWaveScatteringMovie:MissingOutputVariable", "The 'wave-vortex' group does not contain '%s'.",variableName)
    end
end
forcing = selectScatteringForcing(wvt,options.forcingName);
allTime = reshape(group.readVariables("t"),[],1);
if isempty(options.iTime)
    timeIndices = reshape(1:numel(allTime),[],1);
else
    if ~isvector(options.iTime)
        error("GaussianRidgeWaveScatteringMovie:InvalidTimeIndices", "iTime must be a vector of output indices.")
    end
    timeIndices = reshape(options.iTime,[],1);
end
if any(timeIndices > numel(allTime)) || any(diff(timeIndices) <= 0)
    error("GaussianRidgeWaveScatteringMovie:InvalidTimeIndices", "iTime must contain unique increasing indices within the saved time axis.")
end
time = allTime(timeIndices);
if numel(time) < 2 || any(diff(time) <= 0)
    error("GaussianRidgeWaveScatteringMovie:InvalidTimeAxis", "At least two strictly increasing output times are required.")
end

numberOfFrames = numel(time);
Ap = complex(zeros([size(wvt.Ap) numberOfFrames]));
Am = complex(zeros([size(wvt.Am) numberOfFrames]));
A0 = complex(zeros([size(wvt.A0) numberOfFrames]));
for iFrame = 1:numberOfFrames
    [Ap(:,:,iFrame),Am(:,:,iFrame),A0(:,:,iFrame)] = group.readVariablesAtIndexAlongDimension("t",timeIndices(iFrame),"Ap","Am","A0");
end
clear fileCleanup

diagnostics = gaussianRidgeScatteringDiagnostics(wvt,forcing,Ap,Am,A0,time);
[dominantFrequency,firstModeWavelength,groupVelocity] = dominantWaveScales(wvt,Ap(:,:,1),Am(:,:,1));
wavePeriod = 2*pi/dominantFrequency;
travelDistance = groupVelocity*(time-time(1));
[section,sectionZ,energyProfile,fieldUnits,fieldScaleFactor,wetMask] = reconstructMovieFields(wvt,forcing,Ap,Am,A0,time,options.field);
colorLimit = max(abs(section),[],"all","omitmissing");
if colorLimit == 0
    error("GaussianRidgeWaveScatteringMovie:ZeroField", "The selected field '%s' vanishes at every requested frame.",options.field)
end
profileMean = reshape(mean(energyProfile,1),[],1);
energyProfileConsistencyError = max(abs(profileMean-diagnostics.firstOrderEnergy)./max(abs(diagnostics.firstOrderEnergy),eps));
if energyProfileConsistencyError > 1e-9
    error("GaussianRidgeWaveScatteringMovie:InconsistentEnergyProfile", "The horizontal mean of the first-order energy profile differs from the modal diagnostic by %.3g relative error.",energyProfileConsistencyError)
end
posterFrameIndex = find(diagnostics.leftwardFirstModeEnergy+diagnostics.higherModeEnergy == max(diagnostics.leftwardFirstModeEnergy+diagnostics.higherModeEnergy),1,"first");

figureHandle = figure(Name="Gaussian-ridge wave-scattering movie",Color="w",Position=[100 100 1440 800]);
figureCleanup = onCleanup(@()closeValidFigure(figureHandle));
layout = tiledlayout(figureHandle,2,2,TileSpacing="compact",Padding="compact");
layoutTitle = title(layout,"");
sectionAxes = nexttile(layout,[1 2]);
sectionImage = imagesc(sectionAxes,wvt.x/1e3,sectionZ/1e3,section(:,:,1).'*fieldScaleFactor);
sectionImage.AlphaData = wetMask.';
axis(sectionAxes,"xy")
hold(sectionAxes,"on")
bottom = -wvt.Lz+forcing.topographicHeight(:,1);
terrainPatch = fill(sectionAxes,[wvt.x(:); flipud(wvt.x(:))]/1e3,[bottom; -wvt.Lz*ones(size(bottom))]/1e3,[0.32 0.24 0.16],EdgeColor="none");
terrainLine = plot(sectionAxes,wvt.x/1e3,bottom/1e3,Color=[0.12 0.08 0.04],LineWidth=1.2);
clim(sectionAxes,fieldScaleFactor*colorLimit*[-1 1])
colormap(sectionAxes,divergingColormap())
sectionColorbar = colorbar(sectionAxes);
sectionColorbar.Label.String = fieldUnits;
xlabel(sectionAxes,"x (km)")
ylabel(sectionAxes,"z (km)")
pbaspect(sectionAxes,[wvt.Lx options.verticalExaggeration*wvt.Lz 1])
title(sectionAxes,sprintf("First-order mean-depth reconstruction (%g× vertical exaggeration)",options.verticalExaggeration))

energyAxes = nexttile(layout,3);
energyLine = plot(energyAxes,wvt.x/1e3,energyProfile(:,1),LineWidth=1.7);
grid(energyAxes,"on")
xlabel(energyAxes,"x (km)")
ylabel(energyAxes,"first-order depth-integrated wave energy")
title(energyAxes,"Instantaneous packet energy")
xlim(energyAxes,[min(wvt.x) max(wvt.x)]/1e3)
energyLimits = [min(0,1.05*min(energyProfile,[],"all")) max(0,1.05*max(energyProfile,[],"all"))];
ylim(energyAxes,energyLimits)

historyAxes = nexttile(layout,4);
normalization = diagnostics.initialEnergy;
yyaxis(historyAxes,"left")
rightwardLine = plot(historyAxes,travelDistance/firstModeWavelength,diagnostics.rightwardFirstModeEnergy/normalization,LineWidth=1.5);
hold(historyAxes,"on")
ylabel(historyAxes,"rightward mode 1 / initial energy")
yyaxis(historyAxes,"right")
reflectedLine = plot(historyAxes,travelDistance/firstModeWavelength,diagnostics.leftwardFirstModeEnergy/normalization,LineWidth=1.5);
higherLine = plot(historyAxes,travelDistance/firstModeWavelength,diagnostics.higherModeEnergy/normalization,LineWidth=1.5);
ylabel(historyAxes,"scattered energy / initial energy")
yyaxis(historyAxes,"left")
timeMarker = xline(historyAxes,travelDistance(1)/firstModeWavelength,"k--",HandleVisibility="off");
grid(historyAxes,"on")
xlabel(historyAxes,"packet travel distance / mode-1 wavelength")
legend(historyAxes,[rightwardLine reflectedLine higherLine],["rightward mode 1" "reflected mode 1" "higher modes"],Location="best")
title(historyAxes,"Modal scattering")

video = VideoWriter(char(videoPath),"MPEG-4");
video.FrameRate = options.frameRate;
video.Quality = options.quality;
open(video)
videoCleanup = onCleanup(@()closeVideoIfOpen(video));
for iFrame = 1:numberOfFrames
    sectionImage.CData = section(:,:,iFrame).'*fieldScaleFactor;
    energyLine.YData = energyProfile(:,iFrame);
    timeMarker.Value = travelDistance(iFrame)/firstModeWavelength;
    layoutTitle.String = sprintf("Gaussian-ridge wave scattering: (t-t_1)/T = %.2f, packet distance = %.2f wavelengths",(time(iFrame)-time(1))/wavePeriod,travelDistance(iFrame)/firstModeWavelength);
    uistack(terrainPatch,"top")
    uistack(terrainLine,"top")
    drawnow
    writeVideo(video,getframe(figureHandle))
    if iFrame == posterFrameIndex
        exportgraphics(figureHandle,posterPath,Resolution=180)
    end
end
close(video)
clear videoCleanup

reader = VideoReader(char(videoPath));
readableFrames = 0;
frameDimensions = [];
while hasFrame(reader)
    frame = readFrame(reader);
    readableFrames = readableFrames+1;
    if isempty(frameDimensions)
        frameDimensions = size(frame);
    elseif ~isequal(size(frame),frameDimensions)
        error("GaussianRidgeWaveScatteringMovie:InconsistentFrameSize", "The encoded movie contains inconsistent frame dimensions.")
    end
end
if readableFrames ~= numberOfFrames
    error("GaussianRidgeWaveScatteringMovie:IncompleteVideo", "Expected %d readable frames but found %d.",numberOfFrames,readableFrames)
end
if abs(reader.Duration-numberOfFrames/options.frameRate) > 1/options.frameRate
    error("GaussianRidgeWaveScatteringMovie:InvalidDuration", "The encoded movie duration does not match the requested frame count and frame rate.")
end
videoInformation = dir(videoPath);
if isempty(videoInformation) || videoInformation.bytes == 0 || ~isfile(posterPath)
    error("GaussianRidgeWaveScatteringMovie:EmptyOutput", "The movie or poster output is missing or empty.")
end

movie = struct(inputPath=inputPath,videoPath=videoPath,posterPath=posterPath,forcingName=string(forcing.name),field=options.field,time=time,timeIndices=timeIndices,frameRate=options.frameRate,numberOfFrames=numberOfFrames,readableFrames=readableFrames,frameDimensions=frameDimensions,duration=reader.Duration,colorLimits=colorLimit*[-1 1],fieldScaleFactor=fieldScaleFactor,fieldUnits=fieldUnits,verticalExaggeration=options.verticalExaggeration,sectionZ=sectionZ,wetMask=wetMask,numberOfMaskedGridCells=nnz(~wetMask),energyProfile=energyProfile,energyProfileConsistencyError=energyProfileConsistencyError,dominantFrequency=dominantFrequency,wavePeriod=wavePeriod,firstModeWavelength=firstModeWavelength,groupVelocity=groupVelocity,posterFrameIndex=posterFrameIndex);
clear figureCleanup
end

function forcing = selectScatteringForcing(wvt,forcingName)
forcingNames = string(wvt.forcingNames());
isScattering = false(size(forcingNames));
for iForcing = 1:numel(forcingNames)
    isScattering(iForcing) = isa(wvt.forcingWithName(char(forcingNames(iForcing))),"WVBottomWaveScatteringForcing");
end
if strlength(forcingName) == 0
    matchingNames = forcingNames(isScattering);
    if numel(matchingNames) ~= 1
        error("GaussianRidgeWaveScatteringMovie:AmbiguousForcing", "The input must contain exactly one WVBottomWaveScatteringForcing when forcingName is omitted; found %d.",numel(matchingNames))
    end
    forcingName = matchingNames;
elseif ~any(forcingNames == forcingName) || ~isa(wvt.forcingWithName(char(forcingName)),"WVBottomWaveScatteringForcing")
    error("GaussianRidgeWaveScatteringMovie:InvalidForcing", "The input does not contain a WVBottomWaveScatteringForcing named '%s'.",forcingName)
end
forcing = wvt.forcingWithName(char(forcingName));
end

function [frequency,wavelength,groupVelocity] = dominantWaveScales(wvt,Ap,Am)
modalEnergy = wvt.Apm_TE_factor.*(abs(Ap).^2+abs(Am).^2);
[maximumEnergy,index] = max(modalEnergy,[],"all","linear");
if maximumEnergy == 0 || wvt.Kh(index) == 0
    error("GaussianRidgeWaveScatteringMovie:MissingIncidentWave", "The first selected record does not contain a propagating wave.")
end
frequency = wvt.Omega(index);
wavelength = 2*pi/wvt.Kh(index);
N2 = reshape(wvt.N2,[],1);
if max(N2)-min(N2) > 100*eps(max(N2))
    error("GaussianRidgeWaveScatteringMovie:NonconstantStratification", "The Gaussian-ridge movie currently requires constant stratification to diagnose packet travel distance.")
end
horizontalWavenumber = wvt.Kh(index);
verticalWavenumber = wvt.J(index)*pi/wvt.Lz;
groupVelocity = horizontalWavenumber*verticalWavenumber^2*(mean(N2)-wvt.f^2)/(frequency*(horizontalWavenumber^2+verticalWavenumber^2)^2);
end

function [section,sectionZ,energyProfile,units,scaleFactor,wetMask] = reconstructMovieFields(wvt,forcing,Ap,Am,A0,time,fieldName)
numberOfFrames = numel(time);
nativeSection = zeros(wvt.Nx,wvt.Nz,numberOfFrames);
energyProfile = zeros(wvt.Nx,numberOfFrames);
[~,bottomIndex] = min(wvt.z);
for iFrame = 1:numberOfFrames
    [u,v,w,eta] = wvt.transformWaveVortexToUVWEta(Ap(:,:,iFrame),Am(:,:,iFrame),A0(:,:,iFrame),time(iFrame));
    switch fieldName
        case "u"
            field = u;
        case "w"
            field = w;
        case "eta"
            field = eta;
    end
    nativeSection(:,:,iFrame) = squeeze(field(:,1,:));
    energyDensity = 0.5*(u.^2+v.^2+w.^2+reshape(wvt.N2,1,1,[]).*eta.^2);
    referenceEnergyProfile = squeeze(mean(sum(energyDensity.*reshape(wvt.z_int,1,1,[]),3),2));
    bottomCorrectionProfile = 0.5*mean(forcing.topographicHeight.*(u(:,:,bottomIndex).^2+v(:,:,bottomIndex).^2),2);
    energyProfile(:,iFrame) = referenceEnergyProfile-bottomCorrectionProfile;
end
bottom = -wvt.Lz+forcing.topographicHeight(:,1);
sectionZ = linspace(min(wvt.z),max(wvt.z),16*(wvt.Nz-1)+1);
section = zeros(wvt.Nx,numel(sectionZ),numberOfFrames);
for iFrame = 1:numberOfFrames
    section(:,:,iFrame) = interp1(reshape(wvt.z,[],1),nativeSection(:,:,iFrame).',sectionZ,"linear").';
end
wetMask = reshape(sectionZ,1,[]) >= bottom;
section(repmat(~wetMask,1,1,numberOfFrames)) = NaN;
switch fieldName
    case "u"
        units = "u (cm s^{-1})";
        scaleFactor = 1e2;
    case "w"
        units = "w (mm s^{-1})";
        scaleFactor = 1e3;
    case "eta"
        units = "\eta (m)";
        scaleFactor = 1;
end
end

function prepareOutputPath(path,shouldOverwrite)
parent = fileparts(path);
if strlength(parent) > 0 && ~isfolder(parent)
    error("GaussianRidgeWaveScatteringMovie:MissingOutputDirectory", "The output directory '%s' does not exist.",parent)
end
if isfile(path)
    if ~shouldOverwrite
        error("GaussianRidgeWaveScatteringMovie:OutputExists", "The output '%s' already exists. Set shouldOverwriteExisting=true to replace it.",path)
    end
    delete(path)
end
end

function closeVideoIfOpen(video)
try
    close(video)
catch
end
end

function closeValidFigure(figureHandle)
if isgraphics(figureHandle)
    close(figureHandle)
end
end

function colors = divergingColormap()
colors = interp1([0 0.5 1],[0.15 0.30 0.72; 0.98 0.98 0.98; 0.70 0.09 0.17],linspace(0,1,256));
end
