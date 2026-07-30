function run_masked_case(caseName)
arguments
    caseName (1,1) string {mustBeMember(caseName,["eddy" "control"])}
end

scriptDirectory = fileparts(mfilename("fullpath"));
runDirectory = fileparts(scriptDirectory);
repositoryRoot = fileparts(fileparts(fileparts(scriptDirectory)));
workspaceRoot = fileparts(repositoryRoot);
addpath(repositoryRoot)
addpath(fullfile(repositoryRoot,"Examples"))
addpath(fullfile(workspaceRoot,"wave-vortex-model"))
addpath(fullfile(workspaceRoot,"wave-vortex-model-diagnostics"))

if caseName == "eddy"
    includeEddy = true;
    caseFilename = "eddy-tide-topographic-eddy-generation-only-Lxy500km-Nxy256-Ein10-masked-day25-hrms100m-lmin6km-seed2023.nc";
    referenceFilename = "eddy-tide-topographic-eddy-generation-only-Lxy500km-Nxy256-Ein10-hrms100m-lmin6km-seed2023.nc";
else
    includeEddy = false;
    caseFilename = "eddy-tide-topographic-no-eddy-generation-only-Lxy500km-Nxy256-Ein10-masked-day25-hrms100m-lmin6km-seed2023.nc";
    referenceFilename = "eddy-tide-topographic-no-eddy-generation-only-Lxy500km-Nxy256-Ein10-hrms100m-lmin6km-seed2023.nc";
end

outputPath = fullfile(runDirectory,caseFilename);
referencePath = fullfile(repositoryRoot,"output",referenceFilename);
successPath = fullfile(runDirectory,caseName+"-success.mat");
failurePath = fullfile(runDirectory,"logs",caseName+"-failure-"+string(datetime("now","Format","yyyyMMdd-HHmmss"))+".mat");
if isfile(successPath)
    error("MaskedEin10Runner:ExistingOutput", "Refusing to replace completed production output for the %s case.",caseName)
end

try
    usableBytes = java.io.File(runDirectory).getUsableSpace();
    if usableBytes < 20*2^30
        error("MaskedEin10Runner:StorageReserve", "Only %.2f GiB remains before the %s case; at least 20 GiB is required.",double(usableBytes)/2^30,caseName)
    end

    preflightPath = fullfile(runDirectory,caseName+"-preflight.mat");
    if isfile(outputPath)
        verifyPartialTimeAxis(outputPath);
        model = WVModel.modelFromFile(char(outputPath));
        verifyRestartMask(model.wvt);
        loadedPreflight = load(preflightPath,"preflight");
        preflight = loadedPreflight.preflight;
    else
        [wvt,generation,adaptiveDamping] = buildCase(includeEddy);
        preflight = verifyInitialConfiguration(wvt,generation,adaptiveDamping,referencePath,caseName);
        save(preflightPath,"preflight","-v7.3")
        model = WVModel(wvt);
        model.createNetCDFFileForModelOutput(char(outputPath),outputInterval=6*3600,shouldOverwriteExisting=false);
    end
    outputCleanup = onCleanup(@()model.closeNetCDFFile());
    model.integrateToTime(25*86400,shouldShowIntegrationDiagnostics=true,callback=@(~)[]);
    clear outputCleanup

    postflight = verifyCompletedFile(outputPath);
    save(successPath,"caseName","outputPath","preflight","postflight","-v7.3")
catch exception
    failure = struct(identifier=string(exception.identifier),message=string(exception.message),report=string(getReport(exception,"extended","hyperlinks","off")));
    save(failurePath,"caseName","outputPath","failure","-v7.3")
    rethrow(exception)
end
end

function [wvt,generation,adaptiveDamping] = buildCase(includeEddy)
Nxy = 256;
Lxy = 500e3;
Lz = 2000;
latitude = 45;
N0 = sqrt(2e-5);
N2 = @(z)N0*N0*ones(size(z));
Nz = WVStratification.verticalResolutionForHorizontalResolution(Lxy,Lz,Nxy,N2=N2,latitude=latitude);
wvt = WVTransformConstantStratification([Lxy Lxy Lz],[Nxy Nxy Nz],N0=N0,latitude=latitude,isHydrostatic=false,shouldAntialias=true);
wvt.addForcing(WVAdaptiveDamping(wvt));
adaptiveDamping = wvt.forcingWithName("adaptive damping");
wvt.removeAll;

if includeEddy
    x0 = wvt.Lx/2;
    y0 = wvt.Ly/2;
    eddyRadius = 80e3;
    eddyDepth = 300;
    eddySpeed = 0.10;
    verticalStructure = @(z)exp(-(z/eddyDepth/sqrt(2)).^2);
    horizontalStructure = @(x,y)exp(-((x-x0)/eddyRadius).^2-((y-y0)/eddyRadius).^2);
    streamfunction = @(x,y,z)eddySpeed*(eddyRadius/sqrt(2))*exp(1/2)*verticalStructure(z).*(horizontalStructure(x,y)-pi*eddyRadius^2/(wvt.Lx*wvt.Ly));
    wvt.addGeostrophicStreamfunction(streamfunction);
    wvt.A0(wvt.Kh > adaptiveDamping.k_damp) = 0;
end

tidalPeriod = 12.420602*3600;
[~,topographicHeight] = WVBottomWaveGenerationForcing.goffAbyssalHillTopography(wvt,rmsHeight=100,cornerWavenumber=1e-4,minimumWavelength=6e3,randomSeed=2023);
generation = WVBottomWaveGenerationForcing(wvt,topographicHeight=topographicHeight,barotropicVelocityAmplitude=[0.05/sqrt(10); 0],frequency=2*pi/tidalPeriod,rampDuration=tidalPeriod,startTime=0,name="bottom wave generation");
wvt.addForcing(generation);
end

function preflight = verifyInitialConfiguration(wvt,generation,adaptiveDamping,referencePath,caseName)
assert(wvt.Nx == 256 && wvt.Ny == 256 && wvt.Nz == 45 && wvt.Nj == 29)
assert(wvt.Lx == 500e3 && wvt.Ly == 500e3 && wvt.Lz == 2000)
assert(generation.shouldAvoidAdaptiveDamping)
assert(isinf(generation.maximumForcedHorizontalWavenumber))
assert(isinf(generation.maximumForcedVerticalMode))

[mask,components] = generation.spectralGenerationMask();
expectedMask = (logical(wvt.waveComponent.maskAp) | logical(wvt.waveComponent.maskAm)) & adaptiveDamping.damp == 0;
assert(isequal(mask,expectedMask))
assert(isequal(components.adaptiveDamping,adaptiveDamping.damp == 0))

wvt.t = generation.rampDuration;
[Fp,Fm] = generation.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
dampingRegion = adaptiveDamping.damp ~= 0;
assert(~any(Fp(dampingRegion),"all") && ~any(Fm(dampingRegion),"all"))
assert(norm([Fp(~dampingRegion); Fm(~dampingRegion)]) > 0)
wvt.t = 0;

    [reference,referenceFile] = WVTransform.waveVortexTransformFromFile(char(referencePath),iTime=1,shouldReadOnly=true);
    referenceCleanup = onCleanup(@()referenceFile.close());
    referenceGeneration = reference.forcingWithName("bottom wave generation");
    assert(isequal(generation.topographicHeight,referenceGeneration.topographicHeight))
    wvt.Ap = reference.Ap;
    wvt.Am = reference.Am;
    wvt.A0 = reference.A0;
    assert(isequal(wvt.Ap,reference.Ap) && isequal(wvt.Am,reference.Am) && isequal(wvt.A0,reference.A0))
assert(isequal([wvt.Lx wvt.Ly wvt.Lz wvt.Nx wvt.Ny wvt.Nz],[reference.Lx reference.Ly reference.Lz reference.Nx reference.Ny reference.Nz]))
clear referenceCleanup

preflight = struct( ...
    caseName=caseName, ...
    maskCount=nnz(mask), ...
    waveValidCount=nnz(logical(wvt.waveComponent.maskAp) | logical(wvt.waveComponent.maskAm)), ...
    dampingSupportCount=nnz(dampingRegion), ...
    forcedUndampedCoefficientNorm=norm([Fp(~dampingRegion); Fm(~dampingRegion)]), ...
    forcedDampedCoefficientNorm=norm([Fp(dampingRegion); Fm(dampingRegion)]), ...
    terrainMatchesReference=true, ...
    initialStateMatchesReference=true);
end

function verifyPartialTimeAxis(outputPath)
[~,ncfile] = WVTransform.waveVortexTransformFromFile(char(outputPath),iTime=Inf,shouldReadOnly=true);
fileCleanup = onCleanup(@()ncfile.close());
time = ncfile.readVariables("wave-vortex/t");
expectedTime = (0:6*3600:25*86400).';
assert(~isempty(time) && numel(time) <= numel(expectedTime))
assert(isequal(time(:),expectedTime(1:numel(time))))
clear fileCleanup
end

function verifyRestartMask(wvt)
generation = wvt.forcingWithName("bottom wave generation");
adaptiveDamping = wvt.forcingWithName("adaptive damping");
assert(generation.shouldAvoidAdaptiveDamping)
[mask,components] = generation.spectralGenerationMask();
expectedMask = (logical(wvt.waveComponent.maskAp) | logical(wvt.waveComponent.maskAm)) & adaptiveDamping.damp == 0;
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

function postflight = verifyCompletedFile(outputPath)
[restored,ncfile] = WVTransform.waveVortexTransformFromFile(char(outputPath),iTime=Inf,shouldReadOnly=true);
fileCleanup = onCleanup(@()ncfile.close());
time = ncfile.readVariables("wave-vortex/t");
expectedTime = (0:6*3600:25*86400).';
assert(isequal(time(:),expectedTime))
assert(all(isfinite([restored.Ap(:); restored.Am(:); restored.A0(:)])))
generation = restored.forcingWithName("bottom wave generation");
verifyRestartMask(restored)
[mask,~] = generation.spectralGenerationMask();
postflight = struct(recordCount=numel(time),finalTime=time(end),finalWaveEnergy=restored.waveEnergy,finalGeostrophicEnergy=restored.geostrophicEnergy,maskCount=nnz(mask));
clear fileCleanup
end
