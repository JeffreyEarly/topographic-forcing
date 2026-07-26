function benchmark = BottomWaveGenerationPerformanceBenchmark(options)
% Profile construction and runtime application at several resolutions.
%
% Timing results are intended for comparisons on one machine rather than
% as platform-independent acceptance thresholds.
%
% - Topic: Examples
% - Declaration: benchmark = BottomWaveGenerationPerformanceBenchmark(options)
% - Parameter options.resolutions: one transform resolution per row
% - Parameter options.numberOfApplications: minimum forcing applications used by `timeit`
% - Returns benchmark: table of construction, storage, and application measurements
arguments
    options.resolutions (:,3) double {mustBeInteger,mustBePositive} = [16 16 9; 32 32 17; 64 64 33]
    options.numberOfApplications (1,1) double {mustBeInteger,mustBePositive} = 20
end

domainSize = [240e3 240e3 2e3];
N2Function = @(z) 2e-5*exp(z/4000);
numberOfResolutions = size(options.resolutions,1);
constructionSeconds = zeros(numberOfResolutions,1);
applicationSeconds = zeros(numberOfResolutions,1);
estimatedStoredBytes = zeros(numberOfResolutions,1);
numberOfWaveCoefficients = zeros(numberOfResolutions,1);

for iResolution = 1:numberOfResolutions
    wvt = WVTransformBoussinesq(domainSize,options.resolutions(iResolution,:),N2=N2Function,latitude=45,shouldAntialias=true);
    [x,y] = ndgrid(wvt.x,wvt.y);
    topographicHeight = 40*cos(2*pi*x/wvt.Lx)+25*sin(2*pi*y/wvt.Ly)+15*cos(2*pi*(2*x/wvt.Lx+y/wvt.Ly));
    constructionSeconds(iResolution) = timeit(@()constructForcing(wvt,topographicHeight));
    forcing = constructForcing(wvt,topographicHeight);
    Fp = complex(zeros(size(wvt.Ap)));
    Fm = complex(zeros(size(wvt.Am)));
    F0 = complex(zeros(size(wvt.A0)));
    applicationSeconds(iResolution) = timeit(@()repeatedApplication(forcing,wvt,Fp,Fm,F0,options.numberOfApplications))/options.numberOfApplications;
    numberOfWaveCoefficients(iResolution) = nnz(wvt.waveComponent.maskAp)+nnz(wvt.waveComponent.maskAm);
    estimatedStoredBytes(iResolution) = 4*16*numel(wvt.Ap)+3*8*wvt.Nx*wvt.Ny;
end

benchmark = table(options.resolutions(:,1),options.resolutions(:,2),options.resolutions(:,3),numberOfWaveCoefficients,constructionSeconds,applicationSeconds,estimatedStoredBytes,VariableNames=["Nx","Ny","Nz","WaveCoefficients","ConstructionSeconds","ApplicationSeconds","EstimatedStoredBytes"]);
disp(benchmark)
end

function forcing = constructForcing(wvt,topographicHeight)
forcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=topographicHeight,barotropicVelocityAmplitude=[0.05+0.01i; -0.02]);
end

function checksum = repeatedApplication(forcing,wvt,Fp,Fm,F0,numberOfApplications)
checksum = 0;
for iApplication = 1:numberOfApplications
    [FpResult,FmResult,F0Result] = forcing.addSpectralForcing(wvt,Fp,Fm,F0);
    checksum = checksum+sum(abs(FpResult(:)))+sum(abs(FmResult(:)))+sum(abs(F0Result(:)));
end
end
