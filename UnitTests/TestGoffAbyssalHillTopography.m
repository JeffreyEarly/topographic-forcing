classdef TestGoffAbyssalHillTopography < matlab.unittest.TestCase
    % Verify deterministic Goff-spectrum topography generation.

    properties
        wvt
    end

    methods (TestClassSetup)
        function createTransform(testCase)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(repositoryRoot);
            testCase.wvt = WVTransformHydrostatic([750e3 750e3 2e3],[32 32 7],N2=@(z) 2e-5*ones(size(z)),latitude=45,shouldAntialias=false);
        end
    end

    methods (Test)
        function statisticsShapeAndSignConvention(testCase)
            [virtualDepth,topographicHeight,diagnostics] = WVBottomWaveGenerationForcing.goffAbyssalHillTopography(testCase.wvt,rmsHeight=100,minimumWavelength=60e3,randomSeed=2023);

            testCase.verifySize(virtualDepth,[testCase.wvt.Nx testCase.wvt.Ny])
            testCase.verifySize(topographicHeight,[testCase.wvt.Nx testCase.wvt.Ny])
            testCase.verifyTrue(isreal(topographicHeight))
            testCase.verifyLessThanOrEqual(abs(mean(topographicHeight,"all")),1e-12*diagnostics.rmsHeight)
            testCase.verifyEqual(sqrt(mean(topographicHeight.^2,"all")),100,RelTol=1e-13)
            testCase.verifyEqual(virtualDepth+topographicHeight,testCase.wvt.Lz*ones(size(virtualDepth)),AbsTol=1e-12)
            testCase.verifyGreaterThan(min(virtualDepth,[],"all"),0)

            slopeMagnitude = hypot(testCase.wvt.diffX(topographicHeight),testCase.wvt.diffY(topographicHeight));
            testCase.verifyEqual(diagnostics.rmsSlope,sqrt(mean(slopeMagnitude.^2,"all")),RelTol=1e-13)
            testCase.verifyEqual(diagnostics.maximumSlope,max(slopeMagnitude,[],"all"),RelTol=1e-13)
        end

        function seedIsolationAndRepeatability(testCase)
            rng(99,"twister")
            globalStateBefore = rng;
            [virtualDepth1,topographicHeight1,diagnostics1] = WVBottomWaveGenerationForcing.goffAbyssalHillTopography(testCase.wvt,minimumWavelength=60e3,randomSeed=17);
            globalStateAfter = rng;
            [virtualDepth2,topographicHeight2,diagnostics2] = WVBottomWaveGenerationForcing.goffAbyssalHillTopography(testCase.wvt,minimumWavelength=60e3,randomSeed=17);
            [~,topographicHeight3] = WVBottomWaveGenerationForcing.goffAbyssalHillTopography(testCase.wvt,minimumWavelength=60e3,randomSeed=18);

            testCase.verifyEqual(globalStateAfter,globalStateBefore)
            testCase.verifyEqual(virtualDepth2,virtualDepth1)
            testCase.verifyEqual(topographicHeight2,topographicHeight1)
            testCase.verifyEqual(diagnostics2.fourierCoefficients,diagnostics1.fourierCoefficients)
            testCase.verifyNotEqual(topographicHeight3,topographicHeight1)
            testCase.verifyEqual(sqrt(mean(topographicHeight3.^2,"all")),diagnostics1.rmsHeight,RelTol=1e-13)
        end

        function spectralShapeAndCutoff(testCase)
            [~,~,diagnostics] = WVBottomWaveGenerationForcing.goffAbyssalHillTopography(testCase.wvt,rmsHeight=80,cornerWavenumber=1e-4,minimumWavelength=60e3,randomSeed=7);
            [K,L] = ndgrid(testCase.wvt.k_dft,testCase.wvt.l_dft);
            Kh = hypot(K,L);
            realizedPowerSpectrum = testCase.wvt.Lx*testCase.wvt.Ly*abs(diagnostics.fourierCoefficients).^2;
            targetPowerSpectrum = 4*pi*80^2/(1e-4)^2.*(1+(Kh/1e-4).^2).^(-2);
            retained = Kh > 0 & Kh <= diagnostics.cutoffWavenumber;
            excluded = Kh > diagnostics.cutoffWavenumber;
            spectralRatio = realizedPowerSpectrum(retained)./targetPowerSpectrum(retained);

            testCase.verifyLessThan(std(spectralRatio)/mean(spectralRatio),1e-10)
            testCase.verifyLessThanOrEqual(max(abs(diagnostics.fourierCoefficients(excluded))),1e-12*max(abs(diagnostics.fourierCoefficients),[],"all"))
            testCase.verifyGreaterThan(nnz(diagnostics.radialModeCount),0)
            testCase.verifyEqual(diagnostics.minimumWavelength,60e3)
            testCase.verifyEqual(diagnostics.randomSeed,7)
        end

        function invalidConfigurations(testCase)
            minimumResolvedWavelength = 2*pi/(pi/testCase.wvt.effectiveHorizontalGridResolution());
            testCase.verifyError(@()WVBottomWaveGenerationForcing.goffAbyssalHillTopography(testCase.wvt,minimumWavelength=0.9*minimumResolvedWavelength),"WVBottomWaveGenerationForcing:UnresolvedCutoff")
            testCase.verifyError(@()WVBottomWaveGenerationForcing.goffAbyssalHillTopography(testCase.wvt,minimumWavelength=2*testCase.wvt.Lx),"WVBottomWaveGenerationForcing:EmptySpectrum")
            testCase.verifyError(@()WVBottomWaveGenerationForcing.goffAbyssalHillTopography(testCase.wvt,rmsHeight=2*testCase.wvt.Lz,minimumWavelength=60e3),"WVBottomWaveGenerationForcing:NonpositiveVirtualDepth")
            testCase.verifyError(@()WVBottomWaveGenerationForcing.goffAbyssalHillTopography(testCase.wvt,minimumWavelength=60e3,randomSeed=double(intmax("uint32"))+1),"WVBottomWaveGenerationForcing:InvalidRandomSeed")
        end
    end
end
