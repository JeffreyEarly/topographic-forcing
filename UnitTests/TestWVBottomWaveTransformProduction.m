classdef TestWVBottomWaveTransformProduction < matlab.unittest.TestCase
    % Verify resolution, persistence, and integration on newly supported transforms.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(repositoryRoot)
        end
    end

    methods (Test)
        function resolutionAndRestartPreserveBothForcings(testCase)
            for transformType = ["constantNonhydrostatic" "hydrostatic"]
                source = TestWVBottomWaveTransformCompatibility.createTransform(transformType,[8 6 5],false);
                terrain = TestWVBottomWaveTransformCompatibility.topography(source);
                generation = WVBottomWaveGenerationForcing(source,topographicHeight=terrain,barotropicVelocityAmplitude=[0.04+0.01i; -0.02],frequency=1.31e-4,rampDuration=200,startTime=50,name="production generation");
                scattering = WVBottomWaveScatteringForcing(source,topographicHeight=terrain,name="production scattering");
                source.removeAllForcing();
                source.addForcing(generation);
                source.addForcing(scattering);
                TestWVBottomWaveTransformProduction.setWaveState(source)
                source.t = 317;

                [originalGenerationFp,originalGenerationFm,originalGenerationF0] = generation.addSpectralForcing(source,zeros(size(source.Ap)),zeros(size(source.Am)),zeros(size(source.A0)));
                [originalScatteringFp,originalScatteringFm,originalScatteringF0] = scattering.addSpectralForcing(source,zeros(size(source.Ap)),zeros(size(source.Am)),zeros(size(source.A0)));

                target = source.waveVortexTransformWithResolution([12 10 7]);
                convertedGeneration = target.forcingWithName("production generation");
                convertedScattering = target.forcingWithName("production scattering");
                expectedTerrain = TestWVBottomWaveTransformCompatibility.topography(target);
                testCase.verifyClass(convertedGeneration,"WVBottomWaveGenerationForcing")
                testCase.verifyClass(convertedScattering,"WVBottomWaveScatteringForcing")
                testCase.verifyEqual(convertedGeneration.topographicHeight,expectedTerrain,AbsTol=1e-12)
                testCase.verifyEqual(convertedScattering.topographicHeight,expectedTerrain,AbsTol=1e-12)
                [generationFp,generationFm,generationF0] = convertedGeneration.addSpectralForcing(target,zeros(size(target.Ap)),zeros(size(target.Am)),zeros(size(target.A0)));
                [scatteringFp,scatteringFm,scatteringF0] = convertedScattering.addSpectralForcing(target,zeros(size(target.Ap)),zeros(size(target.Am)),zeros(size(target.A0)));
                testCase.verifyTrue(all(isfinite([generationFp(:);generationFm(:);scatteringFp(:);scatteringFm(:)])))
                testCase.verifyEqual(generationF0,zeros(size(target.A0)))
                testCase.verifyEqual(scatteringF0,zeros(size(target.A0)))

                path = string(tempname)+".nc";
                fileCleanup = onCleanup(@()TestWVBottomWaveTransformProduction.deleteFile(path));
                ncfile = source.writeToFile(path,shouldOverwriteExisting=true);
                ncfile.close();
                [restored,restoredFile] = WVTransform.waveVortexTransformFromFile(path);
                restoredFile.close();
                restoredGeneration = restored.forcingWithName("production generation");
                restoredScattering = restored.forcingWithName("production scattering");
                [restoredGenerationFp,restoredGenerationFm,restoredGenerationF0] = restoredGeneration.addSpectralForcing(restored,zeros(size(restored.Ap)),zeros(size(restored.Am)),zeros(size(restored.A0)));
                [restoredScatteringFp,restoredScatteringFm,restoredScatteringF0] = restoredScattering.addSpectralForcing(restored,zeros(size(restored.Ap)),zeros(size(restored.Am)),zeros(size(restored.A0)));
                testCase.verifyClass(restored,class(source))
                testCase.verifyEqual(restoredGeneration.topographicHeight,terrain)
                testCase.verifyEqual(restoredScattering.topographicHeight,terrain)
                testCase.verifyLessThanOrEqual(testCase.relativeError(restoredGenerationFp,originalGenerationFp),1e-12)
                testCase.verifyLessThanOrEqual(testCase.relativeError(restoredGenerationFm,originalGenerationFm),1e-12)
                testCase.verifyEqual(restoredGenerationF0,originalGenerationF0)
                testCase.verifyLessThanOrEqual(testCase.relativeError(restoredScatteringFp,originalScatteringFp),1e-12)
                testCase.verifyLessThanOrEqual(testCase.relativeError(restoredScatteringFm,originalScatteringFm),1e-12)
                testCase.verifyEqual(restoredScatteringF0,originalScatteringF0)
                clear fileCleanup
            end
        end

        function explicitAntialiasConversionPreservesBothForcings(testCase)
            source = TestWVBottomWaveTransformCompatibility.createTransform("constantNonhydrostatic",[8 6 7],true);
            terrain = TestWVBottomWaveTransformCompatibility.topography(source);
            source.removeAllForcing();
            source.addForcing(WVBottomWaveGenerationForcing(source,topographicHeight=terrain,barotropicVelocityAmplitude=[0.04; -0.01],name="explicit generation"));
            source.addForcing(WVBottomWaveScatteringForcing(source,topographicHeight=terrain,name="explicit scattering"));
            TestWVBottomWaveTransformProduction.setWaveState(source)
            source.t = 211;

            explicit = source.waveVortexTransformWithExplicitAntialiasing();
            generation = explicit.forcingWithName("explicit generation");
            scattering = explicit.forcingWithName("explicit scattering");
            testCase.verifyFalse(explicit.shouldAntialias)
            testCase.verifyEqual(explicit.N0,source.N0)
            testCase.verifyEqual(explicit.isHydrostatic,source.isHydrostatic)
            testCase.verifyClass(generation,"WVBottomWaveGenerationForcing")
            testCase.verifyClass(scattering,"WVBottomWaveScatteringForcing")
            [generationFp,generationFm,generationF0] = generation.addSpectralForcing(explicit,zeros(size(explicit.Ap)),zeros(size(explicit.Am)),zeros(size(explicit.A0)));
            [scatteringFp,scatteringFm,scatteringF0] = scattering.addSpectralForcing(explicit,zeros(size(explicit.Ap)),zeros(size(explicit.Am)),zeros(size(explicit.A0)));
            testCase.verifyTrue(all(isfinite([generationFp(:);generationFm(:);scatteringFp(:);scatteringFm(:)])))
            testCase.verifyEqual(generationF0,zeros(size(explicit.A0)))
            testCase.verifyEqual(scatteringF0,zeros(size(explicit.A0)))
        end

        function constantOutputSupportsBottomDisplacementDiagnostic(testCase)
            wvt = TestWVBottomWaveTransformCompatibility.createTransform("constantNonhydrostatic",[8 6 5],false);
            terrain = TestWVBottomWaveTransformCompatibility.topography(wvt);
            forcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=terrain,name="constant diagnostic scattering");
            wvt.removeAllForcing();
            wvt.addForcing(forcing);
            TestWVBottomWaveTransformProduction.setWaveState(wvt)

            path = string(tempname)+".nc";
            fileCleanup = onCleanup(@()TestWVBottomWaveTransformProduction.deleteFile(path));
            model = WVModel(wvt);
            outputFile = model.createNetCDFFileForModelOutput(path,outputInterval=30,shouldOverwriteExisting=true);
            time = outputFile.outputTimesForIntegrationPeriod(0,60);
            for iTime = 1:numel(time)
                wvt.t = time(iTime);
                outputFile.writeTimeStepToOutputFile(time(iTime));
            end
            outputFile.closeNetCDFFile();

            diagnostics = WVBottomWaveScatteringForcing.bottomDisplacementFromFile(path,forcingName="constant diagnostic scattering");
            testCase.verifyEqual(diagnostics.time,time(:),AbsTol=10*eps(time(end)))
            testCase.verifySize(diagnostics.bottomVelocity,[wvt.Nx wvt.Ny numel(time)])
            testCase.verifySize(diagnostics.bottomDisplacement,[wvt.Nx wvt.Ny numel(time)])
            testCase.verifyTrue(all(isfinite([diagnostics.bottomVelocity(:);diagnostics.bottomDisplacement(:)])))
            testCase.verifyEqual(diagnostics.bottomDisplacement(:,:,1),zeros(wvt.Nx,wvt.Ny),AbsTol=1e-13)
            clear fileCleanup
        end

        function constantEddyGenerationAndScatteringIntegrateTogether(testCase)
            wvt = WVTransformConstantStratification([200e3 200e3 2e3],[16 16 9],N0=sqrt(2e-5),latitude=45,shouldAntialias=true,isHydrostatic=false);
            wvt.addForcing(WVAdaptiveDamping(wvt));

            eddyRadius = 40e3;
            eddyDepth = 300;
            eddySpeed = 0.05;
            horizontalStructure = @(x,y) exp(-((x-wvt.Lx/2)/eddyRadius).^2-((y-wvt.Ly/2)/eddyRadius).^2);
            verticalStructure = @(z) exp(-(z/eddyDepth/sqrt(2)).^2);
            psi = @(x,y,z) eddySpeed*(eddyRadius/sqrt(2))*exp(1/2)*verticalStructure(z).*(horizontalStructure(x,y)-pi*eddyRadius^2/(wvt.Lx*wvt.Ly));
            wvt.addGeostrophicStreamfunction(psi);

            terrain = TestWVBottomWaveTransformCompatibility.topography(wvt);
            wvt.addForcing(WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.03; 0],rampDuration=0,startTime=0));
            wvt.addForcing(WVBottomWaveScatteringForcing(wvt,topographicHeight=terrain));
            model = WVModel(wvt);
            model.integrateToTime(60,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);

            testCase.verifyTrue(all(isfinite([wvt.Ap(:);wvt.Am(:);wvt.A0(:)])))
            testCase.verifyGreaterThan(norm([wvt.Ap(:);wvt.Am(:)]),0)
            expectedNames = ["nonlinear advection" "adaptive damping" "bottom wave generation" "bottom wave scattering"];
            testCase.verifyEqual(sort(reshape(wvt.forcingNames(),[],1)),sort(reshape(expectedNames,[],1)))
        end
    end

    methods (Static, Access = private)
        function setWaveState(wvt)
            previousRandomState = rng;
            randomStateCleanup = onCleanup(@()rng(previousRandomState));
            rng(78312,"twister")
            wvt.Ap = 1e-3*(randn(size(wvt.Ap))+1i*randn(size(wvt.Ap))).*wvt.waveComponent.maskAp;
            wvt.Am = 1e-3*(randn(size(wvt.Am))+1i*randn(size(wvt.Am))).*wvt.waveComponent.maskAm;
            wvt.A0(:) = 0;
            clear randomStateCleanup
        end

        function errorValue = relativeError(actual,expected)
            errorValue = norm(actual(:)-expected(:))/max(norm(expected(:)),eps);
        end

        function deleteFile(path)
            if isfile(path)
                delete(path)
            end
        end
    end
end
