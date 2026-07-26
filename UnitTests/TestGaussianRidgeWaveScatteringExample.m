classdef TestGaussianRidgeWaveScatteringExample < matlab.unittest.TestCase
    % Verify the Gaussian-ridge scattering example and movie renderer.

    properties
        HalfHeight
        FullHeight
        LooseTolerance
    end

    methods (TestClassSetup)
        function addRepositoryAndRunExamples(testCase)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(repositoryRoot)
            addpath(fullfile(repositoryRoot,"Examples"))
            common = {"resolution",[32 4 7],"numberOfOutputTimes",9,"shouldAntialias",false,"shouldMakeFigures",false};
            testCase.HalfHeight = GaussianRidgeWaveScatteringExample(common{:},heightRatio=0.025,criticality=0.025,relativeTolerance=1e-8);
            testCase.FullHeight = GaussianRidgeWaveScatteringExample(common{:},heightRatio=0.05,criticality=0.05,relativeTolerance=1e-8);
            testCase.LooseTolerance = GaussianRidgeWaveScatteringExample(common{:},heightRatio=0.05,criticality=0.05,relativeTolerance=1e-6);
        end
    end

    methods (Test)
        function geometryFlatControlAndScattering(testCase)
            result = testCase.FullHeight;
            configuration = result.configuration;
            testCase.verifyEqual(max(result.topographicHeight,[],"all"),configuration.ridgeHeight,RelTol=1e-13)
            testCase.verifyEqual(configuration.heightRatio,0.05,RelTol=1e-13)
            testCase.verifyEqual(configuration.ridgeWidth,configuration.ridgeHeight/(configuration.requestedCriticality*configuration.mu*sqrt(exp(1))),RelTol=1e-13)
            testCase.verifyGreaterThan(configuration.pointsAcrossFWHM,3.5)
            testCase.verifyGreaterThan(configuration.discreteCriticality,0)
            testCase.verifyLessThanOrEqual(configuration.discreteCriticality,configuration.requestedCriticality*(1+1e-12))
            testCase.verifyEqual(result.flatReferenceDiagnostics.scatteredCoefficientNorm,zeros(size(result.time)),AbsTol=1e-14)
            testCase.verifyEqual(result.flatReferenceDiagnostics.leftwardFirstModeEnergy,zeros(size(result.time)),AbsTol=1e-14)
            testCase.verifyEqual(result.flatReferenceDiagnostics.higherModeEnergy,zeros(size(result.time)),AbsTol=1e-14)
            testCase.verifyGreaterThan(result.diagnostics.leftwardFirstModeEnergy(end)+result.diagnostics.higherModeEnergy(end),0)
            testCase.verifyEqual(result.diagnostics.balancedEnergy,zeros(size(result.time)))
            testCase.verifyEqual(result.diagnostics.qgpvNorm,zeros(size(result.time)))
            testCase.verifyEqual(result.integration.integrator,"adaptive ode78")
        end

        function firstOrderScalingAndAdaptiveConvergence(testCase)
            half = testCase.HalfHeight.diagnostics;
            full = testCase.FullHeight.diagnostics;
            amplitudeRatio = full.scatteredCoefficientNorm(end)/half.scatteredCoefficientNorm(end);
            scatteredEnergyRatio = (full.leftwardFirstModeEnergy(end)+full.higherModeEnergy(end))/(half.leftwardFirstModeEnergy(end)+half.higherModeEnergy(end));
            energyResidualRatio = max(abs(full.relativeFirstOrderEnergyChange))/max(abs(half.relativeFirstOrderEnergyChange));
            testCase.verifyGreaterThanOrEqual(amplitudeRatio,1.8)
            testCase.verifyLessThanOrEqual(amplitudeRatio,2.2)
            testCase.verifyGreaterThanOrEqual(scatteredEnergyRatio,3)
            testCase.verifyLessThanOrEqual(scatteredEnergyRatio,5)
            testCase.verifyGreaterThanOrEqual(energyResidualRatio,3)
            testCase.verifyLessThanOrEqual(energyResidualRatio,5)

            tight = testCase.FullHeight.coefficients;
            loose = testCase.LooseTolerance.coefficients;
            scale = norm([tight.Ap(:); tight.Am(:)]);
            coefficientError = norm([tight.Ap(:)-loose.Ap(:); tight.Am(:)-loose.Am(:)])/scale;
            testCase.verifyLessThanOrEqual(coefficientError,1e-4)
        end

        function optionalOutputPersistenceAndMovie(testCase)
            outputPath = string(tempname)+".nc";
            videoPath = erase(outputPath,".nc")+".mp4";
            posterPath = erase(outputPath,".nc")+".png";
            cleanup = onCleanup(@()TestGaussianRidgeWaveScatteringExample.deleteFiles([outputPath videoPath posterPath]));
            result = GaussianRidgeWaveScatteringExample(resolution=[24 4 5],numberOfOutputTimes=5,relativeTolerance=1e-7,shouldAntialias=false,outputPath=outputPath,shouldOverwriteExisting=true,shouldMakeFigures=false);
            testCase.verifyTrue(isfile(outputPath))
            testCase.verifyEqual(result.outputPath,outputPath)
            testCase.verifyError(@()GaussianRidgeWaveScatteringExample(resolution=[24 4 5],numberOfOutputTimes=5,relativeTolerance=1e-7,shouldAntialias=false,outputPath=outputPath,shouldMakeFigures=false),"GaussianRidgeWaveScatteringExample:OutputExists")

            [wvt,ncfile] = WVTransformBoussinesq.waveVortexTransformFromFile(char(outputPath),iTime=2,shouldReadOnly=true);
            fileCleanup = onCleanup(@()ncfile.close());
            forcing = wvt.forcingWithName("Gaussian-ridge wave scattering");
            testCase.verifyClass(forcing,"WVBottomWaveScatteringForcing")
            group = ncfile.groupWithName("wave-vortex");
            testCase.verifyEqual(reshape(group.readVariables("t"),[],1),result.time,AbsTol=10*eps(result.time(end)))
            incomingF0 = complex(ones(size(wvt.A0)));
            [~,~,returnedF0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),incomingF0);
            testCase.verifyEqual(returnedF0,incomingF0)
            clear fileCleanup

            originalVisibility = get(groot,"defaultFigureVisible");
            visibilityCleanup = onCleanup(@()set(groot,"defaultFigureVisible",originalVisibility));
            set(groot,"defaultFigureVisible","off")
            movie = GaussianRidgeWaveScatteringMovie(outputPath,videoPath=videoPath,posterPath=posterPath,frameRate=5,quality=75,shouldOverwriteExisting=true);
            testCase.verifyEqual(movie.numberOfFrames,5)
            testCase.verifyEqual(movie.readableFrames,5)
            testCase.verifyEqual(movie.field,"u")
            testCase.verifyEqual(movie.verticalExaggeration,20)
            testCase.verifyTrue(isfile(movie.videoPath))
            testCase.verifyTrue(isfile(movie.posterPath))
            testCase.verifyGreaterThan(dir(movie.videoPath).bytes,0)
            testCase.verifyGreaterThan(dir(movie.posterPath).bytes,0)
            testCase.verifyGreaterThan(prod(movie.frameDimensions),0)
            testCase.verifyGreaterThan(movie.colorLimits(2),0)
            expectedWetMask = reshape(movie.sectionZ,1,[]) >= (-wvt.Lz+forcing.topographicHeight(:,1));
            testCase.verifyEqual(movie.wetMask,expectedWetMask)
            testCase.verifyEqual(movie.numberOfMaskedGridCells,nnz(~expectedWetMask))
            testCase.verifyGreaterThan(movie.numberOfMaskedGridCells,0)
            wetFieldMaximum = 0;
            for iTime = 1:numel(result.time)
                [u,~,~,~] = wvt.transformWaveVortexToUVWEta(result.coefficients.Ap(:,:,iTime),result.coefficients.Am(:,:,iTime),result.coefficients.A0(:,:,iTime),result.time(iTime));
                nativeSection = squeeze(u(:,1,:));
                section = interp1(reshape(wvt.z,[],1),nativeSection.',movie.sectionZ,"linear").';
                wetFieldMaximum = max(wetFieldMaximum,max(abs(section(expectedWetMask))));
            end
            testCase.verifyEqual(movie.colorLimits,wetFieldMaximum*[-1 1],RelTol=1e-12)
            testCase.verifyEqual(reshape(mean(movie.energyProfile,1),[],1),result.diagnostics.firstOrderEnergy,RelTol=1e-10)
            testCase.verifyLessThan(movie.energyProfileConsistencyError,1e-10)
            testCase.verifyError(@()GaussianRidgeWaveScatteringMovie(outputPath,videoPath=videoPath,posterPath=posterPath),"GaussianRidgeWaveScatteringMovie:OutputExists")
            clear visibilityCleanup cleanup
        end

        function staticFiguresRender(testCase)
            originalVisibility = get(groot,"defaultFigureVisible");
            visibilityCleanup = onCleanup(@()set(groot,"defaultFigureVisible",originalVisibility));
            set(groot,"defaultFigureVisible","off")
            result = GaussianRidgeWaveScatteringExample(resolution=[32 4 7],numberOfOutputTimes=5,relativeTolerance=1e-7,shouldAntialias=true,shouldMakeFigures=true);
            figureCleanup = onCleanup(@()close(result.figureHandles(isgraphics(result.figureHandles))));
            testCase.verifyNumElements(result.figureHandles,2)
            testCase.verifyTrue(all(isgraphics(result.figureHandles,"figure")))
            testCase.verifyTrue(result.configuration.shouldAntialias)
            for figureHandle = reshape(result.figureHandles,1,[])
                axesHandles = findall(figureHandle,Type="axes");
                testCase.verifyGreaterThan(numel(axesHandles),0)
                lineHandles = findall(figureHandle,Type="line");
                imageHandles = findall(figureHandle,Type="image");
                testCase.verifyGreaterThan(numel(lineHandles)+numel(imageHandles),0)
                for lineHandle = reshape(lineHandles,1,[])
                    testCase.verifyTrue(all(isfinite([lineHandle.XData(:); lineHandle.YData(:)])))
                end
                for imageHandle = reshape(imageHandles,1,[])
                    testCase.verifyTrue(all(isfinite(imageHandle.CData(:))))
                end
            end
            clear figureCleanup visibilityCleanup
        end
    end

    methods (Static, Access = private)
        function deleteFiles(paths)
            for path = reshape(paths,1,[])
                if isfile(path)
                    delete(path)
                end
            end
        end
    end
end
