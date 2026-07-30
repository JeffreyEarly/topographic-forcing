classdef TestEddyTideTopographicForcingSimulation < matlab.unittest.TestCase
    % Verify the paired eddy-tide simulation and diagnostics workflow.

    properties
        temporaryDirectory (1,1) string
        eddyFile (1,1) string
        noEddyFile (1,1) string
        scatteringFile (1,1) string
    end

    methods (TestClassSetup)
        function createPairedSimulationOutput(testCase)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(repositoryRoot)
            addpath(fullfile(repositoryRoot,"Examples"))
            testCase.temporaryDirectory = string(tempname);
            mkdir(testCase.temporaryDirectory)
            [~,~,testCase.eddyFile] = EddyTideTopographicForcingSimulation(Nxy=32,maxT=600,outputInterval=300,outputDirectory=testCase.temporaryDirectory,outputFilename="eddy.nc",minimumWavelength=100e3,shouldOverwriteExisting=true,shouldShowIntegrationDiagnostics=false);
            [~,~,testCase.noEddyFile] = EddyTideTopographicForcingSimulation(Nxy=32,includeEddy=false,maxT=600,outputInterval=300,outputDirectory=testCase.temporaryDirectory,outputFilename="no-eddy.nc",minimumWavelength=100e3,shouldOverwriteExisting=true,shouldShowIntegrationDiagnostics=false);
            [~,~,testCase.scatteringFile] = EddyTideTopographicForcingSimulation(Nxy=32,includeEddy=false,shouldUseScattering=true,maxT=600,outputInterval=300,outputDirectory=testCase.temporaryDirectory,outputFilename="scattering.nc",minimumWavelength=100e3,shouldOverwriteExisting=true,shouldShowIntegrationDiagnostics=false);
        end
    end

    methods (TestClassTeardown)
        function removeTemporaryOutput(testCase)
            close(findall(groot,Type="figure"))
            if isfolder(testCase.temporaryDirectory)
                rmdir(testCase.temporaryDirectory,"s")
            end
        end
    end

    methods (Test)
        function pairedRunsHaveMatchedForcingAndExpectedInitialState(testCase)
            [eddyInitial,eddyInitialFile] = WVTransform.waveVortexTransformFromFile(char(testCase.eddyFile),iTime=1,shouldReadOnly=true);
            eddyInitialCleanup = onCleanup(@()eddyInitialFile.close());
            [noEddyInitial,noEddyInitialFile] = WVTransform.waveVortexTransformFromFile(char(testCase.noEddyFile),iTime=1,shouldReadOnly=true);
            noEddyInitialCleanup = onCleanup(@()noEddyInitialFile.close());

            testCase.verifyEqual(eddyInitial.Ap,complex(zeros(size(eddyInitial.Ap))))
            testCase.verifyEqual(eddyInitial.Am,complex(zeros(size(eddyInitial.Am))))
            testCase.verifyGreaterThan(norm(eddyInitial.A0(:)),0)
            testCase.verifyEqual(noEddyInitial.Ap,complex(zeros(size(noEddyInitial.Ap))))
            testCase.verifyEqual(noEddyInitial.Am,complex(zeros(size(noEddyInitial.Am))))
            testCase.verifyEqual(noEddyInitial.A0,complex(zeros(size(noEddyInitial.A0))))

            forcingNames = ["nonlinear advection" "adaptive damping" "bottom wave generation"];
            for iName = 1:numel(forcingNames)
                testCase.verifyEqual(nnz(string(eddyInitial.forcingNames) == forcingNames(iName)),1)
                testCase.verifyEqual(nnz(string(noEddyInitial.forcingNames) == forcingNames(iName)),1)
            end
            testCase.verifyFalse(any(string(eddyInitial.forcingNames) == "bottom wave scattering"))
            testCase.verifyFalse(any(string(noEddyInitial.forcingNames) == "bottom wave scattering"))

            eddyGeneration = eddyInitial.forcingWithName("bottom wave generation");
            noEddyGeneration = noEddyInitial.forcingWithName("bottom wave generation");
            testCase.verifyClass(eddyGeneration,"WVBottomWaveGenerationForcing")
            testCase.verifyEqual(eddyGeneration.topographicHeight,noEddyGeneration.topographicHeight)
            testCase.verifyEqual(eddyGeneration.barotropicVelocityAmplitude,complex([0.05; 0]))
            testCase.verifyEqual(eddyGeneration.rampDuration,12.420602*3600)
            testCase.verifyEqual(eddyGeneration.startTime,0)
            testCase.verifyTrue(eddyGeneration.shouldAvoidAdaptiveDamping)
            testCase.verifyEqual(eddyGeneration.maximumForcedHorizontalWavenumber,Inf)
            testCase.verifyEqual(eddyGeneration.maximumForcedVerticalMode,Inf)

            adaptiveDamping = eddyInitial.forcingWithName("adaptive damping");
            eddyInitial.t = eddyGeneration.rampDuration;
            [Fp,Fm] = eddyGeneration.addSpectralForcing(eddyInitial,zeros(size(eddyInitial.Ap)),zeros(size(eddyInitial.Am)),zeros(size(eddyInitial.A0)));
            dampingRegion = adaptiveDamping.damp ~= 0;
            testCase.verifyEqual(Fp(dampingRegion),zeros(nnz(dampingRegion),1))
            testCase.verifyEqual(Fm(dampingRegion),zeros(nnz(dampingRegion),1))
            testCase.verifyGreaterThan(norm([Fp(~dampingRegion); Fm(~dampingRegion)]),0)

            clear noEddyInitialCleanup eddyInitialCleanup
        end

        function scatteringIsExplicitAndUsesGenerationTerrain(testCase)
            [wvt,ncfile] = WVTransform.waveVortexTransformFromFile(char(testCase.scatteringFile),iTime=1,shouldReadOnly=true);
            cleanup = onCleanup(@()ncfile.close());
            generation = wvt.forcingWithName("bottom wave generation");
            scattering = wvt.forcingWithName("bottom wave scattering");
            testCase.verifyClass(generation,"WVBottomWaveGenerationForcing")
            testCase.verifyClass(scattering,"WVBottomWaveScatteringForcing")
            testCase.verifyEqual(generation.topographicHeight,scattering.topographicHeight)
            testCase.verifyEqual(nnz(string(wvt.forcingNames) == "bottom wave scattering"),1)
            clear cleanup
        end

        function defaultFilenamesIdentifyScatteringConfiguration(testCase)
            [~,~,generationOnlyFile] = EddyTideTopographicForcingSimulation( ...
                Nxy=16,includeEddy=false,maxT=1,outputInterval=1, ...
                outputDirectory=testCase.temporaryDirectory,minimumWavelength=200e3, ...
                shouldOverwriteExisting=true,shouldShowIntegrationDiagnostics=false);
            [~,~,generatedScatteringFile] = EddyTideTopographicForcingSimulation( ...
                Nxy=16,includeEddy=false,shouldUseScattering=true,maxT=1,outputInterval=1, ...
                outputDirectory=testCase.temporaryDirectory,minimumWavelength=200e3, ...
                shouldOverwriteExisting=true,shouldShowIntegrationDiagnostics=false);
            testCase.verifyTrue(contains(generationOnlyFile,"generation-only"))
            testCase.verifyTrue(contains(generatedScatteringFile,"generation-scattering"))
            testCase.verifyNotEqual(generationOnlyFile,generatedScatteringFile)
        end

        function explicitHorizontalDomainControlsResolutionAndFilename(testCase)
            domainSize = 400e3;
            Nxy = 16;
            [~,wvt,outputFile] = EddyTideTopographicForcingSimulation( ...
                Nxy=Nxy,horizontalDomainSize=domainSize,includeEddy=false,maxT=1,outputInterval=1, ...
                outputDirectory=testCase.temporaryDirectory,minimumWavelength=200e3, ...
                shouldOverwriteExisting=true,shouldShowIntegrationDiagnostics=false);
            N2 = @(z)2e-5*ones(size(z));
            expectedNz = WVStratification.verticalResolutionForHorizontalResolution(domainSize,2000,Nxy,N2=N2,latitude=45);
            testCase.verifyEqual(wvt.Lx,domainSize)
            testCase.verifyEqual(wvt.Ly,domainSize)
            testCase.verifyEqual(wvt.Nz,expectedNz)
            testCase.verifyTrue(contains(outputFile,"Lxy400km"))
        end

        function outputDevelopsFiniteWavesAtRequestedCadence(testCase)
            [eddyFinal,eddyFinalFile] = WVTransform.waveVortexTransformFromFile(char(testCase.eddyFile),iTime=Inf,shouldReadOnly=true);
            eddyFinalCleanup = onCleanup(@()eddyFinalFile.close());
            [noEddyFinal,noEddyFinalFile] = WVTransform.waveVortexTransformFromFile(char(testCase.noEddyFile),iTime=Inf,shouldReadOnly=true);
            noEddyFinalCleanup = onCleanup(@()noEddyFinalFile.close());
            testCase.verifyGreaterThan(eddyFinal.waveEnergy,0)
            testCase.verifyGreaterThan(noEddyFinal.waveEnergy,0)
            testCase.verifyTrue(all(isfinite([eddyFinal.Ap(:); eddyFinal.Am(:); eddyFinal.A0(:)])))
            testCase.verifyTrue(all(isfinite([noEddyFinal.Ap(:); noEddyFinal.Am(:); noEddyFinal.A0(:)])))

            eddyTime = eddyFinalFile.readVariables("wave-vortex/t");
            noEddyTime = noEddyFinalFile.readVariables("wave-vortex/t");
            testCase.verifyEqual(eddyTime(:),[0; 300; 600])
            testCase.verifyEqual(noEddyTime(:),eddyTime(:))
            clear noEddyFinalCleanup eddyFinalCleanup
        end

        function standardRestartRestoresForcingsAndContinues(testCase)
            restartFile = fullfile(testCase.temporaryDirectory,"restart.nc");
            copyfile(testCase.eddyFile,restartFile)
            model = WVModel.modelFromFile(char(restartFile));
            testCase.verifyClass(model.wvt,"WVTransformConstantStratification")
            testCase.verifyClass(model.wvt.forcingWithName("bottom wave generation"),"WVBottomWaveGenerationForcing")
            testCase.verifyFalse(any(string(model.wvt.forcingNames) == "bottom wave scattering"))
            model.integrateToTime(900,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.closeNetCDFFile();
            testCase.verifyEqual(model.t,900)

            [restored,restoredFile] = WVTransform.waveVortexTransformFromFile(char(restartFile),iTime=Inf,shouldReadOnly=true);
            restoredCleanup = onCleanup(@()restoredFile.close());
            testCase.verifyEqual(restored.t,900)
            testCase.verifyClass(restored.forcingWithName("bottom wave generation"),"WVBottomWaveGenerationForcing")
            testCase.verifyFalse(any(string(restored.forcingNames) == "bottom wave scattering"))
            clear restoredCleanup
        end

        function diagnosticsAndEnergyFigureUseCommonEddyNormalization(testCase)
            figurePath = fullfile(testCase.temporaryDirectory,"energy.png");
            [figureHandle,energy,diagnosticsFiles] = AnalyzeEddyTideEnergy(testCase.eddyFile,testCase.noEddyFile,figureVisible="off",exportPath=figurePath,shouldOverwriteExisting=true);
            figureCleanup = onCleanup(@()close(figureHandle));

            testCase.verifyTrue(all(isfile(diagnosticsFiles)))
            testCase.verifyTrue(isfile(figurePath))
            testCase.verifyEqual(energy.eddy.geostrophicNormalized(1),1,AbsTol=10*eps)
            testCase.verifyEqual(energy.noEddy.geostrophicNormalized(1),0,AbsTol=10*eps)
            testCase.verifyEqual(energy.eddy.waveNormalized,energy.eddy.wave/energy.normalization,AbsTol=10*eps)
            testCase.verifyEqual(energy.noEddy.waveNormalized,energy.noEddy.wave/energy.normalization,AbsTol=10*eps)
            testCase.verifyEqual(energy.eddy.time(:),[0; 300; 600])
            testCase.verifyEqual(energy.noEddy.time,energy.eddy.time)
            testCase.verifyEqual(string(energy.figurePath),string(figurePath))
            testCase.verifyNumElements(findall(figureHandle,Type="axes"),2)
            testCase.verifyNumElements(findall(figureHandle,Type="line"),4)

            [secondFigure,secondEnergy,secondDiagnosticsFiles] = AnalyzeEddyTideEnergy(testCase.eddyFile,testCase.noEddyFile,figureVisible="off",shouldExport=false);
            secondFigureCleanup = onCleanup(@()close(secondFigure));
            testCase.verifyEqual(secondDiagnosticsFiles,diagnosticsFiles)
            testCase.verifyEqual(secondEnergy.eddy.wave,energy.eddy.wave)
            testCase.verifyEqual(secondEnergy.noEddy.geostrophic,energy.noEddy.geostrophic)
            clear secondFigureCleanup figureCleanup
        end

        function energyEnstrophyFigureUsesCommonEddyNormalization(testCase)
            figurePath = fullfile(testCase.temporaryDirectory,"energy-enstrophy.png");
            [figureHandle,series,diagnosticsFiles] = AnalyzeEddyTideEnergyEnstrophy( ...
                testCase.eddyFile,testCase.noEddyFile,figureVisible="off", ...
                exportPath=figurePath,shouldOverwriteExisting=true);
            figureCleanup = onCleanup(@()close(figureHandle));

            testCase.verifyTrue(all(isfile(diagnosticsFiles)))
            testCase.verifyTrue(isfile(figurePath))
            testCase.verifyEqual(series.eddy.time,series.control.time)
            testCase.verifyEqual(series.eddy.time(:),[0; 300; 600])
            testCase.verifyTrue(isfinite(series.normalization.energy))
            testCase.verifyGreaterThan(series.normalization.energy,0)
            testCase.verifyTrue(isfinite(series.normalization.enstrophy))
            testCase.verifyGreaterThan(series.normalization.enstrophy,0)
            testCase.verifyEqual(series.normalization.energy,series.eddy.raw.energy.total(1))
            testCase.verifyEqual(series.normalization.enstrophy,series.eddy.raw.apvEnstrophy(1))
            testCase.verifyEqual(series.eddy.normalized.energy.total(1),1,AbsTol=10*eps)
            testCase.verifyEqual(series.eddy.normalized.apvEnstrophy(1),1,AbsTol=10*eps)
            testCase.verifyEqual(series.control.normalized.energy.total, ...
                series.control.raw.energy.total/series.normalization.energy)
            testCase.verifyEqual(series.control.normalized.apvEnstrophy, ...
                series.control.raw.apvEnstrophy/series.normalization.enstrophy)
            eddyNormalizedEnergy = [ ...
                series.eddy.normalized.energy.total; ...
                series.eddy.normalized.energy.wave; ...
                series.eddy.normalized.energy.geostrophic; ...
                series.eddy.normalized.energy.geostrophicKinetic; ...
                series.eddy.normalized.energy.geostrophicPotential];
            controlNormalizedEnergy = [ ...
                series.control.normalized.energy.total; ...
                series.control.normalized.energy.wave; ...
                series.control.normalized.energy.geostrophic; ...
                series.control.normalized.energy.geostrophicKinetic; ...
                series.control.normalized.energy.geostrophicPotential];
            testCase.verifyTrue(all(isfinite(eddyNormalizedEnergy)))
            testCase.verifyTrue(all(isfinite(controlNormalizedEnergy)))
            testCase.verifyTrue(all(isfinite(series.eddy.normalized.apvEnstrophy)))
            testCase.verifyTrue(all(isfinite(series.control.normalized.apvEnstrophy)))
            testCase.verifyEqual(string(series.figurePath),string(figurePath))
            testCase.verifyNumElements(findall(figureHandle,Type="axes"),2)

            reopenFiles = [testCase.eddyFile; testCase.noEddyFile; diagnosticsFiles];
            for iFile = 1:numel(reopenFiles)
                writableFile = NetCDFFile(char(reopenFiles(iFile)),shouldReadOnly=false);
                writableCleanup = onCleanup(@()writableFile.close());
                clear writableCleanup writableFile
            end

            testCase.verifyError(@()AnalyzeEddyTideEnergyEnstrophy( ...
                testCase.eddyFile,testCase.noEddyFile,figureVisible="off", ...
                exportPath=figurePath), ...
                "AnalyzeEddyTideEnergyEnstrophy:ExportFileExists")

            [secondFigure,secondSeries,secondDiagnosticsFiles] = AnalyzeEddyTideEnergyEnstrophy( ...
                testCase.eddyFile,testCase.noEddyFile,figureVisible="off",shouldExport=false);
            secondFigureCleanup = onCleanup(@()close(secondFigure));
            testCase.verifyEqual(secondDiagnosticsFiles,diagnosticsFiles)
            testCase.verifyEqual(secondSeries.eddy.raw.energy.total,series.eddy.raw.energy.total)
            testCase.verifyEqual(secondSeries.control.raw.apvEnstrophy,series.control.raw.apvEnstrophy)
            clear secondFigureCleanup figureCleanup
        end

        function basicFigureAnalysisExportsPairedDiagnostics(testCase)
            exportDirectory = fullfile(testCase.temporaryDirectory,"basic-figures");
            [figures,summary] = AnalyzeEddyTideBasicFigures(testCase.eddyFile,testCase.noEddyFile, ...
                figureVisible="off",exportDirectory=exportDirectory,exportPrefix="short-pair",shouldOverwriteExisting=true);
            cleanup = onCleanup(@()close(figures));
            testCase.verifyNumElements(figures,4)
            testCase.verifyTrue(all(isfile(summary.figurePaths)))
            testCase.verifyTrue(isfile(summary.summaryPath))
            testCase.verifyEqual(summary.eddy.time,summary.noEddy.time)
            testCase.verifyEqual(summary.energy.eddy.geostrophicNormalized(1),1,AbsTol=10*eps)
            testCase.verifyEqual(size(summary.eddy.vorticity.waveZetaOverF),[32 32])
            testCase.verifyEqual(size(summary.noEddy.vorticity.geostrophicZetaOverF),[32 32])
            testCase.verifyTrue(all(isfinite(summary.eddy.spectrum.horizontalFraction)))
            testCase.verifyTrue(all(isfinite(summary.noEddy.spectrum.modeFraction)))
            testCase.verifyGreaterThan(summary.comparison.waveVorticityColorLimit,0)
            testCase.verifyGreaterThan(summary.comparison.geostrophicVorticityColorLimit,0)
            testCase.verifyLessThanOrEqual(max(abs(summary.eddy.vorticity.waveZetaOverF),[],"all"),summary.comparison.waveVorticityColorLimit)
            testCase.verifyLessThanOrEqual(max(abs(summary.noEddy.vorticity.geostrophicZetaOverF),[],"all"),summary.comparison.geostrophicVorticityColorLimit)
            clear cleanup
        end

        function adaptiveDampingFigureClosesExactSpatialWork(testCase)
            figurePath = fullfile(testCase.temporaryDirectory,"adaptive-damping.png");
            timeRangeDays = [0 600]/86400;
            [figureHandle,analysis,diagnosticsFile] = AnalyzeEddyTideAdaptiveDamping( ...
                testCase.eddyFile,timeRangeDays=timeRangeDays,radialBinCount=8, ...
                figureVisible="off",exportPath=figurePath,shouldOverwriteExisting=true);
            figureCleanup = onCleanup(@()close(figureHandle));

            testCase.verifyTrue(isfile(figurePath))
            testCase.verifyTrue(isfile(diagnosticsFile))
            testCase.verifyEqual(analysis.time(:),[0; 300; 600])
            testCase.verifyEqual(analysis.timeRangeDays,timeRangeDays)
            testCase.verifyEqual(sum(analysis.timeWeights),1,AbsTol=10*eps)
            testCase.verifyEqual(analysis.eddyOrientationFactor,-1)
            testCase.verifySize(analysis.azimuthalVelocity,[8 numel(analysis.depth)])
            testCase.verifySize(analysis.adaptiveDampingEnergyRemovalRate,[8 numel(analysis.depth)])
            testCase.verifyTrue(all(isfinite(analysis.azimuthalVelocity),"all"))
            testCase.verifyTrue(all(isfinite(analysis.adaptiveDampingEnergyRemovalRate),"all"))
            testCase.verifyTrue(all(analysis.radialCounts > 0))
            testCase.verifyLessThanOrEqual(analysis.validation.maximumRelativeFluxError,1e-10)
            testCase.verifyLessThanOrEqual(max(analysis.validation.absoluteFluxError),1e-14)
            testCase.verifyLessThanOrEqual(max(analysis.validation.diagnosticWork),1e-14)
            testCase.verifyEqual(string(analysis.figurePath),string(figurePath))
            testCase.verifyNumElements(findall(figureHandle,Type="axes"),2)
            testCase.verifyNumElements(findall(figureHandle,Type="colorbar"),2)

            reopenFiles = [testCase.eddyFile; diagnosticsFile];
            for iFile = 1:numel(reopenFiles)
                writableFile = NetCDFFile(char(reopenFiles(iFile)),shouldReadOnly=false);
                writableCleanup = onCleanup(@()writableFile.close());
                clear writableCleanup writableFile
            end

            testCase.verifyError(@()AnalyzeEddyTideAdaptiveDamping( ...
                testCase.eddyFile,timeRangeDays=timeRangeDays,radialBinCount=8, ...
                figureVisible="off",exportPath=figurePath), ...
                "AnalyzeEddyTideAdaptiveDamping:ExportFileExists")

            [secondFigure,secondAnalysis,secondDiagnosticsFile] = AnalyzeEddyTideAdaptiveDamping( ...
                testCase.eddyFile,timeRangeDays=timeRangeDays,radialBinCount=8, ...
                figureVisible="off",exportPath=figurePath,shouldOverwriteExisting=true);
            secondFigureCleanup = onCleanup(@()close(secondFigure));
            testCase.verifyEqual(secondDiagnosticsFile,diagnosticsFile)
            testCase.verifyEqual(secondAnalysis.azimuthalVelocity,analysis.azimuthalVelocity)
            clear secondFigureCleanup figureCleanup
        end

        function forcingBudgetsAndTriadsUseStandardDiagnostics(testCase)
            exportDirectory = fullfile(testCase.temporaryDirectory,"budgets");
            [figures,budget,diagnosticsFiles] = AnalyzeEddyTideBudgets( ...
                testCase.eddyFile,testCase.noEddyFile,figureVisible="off", ...
                exportDirectory=exportDirectory,shouldOverwriteExisting=true);
            cleanup = onCleanup(@()close(figures));
            testCase.verifyNumElements(figures,3)
            testCase.verifyTrue(all(isfile(diagnosticsFiles)))
            testCase.verifyTrue(all(isfile(budget.figurePaths)))
            testCase.verifyEqual(budget.eddy.time,budget.control.time)
            testCase.verifyEqual(budget.excess.state.wave, ...
                budget.eddy.state.wave-budget.control.state.wave)
            testCase.verifyTrue(isfield(budget.eddy.closure,"quadraticEnstrophy"))
            testCase.verifyTrue(isfield(budget.eddy.closure,"exactEnstrophy"))

            forcingNames = string({budget.control.forcing.energy.fancyName});
            generationIndex = forcingNames == "bottom wave generation";
            testCase.verifyEqual(nnz(generationIndex),1)
            testCase.verifyFalse(any(forcingNames == "bottom wave scattering"))
            testCase.verifyEqual(budget.control.forcing.energy(generationIndex).te_g, ...
                zeros(size(budget.control.time)))
            testCase.verifyEqual(budget.control.forcing.enstrophy(generationIndex).Z0, ...
                zeros(size(budget.control.time)))
            testCase.verifyNumElements(budget.control.triads.energy,4)
            testCase.verifyNumElements(budget.control.triads.enstrophy,4)
            clear cleanup
        end

        function scatteringBenchmarkRecoversPerturbationExponents(testCase)
            benchmark = GoffScatteringPerturbationBenchmark(testCase.noEddyFile, ...
                rmsHeights=[20 10 5],timeIndices=2:3,referenceTime=300, ...
                numberOfIntegrationPeriods=0.002,shouldMakeFigure=false);
            testCase.verifyEqual(benchmark.workExponent,1,AbsTol=1e-10)
            testCase.verifyEqual(benchmark.residualExponent,2,AbsTol=1e-10)
            testCase.verifyLessThan(benchmark.strictIdentityError,1e-10)
            testCase.verifyTrue(all(isfinite(benchmark.fractionalEnergyDrift)))
            testCase.verifyEqual(size(benchmark.flatWork),[2 3])
        end

        function invalidRampAndExportCollisionFailClearly(testCase)
            testCase.verifyError(@()EddyTideTopographicForcingSimulation(Nxy=32,maxT=600,outputDirectory=testCase.temporaryDirectory,outputFilename="invalid-ramp.nc",minimumWavelength=100e3,rampDuration=-1),"EddyTideTopographicForcingSimulation:InvalidRampDuration")
            testCase.verifyError(@()EddyTideTopographicForcingSimulation(Nxy=32,horizontalDomainSize=0,maxT=600,outputDirectory=testCase.temporaryDirectory,outputFilename="invalid-domain.nc",minimumWavelength=100e3),"EddyTideTopographicForcingSimulation:InvalidHorizontalDomainSize")
            testCase.verifyError(@()EddyTideTopographicForcingSimulation(Nxy=32,horizontalDomainSize=Inf,maxT=600,outputDirectory=testCase.temporaryDirectory,outputFilename="invalid-domain.nc",minimumWavelength=100e3),"EddyTideTopographicForcingSimulation:InvalidHorizontalDomainSize")
            testCase.verifyError(@()EddyTideTopographicForcingSimulation(Nxy=32,maxT=600,outputDirectory=testCase.temporaryDirectory,outputFilename="eddy.nc",minimumWavelength=100e3),"EddyTideTopographicForcingSimulation:OutputFileExists")

            figurePath = fullfile(testCase.temporaryDirectory,"collision.png");
            firstFigure = AnalyzeEddyTideEnergy(testCase.eddyFile,testCase.noEddyFile,figureVisible="off",exportPath=figurePath,shouldOverwriteExisting=true);
            firstFigureCleanup = onCleanup(@()close(firstFigure));
            testCase.verifyError(@()AnalyzeEddyTideEnergy(testCase.eddyFile,testCase.noEddyFile,figureVisible="off",exportPath=figurePath),"AnalyzeEddyTideEnergy:ExportFileExists")
            clear firstFigureCleanup
        end
    end
end
