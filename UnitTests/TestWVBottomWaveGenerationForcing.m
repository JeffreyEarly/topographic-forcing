classdef TestWVBottomWaveGenerationForcing < matlab.unittest.TestCase
    % Verify prescribed bottom wave generation and model integration.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(repositoryRoot);
        end
    end

    methods (Test)
        function constructorValidatesScientificContract(testCase)
            for shouldAntialias = [false true]
                wvt = TestWVBottomWaveGenerationForcing.createTransform(shouldAntialias);
                terrains = {
                    zeros(wvt.Nx,wvt.Ny)
                    100*ones(wvt.Nx,wvt.Ny)
                    TestWVBottomWaveGenerationForcing.sinusoidalTopography(wvt,50)
                    TestWVBottomWaveGenerationForcing.bandLimitedTopography(wvt)
                    };
                for iTerrain = 1:numel(terrains)
                    forcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrains{iTerrain},barotropicVelocityAmplitude=[0.05; 0]);
                    testCase.verifyEqual(forcing.topographicHeight,terrains{iTerrain})
                end
            end

            wvt = TestWVBottomWaveGenerationForcing.createTransform(false);
            terrain = zeros(wvt.Nx,wvt.Ny);
            forcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05+0.01i; -0.02i],name="terrain waves");
            testCase.verifyEqual(forcing.name,"terrain waves")
            testCase.verifyEqual(forcing.frequency,2*pi/(12.4206012*3600),"RelTol",10*eps)
            testCase.verifyEqual(forcing.rampDuration,0)
            testCase.verifyEqual(forcing.startTime,wvt.t)

            testCase.verifyError(@()WVBottomWaveGenerationForcing(wvt,topographicHeight=zeros(wvt.Nx-1,wvt.Ny),barotropicVelocityAmplitude=[0.05; 0]),"WVBottomWaveGenerationForcing:InvalidTopographicHeightSize")
            testCase.verifyError(@()WVBottomWaveGenerationForcing(wvt,topographicHeight=1i*ones(wvt.Nx,wvt.Ny),barotropicVelocityAmplitude=[0.05; 0]),"WVBottomWaveGenerationForcing:InvalidTopographicHeight")
            invalidTerrain = terrain;
            invalidTerrain(1) = NaN;
            testCase.verifyError(@()WVBottomWaveGenerationForcing(wvt,topographicHeight=invalidTerrain,barotropicVelocityAmplitude=[0.05; 0]),"WVBottomWaveGenerationForcing:InvalidTopographicHeight")
            testCase.verifyError(@()WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05 0]),"WVBottomWaveGenerationForcing:InvalidBarotropicVelocityAmplitude")
            testCase.verifyError(@()WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[Inf; 0]),"WVBottomWaveGenerationForcing:InvalidBarotropicVelocityAmplitude")
            testCase.verifyError(@()WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],frequency=0),"WVBottomWaveGenerationForcing:InvalidFrequency")
            testCase.verifyError(@()WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],rampDuration=-1),"WVBottomWaveGenerationForcing:InvalidRampDuration")
            testCase.verifyError(@()WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],startTime=NaN),"WVBottomWaveGenerationForcing:InvalidStartTime")
            testCase.verifyError(@()WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],name=""),"WVBottomWaveGenerationForcing:InvalidName")

            hydrostatic = WVTransformHydrostatic([wvt.Lx wvt.Ly wvt.Lz],[wvt.Nx wvt.Ny wvt.Nz],N2=@(z)2e-5*ones(size(z)),latitude=45,shouldAntialias=false);
            testCase.verifyError(@()WVBottomWaveGenerationForcing(hydrostatic,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0]),"WVBottomWaveGenerationForcing:UnsupportedTransform")
            variableN2 = WVTransformBoussinesq([wvt.Lx wvt.Ly wvt.Lz],[wvt.Nx wvt.Ny wvt.Nz],N2=@(z)2e-5*exp(z/4000),latitude=45,shouldAntialias=false);
            variableForcing = WVBottomWaveGenerationForcing(variableN2,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0]);
            testCase.verifyClass(variableForcing,"WVBottomWaveGenerationForcing")
        end

        function prescribedVelocityAndRampAreCorrect(testCase)
            wvt = TestWVBottomWaveGenerationForcing.createTransform(true);
            terrain = TestWVBottomWaveGenerationForcing.bandLimitedTopography(wvt);
            amplitude = [0.08+0.02i; -0.03+0.04i];
            frequency = 1.4e-4;
            forcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=amplitude,frequency=frequency,rampDuration=200,startTime=100);

            testCase.verifyEqual(forcing.barotropicVelocityAtTime(99),zeros(2,1))
            testCase.verifyEqual(forcing.barotropicVelocityAtTime(100),zeros(2,1))
            expectedMidRamp = 0.5*real(amplitude*exp(-1i*frequency*100));
            testCase.verifyEqual(forcing.barotropicVelocityAtTime(200),expectedMidRamp,"AbsTol",10*eps)
            expectedAfterRamp = real(amplitude*exp(-1i*frequency*350));
            testCase.verifyEqual(forcing.barotropicVelocityAtTime(450),expectedAfterRamp,"AbsTol",10*eps)

            velocity = forcing.barotropicVelocityAtTime(200);
            expectedBottom = velocity(1)*wvt.diffX(terrain)+velocity(2)*wvt.diffY(terrain);
            testCase.verifyEqual(forcing.bottomVelocityAtTime(200),expectedBottom,"AbsTol",1e-14)

            immediate = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=amplitude,startTime=100);
            testCase.verifyEqual(immediate.barotropicVelocityAtTime(100),real(amplitude))
        end

        function directOracleMatchesFourierSelectionAndKernel(testCase)
            for shouldAntialias = [false true]
                wvt = TestWVBottomWaveGenerationForcing.createTransform(shouldAntialias);
                terrains = {
                    TestWVBottomWaveGenerationForcing.sinusoidalTopography(wvt,50)
                    TestWVBottomWaveGenerationForcing.bandLimitedTopography(wvt)
                    };
                for iTerrain = 1:numel(terrains)
                    amplitude = [0.05+0.01i; -0.02+0.015i];
                    forcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrains{iTerrain},barotropicVelocityAmplitude=amplitude,frequency=1.405e-4);
                    for t = [0 813]
                        wvt.t = t;
                        velocity = forcing.barotropicVelocityAtTime(t);
                        [directFp,directFm,fourierFp,fourierFm,gBottom,gBottomFourier] = referenceBoundaryProjection(wvt,terrains{iTerrain},velocity);
                        [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
                        testCase.verifyLessThanOrEqual(TestWVBottomWaveGenerationForcing.relativeError(directFp,fourierFp),5e-12)
                        testCase.verifyLessThanOrEqual(TestWVBottomWaveGenerationForcing.relativeError(directFm,fourierFm),5e-12)
                        testCase.verifyLessThanOrEqual(TestWVBottomWaveGenerationForcing.relativeError(Fp,directFp),5e-12)
                        testCase.verifyLessThanOrEqual(TestWVBottomWaveGenerationForcing.relativeError(Fm,directFm),5e-12)
                        testCase.verifyEqual(F0,zeros(size(wvt.A0)))

                        reconstructed = wvt.transformToSpatialDomainWithFourier(repmat(gBottomFourier,wvt.Nz,1));
                        testCase.verifyLessThanOrEqual(norm(reconstructed(:,:,1)-gBottom,"fro")/max(norm(gBottom,"fro"),eps),5e-12)
                        testCase.verifyLessThanOrEqual(max(abs(imag(reconstructed)),[],"all"),1e-12*max(1,max(abs(gBottom),[],"all")))
                    end
                end
            end
        end

        function kernelSelectionLimitsAndScaling(testCase)
            wvt = TestWVBottomWaveGenerationForcing.createTransform(true);
            terrain = TestWVBottomWaveGenerationForcing.sinusoidalTopography(wvt,50);
            xForcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0]);
            [Fp,Fm,F0] = xForcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            expectedHorizontal = abs(abs(wvt.K)-2*pi/wvt.Lx) < 100*eps(2*pi/wvt.Lx) & abs(wvt.L) < eps;
            testCase.verifyGreaterThan(norm([Fp(:); Fm(:)]),0)
            testCase.verifyLessThanOrEqual(norm(Fp(~expectedHorizontal)),1e-18)
            testCase.verifyLessThanOrEqual(norm(Fm(~expectedHorizontal)),1e-18)
            testCase.verifyTrue(all(Fp(~wvt.waveComponent.maskAp) == 0))
            testCase.verifyTrue(all(Fm(~wvt.waveComponent.maskAm) == 0))
            testCase.verifyEqual(F0,zeros(size(wvt.A0)))

            transverse = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0; 0.05]);
            [FpTransverse,FmTransverse] = transverse.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyLessThanOrEqual(norm([FpTransverse(:); FmTransverse(:)]),1e-15)

            flat = WVBottomWaveGenerationForcing(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny),barotropicVelocityAmplitude=[0.05; 0]);
            uniform = WVBottomWaveGenerationForcing(wvt,topographicHeight=50*ones(wvt.Nx,wvt.Ny),barotropicVelocityAmplitude=[0.05; 0]);
            zeroCurrent = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=zeros(2,1));
            TestWVBottomWaveGenerationForcing.verifyZeroForcing(testCase,wvt,flat)
            TestWVBottomWaveGenerationForcing.verifyZeroForcing(testCase,wvt,uniform)
            TestWVBottomWaveGenerationForcing.verifyZeroForcing(testCase,wvt,zeroCurrent)

            doubledTerrain = WVBottomWaveGenerationForcing(wvt,topographicHeight=2*terrain,barotropicVelocityAmplitude=[0.05; 0]);
            doubledCurrent = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.10; 0]);
            [FpTerrain,FmTerrain] = doubledTerrain.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            [FpCurrent,FmCurrent] = doubledCurrent.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyEqual(FpTerrain,2*Fp,"RelTol",1e-12,"AbsTol",1e-15)
            testCase.verifyEqual(FmTerrain,2*Fm,"RelTol",1e-12,"AbsTol",1e-15)
            testCase.verifyEqual(FpCurrent,2*Fp,"RelTol",1e-12,"AbsTol",1e-15)
            testCase.verifyEqual(FmCurrent,2*Fm,"RelTol",1e-12,"AbsTol",1e-15)
        end

        function antialiasLayoutsAgreeOnCommonModes(testCase)
            wvtFull = TestWVBottomWaveGenerationForcing.createTransform(false);
            wvtDealiased = TestWVBottomWaveGenerationForcing.createTransform(true);
            terrainFull = TestWVBottomWaveGenerationForcing.bandLimitedTopography(wvtFull);
            terrainDealiased = TestWVBottomWaveGenerationForcing.bandLimitedTopography(wvtDealiased);
            amplitude = [0.05+0.01i; -0.02];
            forcingFull = WVBottomWaveGenerationForcing(wvtFull,topographicHeight=terrainFull,barotropicVelocityAmplitude=amplitude);
            forcingDealiased = WVBottomWaveGenerationForcing(wvtDealiased,topographicHeight=terrainDealiased,barotropicVelocityAmplitude=amplitude);
            [FpFull,FmFull] = forcingFull.addSpectralForcing(wvtFull,zeros(size(wvtFull.Ap)),zeros(size(wvtFull.Am)),zeros(size(wvtFull.A0)));
            [FpDealiased,FmDealiased] = forcingDealiased.addSpectralForcing(wvtDealiased,zeros(size(wvtDealiased.Ap)),zeros(size(wvtDealiased.Am)),zeros(size(wvtDealiased.A0)));

            for index = reshape(find(wvtDealiased.waveComponent.maskAp),1,[])
                match = find(wvtFull.waveComponent.maskAp & wvtFull.K == wvtDealiased.K(index) & wvtFull.L == wvtDealiased.L(index) & wvtFull.J == wvtDealiased.J(index));
                testCase.assertNumElements(match,1)
                testCase.verifyEqual(FpDealiased(index),FpFull(match),"RelTol",1e-12,"AbsTol",1e-15)
                testCase.verifyEqual(FmDealiased(index),FmFull(match),"RelTol",1e-12,"AbsTol",1e-15)
            end
        end

        function callbackAccumulatesWithoutMutatingState(testCase)
            wvt = TestWVBottomWaveGenerationForcing.createTransform(false);
            terrain = TestWVBottomWaveGenerationForcing.bandLimitedTopography(wvt);
            forcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05+0.01i; -0.02i]);
            wvt.t = 647;
            wvt.Ap(2,2) = 0.3+0.2i;
            wvt.Am(2,2) = -0.1+0.4i;
            wvt.A0(2,2) = 0.2-0.3i;
            tBefore = wvt.t;
            ApBefore = wvt.Ap;
            AmBefore = wvt.Am;
            A0Before = wvt.A0;
            [contributionFp,contributionFm] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));

            baseFp = complex(reshape(1:numel(wvt.Ap),size(wvt.Ap)),reshape(numel(wvt.Ap):-1:1,size(wvt.Ap)));
            baseFm = -2*baseFp;
            baseF0 = 3*baseFp;
            [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,baseFp,baseFm,baseF0);
            testCase.verifyEqual(Fp,baseFp+contributionFp)
            testCase.verifyEqual(Fm,baseFm+contributionFm)
            testCase.verifyEqual(F0,baseF0)
            testCase.verifyEqual(wvt.t,tBefore)
            testCase.verifyEqual(wvt.Ap,ApBefore)
            testCase.verifyEqual(wvt.Am,AmBefore)
            testCase.verifyEqual(wvt.A0,A0Before)

            otherTransform = TestWVBottomWaveGenerationForcing.createTransform(false);
            testCase.verifyError(@()forcing.addSpectralForcing(otherTransform,baseFp,baseFm,baseF0),"WVBottomWaveGenerationForcing:TransformMismatch")
            convertedForcing = forcing.forcingWithResolutionOfTransform(otherTransform);
            testCase.verifyClass(convertedForcing,"WVBottomWaveGenerationForcing")
        end

        function noPVAndEnergyWorkIdentities(testCase)
            previousRandomState = rng;
            randomStateCleanup = onCleanup(@()rng(previousRandomState));
            rng(67291,"twister")
            for shouldAntialias = [false true]
                wvt = TestWVBottomWaveGenerationForcing.createTransform(shouldAntialias);
                terrain = TestWVBottomWaveGenerationForcing.sinusoidalTopography(wvt,50);
                amplitude = [0.05+0.01i; -0.02+0.015i];
                forcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=amplitude,frequency=1.405e-4);
                wvt.t = 813;
                [forcingFp,forcingFm] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
                randomAp = (randn(size(wvt.Ap))+1i*randn(size(wvt.Ap))).*wvt.waveComponent.maskAp;
                randomAm = (randn(size(wvt.Am))+1i*randn(size(wvt.Am))).*wvt.waveComponent.maskAm;
                states = {
                    1e4*forcingFp,zeros(size(wvt.Am))
                    zeros(size(wvt.Ap)),1e4*forcingFm
                    1e-3*randomAp+1e4*forcingFp,1e-3*randomAm+1e4*forcingFm
                    };
                for iState = 1:size(states,1)
                    wvt.Ap = states{iState,1};
                    wvt.Am = states{iState,2};
                    wvt.A0(:) = 0;
                    diagnostics = TestWVBottomWaveGenerationForcing.sourceDiagnostics(wvt,forcing);
                    testCase.verifyEqual(diagnostics.F0,zeros(size(wvt.A0)))
                    testCase.verifyEqual(diagnostics.modalQGPVSource,zeros(size(wvt.A0)))
                    testCase.verifyLessThanOrEqual(diagnostics.normalizedQGPV,1e-10)
                    testCase.verifyLessThanOrEqual(diagnostics.powerError,5e-12)
                end

                wvt.Ap = states{3,1};
                wvt.Am = states{3,2};
                baseDiagnostics = TestWVBottomWaveGenerationForcing.sourceDiagnostics(wvt,forcing);
                doubledTerrain = WVBottomWaveGenerationForcing(wvt,topographicHeight=2*terrain,barotropicVelocityAmplitude=amplitude,frequency=forcing.frequency,startTime=forcing.startTime);
                doubledCurrent = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=2*amplitude,frequency=forcing.frequency,startTime=forcing.startTime);
                terrainDiagnostics = TestWVBottomWaveGenerationForcing.sourceDiagnostics(wvt,doubledTerrain);
                currentDiagnostics = TestWVBottomWaveGenerationForcing.sourceDiagnostics(wvt,doubledCurrent);
                testCase.verifyEqual(terrainDiagnostics.modalPower,2*baseDiagnostics.modalPower,RelTol=1e-12,AbsTol=1e-15)
                testCase.verifyEqual(terrainDiagnostics.bottomPower,2*baseDiagnostics.bottomPower,RelTol=1e-12,AbsTol=1e-15)
                testCase.verifyEqual(currentDiagnostics.modalPower,2*baseDiagnostics.modalPower,RelTol=1e-12,AbsTol=1e-15)
                testCase.verifyEqual(currentDiagnostics.bottomPower,2*baseDiagnostics.bottomPower,RelTol=1e-12,AbsTol=1e-15)

                flat = WVBottomWaveGenerationForcing(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny),barotropicVelocityAmplitude=amplitude,frequency=forcing.frequency,startTime=forcing.startTime);
                zeroCurrent = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=zeros(2,1),frequency=forcing.frequency,startTime=forcing.startTime);
                TestWVBottomWaveGenerationForcing.verifyZeroDiagnostics(testCase,wvt,flat)
                TestWVBottomWaveGenerationForcing.verifyZeroDiagnostics(testCase,wvt,zeroCurrent)
            end
            clear randomStateCleanup

            wvt = TestWVBottomWaveGenerationForcing.createTransform(true);
            terrain = TestWVBottomWaveGenerationForcing.sinusoidalTopography(wvt,50);
            forcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0]);
            wvt.removeAllForcing();
            wvt.addForcing(forcing);
            warningState = warning;
            warningCleanup = onCleanup(@()warning(warningState));
            warning("off","all")
            model = WVModel(wvt);
            model.setupIntegrator(integratorType="adaptive",absTolerance=1e-10,relTolerance=1e-8);
            model.integrateToTime(600,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            evolvedDiagnostics = TestWVBottomWaveGenerationForcing.sourceDiagnostics(wvt,forcing);
            testCase.verifyEqual(evolvedDiagnostics.F0,zeros(size(wvt.A0)))
            testCase.verifyEqual(evolvedDiagnostics.modalQGPVSource,zeros(size(wvt.A0)))
            testCase.verifyLessThanOrEqual(evolvedDiagnostics.normalizedQGPV,1e-10)
            testCase.verifyLessThanOrEqual(evolvedDiagnostics.powerError,5e-12)
            clear warningCleanup
        end

        function exactSinusoidalRidgeBenchmarkAndFigures(testCase)
            originalVisibility = get(groot,"defaultFigureVisible");
            visibilityCleanup = onCleanup(@()set(groot,"defaultFigureVisible",originalVisibility));
            set(groot,"defaultFigureVisible","off")
            benchmark = SinusoidalRidgeWaveGenerationBenchmark(resolution=[8 4 5],shouldAntialias=true,shouldMakeFigures=true);
            figureCleanup = onCleanup(@()close(benchmark.figureHandles(isgraphics(benchmark.figureHandles))));

            testCase.verifyTrue(all(diff(benchmark.coefficientError) < 0))
            testCase.verifyLessThan(benchmark.coefficientError(end),1e-8)
            testCase.verifyLessThanOrEqual(benchmark.referenceEnergyWorkError,1e-12)
            testCase.verifyTrue(all(diff(benchmark.energyWorkError) < 0))
            testCase.verifyLessThan(benchmark.energyWorkError(end),1e-8)
            testCase.verifyLessThanOrEqual(max(benchmark.balancedEnergyFraction),1e-20)
            testCase.verifyLessThanOrEqual(benchmark.normalizedQGPV,1e-10)
            testCase.verifyLessThanOrEqual(benchmark.sourceSupportError,1e-12)
            testCase.verifyEqual(sum(benchmark.plusEnergyByMode)+sum(benchmark.minusEnergyByMode),benchmark.referenceEnergy(end),RelTol=1e-12)
            testCase.verifyNumElements(benchmark.figureHandles,2)
            testCase.verifyTrue(all(isgraphics(benchmark.figureHandles,"figure")))
            for figureHandle = benchmark.figureHandles
                testCase.verifyGreaterThanOrEqual(numel(findall(figureHandle,Type="axes")),4)
                lineHandles = findall(figureHandle,Type="line");
                for lineHandle = reshape(lineHandles,1,[])
                    testCase.verifyTrue(all(isfinite([lineHandle.XData(:); lineHandle.YData(:)])))
                end
                imageHandles = findall(figureHandle,Type="image");
                for imageHandle = reshape(imageHandles,1,[])
                    testCase.verifyTrue(all(isfinite(imageHandle.CData),"all"))
                end
            end
            clear figureCleanup visibilityCleanup
        end

        function adaptiveModelSmoke(testCase)
            wvt = TestWVBottomWaveGenerationForcing.createTransform(true);
            terrain = TestWVBottomWaveGenerationForcing.sinusoidalTopography(wvt,50);
            forcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0]);
            wvt.removeAllForcing();
            wvt.addForcing(forcing);
            [Fp,Fm,F0] = wvt.nonlinearFlux();
            testCase.verifyGreaterThan(norm([Fp(:); Fm(:)]),0)
            testCase.verifyEqual(F0,zeros(size(wvt.A0)))

            warningState = warning;
            warningCleanup = onCleanup(@()warning(warningState));
            warning("off","all")
            model = WVModel(wvt);
            model.setupIntegrator(integratorType="adaptive",absTolerance=1e-10,relTolerance=1e-8);
            model.integrateToTime(300,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            testCase.verifyEqual(wvt.t,300,"AbsTol",1e-10)
            testCase.verifyTrue(all(isfinite([wvt.Ap(:); wvt.Am(:); wvt.A0(:)])))
            testCase.verifyGreaterThan(norm([wvt.Ap(:); wvt.Am(:)]),0)
            testCase.verifyEqual(wvt.A0,zeros(size(wvt.A0)))
            clear warningCleanup
        end
    end

    methods (Static, Access = private)
        function wvt = createTransform(shouldAntialias)
            wvt = WVTransformBoussinesq([4e3 4e3 2e3],[8 4 5],N2=@(z)2e-5*ones(size(z)),latitude=45,shouldAntialias=shouldAntialias);
        end

        function terrain = sinusoidalTopography(wvt,amplitude)
            x = reshape(wvt.x,[],1);
            terrain = amplitude*cos(2*pi*x/wvt.Lx).*ones(1,wvt.Ny);
        end

        function terrain = bandLimitedTopography(wvt)
            x = reshape(wvt.x,[],1);
            y = reshape(wvt.y,1,[]);
            previousRandomState = rng;
            randomStateCleanup = onCleanup(@()rng(previousRandomState));
            rng(48271,"twister")
            modes = [1 0; 0 1; 1 1; 2 1];
            amplitudes = 5+20*rand(size(modes,1),1);
            phases = 2*pi*rand(size(modes,1),1);
            terrain = zeros(wvt.Nx,wvt.Ny);
            for iMode = 1:size(modes,1)
                terrain = terrain+amplitudes(iMode)*cos(2*pi*(modes(iMode,1)*x/wvt.Lx+modes(iMode,2)*y/wvt.Ly)+phases(iMode));
            end
            clear randomStateCleanup
        end

        function errorValue = relativeError(actual,expected)
            errorValue = norm(actual(:)-expected(:))/max(norm(expected(:)),eps);
        end

        function verifyZeroForcing(testCase,wvt,forcing)
            [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyEqual(Fp,zeros(size(wvt.Ap)))
            testCase.verifyEqual(Fm,zeros(size(wvt.Am)))
            testCase.verifyEqual(F0,zeros(size(wvt.A0)))
        end

        function diagnostics = sourceDiagnostics(wvt,forcing)
            [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            [Fu,Fv,~,Feta] = wvt.transformWaveVortexToUVWEta(Fp,Fm,F0,wvt.t);
            vorticityX = wvt.diffX(Fv);
            vorticityY = wvt.diffY(Fu);
            stretching = wvt.f*wvt.diffZG(Feta);
            qgpvSource = vorticityX-vorticityY-stretching;
            qgpvScale = norm(vorticityX(:))+norm(vorticityY(:))+norm(stretching(:));
            normalizedQGPV = norm(qgpvSource(:))/max(qgpvScale,eps);
            modalPower = 2*sum(wvt.Apm_TE_factor(:).*real(Fp(:).*conj(wvt.Ap(:))+Fm(:).*conj(wvt.Am(:))));
            [~,iBottom] = min(wvt.z);
            bottomPower = mean((wvt.p(:,:,iBottom)/wvt.rho0).*forcing.bottomVelocityAtTime(wvt.t),"all");
            powerError = abs(modalPower-bottomPower)/max([abs(modalPower) abs(bottomPower) eps]);
            modalQGPVSource = wvt.A0_QGPV_factor.*F0;
            diagnostics = struct(Fp=Fp,Fm=Fm,F0=F0,modalQGPVSource=modalQGPVSource,qgpvSource=qgpvSource,normalizedQGPV=normalizedQGPV,modalPower=modalPower,bottomPower=bottomPower,powerError=powerError);
        end

        function verifyZeroDiagnostics(testCase,wvt,forcing)
            diagnostics = TestWVBottomWaveGenerationForcing.sourceDiagnostics(wvt,forcing);
            testCase.verifyEqual(diagnostics.Fp,zeros(size(wvt.Ap)))
            testCase.verifyEqual(diagnostics.Fm,zeros(size(wvt.Am)))
            testCase.verifyEqual(diagnostics.F0,zeros(size(wvt.A0)))
            testCase.verifyEqual(diagnostics.modalQGPVSource,zeros(size(wvt.A0)))
            testCase.verifyEqual(diagnostics.qgpvSource,zeros(size(diagnostics.qgpvSource)))
            testCase.verifyEqual(diagnostics.modalPower,0)
            testCase.verifyEqual(diagnostics.bottomPower,0)
        end
    end
end
