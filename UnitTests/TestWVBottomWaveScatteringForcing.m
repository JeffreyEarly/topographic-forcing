classdef TestWVBottomWaveScatteringForcing < matlab.unittest.TestCase
    % Verify autonomous first-order wave scattering and diagnostics.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(repositoryRoot)
            addpath(fullfile(repositoryRoot,"Examples"))
        end
    end

    methods (Test)
        function constructorAndEndpointReconstruction(testCase)
            for shouldAntialias = [false true]
                wvt = TestWVBottomWaveScatteringForcing.createTransform([8 6 5],shouldAntialias,false);
                terrain = TestWVBottomWaveScatteringForcing.topography(wvt,30);
                forcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=terrain,name="scattering");
                testCase.verifyEqual(forcing.topographicHeight,terrain)
                testCase.verifyEqual(string(forcing.name),"scattering")

                TestWVBottomWaveScatteringForcing.setRandomWaveState(wvt,1937)
                wvt.A0 = complex(ones(size(wvt.A0)));
                wvt.t = 317;
                [gBottom,bottomFields] = forcing.bottomVelocityFromWaveState(wvt);
                [uReference,vReference,dWdzReference,gReference] = TestWVBottomWaveScatteringForcing.fullEndpointReference(wvt,terrain);
                testCase.verifyLessThanOrEqual(TestWVBottomWaveScatteringForcing.relativeError(bottomFields.u,uReference),1e-12)
                testCase.verifyLessThanOrEqual(TestWVBottomWaveScatteringForcing.relativeError(bottomFields.v,vReference),1e-12)
                testCase.verifyLessThanOrEqual(TestWVBottomWaveScatteringForcing.relativeError(bottomFields.dWdz,dWdzReference),1e-12)
                testCase.verifyLessThanOrEqual(TestWVBottomWaveScatteringForcing.relativeError(gBottom,gReference),1e-12)
            end

            wvt = TestWVBottomWaveScatteringForcing.createTransform([8 6 5],false,false);
            terrain = zeros(wvt.Nx,wvt.Ny);
            testCase.verifyError(@()WVBottomWaveScatteringForcing(wvt,topographicHeight=zeros(wvt.Nx-1,wvt.Ny)),"WVBottomWaveScatteringForcing:InvalidTopographicHeightSize")
            testCase.verifyError(@()WVBottomWaveScatteringForcing(wvt,topographicHeight=1i*ones(wvt.Nx,wvt.Ny)),"WVBottomWaveScatteringForcing:InvalidTopographicHeight")
            terrain(1) = NaN;
            testCase.verifyError(@()WVBottomWaveScatteringForcing(wvt,topographicHeight=terrain),"WVBottomWaveScatteringForcing:InvalidTopographicHeight")
            terrain(1) = 0;
            testCase.verifyError(@()WVBottomWaveScatteringForcing(wvt,topographicHeight=terrain,name=""),"WVBottomWaveScatteringForcing:InvalidName")
            hydrostatic = WVTransformHydrostatic([wvt.Lx wvt.Ly wvt.Lz],[wvt.Nx wvt.Ny wvt.Nz],N2=@(z)2e-5*ones(size(z)),latitude=45,shouldAntialias=false);
            testCase.verifyError(@()WVBottomWaveScatteringForcing(hydrostatic,topographicHeight=terrain),"WVBottomWaveScatteringForcing:UnsupportedTransform")
        end

        function directProjectionAndNoPV(testCase)
            for variableStratification = [false true]
                wvt = TestWVBottomWaveScatteringForcing.createTransform([8 6 5],false,variableStratification);
                terrain = TestWVBottomWaveScatteringForcing.topography(wvt,30);
                forcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=terrain);
                TestWVBottomWaveScatteringForcing.setRandomWaveState(wvt,913+variableStratification)
                wvt.t = 619;
                [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
                [gBottom,~] = forcing.bottomVelocityFromWaveState(wvt);
                [directFp,directFm] = TestWVBottomWaveScatteringForcing.directProjection(wvt,gBottom);
                testCase.verifyLessThanOrEqual(TestWVBottomWaveScatteringForcing.relativeError(Fp,directFp),1e-10)
                testCase.verifyLessThanOrEqual(TestWVBottomWaveScatteringForcing.relativeError(Fm,directFm),1e-10)
                testCase.verifyEqual(F0,zeros(size(wvt.A0)))

                [Fu,Fv,~,Feta] = wvt.transformWaveVortexToUVWEta(Fp,Fm,F0,wvt.t);
                qgpvTendency = wvt.diffX(Fv)-wvt.diffY(Fu)-wvt.f*wvt.diffZG(Feta);
                testCase.verifyLessThanOrEqual(norm(qgpvTendency(:)),1e-12)
            end
        end

        function sinusoidalSidebandsAndScaling(testCase)
            wvt = TestWVBottomWaveScatteringForcing.createTransform([12 4 5],false,false);
            [x,~] = ndgrid(wvt.x,wvt.y);
            terrain = 30*cos(2*pi*x/wvt.Lx);
            forcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=terrain);
            wvt.Ap(:) = 0;
            wvt.Am(:) = 0;
            wvt.A0(:) = 0;
            incident = TestWVBottomWaveScatteringForcing.modeIndex(wvt,2,0,1,"plus");
            wvt.Ap(incident) = 0.01+0.004i;
            wvt.t = 0;
            [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));

            expectedSupport = (abs(wvt.K-2*pi/wvt.Lx) < 100*eps(2*pi/wvt.Lx) | abs(wvt.K-6*pi/wvt.Lx) < 100*eps(6*pi/wvt.Lx)) & abs(wvt.L) < eps;
            totalNorm = norm([Fp(:); Fm(:)]);
            testCase.verifyGreaterThan(totalNorm,0)
            testCase.verifyLessThanOrEqual(norm([Fp(~expectedSupport); Fm(~expectedSupport)])/totalNorm,1e-12)
            testCase.verifyEqual(F0,zeros(size(wvt.A0)))

            doubled = WVBottomWaveScatteringForcing(wvt,topographicHeight=2*terrain);
            [Fp2,Fm2] = doubled.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyEqual(Fp2,2*Fp,RelTol=1e-12,AbsTol=1e-15)
            testCase.verifyEqual(Fm2,2*Fm,RelTol=1e-12,AbsTol=1e-15)

            wvt.Ap = Fp;
            wvt.Am = Fm;
            [secondFp,secondFm] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            wvt.Ap = Fp2;
            wvt.Am = Fm2;
            [secondFp2,secondFm2] = doubled.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyEqual(secondFp2,4*secondFp,RelTol=1e-12,AbsTol=1e-15)
            testCase.verifyEqual(secondFm2,4*secondFm,RelTol=1e-12,AbsTol=1e-15)
        end

        function energyIdentitiesAndQuadraticResidual(testCase)
            wvt = TestWVBottomWaveScatteringForcing.createTransform([12 4 9],false,false);
            wvt.Ap(:) = 0;
            wvt.Am(:) = 0;
            wvt.A0(:) = 0;
            plusIndices = find(wvt.J <= 2 & abs(wvt.K) <= 6*pi/wvt.Lx & abs(wvt.L) < eps & wvt.waveComponent.maskAp);
            minusIndices = find(wvt.J <= 2 & abs(wvt.K) <= 6*pi/wvt.Lx & abs(wvt.L) < eps & wvt.waveComponent.maskAm);
            wvt.Ap(plusIndices) = 0.01*reshape(1:numel(plusIndices),[],1)+0.003i;
            wvt.Am(minusIndices) = 0.007*reshape(1:numel(minusIndices),[],1)-0.002i;
            wvt.t = 431;
            [x,~] = ndgrid(wvt.x,wvt.y);
            baseTerrain = cos(2*pi*x/wvt.Lx);
            scales = [5 10 20 40];
            residual = zeros(size(scales));
            strictError = zeros(size(scales));
            for iScale = 1:numel(scales)
                terrain = scales(iScale)*baseTerrain;
                forcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=terrain);
                [Fp,Fm] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
                flatEnergyTendency = 2*sum(wvt.Apm_TE_factor(:).*real(Fp(:).*conj(wvt.Ap(:))+Fm(:).*conj(wvt.Am(:))));
                phase = exp(wvt.iOmega*(wvt.t-wvt.t0));
                Apt = wvt.Ap.*phase;
                Amt = wvt.Am.*conj(phase);
                dAptFree = wvt.iOmega.*Apt;
                dAmtFree = -wvt.iOmega.*Amt;
                dAptTerrain = Fp.*phase;
                dAmtTerrain = Fm.*conj(phase);
                [~,bottomFields] = forcing.bottomVelocityFromWaveState(wvt);
                [uFree,vFree] = TestWVBottomWaveScatteringForcing.horizontalEndpointTendency(wvt,dAptFree,dAmtFree);
                [uTerrain,vTerrain] = TestWVBottomWaveScatteringForcing.horizontalEndpointTendency(wvt,dAptTerrain,dAmtTerrain);
                bottomFreeTendency = mean(terrain.*(bottomFields.u.*uFree+bottomFields.v.*vFree),"all");
                bottomTerrainTendency = mean(terrain.*(bottomFields.u.*uTerrain+bottomFields.v.*vTerrain),"all");
                strictError(iScale) = abs(flatEnergyTendency-bottomFreeTendency)/max([abs(flatEnergyTendency) abs(bottomFreeTendency) eps]);
                residual(iScale) = abs(flatEnergyTendency-bottomFreeTendency-bottomTerrainTendency);
            end
            testCase.verifyLessThanOrEqual(max(strictError),1e-10)
            exponent = polyfit(log(scales),log(residual),1);
            testCase.verifyGreaterThanOrEqual(exponent(1),1.8)
            testCase.verifyLessThanOrEqual(exponent(1),2.2)
        end

        function persistenceResolutionAndNetCDFDiagnostic(testCase)
            sourceTransform = TestWVBottomWaveScatteringForcing.createTransform([8 6 5],false,true);
            terrain = TestWVBottomWaveScatteringForcing.topography(sourceTransform,30);
            forcing = WVBottomWaveScatteringForcing(sourceTransform,topographicHeight=terrain,name="diagnostic scattering");
            targetTransform = TestWVBottomWaveScatteringForcing.createTransform([12 10 7],false,true);
            converted = forcing.forcingWithResolutionOfTransform(targetTransform);
            expectedTerrain = TestWVBottomWaveScatteringForcing.topography(targetTransform,30);
            testCase.verifyEqual(converted.topographicHeight,expectedTerrain,AbsTol=1e-12)

            path = string(tempname)+".nc";
            cleanup = onCleanup(@()TestWVBottomWaveScatteringForcing.deleteFile(path));
            wvt = TestWVBottomWaveScatteringForcing.createTransform([8 6 5],false,false);
            [x,~] = ndgrid(wvt.x,wvt.y);
            terrain = 30*cos(2*pi*x/wvt.Lx);
            forcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=terrain,name="file scattering");
            wvt.removeAllForcing();
            wvt.addForcing(forcing);
            wvt.Ap(:) = 0;
            wvt.Am(:) = 0;
            wvt.A0(:) = 0;
            incident = TestWVBottomWaveScatteringForcing.modeIndex(wvt,2,0,1,"plus");
            wvt.Ap(incident) = 0.01;
            omega = wvt.Omega(incident);
            numberOfTimes = 97;
            time = linspace(0,2*pi/omega,numberOfTimes);
            model = WVModel(wvt);
            outputFile = model.createNetCDFFileForModelOutput(path,outputInterval=time(2)-time(1),shouldOverwriteExisting=true);
            time = outputFile.outputTimesForIntegrationPeriod(time(1),time(end));
            numberOfTimes = numel(time);
            for iTime = 1:numberOfTimes
                wvt.t = time(iTime);
                outputFile.writeTimeStepToOutputFile(time(iTime));
            end
            outputFile.closeNetCDFFile();

            wvt.t = 0;
            gCos = forcing.bottomVelocityFromWaveState(wvt);
            wvt.t = pi/(2*omega);
            gSin = forcing.bottomVelocityFromWaveState(wvt);
            strides = [4 2 1];
            quadratureError = zeros(size(strides));
            for iStride = 1:numel(strides)
                selectedIndices = reshape(1:strides(iStride):numberOfTimes,[],1);
                diagnostics = WVBottomWaveScatteringForcing.bottomDisplacementFromFile(path,forcingName="file scattering",iTime=selectedIndices);
                expectedDisplacement = zeros(size(diagnostics.bottomDisplacement));
                for iTime = 1:numel(diagnostics.time)
                    expectedDisplacement(:,:,iTime) = gCos*sin(omega*diagnostics.time(iTime))/omega+gSin*(1-cos(omega*diagnostics.time(iTime)))/omega;
                end
                quadratureError(iStride) = TestWVBottomWaveScatteringForcing.relativeError(diagnostics.bottomDisplacement,expectedDisplacement);
            end
            testCase.verifyEqual(diagnostics.time,time(:),AbsTol=10*eps(time(end)))
            testCase.verifyEqual(diagnostics.forcingName,"file scattering")
            testCase.verifyEqual(diagnostics.quadratureMethod,"trapezoidal")
            quadratureExponent = polyfit(log(strides),log(quadratureError),1);
            testCase.verifyGreaterThanOrEqual(quadratureExponent(1),1.8)
            testCase.verifyLessThanOrEqual(quadratureExponent(1),2.2)
            testCase.verifyLessThanOrEqual(quadratureError(end),5e-4)
            testCase.verifyLessThanOrEqual(max(abs(imag(diagnostics.bottomDisplacement)),[],"all"),1e-13)
            initialBottomDisplacement = 0.2*cos(2*pi*x/wvt.Lx);
            initializedDiagnostics = WVBottomWaveScatteringForcing.bottomDisplacementFromFile(path,forcingName="file scattering",iTime=1:4:numberOfTimes,initialBottomDisplacement=initialBottomDisplacement);
            testCase.verifyEqual(initializedDiagnostics.bottomDisplacement(:,:,1),initialBottomDisplacement,AbsTol=1e-13)

            [restoredTransform,restoredFile] = WVTransformBoussinesq.waveVortexTransformFromFile(char(path),iTime=17,shouldReadOnly=true);
            restoredCleanup = onCleanup(@()restoredFile.close());
            restoredForcing = restoredTransform.forcingWithName("file scattering");
            testCase.verifyClass(restoredForcing,"WVBottomWaveScatteringForcing")
            testCase.verifyEqual(restoredForcing.topographicHeight,terrain)
            [restoredFp,restoredFm,restoredF0] = restoredForcing.addSpectralForcing(restoredTransform,zeros(size(restoredTransform.Ap)),zeros(size(restoredTransform.Am)),zeros(size(restoredTransform.A0)));
            wvt.Ap = restoredTransform.Ap;
            wvt.Am = restoredTransform.Am;
            wvt.A0 = restoredTransform.A0;
            wvt.t = restoredTransform.t;
            [originalFp,originalFm,originalF0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyLessThanOrEqual(TestWVBottomWaveScatteringForcing.relativeError(restoredFp,originalFp),1e-12)
            testCase.verifyLessThanOrEqual(TestWVBottomWaveScatteringForcing.relativeError(restoredFm,originalFm),1e-12)
            testCase.verifyEqual(restoredF0,originalF0)
            clear restoredCleanup cleanup
        end

        function explicitAntialiasAndRestartContinuation(testCase)
            implicitTransform = TestWVBottomWaveScatteringForcing.createTransform([8 6 5],true,true);
            implicitForcing = WVBottomWaveScatteringForcing(implicitTransform,topographicHeight=TestWVBottomWaveScatteringForcing.topography(implicitTransform,30),name="antialias scattering");
            implicitTransform.removeAllForcing();
            implicitTransform.addForcing(implicitForcing);
            explicitTransform = implicitTransform.waveVortexTransformWithExplicitAntialiasing();
            explicitForcing = explicitTransform.forcingWithName("antialias scattering");
            testCase.verifyClass(explicitForcing,"WVBottomWaveScatteringForcing")
            TestWVBottomWaveScatteringForcing.setRandomWaveState(explicitTransform,419)
            explicitTransform.t = 213;
            [Fp,Fm,F0] = explicitForcing.addSpectralForcing(explicitTransform,zeros(size(explicitTransform.Ap)),zeros(size(explicitTransform.Am)),zeros(size(explicitTransform.A0)));
            [gBottom,~] = explicitForcing.bottomVelocityFromWaveState(explicitTransform);
            [directFp,directFm] = TestWVBottomWaveScatteringForcing.directProjection(explicitTransform,gBottom);
            testCase.verifyLessThanOrEqual(TestWVBottomWaveScatteringForcing.relativeError(Fp,directFp),1e-10)
            testCase.verifyLessThanOrEqual(TestWVBottomWaveScatteringForcing.relativeError(Fm,directFm),1e-10)
            testCase.verifyEqual(F0,zeros(size(explicitTransform.A0)))

            restartPath = string(tempname)+".nc";
            restartCleanup = onCleanup(@()TestWVBottomWaveScatteringForcing.deleteFile(restartPath));
            [restartModel,checkpointTime,finalTime] = TestWVBottomWaveScatteringForcing.modelForRestart();
            [controlModel,~,~] = TestWVBottomWaveScatteringForcing.modelForRestart();
            restartModel.createNetCDFFileForModelOutput(restartPath,outputInterval=checkpointTime/2,shouldOverwriteExisting=true);
            restartModel.integrateToTime(checkpointTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            restartModel.closeNetCDFFile();
            controlModel.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            resumedModel = WVModel.modelFromFile(restartPath);
            resumedModel.setupIntegrator(integratorType="adaptive",absTolerance=1e-12,relTolerance=1e-10);
            resumedModel.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            resumedModel.closeNetCDFFile();
            testCase.verifyLessThanOrEqual(max(abs(resumedModel.wvt.Ap-controlModel.wvt.Ap),[],"all"),1e-10)
            testCase.verifyLessThanOrEqual(max(abs(resumedModel.wvt.Am-controlModel.wvt.Am),[],"all"),1e-10)
            testCase.verifyEqual(resumedModel.wvt.A0,controlModel.wvt.A0)
            testCase.verifyClass(resumedModel.wvt.forcingWithName("restart scattering"),"WVBottomWaveScatteringForcing")
            clear restartCleanup
        end

        function scientificExamples(testCase)
            originalVisibility = get(groot,"defaultFigureVisible");
            visibilityCleanup = onCleanup(@()set(groot,"defaultFigureVisible",originalVisibility));
            set(groot,"defaultFigureVisible","off")
            uniform = UniformDepthWaveScatteringBenchmark(resolution=[8 4 9],shouldMakeFigures=true);
            uniformFigureCleanup = onCleanup(@()close(uniform.figureHandles(isgraphics(uniform.figureHandles))));
            testCase.verifyGreaterThanOrEqual(uniform.convergenceExponent,1.8)
            testCase.verifyLessThanOrEqual(uniform.convergenceExponent,2.2)
            testCase.verifyLessThanOrEqual(uniform.coefficientError,1e-8)
            testCase.verifyEqual(uniform.balancedEnergy,0)
            testCase.verifyLessThanOrEqual(uniform.qgpvTendencyNorm,1e-12)
            testCase.verifyTrue(all(isgraphics(uniform.figureHandles,"figure")))

            ridge = SinusoidalRidgeWaveScatteringExample(resolution=[12 4 7],numberOfWavePeriods=1,numberOfOutputTimes=33,relativeTolerance=1e-8,shouldAntialias=false,shouldMakeFigures=true);
            ridgeFigureCleanup = onCleanup(@()close(ridge.figureHandles(isgraphics(ridge.figureHandles))));
            testCase.verifyTrue(all(isfinite(ridge.waveEnergy)))
            testCase.verifyGreaterThan(min(max(ridge.sidebandEnergy,[],1)),0)
            testCase.verifyEqual(ridge.balancedEnergy,zeros(size(ridge.balancedEnergy)))
            testCase.verifyLessThanOrEqual(max(ridge.qgpvNorm),1e-12)
            testCase.verifyGreaterThan(max(abs(ridge.bottomDisplacement),[],"all"),0)
            testCase.verifyEqual(ridge.bottomQuadratureMethod,"trapezoidal")
            testCase.verifyNumElements(ridge.figureHandles,2)
            testCase.verifyTrue(all(isgraphics(ridge.figureHandles,"figure")))
            clear ridgeFigureCleanup uniformFigureCleanup visibilityCleanup
        end
    end

    methods (Static, Access = private)
        function wvt = createTransform(resolution,shouldAntialias,variableStratification)
            if variableStratification
                N2 = @(z)2e-5*exp(z/4000);
            else
                N2 = @(z)2e-5*ones(size(z));
            end
            wvt = WVTransformBoussinesq([20e3 20e3 2e3],resolution,N2=N2,latitude=45,shouldAntialias=shouldAntialias);
            wvt.t0 = 0;
            wvt.t = 0;
        end

        function terrain = topography(wvt,amplitude)
            [x,y] = ndgrid(wvt.x,wvt.y);
            terrain = amplitude*(cos(2*pi*x/wvt.Lx)+0.4*sin(2*pi*y/wvt.Ly)+0.2*cos(2*pi*(2*x/wvt.Lx+y/wvt.Ly)));
        end

        function setRandomWaveState(wvt,seed)
            previousState = rng;
            cleanup = onCleanup(@()rng(previousState));
            rng(seed,"twister")
            wvt.Ap = 1e-3*(randn(size(wvt.Ap))+1i*randn(size(wvt.Ap))).*wvt.waveComponent.maskAp;
            wvt.Am = 1e-3*(randn(size(wvt.Am))+1i*randn(size(wvt.Am))).*wvt.waveComponent.maskAm;
            wvt.A0(:) = 0;
            clear cleanup
        end

        function [u,v,dWdz,gBottom] = fullEndpointReference(wvt,terrain)
            waveU = wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*wvt.Apt+wvt.UAm.*wvt.Amt);
            waveV = wvt.transformToSpatialDomainWithF(Apm=wvt.VAp.*wvt.Apt+wvt.VAm.*wvt.Amt);
            waveW = wvt.transformToSpatialDomainWithG(Apm=wvt.WAp.*wvt.Apt+wvt.WAm.*wvt.Amt);
            waveWz = wvt.diffZG(waveW);
            [~,iBottom] = min(wvt.z);
            u = waveU(:,:,iBottom);
            v = waveV(:,:,iBottom);
            dWdz = waveWz(:,:,iBottom);
            gBottom = u.*wvt.diffX(terrain)+v.*wvt.diffY(terrain)-terrain.*dWdz;
        end

        function [directFp,directFm] = directProjection(wvt,gBottom)
            directFp = complex(zeros(size(wvt.Ap)));
            directFm = complex(zeros(size(wvt.Am)));
            [~,iBottom] = min(wvt.z);
            x = reshape(wvt.x,[],1);
            y = reshape(wvt.y,1,[]);
            for index = reshape(find(wvt.waveComponent.maskAp),1,[])
                [~,iHorizontal] = ind2sub(size(wvt.Ap),index);
                coefficient = complex(zeros(size(wvt.Ap)));
                coefficient(index) = wvt.NAp(index);
                pressureFourier = wvt.g*wvt.transformToSpatialDomainWithFw(coefficient);
                pressureBottom = pressureFourier(iBottom,iHorizontal)*wvt.phase(index);
                pressurePlane = pressureBottom*exp(1i*(wvt.K(index)*x+wvt.L(index)*y));
                directFp(index) = mean(conj(pressurePlane).*gBottom,"all")/wvt.Apm_TE_factor(index);
            end
            for index = reshape(find(wvt.waveComponent.maskAm),1,[])
                [~,iHorizontal] = ind2sub(size(wvt.Am),index);
                coefficient = complex(zeros(size(wvt.Am)));
                coefficient(index) = wvt.NAm(index);
                pressureFourier = wvt.g*wvt.transformToSpatialDomainWithFw(coefficient);
                pressureBottom = pressureFourier(iBottom,iHorizontal)*wvt.conjPhase(index);
                pressurePlane = pressureBottom*exp(1i*(wvt.K(index)*x+wvt.L(index)*y));
                directFm(index) = mean(conj(pressurePlane).*gBottom,"all")/wvt.Apm_TE_factor(index);
            end
        end

        function [u,v] = horizontalEndpointTendency(wvt,dApt,dAmt)
            uFull = wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*dApt+wvt.UAm.*dAmt);
            vFull = wvt.transformToSpatialDomainWithF(Apm=wvt.VAp.*dApt+wvt.VAm.*dAmt);
            [~,iBottom] = min(wvt.z);
            u = uFull(:,:,iBottom);
            v = vFull(:,:,iBottom);
        end

        function [model,checkpointTime,finalTime] = modelForRestart()
            wvt = TestWVBottomWaveScatteringForcing.createTransform([8 6 5],false,false);
            wvt.initWithWaveModes(kMode=2,lMode=0,j=1,phi=0,u=0.01,sign=1);
            forcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=TestWVBottomWaveScatteringForcing.topography(wvt,20),name="restart scattering");
            wvt.removeAllForcing();
            wvt.addForcing(forcing);
            model = WVModel(wvt);
            model.setupIntegrator(integratorType="adaptive",absTolerance=1e-12,relTolerance=1e-10);
            checkpointTime = 300;
            finalTime = 600;
        end

        function index = modeIndex(wvt,kMode,lMode,j,branch)
            mask = wvt.J == j & abs(wvt.K-2*pi*kMode/wvt.Lx) < 100*eps(max(1,abs(2*pi*kMode/wvt.Lx))) & abs(wvt.L-2*pi*lMode/wvt.Ly) < 100*eps(max(1,abs(2*pi*lMode/wvt.Ly)));
            if branch == "plus"
                mask = mask & wvt.waveComponent.maskAp;
            else
                mask = mask & wvt.waveComponent.maskAm;
            end
            index = find(mask,1);
            if isempty(index)
                error("TestWVBottomWaveScatteringForcing:MissingMode", "Unable to locate the requested wave mode.")
            end
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
