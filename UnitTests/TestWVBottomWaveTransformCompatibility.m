classdef TestWVBottomWaveTransformCompatibility < matlab.unittest.TestCase
    % Verify generation and scattering across supported transform geometries.

    properties
        wvt
    end

    properties (ClassSetupParameter)
        transformType = struct( ...
            boussinesq="boussinesq", ...
            constantNonhydrostatic="constantNonhydrostatic", ...
            constantHydrostatic="constantHydrostatic", ...
            hydrostatic="hydrostatic")
        shouldAntialias = struct(full=false,dealiased=true)
    end

    methods (TestClassSetup)
        function addRepositoryToPathAndCreateTransform(testCase,transformType,shouldAntialias)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(repositoryRoot)
            testCase.wvt = TestWVBottomWaveTransformCompatibility.createTransform(transformType,[8 6 5],shouldAntialias);
        end
    end

    methods (Test)
        function generationMatchesDirectPressureProjection(testCase)
            wvt = testCase.wvt;
            terrain = testCase.topography(wvt);
            forcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05+0.01i; -0.02+0.015i],frequency=1.405e-4);
            for t = [0 813]
                wvt.t = t;
                velocity = forcing.barotropicVelocityAtTime(t);
                [directFp,directFm] = referenceBoundaryProjection(wvt,terrain,velocity);
                [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
                testCase.verifyLessThanOrEqual(testCase.relativeError(Fp,directFp),1e-10)
                testCase.verifyLessThanOrEqual(testCase.relativeError(Fm,directFm),1e-10)
                testCase.verifyEqual(F0,zeros(size(wvt.A0)))
            end
        end

        function scatteringMatchesNativeEndpointFieldsAndProjection(testCase)
            wvt = testCase.wvt;
            terrain = testCase.topography(wvt);
            forcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=terrain);
            previousRandomState = rng;
            randomStateCleanup = onCleanup(@()rng(previousRandomState));
            rng(82731,"twister")
            wvt.Ap = 1e-3*(randn(size(wvt.Ap))+1i*randn(size(wvt.Ap))).*wvt.waveComponent.maskAp;
            wvt.Am = 1e-3*(randn(size(wvt.Am))+1i*randn(size(wvt.Am))).*wvt.waveComponent.maskAm;
            wvt.A0 = complex(ones(size(wvt.A0)));
            wvt.t = 619;

            [gBottom,bottomFields] = forcing.bottomVelocityFromWaveState(wvt);
            [uReference,vReference,dWdzReference,gReference] = testCase.nativeEndpointReference(wvt,terrain);
            testCase.verifyLessThanOrEqual(testCase.relativeError(bottomFields.u,uReference),1e-12)
            testCase.verifyLessThanOrEqual(testCase.relativeError(bottomFields.v,vReference),1e-12)
            testCase.verifyLessThanOrEqual(testCase.relativeError(bottomFields.dWdz,dWdzReference),1e-12)
            testCase.verifyLessThanOrEqual(testCase.relativeError(gBottom,gReference),1e-12)

            [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            [directFp,directFm] = testCase.directProjection(wvt,gBottom);
            testCase.verifyLessThanOrEqual(testCase.relativeError(Fp,directFp),1e-10)
            testCase.verifyLessThanOrEqual(testCase.relativeError(Fm,directFm),1e-10)
            testCase.verifyEqual(F0,zeros(size(wvt.A0)))
            clear randomStateCleanup
        end
    end

    methods (Static)
        function wvt = createTransform(transformType,resolution,shouldAntialias)
            Lxyz = [40e3 30e3 2e3];
            switch transformType
                case "boussinesq"
                    wvt = WVTransformBoussinesq(Lxyz,resolution,N2=@(z) 2e-5*exp(z/4000),latitude=45,shouldAntialias=shouldAntialias);
                case "constantNonhydrostatic"
                    wvt = WVTransformConstantStratification(Lxyz,resolution,N0=sqrt(2e-5),latitude=45,shouldAntialias=shouldAntialias,isHydrostatic=false);
                case "constantHydrostatic"
                    wvt = WVTransformConstantStratification(Lxyz,resolution,N0=sqrt(2e-5),latitude=45,shouldAntialias=shouldAntialias,isHydrostatic=true);
                case "hydrostatic"
                    wvt = WVTransformHydrostatic(Lxyz,resolution,N2=@(z) 2e-5*exp(z/4000),latitude=45,shouldAntialias=shouldAntialias);
                otherwise
                    error("TestWVBottomWaveTransformCompatibility:UnknownTransform", "Unknown transform type '%s'.", transformType)
            end
        end

        function terrain = topography(wvt)
            [x,y] = ndgrid(wvt.x,wvt.y);
            terrain = 30*cos(2*pi*x/wvt.Lx)+17*sin(2*pi*y/wvt.Ly)+11*cos(2*pi*(2*x/wvt.Lx+y/wvt.Ly));
        end
    end

    methods (Static, Access = private)
        function [uBottom,vBottom,dWdzBottom,gBottom] = nativeEndpointReference(wvt,terrain)
            phase = exp(wvt.iOmega*(wvt.t-wvt.t0));
            Apt = wvt.waveComponent.maskAp.*wvt.Ap.*phase;
            Amt = wvt.waveComponent.maskAm.*wvt.Am.*conj(phase);
            u = wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*Apt+wvt.UAm.*Amt);
            v = wvt.transformToSpatialDomainWithF(Apm=wvt.VAp.*Apt+wvt.VAm.*Amt);
            verticalVelocity = wvt.transformToSpatialDomainWithG(Apm=wvt.WAp.*Apt+wvt.WAm.*Amt);
            dWdz = wvt.diffZG(verticalVelocity);
            [~,iBottom] = min(wvt.z);
            uBottom = u(:,:,iBottom);
            vBottom = v(:,:,iBottom);
            dWdzBottom = dWdz(:,:,iBottom);
            gBottom = uBottom.*wvt.diffX(terrain)+vBottom.*wvt.diffY(terrain)-terrain.*dWdzBottom;
        end

        function [directFp,directFm] = directProjection(wvt,gBottom)
            gBottomTransform = wvt.transformFromSpatialDomainWithFourier(repmat(gBottom,1,1,wvt.Nz));
            gBottomFourier = gBottomTransform(1,:);
            directFp = complex(zeros(size(wvt.Ap)));
            directFm = complex(zeros(size(wvt.Am)));
            [~,iBottom] = min(wvt.z);

            for index = reshape(find(wvt.waveComponent.maskAp),1,[])
                [~,iHorizontal] = ind2sub(size(wvt.Ap),index);
                coefficient = complex(zeros(size(wvt.Ap)));
                coefficient(index) = wvt.NAp(index);
                pressure = wvt.g*wvt.transformToSpatialDomainWithF(Apm=coefficient);
                pressureFourier = wvt.transformFromSpatialDomainWithFourier(pressure);
                pressureBottom = pressureFourier(iBottom,iHorizontal)*wvt.phase(index);
                directFp(index) = conj(pressureBottom)*gBottomFourier(iHorizontal)/wvt.Apm_TE_factor(index);
            end

            for index = reshape(find(wvt.waveComponent.maskAm),1,[])
                [~,iHorizontal] = ind2sub(size(wvt.Am),index);
                coefficient = complex(zeros(size(wvt.Am)));
                coefficient(index) = wvt.NAm(index);
                pressure = wvt.g*wvt.transformToSpatialDomainWithF(Apm=coefficient);
                pressureFourier = wvt.transformFromSpatialDomainWithFourier(pressure);
                pressureBottom = pressureFourier(iBottom,iHorizontal)*wvt.conjPhase(index);
                directFm(index) = conj(pressureBottom)*gBottomFourier(iHorizontal)/wvt.Apm_TE_factor(index);
            end
        end

        function errorValue = relativeError(actual,expected)
            errorValue = norm(actual(:)-expected(:))/max(norm(expected(:)),eps);
        end
    end
end
