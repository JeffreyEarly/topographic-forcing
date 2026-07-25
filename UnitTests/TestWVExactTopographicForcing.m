classdef TestWVExactTopographicForcing < matlab.unittest.TestCase
    % Verify mapped geometry, terrain tendencies, and model integration.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(repositoryRoot);
        end
    end

    methods (Test)
        function testConstructorValidation(testCase)
            wvt = TestWVExactTopographicForcing.createTransform();
            flat = zeros(wvt.Nx,wvt.Ny);
            uniform = 100*ones(wvt.Nx,wvt.Ny);
            sinusoidal = TestWVExactTopographicForcing.sinusoidalTopography(wvt,50);

            flatForce = WVExactTopographicForcing(wvt,topographicHeight=flat);
            uniformForce = WVExactTopographicForcing(wvt,topographicHeight=uniform);
            sinusoidalForce = WVExactTopographicForcing(wvt,topographicHeight=sinusoidal,name="terrain");
            testCase.verifyEqual(flatForce.gamma,ones(wvt.Nx,wvt.Ny));
            testCase.verifyEqual(uniformForce.gamma,0.95*ones(wvt.Nx,wvt.Ny),AbsTol=10*eps);
            testCase.verifyEqual(sinusoidalForce.name,"terrain");

            testCase.verifyError(@() WVExactTopographicForcing(wvt,topographicHeight=zeros(wvt.Nx-1,wvt.Ny)), ...
                "WVExactTopographicForcing:InvalidTopographicHeightSize");
            testCase.verifyError(@() WVExactTopographicForcing(wvt,topographicHeight=1i*ones(wvt.Nx,wvt.Ny)), ...
                "WVExactTopographicForcing:InvalidTopographicHeight");
            invalid = flat;
            invalid(1) = NaN;
            testCase.verifyError(@() WVExactTopographicForcing(wvt,topographicHeight=invalid), ...
                "WVExactTopographicForcing:InvalidTopographicHeight");
            testCase.verifyError(@() WVExactTopographicForcing(wvt,topographicHeight=wvt.Lz*ones(wvt.Nx,wvt.Ny)), ...
                "WVExactTopographicForcing:NonpositiveMappedDepth");

            hydrostatic = WVTransformHydrostatic([wvt.Lx wvt.Ly wvt.Lz],[wvt.Nx wvt.Ny wvt.Nz], ...
                N2=@(z) 2e-5*ones(size(z)),latitude=45,shouldAntialias=false);
            testCase.verifyError(@() WVExactTopographicForcing(hydrostatic,topographicHeight=flat), ...
                "WVExactTopographicForcing:UnsupportedTransform");

            variableN2 = WVTransformBoussinesq([wvt.Lx wvt.Ly wvt.Lz],[wvt.Nx wvt.Ny wvt.Nz], ...
                N2=@(z) 2e-5*exp(z/4000),latitude=45,shouldAntialias=false);
            testCase.verifyError(@() WVExactTopographicForcing(variableN2,topographicHeight=flat), ...
                "WVExactTopographicForcing:NonconstantStratificationUnsupported");

            antialiased = WVTransformBoussinesq([wvt.Lx wvt.Ly wvt.Lz],[wvt.Nx wvt.Ny wvt.Nz], ...
                N2=@(z) 2e-5*ones(size(z)),latitude=45,shouldAntialias=true);
            antialiasedTerrain = zeros(antialiased.Nx,antialiased.Ny);
            testCase.verifyError(@() WVExactTopographicForcing(antialiased,topographicHeight=antialiasedTerrain), ...
                "WVExactTopographicForcing:AntialiasingUnsupported");
        end

        function testMappedGeometry(testCase)
            wvt = TestWVExactTopographicForcing.createTransform();
            wvt.initWithWaveModes(kMode=1,lMode=1,j=1,phi=0.3,u=0.01,sign=1);

            flatForce = WVExactTopographicForcing(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
            [u,v,w] = flatForce.physicalVelocity();
            testCase.verifyEqual(u,wvt.u,AbsTol=1e-14);
            testCase.verifyEqual(v,wvt.v,AbsTol=1e-14);
            testCase.verifyEqual(w,wvt.w,AbsTol=1e-14);
            expectedZ = ones(wvt.Nx,wvt.Ny).*reshape(wvt.z,1,1,[]);
            testCase.verifyEqual(flatForce.physicalZ,expectedZ,AbsTol=1e-14);

            field = reshape(1:prod(wvt.spatialMatrixSize),wvt.spatialMatrixSize);
            testCase.verifyEqual(flatForce.mappedVolumeIntegral(field),wvt.volumeIntegral(field),RelTol=1e-14);
            testCase.verifyError(@() flatForce.mappedVolumeIntegral(zeros(wvt.Nx,wvt.Ny)), ...
                "WVExactTopographicForcing:InvalidSpatialFieldSize");

            h0 = 200;
            uniformForce = WVExactTopographicForcing(wvt,topographicHeight=h0*ones(wvt.Nx,wvt.Ny));
            gamma0 = 1-h0/wvt.Lz;
            [u,v,w] = uniformForce.physicalVelocity();
            testCase.verifyEqual(u,wvt.u/gamma0,RelTol=1e-14,AbsTol=1e-14);
            testCase.verifyEqual(v,wvt.v/gamma0,RelTol=1e-14,AbsTol=1e-14);
            testCase.verifyEqual(w,wvt.w,RelTol=1e-14,AbsTol=1e-14);
            testCase.verifyEqual(uniformForce.physicalZ,gamma0*expectedZ,RelTol=1e-14,AbsTol=1e-14);

            terrainForce = WVExactTopographicForcing(wvt,topographicHeight=TestWVExactTopographicForcing.sinusoidalTopography(wvt,50));
            residual = terrainForce.bottomKinematicResidual();
            velocityScale = max(abs(terrainForce.physicalVelocity()),[],"all");
            testCase.verifyLessThanOrEqual(max(abs(residual),[],"all"),1e-12*max(1,velocityScale));
        end

        function testLinearTerrainKernel(testCase)
            wvt = TestWVExactTopographicForcing.createTransform();
            wvt.initWithWaveModes(kMode=1,lMode=1,j=1,phi=0.4,u=0.02,sign=1);
            wvt.addWaveModes(kMode=2,lMode=0,j=2,phi=1.1,u=0.01,sign=-1);
            wvt.t = 700;

            flatForce = WVExactTopographicForcing(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
            [Fu,Fv,Fw,Feta] = flatForce.linearTerrainTendency();
            testCase.verifyEqual(Fu,zeros(wvt.spatialMatrixSize));
            testCase.verifyEqual(Fv,zeros(wvt.spatialMatrixSize));
            testCase.verifyEqual(Fw,zeros(wvt.spatialMatrixSize));
            testCase.verifyEqual(Feta,zeros(wvt.spatialMatrixSize));

            gamma0 = 0.9;
            uniformForce = WVExactTopographicForcing(wvt,topographicHeight=(1-gamma0)*wvt.Lz*ones(wvt.Nx,wvt.Ny));
            [Fu,Fv,Fw,Feta] = uniformForce.linearTerrainTendency();
            [~,~,p] = wvt.variableWithName('u','v','p');
            expectedFu = -(gamma0-1)*wvt.diffX(p)/wvt.rho0;
            expectedFv = -(gamma0-1)*wvt.diffY(p)/wvt.rho0;
            expectedFw = -(1/gamma0-1)*wvt.diffZF(p)/wvt.rho0;
            testCase.verifyEqual(Fu,expectedFu,RelTol=1e-13,AbsTol=1e-14);
            testCase.verifyEqual(Fv,expectedFv,RelTol=1e-13,AbsTol=1e-14);
            testCase.verifyEqual(Fw,expectedFw,RelTol=1e-13,AbsTol=1e-14);
            testCase.verifyEqual(Feta,zeros(wvt.spatialMatrixSize));

            terrain = TestWVExactTopographicForcing.sinusoidalTopography(wvt,50);
            force = WVExactTopographicForcing(wvt,topographicHeight=terrain);
            tBefore = wvt.t;
            ApBefore = wvt.Ap;
            AmBefore = wvt.Am;
            A0Before = wvt.A0;
            [Fu,Fv,Fw,Feta] = force.linearTerrainTendency();
            [FuExpected,FvExpected,FwExpected,FetaExpected] = TestWVExactTopographicForcing.directTerrainTendency(wvt,terrain);
            testCase.verifyEqual(Fu,FuExpected,RelTol=1e-13,AbsTol=1e-14);
            testCase.verifyEqual(Fv,FvExpected,RelTol=1e-13,AbsTol=1e-14);
            testCase.verifyEqual(Fw,FwExpected,RelTol=1e-13,AbsTol=1e-14);
            testCase.verifyEqual(Feta,FetaExpected,RelTol=1e-13,AbsTol=1e-14);
            testCase.verifyTrue(isreal(Fu) && isreal(Fv) && isreal(Fw) && isreal(Feta));
            testCase.verifyEqual(wvt.t,tBefore);
            testCase.verifyEqual(wvt.Ap,ApBefore);
            testCase.verifyEqual(wvt.Am,AmBefore);
            testCase.verifyEqual(wvt.A0,A0Before);
            [FuRepeat,FvRepeat,FwRepeat,FetaRepeat] = force.linearTerrainTendency();
            testCase.verifyEqual(FuRepeat,Fu);
            testCase.verifyEqual(FvRepeat,Fv);
            testCase.verifyEqual(FwRepeat,Fw);
            testCase.verifyEqual(FetaRepeat,Feta);

            [Fp,Fm] = wvt.transformUVWEtaToWaveVortex(Fu,Fv,Fw,Feta);
            testCase.verifyEqual(Fm(:,1),conj(Fp(:,1)),AbsTol=1e-13);
        end

        function testForcingProjectionAndAccumulation(testCase)
            wvt = TestWVExactTopographicForcing.createTransform();
            wvt.initWithWaveModes(kMode=1,lMode=1,j=1,phi=0.2,u=0.02,sign=1);
            force = WVExactTopographicForcing(wvt,topographicHeight=TestWVExactTopographicForcing.sinusoidalTopography(wvt,50));
            wvt.removeAllForcing();
            wvt.addForcing(force);

            TestWVExactTopographicForcing.verifyProjectionAtTime(testCase,wvt,force,0);
            TestWVExactTopographicForcing.verifyProjectionAtTime(testCase,wvt,force,900);

            base = reshape(linspace(-1,1,prod(wvt.spatialMatrixSize)),wvt.spatialMatrixSize);
            [terrainFu,terrainFv,terrainFw,terrainFeta] = force.linearTerrainTendency();
            [Fu,Fv,Fw,Feta] = force.addNonhydrostaticSpatialForcing(wvt,base,2*base,3*base,4*base);
            testCase.verifyEqual(Fu,base+terrainFu);
            testCase.verifyEqual(Fv,2*base+terrainFv);
            testCase.verifyEqual(Fw,3*base+terrainFw);
            testCase.verifyEqual(Feta,4*base+terrainFeta);

            otherTransform = TestWVExactTopographicForcing.createTransform();
            zeroField = zeros(otherTransform.spatialMatrixSize);
            testCase.verifyError(@() force.addNonhydrostaticSpatialForcing(otherTransform,zeroField,zeroField,zeroField,zeroField), ...
                "WVExactTopographicForcing:TransformMismatch");
            testCase.verifyError(@() force.forcingWithResolutionOfTransform(otherTransform), ...
                "WVExactTopographicForcing:ResolutionChangeUnsupported");
        end

        function testAdaptiveModelSmoke(testCase)
            wvt = TestWVExactTopographicForcing.createTransform();
            wvt.initWithWaveModes(kMode=1,lMode=0,j=1,phi=0,u=0.01,sign=1);
            force = WVExactTopographicForcing(wvt,topographicHeight=TestWVExactTopographicForcing.sinusoidalTopography(wvt,50));
            wvt.removeAllForcing();
            wvt.addForcing(force);
            [Fp,Fm,F0] = wvt.nonlinearFlux();
            testCase.verifyGreaterThan(norm([Fp(:); Fm(:); F0(:)]),0);

            ApInitial = wvt.Ap;
            AmInitial = wvt.Am;
            A0Initial = wvt.A0;
            warningState = warning;
            warningCleanup = onCleanup(@()warning(warningState));
            warning("off","all")
            model = WVModel(wvt);
            model.setupIntegrator(integratorType="adaptive",absTolerance=1e-10,relTolerance=1e-7);
            model.integrateToTime(200,shouldShowIntegrationDiagnostics=false,callback=@(~) []);

            testCase.verifyEqual(wvt.t,200,AbsTol=1e-10);
            testCase.verifyTrue(all(isfinite([wvt.Ap(:); wvt.Am(:); wvt.A0(:)])));
            coefficientChange = norm([wvt.Ap(:)-ApInitial(:); wvt.Am(:)-AmInitial(:); wvt.A0(:)-A0Initial(:)]);
            testCase.verifyGreaterThan(coefficientChange,0);
        end
    end

    methods (Static, Access = private)
        function wvt = createTransform()
            wvt = WVTransformBoussinesq([4e3 4e3 2e3],[8 4 5], ...
                N2=@(z) 2e-5*ones(size(z)),latitude=45,shouldAntialias=false);
        end

        function terrain = sinusoidalTopography(wvt,amplitude)
            x = repmat(reshape(wvt.x,[],1),1,wvt.Ny);
            terrain = amplitude*cos(2*pi*x/wvt.Lx);
        end

        function [Fu,Fv,Fw,Feta] = directTerrainTendency(wvt,terrain)
            gamma = 1-terrain/wvt.Lz;
            xi = reshape(wvt.z,1,1,[]);
            dGammaDx = wvt.diffX(gamma);
            dGammaDy = wvt.diffY(gamma);
            [hatU,hatV,p] = wvt.variableWithName('u','v','p');
            dPdx = wvt.diffX(p);
            dPdy = wvt.diffY(p);
            dPdxi = wvt.diffZF(p);
            Tu = ((gamma-1).*dPdx-xi.*dGammaDx.*dPdxi)/wvt.rho0;
            Tv = ((gamma-1).*dPdy-xi.*dGammaDy.*dPdxi)/wvt.rho0;
            HuL = wvt.f*hatV-dPdx/wvt.rho0-Tu;
            HvL = -wvt.f*hatU-dPdy/wvt.rho0-Tv;
            Tw = xi.*(HuL.*dGammaDx./gamma+HvL.*dGammaDy./gamma) ...
                +(1./gamma-1).*dPdxi/wvt.rho0;
            Teta = -xi.*(hatU.*dGammaDx./gamma+hatV.*dGammaDy./gamma);
            Fu = -Tu;
            Fv = -Tv;
            Fw = -Tw;
            Feta = -Teta;
        end

        function verifyProjectionAtTime(testCase,wvt,force,t)
            wvt.t = t;
            [Fu,Fv,Fw,Feta] = force.linearTerrainTendency();
            [expectedFp,expectedFm,expectedF0] = wvt.transformUVWEtaToWaveVortex(Fu,Fv,Fw,Feta);
            [Fp,Fm,F0] = wvt.nonlinearFlux();
            testCase.verifyEqual(Fp,expectedFp,RelTol=1e-12,AbsTol=1e-14);
            testCase.verifyEqual(Fm,expectedFm,RelTol=1e-12,AbsTol=1e-14);
            testCase.verifyEqual(F0,expectedF0,RelTol=1e-12,AbsTol=1e-14);
        end
    end
end
