classdef TestWVTerrainEnergyBottomInversion < matlab.unittest.TestCase
    % Verify the complete balanced bottom-inversion coordinate.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test)
        function constantStratificationMatchesAnalyticSolution(testCase)
            N0 = 5e-3;
            wvt = TestWVTerrainEnergyBottomInversion.createTransform(@(z)N0^2+0*z,7);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
            positive = find(hypot(problem.horizontalLayout.k,problem.horizontalLayout.l) > 0);
            for iK = reshape(positive,1,[])
                profile = problem.bottomInversionProfiles{iK};
                mu = profile.kappa*N0/abs(wvt.f);
                eta = sinh(-mu*wvt.z(:))/sinh(mu*wvt.Lz);
                etaXi = -mu*cosh(-mu*wvt.z(:))/sinh(mu*wvt.Lz);
                psi = -(wvt.f/profile.kappa^2)*etaXi;
                testCase.verifyEqual(profile.eta,eta,"RelTol",1e-11,"AbsTol",2e-13)
                testCase.verifyEqual(profile.etaXi,etaXi,"RelTol",1e-11,"AbsTol",2e-13)
                testCase.verifyEqual(profile.psi,psi,"RelTol",1e-11,"AbsTol",2e-13)
            end
            diagnostics = problem.constructionDiagnostics.bottomInversion;
            testCase.verifyLessThanOrEqual(diagnostics.maximumBottomValueError,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.maximumSurfaceValueError,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.maximumInversionResidual,1e-11)
        end

        function variableStratificationResidualAndRefinement(testCase)
            N0 = 5e-3;
            N2 = @(z)N0^2*exp(2*z/1300);
            residual = zeros(3,1);
            for i = 1:3
                Nz = [7 9 13];
                wvt = TestWVTerrainEnergyBottomInversion.createTransform(N2,Nz(i));
                problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
                residual(i) = problem.constructionDiagnostics.bottomInversion.maximumInversionResidual;
            end
            testCase.verifyLessThanOrEqual(max(residual),2e-10)
            testCase.verifyLessThanOrEqual(problem.constructionDiagnostics.bottomInversion.maximumAPVResidual,1e-13)
        end

        function linkedFieldsAreBalancedAndZeroAPV(testCase)
            wvt = TestWVTerrainEnergyBottomInversion.createTransform(@(z)2e-5+0*z,7);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
            iK = find(problem.horizontalLayout.kMode == 1 & problem.horizontalLayout.lMode == 1,1);
            rows = find(problem.stateLayout.horizontalIndex == iK);
            localBottom = find(problem.stateLayout.component(rows) == "etaB",1);
            block = problem.basisBlocks{iK};
            profile = problem.bottomInversionProfiles{iK};
            k = problem.horizontalLayout.k(iK);
            l = problem.horizontalLayout.l(iK);
            rho0 = wvt.rho0;

            horizontalScale = max(abs(wvt.f)*norm([block.uHat(:,localBottom);block.vHat(:,localBottom)]),realmin);
            horizontalResidual = norm([-wvt.f*block.vHat(:,localBottom)+(1i*k/rho0)*profile.pressure; ...
                wvt.f*block.uHat(:,localBottom)+(1i*l/rho0)*profile.pressure]);
            verticalScale = max(norm(wvt.N2(:).*block.etaHat(:,localBottom)),realmin);
            verticalResidual = norm(profile.psiXi*wvt.f+wvt.N2(:).*block.etaHat(:,localBottom));
            q = 1i*k*block.vHat(:,localBottom)-1i*l*block.uHat(:,localBottom)-wvt.f*block.etaHatXi(:,localBottom);

            testCase.verifyLessThanOrEqual(horizontalResidual/horizontalScale,2e-12)
            testCase.verifyLessThanOrEqual(verticalResidual/verticalScale,2e-12)
            testCase.verifyEqual(block.wHat(:,localBottom),zeros(wvt.Nz,1),"AbsTol",0)
            testCase.verifyLessThanOrEqual(norm(q)/max(norm(1i*k*block.vHat(:,localBottom))+norm(1i*l*block.uHat(:,localBottom)),realmin),1e-13)

            flat = problem.flatModeBlocks{iK};
            testCase.verifyEqual(flat.J(:,localBottom),zeros(size(flat.J,1),1),"AbsTol",0)
            testCase.verifyEqual(flat.J(localBottom,:),zeros(1,size(flat.J,2)),"AbsTol",0)
            testCase.verifyEqual(flat.Q(:,localBottom),zeros(wvt.Nz,1),"AbsTol",0)
            testCase.verifyGreaterThan(real(flat.E(localBottom,localBottom)),0)
        end

        function zeroWavenumberCoordinateIsMeanDensityState(testCase)
            wvt = TestWVTerrainEnergyBottomInversion.createTransform(@(z)2e-5+0*z,7);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
            iK = find(problem.horizontalLayout.kMode == 0 & problem.horizontalLayout.lMode == 0,1);
            rows = find(problem.stateLayout.horizontalIndex == iK);
            localBottom = find(problem.stateLayout.component(rows) == "etaB",1);
            block = problem.basisBlocks{iK};
            testCase.verifyEqual(block.uHat(:,localBottom),zeros(wvt.Nz,1),"AbsTol",0)
            testCase.verifyEqual(block.vHat(:,localBottom),zeros(wvt.Nz,1),"AbsTol",0)
            testCase.verifyEqual(block.etaHat([1 end],localBottom),[1;0],"AbsTol",0)
            testCase.verifyGreaterThan(norm(-wvt.f*block.etaHatXi(:,localBottom)),0)
        end

        function zeroCoriolisIsRejected(testCase)
            wvt = WVTransformBoussinesq([20e3 20e3 1000],[4 4 7], ...
                N2Function=@(z)2e-5+0*z,latitude=0,shouldAntialias=false);
            testCase.verifyError(@()WVTerrainEnergyGalerkin.fromTopography(wvt, ...
                topographicHeight=zeros(wvt.Nx,wvt.Ny)),"WVTerrainEnergyGalerkin:ZeroCoriolis")
        end
    end

    methods (Static)
        function wvt = createTransform(N2,Nz)
            wvt = WVTransformBoussinesq([20e3 20e3 1000],[4 4 Nz], ...
                N2Function=N2,latitude=30,shouldAntialias=false);
        end
    end
end
