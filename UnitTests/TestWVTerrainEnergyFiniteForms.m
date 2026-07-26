classdef TestWVTerrainEnergyFiniteForms < matlab.unittest.TestCase
    % Verify the dense finite-terrain energy, exchange, and APV forms.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test)
        function flatLimitAndStructure(testCase)
            for shouldAntialias = [false true]
                wvt = TestWVTerrainEnergyFiniteForms.createTransform(@(z)2e-5+0*z,shouldAntialias,5);
                problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
                forms = problem.finiteTerrainForms;
                diagnostics = problem.constructionDiagnostics.finiteTerrain;
                testCase.verifyLessThanOrEqual(diagnostics.hermitianDefect,1e-13)
                testCase.verifyLessThanOrEqual(diagnostics.skewHermitianDefect,1e-13)
                testCase.verifyGreaterThan(diagnostics.scaledEnergyRcond,1e-12)
                testCase.verifyLessThanOrEqual(diagnostics.flatEnergyRelativeResidual,1e-12)
                testCase.verifyLessThanOrEqual(diagnostics.flatExchangeRelativeResidual,1e-12)
                testCase.verifyLessThanOrEqual(diagnostics.oversamplingAdjointDefect,1e-12)
                testCase.verifyLessThanOrEqual(norm(forms.energyMatrix-forms.flatEnergyMatrix,"fro")/norm(forms.flatEnergyMatrix,"fro"),1e-12)
                testCase.verifyLessThanOrEqual(norm(forms.exchangeMatrix-forms.flatExchangeMatrix,"fro")/norm(forms.flatExchangeMatrix,"fro"),1e-12)
            end
        end

        function sinusoidalTerrainPreservesDiscreteEnergyStructure(testCase)
            N0 = 5e-3;
            N2 = @(z)N0^2*exp(2*z/1400);
            wvt = TestWVTerrainEnergyFiniteForms.createTransform(N2,false,5);
            h = 60*cos(2*pi*wvt.x/wvt.Lx).*ones(1,wvt.Ny);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h);
            forms = problem.finiteTerrainForms;
            diagnostics = problem.constructionDiagnostics.finiteTerrain;
            testCase.verifyLessThanOrEqual(diagnostics.hermitianDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.skewHermitianDefect,1e-12)
            testCase.verifyGreaterThan(diagnostics.scaledEnergyRcond,1e-12)
            testCase.verifyGreaterThan(norm(forms.energyMatrix-forms.flatEnergyMatrix,"fro"),0)
            testCase.verifyGreaterThan(norm(forms.exchangeMatrix-forms.flatExchangeMatrix,"fro"),0)
            a = complex((1:height(problem.stateLayout))',mod((1:height(problem.stateLayout))',11)-5);
            tendencyScale = max(norm(forms.exchangeMatrix,"fro")*norm(a)^2,realmin);
            testCase.verifyLessThanOrEqual(abs(real(a'*forms.exchangeMatrix*a))/tendencyScale,1e-13)
            nXY = prod(forms.oversampledSize(1:2));
            expectedN2 = N2(repmat(forms.oversampledGamma(:),wvt.Nz,1).*kron(wvt.z(:),ones(nXY,1)));
            testCase.verifyEqual(forms.oversampledN2(:),expectedN2,"RelTol",1e-13)
        end

        function repeatedConstructionIsDeterministic(testCase)
            wvt = TestWVTerrainEnergyFiniteForms.createTransform(@(z)2e-5+0*z,true,5);
            h = 40*cos(2*pi*wvt.x/wvt.Lx).*cos(2*pi*wvt.y.'/wvt.Ly);
            first = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h,horizontalOversamplingFactor=2);
            second = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h,horizontalOversamplingFactor=2);
            testCase.verifyEqual(first.finiteTerrainForms.energyMatrix,second.finiteTerrainForms.energyMatrix,"RelTol",1e-13)
            testCase.verifyEqual(first.finiteTerrainForms.exchangeMatrix,second.finiteTerrainForms.exchangeMatrix,"RelTol",1e-13)
            testCase.verifyEqual(first.finiteTerrainForms.apvMatrix,second.finiteTerrainForms.apvMatrix,"RelTol",1e-13)
        end

        function oversamplingConvergesForNonlinearTerrainWeights(testCase)
            wvt = TestWVTerrainEnergyFiniteForms.createTransform(@(z)2e-5+0*z,false,5);
            h = 140*cos(2*pi*wvt.x/wvt.Lx).*ones(1,wvt.Ny);
            problem1 = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h,horizontalOversamplingFactor=1);
            problem2 = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h,horizontalOversamplingFactor=2);
            problem3 = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h,horizontalOversamplingFactor=3);
            difference12 = norm(problem1.finiteTerrainForms.energyMatrix-problem2.finiteTerrainForms.energyMatrix,"fro");
            difference23 = norm(problem2.finiteTerrainForms.energyMatrix-problem3.finiteTerrainForms.energyMatrix,"fro");
            testCase.verifyLessThan(difference23,difference12)
        end

        function uniformDepthMatchesIndependentPhysicalDepth(testCase)
            N0 = 5e-3;
            gamma = 0.82;
            wvt = TestWVTerrainEnergyFiniteForms.createTransform(@(z)N0^2+0*z,false,7);
            h = (1-gamma)*wvt.Lz*ones(wvt.Nx,wvt.Ny);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h);
            iK = find(problem.horizontalLayout.kMode == 1 & problem.horizontalLayout.lMode == 0,1);
            rows = find(problem.stateLayout.horizontalIndex == iK);
            E = problem.finiteTerrainForms.energyMatrix(rows,rows);
            J = problem.finiteTerrainForms.exchangeMatrix(rows,rows);
            [C,frequency] = eig(1i*J,E,"vector");
            [frequency,order] = sort(real(frequency));
            C = C(:,order);
            isPositive = frequency > 1e-10;

            H = gamma*wvt.Lz;
            physical = WVTransformBoussinesq([wvt.Lx wvt.Ly H],[wvt.Nx wvt.Ny wvt.Nz], ...
                N2Function=@(z)N0^2+0*z,latitude=wvt.latitude,shouldAntialias=false,z=gamma*wvt.z);
            iNative = find(physical.kMode_wv == 1 & physical.lMode_wv == 0,1);
            expectedFrequency = sort(physical.Omega(physical.j > 0,iNative));
            testCase.verifyEqual(frequency(isPositive),expectedFrequency,"RelTol",2e-11,"AbsTol",1e-13)

            [~,iMode] = min(abs(frequency-expectedFrequency(1)));
            c = C(:,iMode)/sqrt(real(C(:,iMode)'*E*C(:,iMode)));
            basis = problem.basisBlocks{iK};
            computed = struct("u",basis.uHat*c/gamma,"v",basis.vHat*c/gamma,"w",basis.wHat*c,"eta",basis.etaHat*c);
            internal = find(physical.j > 0);
            [~,jOrder] = sort(physical.Omega(internal,iNative));
            j = internal(jOrder(1));
            F = physical.FwInvMatrix(1,0);
            G = physical.GwInvMatrix(1,0);
            exact = struct("u",F(:,j)*physical.UAm(j,iNative),"v",F(:,j)*physical.VAm(j,iNative), ...
                "w",G(:,j)*physical.WAm(j,iNative),"eta",G(:,j)*physical.NAm(j,iNative));
            overlap = TestWVTerrainEnergyFiniteForms.energyProduct(computed,exact,physical.z_int,physical.N2,physical.rho0);
            computedNorm = TestWVTerrainEnergyFiniteForms.energyProduct(computed,computed,physical.z_int,physical.N2,physical.rho0);
            exactNorm = TestWVTerrainEnergyFiniteForms.energyProduct(exact,exact,physical.z_int,physical.N2,physical.rho0);
            eigenfunctionError = sqrt(max(0,1-abs(overlap)^2/real(computedNorm*exactNorm)));
            testCase.verifyLessThanOrEqual(eigenfunctionError,1e-9)
        end
    end

    methods (Static)
        function wvt = createTransform(N2,shouldAntialias,Nz)
            wvt = WVTransformBoussinesq([20e3 16e3 1000],[4 4 Nz],N2Function=N2,latitude=30,shouldAntialias=shouldAntialias);
        end

        function value = energyProduct(first,second,zWeight,N2,rho0)
            value = rho0*(first.u'*(zWeight.*second.u)+first.v'*(zWeight.*second.v) ...
                +first.w'*(zWeight.*second.w)+first.eta'*(zWeight.*N2.*second.eta));
        end
    end
end
