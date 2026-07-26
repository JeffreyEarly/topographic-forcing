classdef TestWVTerrainEnergyFlatOracle < matlab.unittest.TestCase
    % Verify recovery of the complete flat nonhydrostatic problem.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test)
        function constantStratificationDispersionAndStructure(testCase)
            for shouldAntialias = [false true]
                N0 = 5e-3;
                wvt = TestWVTerrainEnergyFlatOracle.createTransform(@(z)N0^2+0*z,shouldAntialias,7);
                problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
                testCase.verifyLessThanOrEqual(problem.constructionDiagnostics.maximumHermitianDefect,1e-13)
                testCase.verifyLessThanOrEqual(problem.constructionDiagnostics.maximumSkewHermitianDefect,1e-13)
                testCase.verifyGreaterThan(problem.constructionDiagnostics.minimumScaledEnergyRcond,1e-12)
                testCase.verifyLessThanOrEqual(problem.constructionDiagnostics.maximumEigenResidual,1e-12)
                for iK = 1:numel(problem.flatModeBlocks)
                    block = problem.flatModeBlocks{iK};
                    testCase.verifyLessThanOrEqual(norm(block.eigenvectors'*block.E*block.eigenvectors-eye(size(block.E)),"fro"),1e-11)
                    state = randn(size(block.E,1),1)+1i*randn(size(block.E,1),1);
                    tendencyScale = max(norm(block.J,"fro")*norm(state)^2,realmin);
                    testCase.verifyLessThanOrEqual(abs(real(state'*block.J*state))/tendencyScale,1e-13)
                end

                retained = find(problem.horizontalLayout.isPrimary & (problem.horizontalLayout.kMode ~= 0 | problem.horizontalLayout.lMode ~= 0));
                for iK = reshape(retained,1,[])
                    block = problem.flatModeBlocks{iK};
                    positive = block.frequency(block.frequency > 1e-10);
                    j = problem.verticalModeIndices(problem.verticalModeIndices > 0);
                    m = j*pi/wvt.Lz;
                    kappa = hypot(problem.horizontalLayout.k(iK),problem.horizontalLayout.l(iK));
                    expected = sort(sqrt((N0^2*kappa^2+wvt.f^2*m.^2)./(kappa^2+m.^2)));
                    testCase.verifyEqual(positive,expected,"RelTol",2e-11,"AbsTol",1e-13)
                    TestWVTerrainEnergyFlatOracle.verifyWaveAndNullSpaces(testCase,problem,iK)
                end
            end
        end

        function zeroWavenumberSubspaces(testCase)
            wvt = TestWVTerrainEnergyFlatOracle.createTransform(@(z)2e-5+0*z,false,7);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
            block = problem.flatModesForHorizontalMode(0,0);
            tolerance = 1e-10*max(abs(block.frequency));
            testCase.verifyEqual(nnz(block.frequency < -tolerance),wvt.Nj)
            testCase.verifyEqual(nnz(block.frequency > tolerance),wvt.Nj)
            testCase.verifyEqual(nnz(abs(block.frequency) <= tolerance),wvt.Nj)
            testCase.verifyEqual(block.frequency(1:wvt.Nj),-wvt.f*ones(wvt.Nj,1),"RelTol",1e-12)
            testCase.verifyEqual(block.frequency(end-wvt.Nj+1:end),wvt.f*ones(wvt.Nj,1),"RelTol",1e-12)
        end

        function variableStratificationMatchesBoussinesqTransform(testCase)
            N0 = 5e-3;
            N2 = @(z)N0^2*exp(2*z/1300);
            frequencyError = zeros(3,1);
            resolutions = [7 9 13];
            for iResolution = 1:numel(resolutions)
                wvt = TestWVTerrainEnergyFlatOracle.createTransform(N2,false,resolutions(iResolution));
                problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
                block = problem.flatModesForHorizontalMode(1,0);
                positive = block.frequency(block.frequency > 1e-10);
                iNative = find(wvt.kMode_wv == 1 & wvt.lMode_wv == 0,1);
                expected = sort(wvt.Omega(wvt.j > 0,iNative));
                nCompare = min(3,numel(positive));
                frequencyError(iResolution) = norm(positive(1:nCompare)-expected(1:nCompare))/norm(expected(1:nCompare));
            end
            testCase.verifyTrue(all(diff(frequencyError) < 0))
            testCase.verifyLessThanOrEqual(frequencyError(end),2e-4)
            eigenfunctionError = TestWVTerrainEnergyFlatOracle.waveEigenfunctionErrors(problem,wvt,1,0);
            testCase.verifyLessThanOrEqual(max(eigenfunctionError(end-2:end)),5e-4)
            testCase.verifySize(problem.constructionDiagnostics.nonhydrostaticIndicator,[wvt.Nj-1 height(problem.horizontalLayout)])
            TestWVTerrainEnergyFlatOracle.verifyWaveAndNullSpaces(testCase,problem,find(problem.horizontalLayout.kMode == 1 & problem.horizontalLayout.lMode == 0,1))
        end

        function finiteTerrainAPVMapHasCorrectFlatLimit(testCase)
            wvt = TestWVTerrainEnergyFlatOracle.createTransform(@(z)2e-5+0*z,false,7);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
            iK = find(problem.horizontalLayout.kMode == 1 & problem.horizontalLayout.lMode == 0,1);
            block = problem.flatModeBlocks{iK};
            [~,iMode] = max(block.frequency);
            a = zeros(height(problem.stateLayout),1);
            rows = problem.stateLayout.horizontalIndex == iK;
            a(rows) = block.eigenvectors(:,iMode);
            partnerK = find(problem.horizontalLayout.kMode == -1 & problem.horizontalLayout.lMode == 0,1);
            partnerRows = problem.stateLayout.horizontalIndex == partnerK;
            a(partnerRows) = conj(a(problem.conjugateCoordinateIndex(partnerRows)));
            q = problem.evaluateAPV(a);
            fields = problem.reconstructState(a,outputDomain="spatial");
            scale = max(wvt.f*max(abs(fields.etaHatXi),[],"all"),hypot(problem.horizontalLayout.k(iK),problem.horizontalLayout.l(iK))*max(hypot(fields.uHat,fields.vHat),[],"all"));
            testCase.verifyLessThanOrEqual(max(abs(q),[],"all")/max(scale,realmin),2e-10)
        end
    end

    methods (Static)
        function wvt = createTransform(N2,shouldAntialias,Nz)
            wvt = WVTransformBoussinesq([20e3 20e3 1000],[4 4 Nz],N2Function=N2,latitude=30,shouldAntialias=shouldAntialias);
        end

        function verifyWaveAndNullSpaces(testCase,problem,iK)
            block = problem.flatModeBlocks{iK};
            Nj = numel(problem.verticalModeIndices);
            tolerance = 1e-10*max(abs(block.frequency));
            testCase.verifyEqual(nnz(block.frequency < -tolerance),Nj-1)
            testCase.verifyEqual(nnz(block.frequency > tolerance),Nj-1)
            testCase.verifyEqual(nnz(abs(block.frequency) <= tolerance),Nj+1)

            rows = problem.stateLayout.horizontalIndex == iK;
            bottomLocal = problem.stateLayout.component(rows) == "etaB";
            waveModes = find(abs(block.frequency) > tolerance);
            for iMode = reshape(waveModes,1,[])
                c = block.eigenvectors(:,iMode);
                basis = problem.basisBlocks{iK};
                q = block.Q*c;
                qScale = norm(1i*problem.horizontalLayout.k(iK)*basis.vHat*c)+norm(1i*problem.horizontalLayout.l(iK)*basis.uHat*c)+norm(problem.originatingTransform.f*basis.etaHatXi*c);
                testCase.verifyLessThanOrEqual(norm(q)/max(qScale,realmin),2e-10)
                eta = basis.etaHat*c;
                testCase.verifyLessThanOrEqual(abs(c(bottomLocal))/max(norm(eta),realmin),2e-10)
            end
        end

        function error = waveEigenfunctionErrors(problem,wvt,kMode,lMode)
            iK = find(problem.horizontalLayout.kMode == kMode & problem.horizontalLayout.lMode == lMode,1);
            block = problem.flatModeBlocks{iK};
            positive = find(block.frequency > 1e-10);
            basis = problem.basisBlocks{iK};
            iNative = find(wvt.kMode_wv == kMode & wvt.lMode_wv == lMode,1);
            F = wvt.FwInvMatrix(kMode,lMode);
            G = wvt.GwInvMatrix(kMode,lMode);
            internal = find(wvt.j > 0);
            [~,order] = sort(wvt.Omega(internal,iNative));
            zWeight = wvt.z_int(:);
            N2Weight = zWeight.*wvt.N2(:);
            error = zeros(numel(positive),1);
            for iMode = 1:numel(positive)
                c = block.eigenvectors(:,positive(iMode));
                computed = struct("u",basis.uHat*c,"v",basis.vHat*c,"w",basis.wHat*c,"eta",basis.etaHat*c);
                j = internal(order(iMode));
                exact = struct("u",F(:,j)*wvt.UAm(j,iNative),"v",F(:,j)*wvt.VAm(j,iNative), ...
                    "w",G(:,j)*wvt.WAm(j,iNative),"eta",G(:,j)*wvt.NAm(j,iNative));
                overlap = abs(TestWVTerrainEnergyFlatOracle.energyProduct(computed,exact,zWeight,N2Weight,wvt.rho0));
                overlap = overlap/sqrt(real(TestWVTerrainEnergyFlatOracle.energyProduct(computed,computed,zWeight,N2Weight,wvt.rho0)*TestWVTerrainEnergyFlatOracle.energyProduct(exact,exact,zWeight,N2Weight,wvt.rho0)));
                error(iMode) = sqrt(max(0,1-overlap^2));
            end
        end

        function value = energyProduct(first,second,zWeight,N2Weight,rho0)
            value = rho0*(first.u'*(zWeight.*second.u)+first.v'*(zWeight.*second.v)+first.w'*(zWeight.*second.w)+first.eta'*(N2Weight.*second.eta));
        end
    end
end
