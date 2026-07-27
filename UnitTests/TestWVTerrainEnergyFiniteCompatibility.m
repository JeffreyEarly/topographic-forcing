classdef TestWVTerrainEnergyFiniteCompatibility < matlab.unittest.TestCase
    % Verify the unmodified finite-terrain compatibility gate.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test)
        function flatAndUniformDepthAreExactReferences(testCase)
            for shouldAntialias = [false true]
                wvt = TestWVTerrainEnergyFiniteCompatibility.createTransform([4 4 5],shouldAntialias,@(z)2e-5+0*z);
                terrains = {zeros(wvt.Nx,wvt.Ny),100*ones(wvt.Nx,wvt.Ny)};
                for iTerrain = 1:numel(terrains)
                    problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=terrains{iTerrain});
                    audit = problem.auditFiniteTerrainCompatibility();
                    testCase.verifyEqual(audit.status,"compatible")
                    testCase.verifyTrue(audit.isCompatible)
                    testCase.verifyEmpty(audit.failedConditions)
                    TestWVTerrainEnergyFiniteCompatibility.verifyExactAudit(testCase,audit)
                    testCase.verifyLessThanOrEqual(norm(audit.bottomTendencyMatrix,"fro"),1e-15)
                    testCase.verifyLessThanOrEqual(norm(audit.bottomValueMatrix*audit.rawGenerator,"fro"),1e-12*max(norm(audit.rawGenerator,"fro"),1))

                    if iTerrain == 1
                        expected = TestWVTerrainEnergyFiniteCompatibility.globalFlatGenerator(problem);
                        testCase.verifyEqual(audit.rawGenerator,expected,"RelTol",1e-12,"AbsTol",1e-14)
                    end
                end
            end
        end

        function directBottomQuadratureConvergesToFlatOracle(testCase)
            resolution = [5 7 9 13];
            energyError = zeros(size(resolution));
            exchangeError = zeros(size(resolution));
            for iResolution = 1:numel(resolution)
                wvt = TestWVTerrainEnergyFiniteCompatibility.createTransform([4 4 resolution(iResolution)],false,@(z)2e-5+0*z);
                problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
                diagnostics = problem.constructionDiagnostics.finiteTerrain;
                energyError(iResolution) = diagnostics.quadratureEnergyRelativeResidual;
                exchangeError(iResolution) = diagnostics.quadratureExchangeRelativeResidual;
                testCase.verifySize(problem.finiteTerrainForms.quadratureEnergyMatrix,size(problem.finiteTerrainForms.energyMatrix))
                testCase.verifySize(problem.finiteTerrainForms.quadratureExchangeMatrix,size(problem.finiteTerrainForms.exchangeMatrix))
            end
            testCase.verifyTrue(all(diff(energyError) < 0))
            testCase.verifyTrue(all(diff(exchangeError) < 0))
            testCase.verifyLessThan(energyError(end),energyError(1)/5)
            testCase.verifyLessThan(exchangeError(end),exchangeError(1)/4)
        end

        function sinusoidalTerrainHasResolvedIncompatibility(testCase)
            problem = TestWVTerrainEnergyFiniteCompatibility.sinusoidalProblem([4 4 5],2,20,@(z)2e-5+0*z);
            formsBefore = problem.finiteTerrainForms;
            first = problem.auditFiniteTerrainCompatibility();
            second = problem.auditFiniteTerrainCompatibility();

            testCase.verifyEqual(first.status,"incompatible")
            testCase.verifyFalse(first.isCompatible)
            testCase.verifyTrue(all(ismember(["apv-compatibility";"bottom-evolution-compatibility";"quadratic-enstrophy-compatibility"],first.failedConditions)))
            testCase.verifyLessThanOrEqual(first.diagnostics.energyDefect,1e-12)
            testCase.verifyLessThanOrEqual(first.diagnostics.conjugacyDefect,1e-12)
            testCase.verifyGreaterThan(first.diagnostics.apvDefect,1e-10)
            testCase.verifyGreaterThan(first.diagnostics.bottomDefect,1e-10)
            testCase.verifyGreaterThan(first.diagnostics.enstrophyDefect,1e-10)
            testCase.verifyLessThanOrEqual(first.diagnostics.rawHermitianDefect,1e-12)
            testCase.verifyLessThanOrEqual(first.diagnostics.rawSkewHermitianDefect,1e-12)
            testCase.verifyGreaterThan(first.diagnostics.scaledEnergyRcond,1e-12)
            testCase.verifyEqual(first.diagnostics,second.diagnostics)
            testCase.verifyEqual(problem.finiteTerrainForms,formsBefore)
        end

        function diagnosticsResolveSubspacesRowsRanksAndRandomStates(testCase)
            problem = TestWVTerrainEnergyFiniteCompatibility.sinusoidalProblem([4 4 5],2,20,@(z)2e-5+0*z);
            audit = problem.auditFiniteTerrainCompatibility();
            diagnostics = audit.diagnostics;

            testCase.verifyEqual(string({diagnostics.subspace.name}),["waves" "apv-bearing-geostrophic" "bottom-inversion" "mean-density"])
            testCase.verifyEqual(sum([diagnostics.subspace.dimension]),height(problem.stateLayout))
            testCase.verifyTrue(all(isfinite([[diagnostics.subspace.energyDefect];[diagnostics.subspace.apvDefect];[diagnostics.subspace.bottomDefect];[diagnostics.subspace.enstrophyDefect]]),"all"))
            testCase.verifySize(diagnostics.apvResidualByVerticalLevel,[problem.originatingTransform.Nz 1])
            testCase.verifySize(diagnostics.unweightedAPVResidualByVerticalLevel,[problem.originatingTransform.Nz 1])
            testCase.verifySize(diagnostics.bottomResidualByHorizontalMode,[height(problem.horizontalLayout) 1])
            testCase.verifyEqual(norm(diagnostics.apvResidualByVerticalLevel),1,"RelTol",2e-13)
            testCase.verifyEqual(norm(diagnostics.unweightedAPVResidualByVerticalLevel),1,"RelTol",2e-13)
            testCase.verifyGreaterThan(max(diagnostics.bottomResidualByHorizontalMode),0)
            testCase.verifyGreaterThan(diagnostics.generatorRank,diagnostics.apvNullity)
            testCase.verifyFalse(diagnostics.generatorRangeFitsAPVNullspace)
            testCase.verifyLessThanOrEqual(diagnostics.randomStates.complex.energyTendency,1e-12)
            testCase.verifyGreaterThan(diagnostics.randomStates.complex.apvTendency,1e-10)
            testCase.verifyGreaterThan(diagnostics.randomStates.complex.bottomTendency,1e-10)
            testCase.verifyGreaterThan(diagnostics.randomStates.complex.enstrophyTendency,1e-10)
            testCase.verifyLessThanOrEqual(diagnostics.randomStates.physicalConjugacyDefect,1e-13)
            testCase.verifyLessThanOrEqual(diagnostics.rawRestoredGeneratorAgreement,1e-12)
        end

        function incompatibilityPersistsUnderIndependentRefinement(testCase)
            amplitudes = [1 5 20];
            amplitudeMetrics = zeros(numel(amplitudes),3);
            for iCase = 1:numel(amplitudes)
                amplitudeMetrics(iCase,:) = TestWVTerrainEnergyFiniteCompatibility.metrics([4 4 5],2,amplitudes(iCase),@(z)2e-5+0*z);
            end
            testCase.verifyTrue(all(amplitudeMetrics > 1e-10,"all"))

            horizontalResolution = [4 6 8];
            horizontalMetrics = zeros(numel(horizontalResolution),3);
            for iCase = 1:numel(horizontalResolution)
                horizontalMetrics(iCase,:) = TestWVTerrainEnergyFiniteCompatibility.metrics([horizontalResolution(iCase) 4 5],2,20,@(z)2e-5+0*z);
            end
            horizontalFinalChange = abs(horizontalMetrics(end,:)-horizontalMetrics(end-1,:))./horizontalMetrics(end-1,:);
            testCase.verifyLessThan(horizontalFinalChange,0.25)

            verticalResolution = [5 7 9];
            verticalMetrics = zeros(numel(verticalResolution),3);
            for iCase = 1:numel(verticalResolution)
                verticalMetrics(iCase,:) = TestWVTerrainEnergyFiniteCompatibility.metrics([4 4 verticalResolution(iCase)],2,20,@(z)2e-5+0*z);
            end
            testCase.verifyTrue(all(diff(verticalMetrics(:,1)) < 0))
            testCase.verifyTrue(all(diff(verticalMetrics(:,2)) < 0))
            testCase.verifyTrue(all(diff(verticalMetrics(:,3)) < 0))

            oversampling = [1 2 3];
            oversamplingMetrics = zeros(numel(oversampling),3);
            for iCase = 1:numel(oversampling)
                oversamplingMetrics(iCase,:) = TestWVTerrainEnergyFiniteCompatibility.metrics([4 4 5],oversampling(iCase),20,@(z)2e-5+0*z);
            end
            oversamplingFinalChange = abs(oversamplingMetrics(end,:)-oversamplingMetrics(end-1,:))./oversamplingMetrics(end-1,:);
            testCase.verifyLessThan(oversamplingFinalChange,0.01)
            testCase.verifyTrue(all(oversamplingMetrics(end,:) > 1e-10))
        end

        function variableStratificationProducesFiniteStructuralDiagnostics(testCase)
            N0 = 5e-3;
            problem = TestWVTerrainEnergyFiniteCompatibility.sinusoidalProblem([4 4 7],2,20,@(z)N0^2*exp(2*z/1300));
            audit = problem.auditFiniteTerrainCompatibility();
            testCase.verifyEqual(audit.status,"incompatible")
            testCase.verifyLessThanOrEqual(audit.diagnostics.energyDefect,1e-12)
            testCase.verifyLessThanOrEqual(audit.diagnostics.conjugacyDefect,1e-12)
            testCase.verifyTrue(all(isfinite([audit.diagnostics.apvDefect audit.diagnostics.bottomDefect audit.diagnostics.enstrophyDefect])))
        end
    end

    methods (Static)
        function wvt = createTransform(resolution,shouldAntialias,N2)
            wvt = WVTransformBoussinesq([20e3 16e3 1000],resolution,N2Function=N2,latitude=30,shouldAntialias=shouldAntialias);
        end

        function problem = sinusoidalProblem(resolution,oversampling,amplitude,N2)
            wvt = TestWVTerrainEnergyFiniteCompatibility.createTransform(resolution,false,N2);
            h = amplitude*cos(2*pi*wvt.x/wvt.Lx).*ones(1,wvt.Ny);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h,horizontalOversamplingFactor=oversampling);
        end

        function values = metrics(resolution,oversampling,amplitude,N2)
            problem = TestWVTerrainEnergyFiniteCompatibility.sinusoidalProblem(resolution,oversampling,amplitude,N2);
            audit = problem.auditFiniteTerrainCompatibility();
            values = [audit.diagnostics.apvDefect audit.diagnostics.bottomDefect audit.diagnostics.enstrophyDefect];
        end

        function expected = globalFlatGenerator(problem)
            expected = zeros(height(problem.stateLayout));
            for iK = 1:numel(problem.flatModeBlocks)
                rows = find(problem.stateLayout.horizontalIndex == iK);
                block = problem.flatModeBlocks{iK};
                expected(rows,rows) = block.E\block.J;
            end
        end

        function verifyExactAudit(testCase,audit)
            diagnostics = audit.diagnostics;
            testCase.verifyLessThanOrEqual(diagnostics.energyDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.apvDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.bottomDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.enstrophyDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.conjugacyDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.preRestorationEnergyDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.preRestorationAPVDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.preRestorationBottomDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.preRestorationEnstrophyDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.preRestorationConjugacyDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.rawRestoredGeneratorAgreement,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.randomStates.complex.energyTendency,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.randomStates.complex.apvTendency,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.randomStates.complex.bottomTendency,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.randomStates.complex.enstrophyTendency,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.randomStates.physical.energyTendency,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.randomStates.physical.apvTendency,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.randomStates.physical.bottomTendency,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.randomStates.physical.enstrophyTendency,1e-12)
        end
    end
end
