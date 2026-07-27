classdef TestWVTerrainEnergyConstrainedClosure < matlab.unittest.TestCase
    % Verify the dense constrained-closure feasibility audit.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test)
        function flatAndUniformDepthAreAdmissible(testCase)
            for shouldAntialias = [false true]
                wvt = TestWVTerrainEnergyConstrainedClosure.createTransform([4 4 5],shouldAntialias);
                terrains = {zeros(wvt.Nx,wvt.Ny),100*ones(wvt.Nx,wvt.Ny)};
                for iTerrain = 1:numel(terrains)
                    problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=terrains{iTerrain});
                    audit = problem.auditConstrainedClosure();
                    testCase.verifyEqual(audit.status,"admissible")
                    testCase.verifyTrue(audit.constraintSystemIsFeasible)
                    testCase.verifyTrue(audit.waveSpaceIsComplete)
                    testCase.verifyTrue(audit.isScientificallyAdmissible)
                    testCase.verifyEmpty(audit.failedConditions)
                    testCase.verifyEqual(audit.diagnostics.apvNullity,audit.diagnostics.expectedZeroAPVDimension)
                    testCase.verifyLessThanOrEqual(audit.diagnostics.correctionNorm,1e-12)
                    testCase.verifyEqual(audit.constrainedGenerator,audit.rawGenerator,"RelTol",1e-12,"AbsTol",1e-15)
                    testCase.verifyEqual(audit.constrainedExchangeMatrix,audit.energyFactor'*audit.constrainedGenerator*audit.energyFactor,"RelTol",1e-13,"AbsTol",1e-15)
                    testCase.verifyLessThanOrEqual(audit.diagnostics.apvConstraintResidual,1e-12)
                    testCase.verifyLessThanOrEqual(audit.diagnostics.bottomConstraintResidual,1e-12)
                    testCase.verifyLessThanOrEqual(audit.diagnostics.energyConstraintResidual,1e-13)
                    testCase.verifyLessThanOrEqual(audit.diagnostics.constrainedConjugacyDefect,1e-12)
                    testCase.verifyLessThanOrEqual(audit.diagnostics.bottom.phaseGramDefect,1e-12)
                    TestWVTerrainEnergyConstrainedClosure.verifyRandomStateIdentities(testCase,audit)
                end
            end
        end

        function bottomMapsMatchDirectResolvedCalculation(testCase)
            wvt = TestWVTerrainEnergyConstrainedClosure.createTransform([6 4 5],false);
            h = 30*cos(2*pi*wvt.x/wvt.Lx).*ones(1,wvt.Ny);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h);
            audit = problem.auditConstrainedClosure();
            zeroK = find(problem.horizontalLayout.kMode == 0 & problem.horizontalLayout.lMode == 0,1);
            rows = find(problem.stateLayout.horizontalIndex == zeroK);
            velocityRow = rows(find(problem.stateLayout.component(rows) == "Ap",1));
            a = zeros(height(problem.stateLayout),1);
            a(velocityRow) = 1+0.25i;
            b = audit.energyFactor*a;

            fields = problem.reconstructState(a,outputDomain="spatial");
            Nx = problem.finiteTerrainForms.oversampledSize(1);
            Ny = problem.finiteTerrainForms.oversampledSize(2);
            uHat = interpft(interpft(fields.uHat(:,:,1),Nx,1),Ny,2);
            vHat = interpft(interpft(fields.vHat(:,:,1),Nx,1),Ny,2);
            [x,y] = ndgrid((0:Nx-1)'*wvt.Lx/Nx,(0:Ny-1)'*wvt.Ly/Ny);
            hOversampled = 30*cos(2*pi*x/wvt.Lx);
            gamma = 1-hOversampled/wvt.Lz;
            hX = -30*(2*pi/wvt.Lx)*sin(2*pi*x/wvt.Lx);
            gBottom = (uHat./gamma).*hX+(vHat./gamma).*zeros(size(y));
            phase = exp(1i*(x(:)*problem.horizontalLayout.k.'+y(:)*problem.horizontalLayout.l.'));
            expected = phase'*gBottom(:)/(Nx*Ny);
            actual = audit.bottomTendencyMap*b;
            testCase.verifyEqual(actual,expected,"RelTol",2e-12,"AbsTol",2e-14)

            bottomRow = rows(problem.stateLayout.component(rows) == "etaB");
            a(bottomRow) = 2-0.5i;
            b = audit.energyFactor*a;
            actualBottomValue = audit.bottomValueMap*b;
            expectedBottomValue = zeros(height(problem.horizontalLayout),1);
            expectedBottomValue(zeroK) = a(bottomRow);
            testCase.verifyEqual(actualBottomValue,expectedBottomValue,"AbsTol",2e-13)
        end

        function factorOneOversamplingRetainsCompleteAPVNullspace(testCase)
            wvt = TestWVTerrainEnergyConstrainedClosure.createTransform([4 4 5],false);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny),horizontalOversamplingFactor=1);
            audit = problem.auditConstrainedClosure();
            testCase.verifyTrue(audit.isScientificallyAdmissible)
            testCase.verifyEqual(audit.diagnostics.apvNullity,audit.diagnostics.expectedZeroAPVDimension)
            testCase.verifyEqual(size(audit.apvMap,2),height(problem.stateLayout))
            testCase.verifyGreaterThan(size(audit.apvMap,2),size(audit.apvMap,1))
        end

        function sinusoidalTerrainAuditIsDeterministicAndNonmutating(testCase)
            amplitudes = [1 20];
            for amplitude = amplitudes
                wvt = TestWVTerrainEnergyConstrainedClosure.createTransform([4 4 5],false);
                h = amplitude*cos(2*pi*wvt.x/wvt.Lx).*ones(1,wvt.Ny);
                problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h);
                formsBefore = problem.finiteTerrainForms;
                first = problem.auditConstrainedClosure();
                second = problem.auditConstrainedClosure();
                testCase.verifyEqual(first.status,second.status)
                testCase.verifyEqual(first.constraintSystemIsFeasible,second.constraintSystemIsFeasible)
                testCase.verifyEqual(first.diagnostics.apvNullity,second.diagnostics.apvNullity)
                testCase.verifyEqual(first.diagnostics.minimumConstraintResidual,second.diagnostics.minimumConstraintResidual,"RelTol",1e-13)
                testCase.verifyTrue(all(isfinite(first.diagnostics.apvSingularValues)))
                testCase.verifyTrue(all(isfinite(first.diagnostics.constraintSingularValues)))
                testCase.verifyEqual(problem.finiteTerrainForms,formsBefore)
            end
        end

        function incompatibilityDiagnosticsAreDeterministicAcrossResolution(testCase)
            resolutions = {[4 4 5],[6 4 5]};
            for iResolution = 1:numel(resolutions)
                wvt = TestWVTerrainEnergyConstrainedClosure.createTransform(resolutions{iResolution},false);
                h = 20*cos(2*pi*wvt.x/wvt.Lx).*ones(1,wvt.Ny);
                problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h);
                first = problem.auditConstrainedClosure();
                second = problem.auditConstrainedClosure();
                testCase.verifyEqual(first.status,second.status)
                testCase.verifyEqual(first.diagnostics.apvRank,second.diagnostics.apvRank)
                testCase.verifyEqual(first.diagnostics.apvNullity,second.diagnostics.apvNullity)
                testCase.verifyEqual(first.diagnostics.minimumConstraintResidual,second.diagnostics.minimumConstraintResidual,"RelTol",1e-13)
                testCase.verifyTrue(all(isfinite(first.diagnostics.apvSingularValues)))
                testCase.verifyTrue(all(isfinite(first.diagnostics.constraintSingularValues)))
                testCase.verifyEqual(first.isScientificallyAdmissible,second.isScientificallyAdmissible)
            end
        end
    end

    methods (Static)
        function wvt = createTransform(resolution,shouldAntialias)
            wvt = WVTransformBoussinesq([20e3 16e3 1000],resolution,N2Function=@(z)2e-5+0*z,latitude=30,shouldAntialias=shouldAntialias);
        end

        function verifyRandomStateIdentities(testCase,audit)
            b = complex((1:size(audit.rawGenerator,1))',mod((1:size(audit.rawGenerator,1))',7)-3);
            energyScale = max(1,norm(audit.constrainedGenerator,"fro")*norm(b)^2);
            testCase.verifyLessThanOrEqual(abs(real(b'*audit.constrainedGenerator*b))/energyScale,1e-13)
            testCase.verifyLessThanOrEqual(norm(audit.apvMap*audit.constrainedGenerator*b)/max(1,norm(b)),1e-12)
            testCase.verifyLessThanOrEqual(norm((audit.bottomValueMap*audit.constrainedGenerator-audit.bottomTendencyMap)*b)/max(1,norm(b)),1e-12)
            frequency = eig(1i*audit.constrainedGenerator);
            testCase.verifyLessThanOrEqual(max(abs(imag(frequency))),1e-12*max(1,max(abs(real(frequency)))))
        end
    end
end
