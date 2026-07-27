classdef TestWVTerrainEnergyQuadraticInvariantClosure < matlab.unittest.TestCase
    % Verify the quadratic energy-enstrophy closure audit.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test)
        function flatAndUniformDepthReturnTheRawGenerator(testCase)
            for shouldAntialias = [false true]
                wvt = TestWVTerrainEnergyQuadraticInvariantClosure.createTransform([4 4 5],shouldAntialias);
                terrains = {zeros(wvt.Nx,wvt.Ny),100*ones(wvt.Nx,wvt.Ny)};
                for iTerrain = 1:numel(terrains)
                    problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=terrains{iTerrain});
                    audit = problem.auditQuadraticInvariantClosure();
                    testCase.verifyEqual(audit.status,"admissible")
                    testCase.verifyTrue(audit.constraintSystemIsFeasible)
                    testCase.verifyTrue(audit.isScientificallyAdmissible)
                    testCase.verifyEmpty(audit.failedConditions)
                    testCase.verifyTrue(audit.diagnostics.rawSatisfiesConstraints)
                    testCase.verifyEqual(audit.diagnostics.correctionNorm,0,"AbsTol",1e-15)
                    testCase.verifyEqual(audit.constrainedGenerator,audit.rawGenerator,"AbsTol",1e-15)
                    testCase.verifyLessThanOrEqual(audit.diagnostics.energyConstraintResidual,1e-13)
                    testCase.verifyLessThanOrEqual(audit.diagnostics.enstrophyConstraintResidual,1e-12)
                    testCase.verifyLessThanOrEqual(audit.diagnostics.bottomConstraintResidual,1e-12)
                    testCase.verifyLessThanOrEqual(audit.diagnostics.constrainedConjugacyDefect,1e-12)
                    testCase.verifyLessThanOrEqual(audit.diagnostics.maximumFrequencyImaginaryPart,1e-12)
                    TestWVTerrainEnergyQuadraticInvariantClosure.verifyStructuralIdentities(testCase,problem,audit)
                    TestWVTerrainEnergyQuadraticInvariantClosure.verifyRandomStateIdentities(testCase,problem,audit)
                end
            end
        end

        function sinusoidalTerrainAuditIsDeterministicAndNonmutating(testCase)
            wvt = TestWVTerrainEnergyQuadraticInvariantClosure.createTransform([4 4 5],false);
            for amplitude = [1 20]
                h = amplitude*cos(2*pi*wvt.x/wvt.Lx).*ones(1,wvt.Ny);
                problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h);
                formsBefore = problem.finiteTerrainForms;
                first = problem.auditQuadraticInvariantClosure();
                second = problem.auditQuadraticInvariantClosure();
                testCase.verifyEqual(first.status,second.status)
                testCase.verifyEqual(first.constraintSystemIsFeasible,second.constraintSystemIsFeasible)
                testCase.verifyEqual(first.diagnostics.minimumConstraintResidual,second.diagnostics.minimumConstraintResidual,"RelTol",1e-13)
                testCase.verifyEqual(first.diagnostics.eigenvalues,second.diagnostics.eigenvalues,"RelTol",1e-13,"AbsTol",1e-30)
                testCase.verifyLessThanOrEqual(max([first.diagnostics.block.rawCoordinateRoundTripDefect]),1e-12)
                testCase.verifyLessThanOrEqual(max([first.diagnostics.block.normalEquationRelativeResidual]),1e-12)
                testCase.verifyEqual(problem.finiteTerrainForms,formsBefore)
                TestWVTerrainEnergyQuadraticInvariantClosure.verifyStructuralIdentities(testCase,problem,first)
            end
        end

        function historicalAuditResidualsRemainFinite(testCase)
            wvt = TestWVTerrainEnergyQuadraticInvariantClosure.createTransform([4 4 5],false);
            amplitudes = [1 5 20];
            residuals = zeros(size(amplitudes));
            for iAmplitude = 1:numel(amplitudes)
                h = amplitudes(iAmplitude)*cos(2*pi*wvt.x/wvt.Lx).*ones(1,wvt.Ny);
                problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h);
                audit = problem.auditQuadraticInvariantClosure();
                residuals(iAmplitude) = audit.diagnostics.minimumConstraintResidual;
            end
            testCase.verifyTrue(all(isfinite(residuals)))
            testCase.verifyTrue(all(residuals >= 0))
        end

        function auditIsDeterministicAcrossDiscretizations(testCase)
            resolutions = {[4 4 5],[4 4 5],[4 4 7],[6 4 5]};
            oversamplingFactors = [1 3 2 2];
            for iCase = 1:numel(resolutions)
                wvt = TestWVTerrainEnergyQuadraticInvariantClosure.createTransform(resolutions{iCase},false);
                h = 20*cos(2*pi*wvt.x/wvt.Lx).*ones(1,wvt.Ny);
                problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h,horizontalOversamplingFactor=oversamplingFactors(iCase));
                first = problem.auditQuadraticInvariantClosure();
                second = problem.auditQuadraticInvariantClosure();
                testCase.verifyEqual(second.status,first.status)
                testCase.verifyEqual(second.diagnostics.minimumConstraintResidual,first.diagnostics.minimumConstraintResidual,"RelTol",1e-13)
                testCase.verifyEqual(second.diagnostics.eigenvalues,first.diagnostics.eigenvalues,"RelTol",1e-13,"AbsTol",1e-30)
                testCase.verifyTrue(isfinite(first.diagnostics.minimumConstraintRelativeResidual))
                testCase.verifyLessThanOrEqual(first.diagnostics.maximumRealificationDefect,1e-12)
            end
        end
    end

    methods (Static)
        function wvt = createTransform(resolution,shouldAntialias)
            wvt = WVTransformBoussinesq([20e3 16e3 1000],resolution,N2Function=@(z)2e-5+0*z,latitude=30,shouldAntialias=shouldAntialias);
        end

        function verifyStructuralIdentities(testCase,problem,audit)
            n = height(problem.stateLayout);
            testCase.verifyEqual(audit.stateRealification'*audit.stateRealification,eye(n),"AbsTol",2e-13)
            energyMatrix = problem.finiteTerrainForms.energyMatrix;
            energyFactorDefect = norm(audit.energyFactor'*audit.energyFactor-energyMatrix,"fro")/norm(energyMatrix,"fro");
            testCase.verifyLessThanOrEqual(energyFactorDefect,2e-13)
            expectedZ = problem.finiteTerrainForms.apvMatrix'*problem.finiteTerrainForms.apvMatrix;
            testCase.verifyEqual(audit.quadraticEnstrophyMatrix,expectedZ)
            quadraticFormDefect = norm(audit.energyFactor'*audit.enstrophyMetric*audit.energyFactor-expectedZ,"fro")/norm(expectedZ,"fro");
            testCase.verifyLessThanOrEqual(quadraticFormDefect,2e-13)
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximumRealificationDefect,1e-12)
            testCase.verifyLessThanOrEqual(audit.diagnostics.bottom.phaseGramDefect,1e-12)
        end

        function verifyRandomStateIdentities(testCase,problem,audit)
            n = size(audit.constrainedGenerator,1);
            b = complex((1:n)',mod((1:n)',7)-3);
            K = audit.constrainedGenerator;
            G = audit.enstrophyMetric;
            energyScale = max(norm(K,"fro")*norm(b)^2,realmin);
            enstrophyScale = max(norm(G,"fro")*norm(K,"fro")*norm(b)^2,realmin);
            testCase.verifyLessThanOrEqual(abs(real(b'*K*b))/energyScale,1e-13)
            testCase.verifyLessThanOrEqual(abs(real(b'*G*K*b))/enstrophyScale,1e-12)
            bottomScale = max((norm(audit.bottomValueMap,"fro")*norm(K,"fro")+norm(audit.bottomTendencyMap,"fro"))*norm(b),realmin);
            testCase.verifyLessThanOrEqual(norm((audit.bottomValueMap*K-audit.bottomTendencyMap)*b)/bottomScale,1e-12)

            realState = (1:n)';
            a = audit.energyFactor\realState;
            conjugateState = conj(a(problem.conjugateCoordinateIndex));
            testCase.verifyEqual(a,conjugateState,"RelTol",2e-13,"AbsTol",2e-13)
            q = problem.finiteTerrainForms.apvMatrix*a;
            directZ = real(q'*q)/2;
            matrixZ = real(a'*audit.quadraticEnstrophyMatrix*a)/2;
            energyCoordinateZ = real(realState'*G*realState)/2;
            testCase.verifyEqual(matrixZ,directZ,"RelTol",2e-13,"AbsTol",2e-20)
            testCase.verifyEqual(energyCoordinateZ,directZ,"RelTol",2e-13,"AbsTol",2e-20)
        end
    end
end
