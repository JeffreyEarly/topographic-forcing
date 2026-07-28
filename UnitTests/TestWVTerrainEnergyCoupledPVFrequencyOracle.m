classdef TestWVTerrainEnergyCoupledPVFrequencyOracle < matlab.unittest.TestCase
    % Verify the coupled volume--boundary PV frequency oracle.

    properties
        problem
    end

    methods (TestClassSetup)
        function createProblem(testCase)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            testCase.problem = TestWVTerrainEnergyCoupledPVFrequencyOracle.createProblemWithResolution(9,false,@(z)2e-5+0*z);
        end
    end

    methods (Test)
        function directAndEndpointProblemsAgree(testCase)
            audit = testCase.problem.auditCoupledPVFrequencyOracle(polynomialDegree=20);
            testCase.verifyEqual(audit.status,"passed")
            testCase.verifyLessThan(audit.direct.diagnostics.inversionDefect,1e-12)
            testCase.verifyLessThan(audit.direct.diagnostics.greenIdentityDefect,1e-12)
            testCase.verifyLessThan(audit.direct.diagnostics.energyDefect,1e-12)
            testCase.verifyLessThan(audit.direct.diagnostics.pseudoenstrophyDefect,1e-12)
            testCase.verifyLessThan(audit.direct.diagnostics.frequencyImaginaryDefect,1e-11)
            testCase.verifyLessThan(audit.endpoint.diagnostics.resolvedFrequencyDefect,1e-10)
            testCase.verifyLessThan(audit.endpoint.diagnostics.resolvedEigenfunctionDefect,1e-9)
            testCase.verifyEqual(audit.endpoint.diagnostics.numberOfFiniteModes,audit.endpoint.diagnostics.expectedNumberOfFiniteModes)
            testCase.verifyEqual(audit.direct.pseudoenstrophy.classification,"definite")
            testCase.verifyEqual(audit.direct.pseudoenstrophy.inertia(2:3),[0 0])
        end

        function fPlaneLimitContainsActiveZeroVolumeAPVBoundaryMode(testCase)
            audit = testCase.problem.auditCoupledPVFrequencyOracle(volumePVGradient=[0 0],polynomialDegree=20);
            testCase.verifyEqual(audit.status,"passed")
            testCase.verifyFalse(audit.endpoint.isAvailable)
            testCase.verifyFalse(audit.direct.pseudoenstrophy.isAvailable)
            frequencyTolerance = 1e-11*max(abs(audit.direct.frequency));
            nonzero = abs(audit.direct.frequency) > frequencyTolerance;
            testCase.verifyEqual(nnz(nonzero),1)
            testCase.verifyLessThan(max(audit.direct.qNorm(nonzero)),1e-10)
            testCase.verifyGreaterThan(min(audit.direct.bottomParticipation(nonzero)),0.9)
            testCase.verifyLessThan(audit.direct.diagnostics.energyDefect,1e-12)
        end

        function zeroGradientsProduceStationaryPVState(testCase)
            flatBottom = testCase.problem.auditCoupledPVFrequencyOracle(bottomSlope=[0 0],polynomialDegree=20);
            testCase.verifyEqual(flatBottom.status,"passed")
            testCase.verifyFalse(flatBottom.endpoint.isAvailable)
            testCase.verifyFalse(flatBottom.direct.pseudoenstrophy.isAvailable)
            testCase.verifyEqual(flatBottom.direct.generator(end,:),zeros(1,size(flatBottom.direct.generator,2)),"AbsTol",1e-14)

            audit = testCase.problem.auditCoupledPVFrequencyOracle(volumePVGradient=[0 0],bottomSlope=[0 0],polynomialDegree=8);
            testCase.verifyEqual(audit.status,"passed")
            testCase.verifyEqual(audit.direct.generator,zeros(size(audit.direct.generator)),"AbsTol",1e-14)
            testCase.verifyEqual(audit.direct.frequency,zeros(size(audit.direct.frequency)),"AbsTol",1e-14)
        end

        function gradientSignsAndFourierConjugacyAreRetained(testCase)
            positive = testCase.problem.auditCoupledPVFrequencyOracle(polynomialDegree=20);
            negative = testCase.problem.auditCoupledPVFrequencyOracle(horizontalMode=[-1 0],polynomialDegree=20);
            signed = testCase.problem.auditCoupledPVFrequencyOracle(volumePVGradient=[0 -2e-11],polynomialDegree=20);
            testCase.verifyEqual(positive.status,"passed")
            testCase.verifyEqual(negative.status,"passed")
            testCase.verifyEqual(signed.status,"passed")
            testCase.verifyEqual(negative.direct.frequency,-flip(positive.direct.frequency),"RelTol",1e-10,"AbsTol",1e-15)
            testCase.verifyEqual(signed.direct.pseudoenstrophy.classification,"pontryagin")
            testCase.verifyGreaterThan(signed.direct.pseudoenstrophy.inertia(1),0)
            testCase.verifyGreaterThan(signed.direct.pseudoenstrophy.inertia(2),0)
            testCase.verifyLessThan(signed.direct.diagnostics.pseudoenstrophyDefect,1e-12)
        end

        function balancedBasisRepresentationConverges(testCase)
            resolution = [5 7 9];
            defect = zeros(size(resolution));
            for iResolution = 1:numel(resolution)
                currentProblem = TestWVTerrainEnergyCoupledPVFrequencyOracle.createProblemWithResolution(resolution(iResolution),false,@(z)2e-5+0*z);
                audit = currentProblem.auditCoupledPVFrequencyOracle(polynomialDegree=20);
                testCase.verifyEqual(audit.status,"passed")
                testCase.verifyEqual(audit.basis.diagnostics.rank,audit.basis.diagnostics.expectedRank)
                testCase.verifyLessThan(audit.basis.diagnostics.maximumBottomValueDefect,1e-13)
                defect(iResolution) = audit.basis.diagnostics.resolvedReconstructionDefect;
            end
            testCase.verifyTrue(all(diff(defect) < 0))
            testCase.verifyLessThan(defect(end),0.02)
        end

        function arbitraryStratificationAndAntialiasingRetainOracle(testCase)
            N2 = @(z)1.2e-5*exp(z/1800);
            currentProblem = TestWVTerrainEnergyCoupledPVFrequencyOracle.createProblemWithResolution(11,true,N2);
            audit = currentProblem.auditCoupledPVFrequencyOracle(polynomialDegree=24);
            testCase.verifyEqual(audit.status,"passed")
            testCase.verifyLessThan(audit.endpoint.diagnostics.resolvedFrequencyDefect,1e-10)
            testCase.verifyLessThan(audit.endpoint.diagnostics.resolvedEigenfunctionDefect,1e-9)
            testCase.verifyLessThan(audit.direct.diagnostics.energyDefect,1e-12)
        end

        function auditRejectsInvalidInputsWithoutChangingLayout(testCase)
            horizontalBefore = testCase.problem.horizontalLayout;
            stateBefore = testCase.problem.stateLayout;
            testCase.problem.auditCoupledPVFrequencyOracle(polynomialDegree=20);
            testCase.verifyEqual(testCase.problem.horizontalLayout,horizontalBefore)
            testCase.verifyEqual(testCase.problem.stateLayout,stateBefore)
            testCase.verifyError(@()testCase.problem.auditCoupledPVFrequencyOracle(horizontalMode=[0 0]), ...
                "WVTerrainEnergyGalerkin:InvalidCoupledPVHorizontalMode")
            testCase.verifyError(@()testCase.problem.auditCoupledPVFrequencyOracle(horizontalMode=[7 0]), ...
                "WVTerrainEnergyGalerkin:UnknownHorizontalMode")
            testCase.verifyError(@()testCase.problem.auditCoupledPVFrequencyOracle(volumePVGradient=[0 NaN]), ...
                "WVTerrainEnergyGalerkin:InvalidVolumePVGradient")
            testCase.verifyError(@()testCase.problem.auditCoupledPVFrequencyOracle(bottomSlope=[Inf 0]), ...
                "WVTerrainEnergyGalerkin:InvalidBottomSlope")
            testCase.verifyError(@()testCase.problem.auditCoupledPVFrequencyOracle(polynomialDegree=2), ...
                "WVTerrainEnergyGalerkin:InvalidPolynomialDegree")
            testCase.verifyError(@()testCase.problem.auditCoupledPVFrequencyOracle(polynomialDegree=8,quadratureOrder=8), ...
                "WVTerrainEnergyGalerkin:InsufficientQuadratureOrder")
        end
    end

    methods (Static,Access=private)
        function currentProblem = createProblemWithResolution(Nz,shouldAntialias,N2)
            z = linspace(-1200,0,Nz)';
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[4 4 Nz], ...
                N2Function=N2,latitude=45,shouldAntialias=shouldAntialias,z=z);
            currentProblem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(4,4));
        end
    end
end
