classdef TestWVTerrainEnergyGlobalFirstOrderTerrainDressing < matlab.unittest.TestCase
    % Verify the Milestone-9.4 global first-order dressing oracle.

    properties
        problem
        audit
        horizontalLayout
        stateLayout
    end

    methods (TestClassSetup)
        function createGlobalDressingAudit(testCase)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            N0 = sqrt(2e-5);
            z = linspace(-1200,0,5)';
            wvt = WVTransformBoussinesq( ...
                [24e3 20e3 1200],[6 10 5], ...
                N2Function=@(z)N0^2+0*z,latitude=45, ...
                shouldAntialias=false,z=z);
            [~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
                (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
            terrain = 2.5*cos(2*pi*y/wvt.Ly);
            testCase.problem = WVTerrainEnergyGalerkin.fromTopography( ...
                wvt,topographicHeight=terrain);
            testCase.horizontalLayout = testCase.problem.horizontalLayout;
            testCase.stateLayout = testCase.problem.stateLayout;
            testCase.audit = testCase.problem.auditGlobalFirstOrderTerrainDressing;
        end
    end

    methods (Test)
        function outcomeUsesDeclaredClassification(testCase)
            audit = testCase.audit;
            testCase.verifyTrue(ismember(audit.classification, ...
                ["global-dressing-acceleration", ...
                "global-dressing-equivalent", ...
                "global-dressing-seed", ...
                "global-dressing-incompatible"]))
            testCase.verifyEqual(audit.status,audit.classification)
            testCase.verifyEqual(audit.nextScope, ...
                "milestone-10-only-if-separately-authorized")
        end

        function analyticTerrainDerivativesAgreeWithCenteredDifferences(testCase)
            results = [testCase.audit.details(:);testCase.audit.comparison(:)];
            defect = cellfun(@(value)value.tangentAgreementDefect,results);
            testCase.verifyLessThan(max(defect),1e-9)
        end

        function globalBlocksAreCompleteAndFourierCompatible(testCase)
            results = [testCase.audit.details(:);testCase.audit.comparison(:)];
            orthogonality = cellfun(@(value)value.flatEnergyOrthogonalityDefect,results);
            conjugacy = cellfun(@(value)value.conjugacyDefect,results);
            selection = cellfun(@(value)value.fourierSelectionDefect,results);
            testCase.verifyLessThan(max(orthogonality),1e-11)
            testCase.verifyLessThan(max(conjugacy),1e-11)
            testCase.verifyLessThan(max(selection),1e-11)
            testCase.verifyGreaterThan(testCase.audit.primary.zeroFrequency.numberOfStationaryDirections,0)
            testCase.verifyEqual( ...
                testCase.audit.primary.zeroFrequency.numberOfActiveDirections ...
                +testCase.audit.primary.zeroFrequency.numberOfStationaryDirections, ...
                numel(testCase.audit.primary.zeroFrequency.firstOrderFrequency))
        end

        function firstOrderCorrectionClosesWeakAndBottomEquations(testCase)
            results = [testCase.audit.details(:);testCase.audit.comparison(:)];
            correction = cellfun(@(value)value.correctionResidual,results);
            bottom = cellfun(@(value)value.firstOrderBottomDefect,results);
            testCase.verifyLessThan(max(correction),1e-10)
            testCase.verifyLessThan(max(bottom),1e-10)
            testCase.verifyLessThan( ...
                testCase.audit.weakSequenceDiagnostics.maximumGreenIdentityDefect,1e-10)
            testCase.verifyLessThan( ...
                testCase.audit.weakSequenceDiagnostics.maximumStationaryRowDefect,1e-10)
        end

        function dressedResidualsHaveSecondOrderScaling(testCase)
            scaling = testCase.audit.primary.scaling;
            testCase.verifyGreaterThanOrEqual(scaling.minimumDressedWeakOrder,1.8)
            testCase.verifyGreaterThanOrEqual(scaling.minimumDressedBottomOrder,1.8)
            testCase.verifyGreaterThanOrEqual(scaling.weakImprovement,4)
            testCase.verifyGreaterThanOrEqual(scaling.bottomImprovement,4)
        end

        function exactReducedFormsPreservePhysicalEnergyStructure(testCase)
            results = [testCase.audit.details(:);testCase.audit.comparison(:)];
            energy = cellfun(@(value)value.reducedEnergyHermitianDefect,results);
            exchange = cellfun(@(value)value.reducedExchangeSkewHermitianDefect,results);
            testCase.verifyLessThan(max(energy),1e-12)
            testCase.verifyLessThan(max(exchange),1e-12)
            testCase.verifyEqual(testCase.audit.primary.finite.decomposition.completenessDefect,0,AbsTol=1e-10)
            testCase.verifyLessThan(testCase.audit.primary.finite.decomposition.maximumOrthogonalityDefect,1e-10)
            testCase.verifyGreaterThanOrEqual(testCase.audit.primary.finite.internalProjectorImprovement,2)
        end

        function publicLayoutsRemainUnchanged(testCase)
            testCase.verifyEqual(testCase.problem.horizontalLayout,testCase.horizontalLayout)
            testCase.verifyEqual(testCase.problem.stateLayout,testCase.stateLayout)
        end

        function validationRejectsIncompleteInputs(testCase)
            testCase.verifyError(@() ...
                testCase.problem.auditGlobalFirstOrderTerrainDressing( ...
                supportModeBounds=[1 2;1 3], ...
                primitivePolynomialDegrees=[4;6]), ...
                "WVTerrainEnergyGalerkin:InvalidGlobalDressingDegrees")
            testCase.verifyError(@() ...
                testCase.problem.auditGlobalFirstOrderTerrainDressing( ...
                paddingFactors=2), ...
                "WVTerrainEnergyGalerkin:InvalidGlobalDressingPadding")
            testCase.verifyError(@() ...
                testCase.problem.auditGlobalFirstOrderTerrainDressing( ...
                terrainScales=[0.25;0.5;0.75]), ...
                "WVTerrainEnergyGalerkin:InvalidGlobalDressingTerrainScales")
        end
    end
end
