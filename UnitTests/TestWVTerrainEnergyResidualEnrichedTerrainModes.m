classdef TestWVTerrainEnergyResidualEnrichedTerrainModes < matlab.unittest.TestCase
    % Verify the Milestone-10 exact-residual enrichment oracle.

    properties
        problem
        audit
        horizontalLayout
        stateLayout
    end

    methods (TestClassSetup)
        function createResidualEnrichmentAudit(testCase)
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
            testCase.audit = testCase.problem.auditResidualEnrichedTerrainModes;
        end
    end

    methods (Test)
        function outcomeIsResidualEnrichmentAcceleration(testCase)
            audit = testCase.audit;
            testCase.verifyEqual(audit.classification, ...
                "residual-enrichment-acceleration")
            testCase.verifyTrue(audit.isCompatible)
            testCase.verifyEqual(audit.nextScope, ...
                "milestone-11-only-if-separately-authorized")
            testCase.verifyGreaterThanOrEqual( ...
                audit.primary.dimension.compressionFactor,2)
        end

        function exactResidualIterationConverges(testCase)
            history = testCase.audit.primary.history;
            testCase.verifyGreaterThan(numel(history),1)
            testCase.verifyLessThanOrEqual( ...
                testCase.audit.primary.numberOfIterations, ...
                testCase.audit.maximumIterations)
            testCase.verifyGreaterThan( ...
                history(1).maximumInternalResidual, ...
                history(end).maximumInternalResidual)
            testCase.verifyLessThan( ...
                history(end).maximumInternalResidual,1e-10)
            testCase.verifyTrue(history(end).passesMandatoryGates)
        end

        function enrichedBlockMatchesIndependentDenseOracle(testCase)
            primary = testCase.audit.primary;
            testCase.verifyLessThan(primary.internalProjectorDefect,1e-9)
            testCase.verifyLessThan(primary.internalFrequencyDefect,1e-9)
            testCase.verifyLessThan(primary.maximumInternalAPVDefect,1e-8)
            testCase.verifyLessThan(primary.maximumInternalBottomDefect,1e-8)
            testCase.verifyLessThan(primary.maximumInternalStrongResidual,1e-5)
            testCase.verifyLessThan(testCase.audit.paddingDefect,1e-8)
            testCase.verifyLessThan(testCase.audit.nestedExcessDefect,1e-8)
        end

        function physicalEnergyAndStationarySpaceRemainExact(testCase)
            results = [testCase.audit.details(:);testCase.audit.comparison(:)];
            energy = cellfun(@(value)value.energyHermitianDefect,results);
            exchange = cellfun(@(value)value.exchangeSkewHermitianDefect,results);
            stationary = cellfun(@(value)value.stationaryRowDefect,results);
            testCase.verifyLessThan(max(energy),1e-12)
            testCase.verifyLessThan(max(exchange),1e-12)
            testCase.verifyLessThan(max(stationary),1e-10)
            testCase.verifyLessThan( ...
                testCase.audit.primary.stationaryOrthogonalityDefect,1e-10)
        end

        function physicalAndRemainderProjectorsAreComplete(testCase)
            decomposition = testCase.audit.primary.decomposition;
            testCase.verifyLessThan( ...
                decomposition.maximumOrthogonalityDefect,1e-10)
            testCase.verifyLessThan(decomposition.completenessDefect,1e-10)
            testCase.verifyLessThan(decomposition.conjugacyDefect,1e-10)
            dimension = testCase.audit.primary.dimension;
            testCase.verifyEqual(dimension.stationary ...
                +dimension.validatedInternal+dimension.unresolved, ...
                dimension.ambient)
        end

        function activeBottomDirectionsRemainUnresolved(testCase)
            final = testCase.audit.primary;
            testCase.verifyGreaterThan( ...
                numel(final.zeroCandidateModeIndices),0)
            testCase.verifyEqual(final.zeroCandidateStatus,"unresolved")
            testCase.verifyGreaterThanOrEqual( ...
                final.minimumZeroCandidateFrequencyResolutionRatio,0)
            testCase.verifyGreaterThanOrEqual( ...
                final.maximumZeroCandidateBottomParticipation,0)
        end

        function globalDressingSeedRemainsVerified(testCase)
            seed = testCase.audit.seed;
            testCase.verifyEqual(seed.classification,"global-dressing-seed")
            testCase.verifyGreaterThanOrEqual( ...
                seed.primary.scaling.minimumDressedWeakOrder,1.8)
            testCase.verifyGreaterThanOrEqual( ...
                seed.primary.scaling.minimumDressedBottomOrder,1.8)
        end

        function publicLayoutsRemainUnchanged(testCase)
            testCase.verifyEqual(testCase.problem.horizontalLayout, ...
                testCase.horizontalLayout)
            testCase.verifyEqual(testCase.problem.stateLayout, ...
                testCase.stateLayout)
        end

        function validationRejectsIncompleteInputs(testCase)
            testCase.verifyError(@() ...
                testCase.problem.auditResidualEnrichedTerrainModes( ...
                supportModeBounds=[1 2;1 3], ...
                primitivePolynomialDegrees=[4;6]), ...
                "WVTerrainEnergyGalerkin:InvalidResidualEnrichmentDegrees")
            testCase.verifyError(@() ...
                testCase.problem.auditResidualEnrichedTerrainModes( ...
                paddingFactors=2), ...
                "WVTerrainEnergyGalerkin:InvalidResidualEnrichmentPadding")
            testCase.verifyError(@() ...
                testCase.problem.auditResidualEnrichedTerrainModes( ...
                terrainScales=[0.25;0.5;0.75]), ...
                "WVTerrainEnergyGalerkin:InvalidResidualEnrichmentTerrainScales")
        end
    end
end
