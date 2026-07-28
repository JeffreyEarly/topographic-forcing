classdef TestWVTerrainEnergyBoundaryGreenIdentity < matlab.unittest.TestCase
    % Verify the local constant-slope Green identity and branch audit.

    properties
        problem
    end

    methods (TestClassSetup)
        function createProblem(testCase)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            N2 = @(z) 2e-5+0*z;
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[4 4 5], ...
                N2Function=N2,latitude=45,shouldAntialias=false);
            testCase.problem = WVTerrainEnergyGalerkin.fromTopography( ...
                wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
        end
    end

    methods (Test)
        function flatLimitIsPhysicalEnergyBranch(testCase)
            audit = testCase.problem.auditBoundaryGreenIdentity(bottomSlope=[0 0]);
            testCase.verifyTrue(audit.isCompatible)
            testCase.verifyEqual(audit.branch,"P")
            testCase.verifyEqual(audit.continuumBranch,"P")
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.energySkewDefect,1e-11)
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.apvTendencyDefect,1e-11)
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.potentialEnstrophyDefect,1e-11)
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.bottomEvolutionDefect,1e-11)
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.pressureGreenDefect,1e-11)
            testCase.verifyEqual(audit.diagnostics.activeBoundaryModeCount,0)
            testCase.verifyGreaterThan(audit.diagnostics.minimumEnergyEigenvalue,0)
        end

        function finiteSlopeRecordsIncompatibleDiscreteExit(testCase)
            audit = testCase.problem.auditBoundaryGreenIdentity(bottomSlope=[0.01 0]);
            testCase.verifyFalse(audit.isCompatible)
            testCase.verifyEqual(audit.branch,"incompatible")
            testCase.verifyEqual(audit.continuumBranch,"P")
            testCase.verifyTrue(all(ismember( ...
                ["energySkewDefect","apvTendencyDefect","potentialEnstrophyDefect"], ...
                audit.failedConditions)))
            testCase.verifyFalse(ismember("bottomEvolutionDefect",audit.failedConditions))
            testCase.verifyFalse(ismember("pressureGreenDefect",audit.failedConditions))
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.bottomEvolutionDefect,1e-11)
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.pressureGreenDefect,1e-11)
            testCase.verifyGreaterThan(audit.diagnostics.activeBoundaryModeCount,0)
            testCase.verifyGreaterThan(audit.diagnostics.slopeBottomValueQuadraticDefect,1e-3)
            testCase.verifyLessThan(audit.diagnostics.flatBottomValueQuadraticDefect,1e-10)
            testCase.verifyLessThan(audit.diagnostics.tangentConvergence,1e-7)
            testCase.verifyEqual(audit.diagnostics.generalizedBoundaryMetricRank,0)
            testCase.verifyEqual(audit.diagnostics.generalizedBoundaryMetricDefiniteness,"absent")
        end

        function rejectsInvalidSlope(testCase)
            testCase.verifyError(@()testCase.problem.auditBoundaryGreenIdentity( ...
                bottomSlope=[NaN 0]),"WVTerrainEnergyGalerkin:InvalidBottomSlope")
        end
    end
end
