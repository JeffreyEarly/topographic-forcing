classdef TestWVTerrainEnergyBranchPDiscreteOracle < matlab.unittest.TestCase
    % Verify the independent Branch-P polynomial oracle and blocker diagnosis.

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
        function polynomialSpacesAndBottomCoordinateAreIndependent(testCase)
            audit = testCase.problem.auditBranchPDiscreteOracle( ...
                bottomSlope=[0.01 0],polynomialDegree=6);
            geometry = audit.polynomialOracle.geometryDiagnostics;
            testCase.verifyLessThanOrEqual(geometry.endpointResidual,1e-13)
            testCase.verifyLessThanOrEqual(geometry.greenIdentityDefect,1e-13)
            for iK = 1:numel(audit.polynomialOracle.blocks)
                block = audit.polynomialOracle.blocks{iK};
                iBottom = block.bottomCoordinateIndex;
                raw = block.rawReconstruction;
                testCase.verifyEqual(raw.uHat(:,iBottom),zeros(size(raw.uHat,1),1),AbsTol=1e-14)
                testCase.verifyEqual(raw.vHat(:,iBottom),zeros(size(raw.vHat,1),1),AbsTol=1e-14)
                testCase.verifyEqual(raw.wHat(:,iBottom),zeros(size(raw.wHat,1),1),AbsTol=1e-14)
                testCase.verifyEqual(raw.etaHat(:,iBottom), ...
                    audit.polynomialOracle.spaces.chi,AbsTol=1e-14)
            end
        end

        function branchPOracleEstablishesDiscreteAPVBlocker(testCase)
            audit = testCase.problem.auditBranchPDiscreteOracle( ...
                bottomSlope=[0.01 0],polynomialDegree=6);
            diagnostics = audit.polynomialOracle.maximum.required;
            testCase.verifyEqual(audit.status,"scientific-blocker")
            testCase.verifyFalse(audit.isCompatible)
            testCase.verifyEqual(audit.repairCandidate,"none")
            testCase.verifyLessThanOrEqual(diagnostics.weakEvolutionDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.energyDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.bottomDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.pressureWorkDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.constraintTangencyDefect,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.saddleResidual,1e-12)
            testCase.verifyLessThanOrEqual(diagnostics.analyticDerivativeDefect,2e-10)
            testCase.verifyGreaterThan(diagnostics.apvDefect,0.5)
            testCase.verifyGreaterThan(diagnostics.enstrophyDefect,0.3)
            flat = audit.polynomialOracle.maximum.flat;
            testCase.verifyLessThanOrEqual(flat.weakEvolutionDefect,1e-12)
            testCase.verifyLessThanOrEqual(flat.energyDefect,1e-12)
            testCase.verifyLessThanOrEqual(flat.apvDefect,1e-12)
            testCase.verifyLessThanOrEqual(flat.enstrophyDefect,1e-12)
            testCase.verifyLessThanOrEqual(flat.bottomDefect,2e-12)
            testCase.verifyGreaterThan(audit.comparison.maximum.nativeRequired.energySkewDefect,1e-4)
            testCase.verifyGreaterThan(audit.comparison.maximum.nativeRequired.apvTendencyDefect,0.5)
        end

        function blockerPersistsUnderPolynomialRefinement(testCase)
            degree = [4 6 8 10];
            apv = zeros(size(degree));
            enstrophy = zeros(size(degree));
            for i = 1:numel(degree)
                audit = testCase.problem.auditBranchPDiscreteOracle( ...
                    bottomSlope=[0.01 0],polynomialDegree=degree(i));
                diagnostics = audit.polynomialOracle.maximum.required;
                apv(i) = diagnostics.apvDefect;
                enstrophy(i) = diagnostics.enstrophyDefect;
                testCase.verifyLessThanOrEqual(diagnostics.energyDefect,2e-12)
                testCase.verifyLessThanOrEqual(diagnostics.bottomDefect,1e-12)
            end
            testCase.verifyGreaterThan(min(apv),0.5)
            testCase.verifyGreaterThan(min(enstrophy),0.3)
            audit8 = testCase.problem.auditBranchPDiscreteOracle( ...
                bottomSlope=[0.01 0],polynomialDegree=8);
            testCase.verifyLessThanOrEqual( ...
                audit8.polynomialOracle.maximum.flatModeOneDispersionDefect,1e-11)
        end

        function doubledQuadratureDoesNotChangeTheOracle(testCase)
            reference = testCase.problem.auditBranchPDiscreteOracle( ...
                bottomSlope=[0.01 0],polynomialDegree=6);
            doubled = testCase.problem.auditBranchPDiscreteOracle( ...
                bottomSlope=[0.01 0],polynomialDegree=6, ...
                quadratureOrder=2*reference.quadratureOrder);
            referenceDiagnostics = reference.polynomialOracle.maximum.required;
            doubledDiagnostics = doubled.polynomialOracle.maximum.required;
            testCase.verifyLessThanOrEqual(abs( ...
                referenceDiagnostics.apvDefect-doubledDiagnostics.apvDefect),5e-4)
            testCase.verifyLessThanOrEqual(abs( ...
                referenceDiagnostics.enstrophyDefect-doubledDiagnostics.enstrophyDefect),1e-12)
            testCase.verifyLessThanOrEqual(abs( ...
                reference.polynomialOracle.maximum.flatModeOneDispersionDefect ...
                -doubled.polynomialOracle.maximum.flatModeOneDispersionDefect),1e-12)
            testCase.verifyLessThanOrEqual( ...
                doubled.polynomialOracle.geometryDiagnostics.greenIdentityDefect,1e-13)
        end

        function slopeDirectionsPreserveEnergyAndBottomButNotAPV(testCase)
            slopes = [0.01 0;0 0.01;0.006 0.008];
            for slope = slopes'
                audit = testCase.problem.auditBranchPDiscreteOracle( ...
                    bottomSlope=slope',polynomialDegree=4);
                diagnostics = audit.polynomialOracle.maximum.required;
                testCase.verifyLessThanOrEqual(diagnostics.energyDefect,2e-12)
                testCase.verifyLessThanOrEqual(diagnostics.bottomDefect,2e-12)
                testCase.verifyGreaterThan(diagnostics.apvDefect,0.4)
            end
        end

        function finiteSlopeRemainderIsSecondOrder(testCase)
            amplitudes = [0.02 0.01 0.005];
            relativeRemainder = zeros(size(amplitudes));
            for i = 1:numel(amplitudes)
                audit = testCase.problem.auditBranchPDiscreteOracle( ...
                    bottomSlope=[amplitudes(i) 0],polynomialDegree=4);
                relativeRemainder(i) = audit.polynomialOracle.maximum.finiteSlopeLinearizationDefect;
            end
            observedOrder = log2(relativeRemainder(1:end-1)./relativeRemainder(2:end));
            testCase.verifyGreaterThan(min(observedOrder),0.95)
        end

        function auditDoesNotChangePublicLayout(testCase)
            horizontalBefore = testCase.problem.horizontalLayout;
            stateBefore = testCase.problem.stateLayout;
            conjugateBefore = testCase.problem.conjugateCoordinateIndex;
            testCase.problem.auditBranchPDiscreteOracle( ...
                bottomSlope=[0.01 0],polynomialDegree=4);
            testCase.verifyEqual(testCase.problem.horizontalLayout,horizontalBefore)
            testCase.verifyEqual(testCase.problem.stateLayout,stateBefore)
            testCase.verifyEqual(testCase.problem.conjugateCoordinateIndex,conjugateBefore)
        end

        function rejectsInvalidOptions(testCase)
            testCase.verifyError(@()testCase.problem.auditBranchPDiscreteOracle( ...
                bottomSlope=[NaN 0]),"WVTerrainEnergyGalerkin:InvalidBottomSlope")
            testCase.verifyError(@()testCase.problem.auditBranchPDiscreteOracle( ...
                polynomialDegree=1),"WVTerrainEnergyGalerkin:InvalidPolynomialDegree")
            testCase.verifyError(@()testCase.problem.auditBranchPDiscreteOracle( ...
                polynomialDegree=6,quadratureOrder=4), ...
                "WVTerrainEnergyGalerkin:InsufficientQuadratureOrder")
        end
    end
end
