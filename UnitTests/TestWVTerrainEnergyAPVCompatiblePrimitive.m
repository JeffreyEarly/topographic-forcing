classdef TestWVTerrainEnergyAPVCompatiblePrimitive < matlab.unittest.TestCase
    % Verify the local APV-compatible primitive investigation.

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
        function crossSlopeModeEstablishesAnalyticBlocker(testCase)
            audit = testCase.problem.auditAPVCompatiblePrimitive();
            testCase.verifyEqual(audit.status,"mathematical-blocker")
            testCase.verifyFalse(audit.isCompatible)
            testCase.verifyTrue(audit.analyticObstruction.hasCrossSlopeObstruction)
            testCase.verifyEqual(audit.analyticObstruction.sourceFormula, ...
                "q_t=i(k s_y-l s_x)p/(rho_0 D)")
            testCase.verifyLessThan(audit.analyticObstruction.algebraicIdentityDefect,1e-12)
            testCase.verifyLessThan(audit.analyticObstruction.sourceComparisonDefect,0.02)
            testCase.verifyGreaterThan(audit.energyWeak.diagnostics.apvDefect,1e-7)
            testCase.verifyLessThan(audit.energyWeak.diagnostics.energyDefect,1e-12)
            testCase.verifyLessThan(audit.energyWeak.diagnostics.bottomDefect,1e-11)
            testCase.verifyLessThan(audit.energyWeak.diagnostics.weakEvolutionDefect,1e-12)
        end

        function alignedSpecialCaseClosesAllIdentities(testCase)
            audit = testCase.problem.auditAPVCompatiblePrimitive(bottomSlope=[0.01 0]);
            testCase.verifyEqual(audit.status,"compatible-aligned-special-case")
            testCase.verifyTrue(audit.isCompatible)
            testCase.verifyFalse(audit.analyticObstruction.hasCrossSlopeObstruction)
            testCase.verifyEqual(audit.analyticObstruction.sourceComparisonDefect,0)
            testCase.verifyLessThan(audit.energyWeak.diagnostics.energyDefect,1e-12)
            testCase.verifyLessThan(audit.energyWeak.diagnostics.apvDefect,1e-12)
            testCase.verifyLessThan(audit.energyWeak.diagnostics.enstrophyDefect,1e-12)
            testCase.verifyLessThan(audit.energyWeak.diagnostics.bottomDefect,1e-11)
        end

        function vorticityAndExplicitAPVTradeEnergyForAPV(testCase)
            audit = testCase.problem.auditAPVCompatiblePrimitive();
            vorticity = audit.vorticityDivergence.diagnostics;
            testCase.verifyLessThan(vorticity.apvDefect,1e-12)
            testCase.verifyLessThan(vorticity.enstrophyDefect,1e-12)
            testCase.verifyLessThan(vorticity.bottomDefect,1e-11)
            testCase.verifyGreaterThan(vorticity.energyDefect,1e-9)
            testCase.verifyGreaterThan(vorticity.weakEvolutionDefect,1e-7)
            explicit = audit.explicitAPV;
            testCase.verifyEqual(explicit.apvCoordinateRank,explicit.coordinateDimension)
            testCase.verifyEqual(explicit.vorticityDivergenceRank,explicit.coordinateDimension)
            testCase.verifyLessThan(explicit.qRowTendencyDefect,1e-12)
            testCase.verifyGreaterThan(explicit.energyDefect,1e-13)
        end

        function compatiblePolynomialComplexCommutes(testCase)
            degree = [4 6 8 10];
            for n = degree
                audit = testCase.problem.auditAPVCompatiblePrimitive(polynomialDegree=n);
                testCase.verifyLessThan(audit.commutingDiagnostics.productRuleDefect,2e-11)
                testCase.verifyEqual(audit.commutingDiagnostics.continuityRank, ...
                    audit.commutingDiagnostics.expectedContinuityRank)
                testCase.verifyEqual(audit.commutingDiagnostics.admissibleDimension, ...
                    audit.commutingDiagnostics.expectedAdmissibleDimension)
                testCase.verifyLessThan(audit.energyWeak.diagnostics.energyDefect,2e-12)
                testCase.verifyLessThan(audit.energyWeak.diagnostics.bottomDefect,2e-11)
                testCase.verifyGreaterThan(audit.energyWeak.diagnostics.apvDefect,1e-7)
                testCase.verifyLessThan(audit.analyticObstruction.sourceComparisonDefect,0.02)
            end
            audit = testCase.problem.auditAPVCompatiblePrimitive(polynomialDegree=8);
            testCase.verifyLessThan(audit.flatReference.diagnostics.modeOneDispersionDefect,1e-10)
        end

        function sourceChangesSignWithSlopeAndMode(testCase)
            positive = testCase.problem.auditAPVCompatiblePrimitive( ...
                bottomSlope=[0 0.01],horizontalMode=[1 0],polynomialDegree=4);
            negative = testCase.problem.auditAPVCompatiblePrimitive( ...
                bottomSlope=[0 -0.01],horizontalMode=[1 0],polynomialDegree=4);
            oblique = testCase.problem.auditAPVCompatiblePrimitive( ...
                bottomSlope=[0.006 0.008],horizontalMode=[1 1],polynomialDegree=4);
            testCase.verifyEqual(positive.analyticObstruction.crossSlopeWavenumber, ...
                -negative.analyticObstruction.crossSlopeWavenumber,RelTol=1e-14)
            expectedOblique = oblique.horizontalWavenumber(1)*oblique.bottomSlope(2) ...
                -oblique.horizontalWavenumber(2)*oblique.bottomSlope(1);
            testCase.verifyEqual(oblique.analyticObstruction.crossSlopeWavenumber, ...
                expectedOblique,RelTol=1e-14)
            testCase.verifyEqual(oblique.status,"mathematical-blocker")
        end

        function doubledQuadratureLeavesClassificationUnchanged(testCase)
            reference = testCase.problem.auditAPVCompatiblePrimitive(polynomialDegree=6);
            doubled = testCase.problem.auditAPVCompatiblePrimitive( ...
                polynomialDegree=6,quadratureOrder=2*reference.quadratureOrder);
            testCase.verifyEqual(doubled.status,reference.status)
            testCase.verifyLessThan(abs(doubled.energyWeak.diagnostics.apvDefect ...
                -reference.energyWeak.diagnostics.apvDefect),1e-6)
            testCase.verifyLessThan(abs(doubled.analyticObstruction.sourceComparisonDefect ...
                -reference.analyticObstruction.sourceComparisonDefect),1e-3)
        end

        function arbitraryStratificationRetainsAlgebraicClassification(testCase)
            N2 = @(z) 1.2e-5*exp(z/1800);
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[4 4 5], ...
                N2Function=N2,latitude=45,shouldAntialias=false);
            exponentialProblem = WVTerrainEnergyGalerkin.fromTopography( ...
                wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
            audit = exponentialProblem.auditAPVCompatiblePrimitive(polynomialDegree=4);
            testCase.verifyEqual(audit.status,"mathematical-blocker")
            testCase.verifyLessThan(audit.energyWeak.diagnostics.energyDefect,1e-12)
            testCase.verifyLessThan(audit.energyWeak.diagnostics.bottomDefect,1e-11)
            testCase.verifyLessThan(audit.vorticityDivergence.diagnostics.apvDefect,1e-12)
        end

        function auditPreservesPublicLayoutAndRejectsInvalidInputs(testCase)
            horizontalBefore = testCase.problem.horizontalLayout;
            stateBefore = testCase.problem.stateLayout;
            conjugateBefore = testCase.problem.conjugateCoordinateIndex;
            testCase.problem.auditAPVCompatiblePrimitive(polynomialDegree=4);
            testCase.verifyEqual(testCase.problem.horizontalLayout,horizontalBefore)
            testCase.verifyEqual(testCase.problem.stateLayout,stateBefore)
            testCase.verifyEqual(testCase.problem.conjugateCoordinateIndex,conjugateBefore)
            testCase.verifyError(@()testCase.problem.auditAPVCompatiblePrimitive( ...
                bottomSlope=[NaN 0]),"WVTerrainEnergyGalerkin:InvalidBottomSlope")
            testCase.verifyError(@()testCase.problem.auditAPVCompatiblePrimitive( ...
                horizontalMode=[0 0]),"WVTerrainEnergyGalerkin:ZeroHorizontalMode")
            testCase.verifyError(@()testCase.problem.auditAPVCompatiblePrimitive( ...
                horizontalMode=[17 19]),"WVTerrainEnergyGalerkin:HorizontalModeNotRetained")
            testCase.verifyError(@()testCase.problem.auditAPVCompatiblePrimitive( ...
                polynomialDegree=1),"WVTerrainEnergyGalerkin:InvalidPolynomialDegree")
            testCase.verifyError(@()testCase.problem.auditAPVCompatiblePrimitive( ...
                polynomialDegree=6,quadratureOrder=4), ...
                "WVTerrainEnergyGalerkin:InsufficientQuadratureOrder")
        end
    end
end
