classdef TestWVTerrainEnergyDealiasedProjectedPrimitiveTangent < matlab.unittest.TestCase
    % Verify the Milestone-6.7 projected primitive tangent oracle.

    properties
        problem
        audit
    end

    methods (TestClassSetup)
        function createReferenceAudit(testCase)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            N2 = @(z) 2e-5+0*z;
            z = linspace(-1200,0,5)';
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[6 8 5], ...
                N2Function=N2,latitude=45,shouldAntialias=false,z=z);
            [~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
                (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
            testCase.problem = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
                topographicHeight=20*cos(2*pi*y/wvt.Ly));
            testCase.audit = testCase.problem.auditDealiasedProjectedPrimitiveTangent( ...
                trustedModeBounds=[1 0],supportModeBounds=[1 1;1 2;1 3], ...
                polynomialDegrees=[2;3;4],paddingFactors=[2;3]);
        end
    end

    methods (Test)
        function projectionMatchesExactConvolution(testCase)
            summaries = [testCase.audit.diagonalRefinement; ...
                testCase.audit.horizontalRefinement;testCase.audit.verticalRefinement];
            testCase.verifyLessThan(max([summaries.projectionDefect]),1e-12)
            testCase.verifyTrue(testCase.audit.externalSidebands.isPresent)
            testCase.verifyGreaterThan(max(testCase.audit.externalSidebands.maximumDiscardedOperatorNorm),0)
        end

        function primitiveIdentitiesAndPaddingPass(testCase)
            structural = testCase.audit.structural;
            testCase.verifyLessThan(structural.maximumTangentDefect,1e-9)
            testCase.verifyLessThan(structural.maximumWeakDefect,1e-11)
            testCase.verifyLessThan(structural.maximumEnergyDefect,1e-11)
            testCase.verifyLessThan(structural.maximumBottomDefect,1e-11)
            testCase.verifyLessThan(structural.maximumConjugacyDefect,1e-11)
            testCase.verifyLessThan(structural.maximumPaddingDefect,1e-10)
        end

        function classificationReportsRatherThanRepairs(testCase)
            audit = testCase.audit;
            testCase.verifyEqual(audit.status,"nonconvergent-projected")
            testCase.verifyFalse(audit.isCompatible)
            testCase.verifyGreaterThan(audit.diagonalRefinement(end).trustedAPVDefect, ...
                audit.requiredTolerance.convergentProjected)
            testCase.verifyEqual(audit.nextScope, ...
                "milestone-7-only-if-exact-or-convergent-projected")
        end

        function arbitraryStratificationAndAntialiasingAreSupported(testCase)
            N2 = @(z) 1.2e-5*exp(z/1800);
            z = linspace(-900,0,5)';
            wvt = WVTransformBoussinesq([18e3 24e3 900],[10 12 5], ...
                N2Function=N2,latitude=45,shouldAntialias=true,z=z);
            [~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
                (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
            exponentialProblem = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
                topographicHeight=10*cos(2*pi*y/wvt.Ly));
            exponentialAudit = exponentialProblem.auditDealiasedProjectedPrimitiveTangent( ...
                trustedModeBounds=[0 0],supportModeBounds=[0 1;0 2;0 3], ...
                polynomialDegrees=[2;3;4],paddingFactors=[2;3],tangentStep=1e-4);
            testCase.verifyLessThan(exponentialAudit.structural.maximumProjectionDefect,1e-12)
            testCase.verifyLessThan(exponentialAudit.structural.maximumEnergyDefect,1e-11)
            testCase.verifyLessThan(exponentialAudit.structural.maximumBottomDefect,1e-11)
        end

        function auditPreservesStateAndRejectsInvalidBands(testCase)
            horizontalBefore = testCase.problem.horizontalLayout;
            stateBefore = testCase.problem.stateLayout;
            conjugateBefore = testCase.problem.conjugateCoordinateIndex;
            testCase.verifyEqual(testCase.problem.horizontalLayout,horizontalBefore)
            testCase.verifyEqual(testCase.problem.stateLayout,stateBefore)
            testCase.verifyEqual(testCase.problem.conjugateCoordinateIndex,conjugateBefore)

            testCase.verifyError(@()testCase.problem.auditDealiasedProjectedPrimitiveTangent( ...
                trustedModeBounds=[1 0],supportModeBounds=[1 1;1 2], ...
                polynomialDegrees=[2;3]), ...
                "WVTerrainEnergyGalerkin:InvalidProjectedPrimitiveRefinement")
            testCase.verifyError(@()testCase.problem.auditDealiasedProjectedPrimitiveTangent( ...
                trustedModeBounds=[1 1],supportModeBounds=[1 1;1 2;1 3], ...
                polynomialDegrees=[2;3;4]), ...
                "WVTerrainEnergyGalerkin:InsufficientProjectedPrimitiveGuard")
            testCase.verifyError(@()testCase.problem.auditDealiasedProjectedPrimitiveTangent( ...
                trustedModeBounds=[1 0],supportModeBounds=[3 1;3 2;3 3], ...
                polynomialDegrees=[2;3;4]), ...
                "WVTerrainEnergyGalerkin:UnavailableProjectedPrimitiveSupport")

            wvt = testCase.problem.originatingTransform;
            nyquistTerrain = 5*repmat((-1).^(0:wvt.Ny-1),wvt.Nx,1);
            nyquistProblem = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
                topographicHeight=nyquistTerrain);
            testCase.verifyError(@()nyquistProblem.auditDealiasedProjectedPrimitiveTangent( ...
                trustedModeBounds=[1 0],supportModeBounds=[1 1;1 2;1 3], ...
                polynomialDegrees=[2;3;4]), ...
                "WVTerrainEnergyGalerkin:ProjectedPrimitiveTerrainNyquist")
        end
    end
end
