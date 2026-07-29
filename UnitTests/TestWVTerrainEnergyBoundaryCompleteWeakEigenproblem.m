classdef TestWVTerrainEnergyBoundaryCompleteWeakEigenproblem < matlab.unittest.TestCase
    % Verify the Milestone-6.8 boundary-complete weak eigenproblem oracle.

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
            testCase.audit = testCase.problem.auditBoundaryCompleteWeakEigenproblem( ...
                trustedModeBounds=[1 0],supportModeBounds=[1 2], ...
                polynomialDegrees=[2;3;4],paddingFactor=2);
        end
    end

    methods (Test)
        function geostrophicTestSequenceCloses(testCase)
            audit = testCase.audit;
            testCase.verifyEqual(audit.status,"compatible-boundary-complete-weak-oracle")
            testCase.verifyTrue(audit.isCompatible)
            testCase.verifyLessThan(max([audit.refinement.scalarBoundaryDefect]),1e-12)
            testCase.verifyLessThan(audit.refinement(end).trustedStateRepresentationDefect,1e-10)
            testCase.verifyLessThan(audit.finest.flatGreenIdentityDefect,1e-10)
            testCase.verifyLessThan(audit.finest.trustedTangentGreenIdentityDefect,1e-10)
            testCase.verifyLessThan(audit.finest.flatStationaryRowDefect,1e-10)
            testCase.verifyLessThan(audit.finest.trustedTangentStationaryRowDefect,1e-10)
        end

        function weakAPVIsAConsequenceOfPrimitiveRows(testCase)
            finest = testCase.audit.finest;
            testCase.verifyLessThan(finest.trustedWeakAPVDefect,1e-10)
            testCase.verifyLessThan(finest.trustedDirectAPVDefect,1e-10)
            testCase.verifyLessThan(finest.trustedWeakStrongAgreementDefect,1e-10)
            testCase.verifyGreaterThan(finest.scalar.numberOfAPVTests, ...
                finest.scalar.numberOfTrustedTests)
        end

        function tangentModesPreserveTheDerivedStructure(testCase)
            modes = testCase.audit.finest.modeDiagnostics;
            testCase.verifyLessThan(modes.eigenproblemDefect,1e-10)
            testCase.verifyLessThan(max(modes.frequencyImaginaryDefectByScale),1e-10)
            testCase.verifyLessThan(max(modes.orthogonalityDefectByScale),1e-8)
            testCase.verifyLessThan(modes.trustedWaveAPVDefect,1e-10)
            testCase.verifyGreaterThan(min(modes.numberOfBoundaryActiveWaves),0)
            testCase.verifyGreaterThan(min(modes.bottomResidualByScale(1:end-1) ...
                ./modes.bottomResidualByScale(2:end)),1.8)
        end

        function arbitraryStratificationConvergesWithVerticalSupport(testCase)
            N2 = @(z) 1.2e-5*exp(z/1800);
            z = linspace(-900,0,5)';
            wvt = WVTransformBoussinesq([18e3 24e3 900],[10 12 5], ...
                N2Function=N2,latitude=45,shouldAntialias=true,z=z);
            [~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
                (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
            variableProblem = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
                topographicHeight=10*cos(2*pi*y/wvt.Ly));
            variableAudit = variableProblem.auditBoundaryCompleteWeakEigenproblem( ...
                trustedModeBounds=[0 0],supportModeBounds=[0 2], ...
                polynomialDegrees=[2;4;8],paddingFactor=2);
            representation = [variableAudit.refinement.trustedStateRepresentationDefect];
            testCase.verifyLessThan(representation(3),1e-10)
            testCase.verifyGreaterThan(min(representation(1:2)./representation(2:3)),100)
            testCase.verifyLessThan(variableAudit.finest.trustedTangentGreenIdentityDefect,1e-10)
            testCase.verifyLessThan(variableAudit.finest.trustedTangentStationaryRowDefect,1e-10)
        end

        function auditRejectsIncompleteSupportAndPreservesLayouts(testCase)
            horizontalBefore = testCase.problem.horizontalLayout;
            stateBefore = testCase.problem.stateLayout;
            testCase.verifyError(@()testCase.problem.auditBoundaryCompleteWeakEigenproblem( ...
                trustedModeBounds=[1 0],supportModeBounds=[1 0], ...
                polynomialDegrees=[2;3;4]), ...
                "WVTerrainEnergyGalerkin:InsufficientWeakEigenproblemGuard")
            testCase.verifyError(@()testCase.problem.auditBoundaryCompleteWeakEigenproblem( ...
                trustedModeBounds=[1 0],supportModeBounds=[3 2], ...
                polynomialDegrees=[2;3;4]), ...
                "WVTerrainEnergyGalerkin:UnavailableWeakEigenproblemSupport")
            testCase.verifyEqual(testCase.problem.horizontalLayout,horizontalBefore)
            testCase.verifyEqual(testCase.problem.stateLayout,stateBefore)
        end
    end
end
