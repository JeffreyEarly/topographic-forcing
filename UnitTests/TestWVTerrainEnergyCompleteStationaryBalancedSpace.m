classdef TestWVTerrainEnergyCompleteStationaryBalancedSpace < matlab.unittest.TestCase
    % Verify the Milestone-8 complete stationary balanced-space gate.

    properties
        problem
        audit
        flatAudit
        uniformAudit
        variableStratificationAudit
    end

    methods (TestClassSetup)
        function createReferenceAudits(testCase)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            z = linspace(-1200,0,5)';
            N2 = @(z) 2e-5+0*z;
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[6 10 5], ...
                N2Function=N2,latitude=45,shouldAntialias=false,z=z);
            [~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
                (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
            terrain = 20*cos(2*pi*y/wvt.Ly);
            testCase.problem = WVTerrainEnergyGalerkin.fromTopography( ...
                wvt,topographicHeight=terrain);
            testCase.audit = testCase.problem.auditCompleteStationaryBalancedSpace( ...
                trustedModeBounds=[1 1],supportModeBounds=[1 2;1 3;1 4], ...
                primitivePolynomialDegrees=[2;3;4],paddingFactors=[2;3]);

            flatProblem = WVTerrainEnergyGalerkin.fromTopography( ...
                wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
            testCase.flatAudit = flatProblem.auditCompleteStationaryBalancedSpace( ...
                trustedModeBounds=[1 1],supportModeBounds=[1 1;1 1;1 1], ...
                primitivePolynomialDegrees=[2;3;4],paddingFactors=[2;3]);

            uniformProblem = WVTerrainEnergyGalerkin.fromTopography( ...
                wvt,topographicHeight=20*ones(wvt.Nx,wvt.Ny));
            testCase.uniformAudit = uniformProblem.auditCompleteStationaryBalancedSpace( ...
                trustedModeBounds=[1 1],supportModeBounds=[1 1;1 1;1 1], ...
                primitivePolynomialDegrees=[2;3;4],paddingFactors=[2;3]);

            variableN2 = @(z) 1e-5+1e-5*exp(z/1000);
            variableTransform = WVTransformBoussinesq( ...
                [24e3 20e3 1200],[6 10 5],N2Function=variableN2, ...
                latitude=45,shouldAntialias=false,z=z);
            variableProblem = WVTerrainEnergyGalerkin.fromTopography( ...
                variableTransform,topographicHeight=terrain);
            testCase.variableStratificationAudit = ...
                variableProblem.auditCompleteStationaryBalancedSpace( ...
                trustedModeBounds=[1 0], ...
                supportModeBounds=[1 1;1 2;1 3], ...
                primitivePolynomialDegrees=[2;3;5],paddingFactors=[2;3]);
        end
    end

    methods (Test)
        function completeStationarySpacePasses(testCase)
            audit = testCase.audit;
            finest = audit.finest;
            testCase.verifyEqual(audit.status,"complete-stationary-space-oracle")
            testCase.verifyTrue(audit.isCompatible)
            testCase.verifyLessThan(finest.stateRepresentationDefect,1e-10)
            testCase.verifyLessThan(finest.greenIdentityDefect,1e-10)
            testCase.verifyLessThan(finest.stationaryRowDefect,1e-10)
            testCase.verifyLessThan(finest.projectedBottomTangencyDefect,1e-10)
            testCase.verifyLessThan(finest.strongBottomTangencyDefect,1e-10)
            testCase.verifyLessThan(finest.conjugacyDefect,1e-11)
        end

        function bottomCompleteGreenIdentityIsActive(testCase)
            finest = testCase.audit.finest;
            testCase.verifyGreaterThan(finest.greenBoundaryRowNorm,0)
            testCase.verifyGreaterThan(finest.stationaryVolumeAPVRank,0)
            testCase.verifyGreaterThan( ...
                finest.stationaryZeroVolumeAPVCandidateDimension,0)
            testCase.verifyGreaterThan( ...
                finest.zeroVolumeAPVBottomCandidateRank,0)
            testCase.verifyGreaterThan(finest.meanDensityAnomalyDimension,0)
            testCase.verifyLessThan( ...
                finest.baselineBottomInversionMaximumAPVResidual,1e-10)
            testCase.verifyLessThan( ...
                finest.baselineBottomInversionMaximumResidual,1e-10)
        end

        function nonTangentBottomDirectionsRemainDynamical(testCase)
            finest = testCase.audit.finest;
            testCase.verifyGreaterThan(finest.nonstationaryBottomDimension,0)
            testCase.verifyGreaterThan(finest.fullSupportBottomTangencyDefect, ...
                finest.projectedBottomTangencyDefect)
            streamfunction = finest.bottomStreamfunctionByHorizontalMode;
            nonTangentRows = finest.horizontalLayout.kMode ~= 0;
            testCase.verifyLessThan(norm(streamfunction(nonTangentRows,:),"fro"), ...
                1e-10*max(norm(streamfunction,"fro"),1))
        end

        function rankAndPaddingDiagnosticsAreStable(testCase)
            finest = testCase.audit.finest;
            testCase.verifyTrue(finest.scalarRankStable)
            testCase.verifyTrue(finest.stationaryNullityStable)
            testCase.verifyGreaterThan(finest.stationaryNullspaceGap,100)
            testCase.verifyEqual(numel(unique( ...
                testCase.audit.convergence.numberOfStationaryStates)),1)
            testCase.verifyLessThan(testCase.audit.padding.maximumDefect,1e-10)
            testCase.verifyGreaterThan(finest.numberOfWeakNullStates, ...
                finest.numberOfStationaryStates)
        end

        function flatAndUniformDepthLimitsPass(testCase)
            flat = testCase.flatAudit;
            uniform = testCase.uniformAudit;
            testCase.verifyTrue(flat.isCompatible)
            testCase.verifyTrue(uniform.isCompatible)
            testCase.verifyEqual(flat.finest.numberOfStationaryStates, ...
                uniform.finest.numberOfStationaryStates)
            testCase.verifyEqual(flat.finest.nonstationaryBottomDimension,0)
            testCase.verifyEqual(uniform.finest.nonstationaryBottomDimension,0)
            testCase.verifyLessThan(flat.finest.greenIdentityDefect,1e-10)
            testCase.verifyLessThan(uniform.finest.greenIdentityDefect,1e-10)
            testCase.verifyLessThan(flat.finest.stationaryRowDefect,1e-10)
            testCase.verifyLessThan(uniform.finest.stationaryRowDefect,1e-10)
        end

        function variableStratificationPasses(testCase)
            variableAudit = testCase.variableStratificationAudit;
            testCase.verifyTrue(variableAudit.isCompatible)
            testCase.verifyLessThan(variableAudit.finest.stateRepresentationDefect,1e-10)
            testCase.verifyLessThan(variableAudit.finest.greenIdentityDefect,1e-10)
            testCase.verifyLessThan(variableAudit.finest.stationaryRowDefect,1e-10)
            testCase.verifyLessThan(variableAudit.finest.projectedBottomTangencyDefect,1e-10)
            testCase.verifyLessThan(variableAudit.padding.maximumDefect,1e-10)
        end

        function antialiasedLayoutPreservesConjugacy(testCase)
            N2 = @(z) 2e-5+0*z;
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[12 12 5], ...
                N2Function=N2,latitude=45,shouldAntialias=true);
            [~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
                (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
            antialiasedProblem = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
                topographicHeight=20*cos(2*pi*y/wvt.Ly));
            antialiasedAudit = antialiasedProblem.auditCompleteStationaryBalancedSpace( ...
                trustedModeBounds=[1 0], ...
                supportModeBounds=[1 1;1 1;1 1], ...
                primitivePolynomialDegrees=[2;3;4],paddingFactors=[2;3]);
            testCase.verifyLessThan( ...
                antialiasedAudit.finest.stateRepresentationDefect,1e-10)
            testCase.verifyLessThan( ...
                antialiasedAudit.finest.projectedBottomTangencyDefect,1e-10)
            testCase.verifyLessThan(antialiasedAudit.finest.conjugacyDefect,1e-11)
            testCase.verifyLessThan(antialiasedAudit.padding.maximumDefect,1e-10)
        end

        function validationPreservesPublicLayouts(testCase)
            horizontalBefore = testCase.problem.horizontalLayout;
            stateBefore = testCase.problem.stateLayout;
            testCase.verifyError(@() ...
                testCase.problem.auditCompleteStationaryBalancedSpace( ...
                trustedModeBounds=[1 1], ...
                supportModeBounds=[1 2;1 3], ...
                primitivePolynomialDegrees=[2;3],paddingFactors=[2;3]), ...
                "WVTerrainEnergyGalerkin:InvalidStationarySpaceDegrees")
            testCase.verifyError(@() ...
                testCase.problem.auditCompleteStationaryBalancedSpace( ...
                trustedModeBounds=[1 1], ...
                supportModeBounds=[1 1;1 2;1 3], ...
                primitivePolynomialDegrees=[2;3;4],paddingFactors=[2;3]), ...
                "WVTerrainEnergyGalerkin:InsufficientStationaryGuard")
            testCase.verifyEqual(testCase.problem.horizontalLayout,horizontalBefore)
            testCase.verifyEqual(testCase.problem.stateLayout,stateBefore)
        end
    end
end
