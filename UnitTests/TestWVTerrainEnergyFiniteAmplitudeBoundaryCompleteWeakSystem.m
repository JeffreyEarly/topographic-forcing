classdef TestWVTerrainEnergyFiniteAmplitudeBoundaryCompleteWeakSystem < matlab.unittest.TestCase
    % Verify the Milestone-7 finite-amplitude boundary-complete weak gate.

    properties
        problem
        audit
        uniformAudit
        variableStratificationAudit
    end

    methods (TestClassSetup)
        function createReferenceAudit(testCase)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            N2 = @(z) 2e-5+0*z;
            z = linspace(-1200,0,5)';
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[6 10 5], ...
                N2Function=N2,latitude=45,shouldAntialias=false,z=z);
            [~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
                (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
            testCase.problem = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
                topographicHeight=20*cos(2*pi*y/wvt.Ly));
            testCase.audit = testCase.problem.auditFiniteAmplitudeBoundaryCompleteWeakSystem( ...
                trustedModeBounds=[1 0],supportModeBounds=[1 1;1 2;1 4], ...
                scalarPolynomialDegree=2,primitivePolynomialDegrees=[2;3;5], ...
                paddingFactors=[2;3],terrainScales=[0.5;1]);
            uniformProblem = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
                topographicHeight=20*ones(wvt.Nx,wvt.Ny));
            testCase.uniformAudit = uniformProblem.auditFiniteAmplitudeBoundaryCompleteWeakSystem( ...
                trustedModeBounds=[1 0],supportModeBounds=[1 0;1 0;1 0], ...
                scalarPolynomialDegree=2,primitivePolynomialDegrees=[2;3;5], ...
                paddingFactors=[2;3],terrainScales=1);
            variableN2 = @(z) 1e-5+1e-5*exp(z/1000);
            variableTransform = WVTransformBoussinesq([24e3 20e3 1200],[6 10 5], ...
                N2Function=variableN2,latitude=45,shouldAntialias=false,z=z);
            variableProblem = WVTerrainEnergyGalerkin.fromTopography( ...
                variableTransform,topographicHeight=20*cos(2*pi*y/wvt.Ly));
            testCase.variableStratificationAudit = ...
                variableProblem.auditFiniteAmplitudeBoundaryCompleteWeakSystem( ...
                trustedModeBounds=[1 0],supportModeBounds=[1 1;1 2;1 4], ...
                scalarPolynomialDegree=2,primitivePolynomialDegrees=[2;4;8], ...
                paddingFactors=[2;3],terrainScales=1);
        end
    end

    methods (Test)
        function finiteAmplitudeWeakSequencePasses(testCase)
            audit = testCase.audit;
            finest = audit.finest;
            testCase.verifyEqual(audit.status,"compatible-finite-amplitude-weak-oracle")
            testCase.verifyTrue(audit.isCompatible)
            testCase.verifyLessThan(finest.stateRepresentationDefect,1e-10)
            testCase.verifyLessThan(finest.greenIdentityDefect,1e-10)
            testCase.verifyLessThan(finest.stationaryRowDefect,1e-10)
            testCase.verifyLessThan(finest.weakAPVDefect,1e-10)
            testCase.verifyLessThan(finest.strongAPVDefect,1e-10)
            testCase.verifyLessThan(finest.weakStrongAgreementDefect,1e-10)
        end

        function primitiveEnergyBottomAndProjectionPass(testCase)
            audit = testCase.audit;
            finest = audit.finest;
            testCase.verifyLessThan(finest.weakEvolutionDefect,1e-10)
            testCase.verifyLessThan(finest.primitiveGeneratorAgreementDefect,1e-10)
            testCase.verifyLessThan(finest.energyDefect,1e-11)
            testCase.verifyLessThan(finest.bottomDefect,1e-10)
            testCase.verifyLessThan(finest.conjugacyDefect,1e-11)
            testCase.verifyLessThan(finest.enstrophyDefect,1e-10)
            testCase.verifyLessThan(audit.padding.maximumDefect,1e-10)
            testCase.verifyGreaterThan(finest.fullSupportBottomDefect,finest.bottomDefect)
        end

        function strongDiagnosticsConvergeWithSupport(testCase)
            convergence = testCase.audit.convergence;
            testCase.verifyLessThan( ...
                convergence.primitiveGeneratorAgreementDefect(end), ...
                convergence.primitiveGeneratorAgreementDefect(1))
            testCase.verifyLessThan( ...
                convergence.primitiveWeakEvolutionDefect(end), ...
                convergence.primitiveWeakEvolutionDefect(1))
            testCase.verifyLessThan( ...
                convergence.primitiveEnergyDefect(end), ...
                convergence.primitiveEnergyDefect(1))
            testCase.verifyGreaterThan(convergence.bottomDefectReduction(end),100)
        end

        function projectedProductsReportDiscardedSidebands(testCase)
            finest = testCase.audit.finest;
            testCase.verifyLessThan(finest.projectionAdjointDefect,1e-12)
            testCase.verifyLessThan(finest.exactConvolutionDefect,1e-12)
            testCase.verifyGreaterThan( ...
                finest.projection.maximumDiscardedOperatorNorm,0)
            testCase.verifyGreaterThan( ...
                finest.projection.numberOfEdgeHorizontalModes,0)
        end

        function nonlinearRemainderIsSecondOrder(testCase)
            amplitude = testCase.audit.amplitude;
            testCase.verifyGreaterThan(amplitude.energyReduction(1),3.5)
            testCase.verifyGreaterThan(amplitude.exchangeReduction(1),3.5)
            testCase.verifyGreaterThan(amplitude.geostrophicInclusionReduction(1),3.5)
        end

        function flatLimitCloses(testCase)
            flat = testCase.audit.flat;
            testCase.verifyLessThan(flat.greenIdentityDefect,1e-10)
            testCase.verifyLessThan(flat.stationaryRowDefect,1e-10)
            testCase.verifyLessThan(flat.weakAPVDefect,1e-10)
            testCase.verifyLessThan(flat.strongAPVDefect,1e-10)
            testCase.verifyLessThan(flat.energyDefect,1e-11)
            testCase.verifyLessThan(flat.bottomResidualNorm,1e-12)
        end

        function uniformDepthCloses(testCase)
            uniform = testCase.uniformAudit;
            testCase.verifyTrue(uniform.isCompatible)
            testCase.verifyLessThan(uniform.finest.greenIdentityDefect,1e-10)
            testCase.verifyLessThan(uniform.finest.stationaryRowDefect,1e-10)
            testCase.verifyLessThan(uniform.finest.strongAPVDefect,1e-10)
            testCase.verifyLessThan(uniform.finest.energyDefect,1e-11)
            testCase.verifyLessThan(uniform.finest.bottomResidualNorm,1e-12)
        end

        function variableStratificationConverges(testCase)
            variable = testCase.variableStratificationAudit;
            testCase.verifyTrue(variable.isCompatible)
            testCase.verifyLessThan(variable.finest.stateRepresentationDefect,1e-10)
            testCase.verifyLessThan(variable.finest.greenIdentityDefect,1e-10)
            testCase.verifyLessThan(variable.finest.stationaryRowDefect,1e-10)
            testCase.verifyLessThan(variable.finest.strongAPVDefect,1e-10)
            testCase.verifyLessThan(variable.finest.bottomDefect,1e-10)
            testCase.verifyLessThan(variable.finest.energyDefect,1e-11)
            testCase.verifyGreaterThan( ...
                variable.convergence.stateRepresentationDefectReduction(end),100)
        end

        function validationPreservesPublicLayouts(testCase)
            horizontalBefore = testCase.problem.horizontalLayout;
            stateBefore = testCase.problem.stateLayout;
            testCase.verifyError(@()testCase.problem.auditFiniteAmplitudeBoundaryCompleteWeakSystem( ...
                trustedModeBounds=[1 0],supportModeBounds=[1 1;1 2], ...
                scalarPolynomialDegree=2,primitivePolynomialDegrees=[2;3], ...
                paddingFactors=[2;3]), ...
                "WVTerrainEnergyGalerkin:InvalidFiniteAmplitudeDegrees")
            testCase.verifyError(@()testCase.problem.auditFiniteAmplitudeBoundaryCompleteWeakSystem( ...
                trustedModeBounds=[1 0],supportModeBounds=[1 1;1 2;1 4], ...
                scalarPolynomialDegree=2,primitivePolynomialDegrees=[2;3;5], ...
                paddingFactors=[2;3],terrainScales=[0.25;0.5]), ...
                "WVTerrainEnergyGalerkin:InvalidFiniteAmplitudeScales")
            testCase.verifyEqual(testCase.problem.horizontalLayout,horizontalBefore)
            testCase.verifyEqual(testCase.problem.stateLayout,stateBefore)
        end
    end
end
