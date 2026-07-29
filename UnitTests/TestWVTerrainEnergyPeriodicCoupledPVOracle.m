classdef TestWVTerrainEnergyPeriodicCoupledPVOracle < matlab.unittest.TestCase
    % Verify the projected periodic volume--boundary PV oracle.

    properties
        problem
    end

    methods (TestClassSetup)
        function createProblem(testCase)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            testCase.problem = TestWVTerrainEnergyPeriodicCoupledPVOracle.sinusoidalProblem(9,false,2,@(z)2e-5+0*z);
        end
    end

    methods (Test)
        function exactProjectedDynamicsCloseAtRoundoff(testCase)
            audit = testCase.problem.auditPeriodicCoupledPVOracle(polynomialDegree=12);
            testCase.verifyEqual(audit.status,"passed")
            testCase.verifyLessThan(audit.forms.diagnostics.maximumInversionDefect,1e-12)
            testCase.verifyLessThan(audit.forms.diagnostics.maximumGreenIdentityDefect,1e-12)
            testCase.verifyGreaterThan(audit.forms.diagnostics.minimumEnergyEigenvalue,0)
            testCase.verifyLessThan(audit.dynamics.diagnostics.bottomOperatorSkewDefect,1e-12)
            testCase.verifyLessThan(audit.dynamics.diagnostics.energyDefect,1e-12)
            testCase.verifyLessThan(audit.dynamics.diagnostics.apvDefect,1e-12)
            testCase.verifyLessThan(audit.dynamics.diagnostics.enstrophyDefect,1e-12)
            testCase.verifyLessThan(audit.dynamics.diagnostics.bottomDefect,1e-12)
            testCase.verifyLessThan(audit.dynamics.diagnostics.conjugacyDefect,1e-12)
            testCase.verifyLessThan(audit.dynamics.diagnostics.frequencyImaginaryDefect,1e-11)
            testCase.verifyGreaterThan(audit.fourier.numberOfNonzeroFrequencies,0)
            testCase.verifyLessThan(audit.fourier.nonzeroModeVolumeAPVDefect,1e-12)
        end

        function exactConvolutionMatchesPseudospectralJacobian(testCase)
            audit = testCase.problem.auditPeriodicCoupledPVOracle(polynomialDegree=12);
            testCase.verifyLessThan(audit.pseudospectral.diagnostics.exactConvolutionDefect,1e-11)
            testCase.verifyLessThan(audit.pseudospectral.diagnostics.meanTendencyDefect,1e-11)
            testCase.verifyLessThan(audit.pseudospectral.phaseGramDefect,1e-12)
            testCase.verifyLessThan(audit.fourier.couplingLeakage,1e-12)
            testCase.verifyEqual(sortrows(audit.fourier.terrainModes),[0 -1;0 1])

            iIn = find(audit.layout.kMode == 1 & audit.layout.lMode == 0,1);
            iOut = find(audit.layout.kMode == 1 & audit.layout.lMode == 1,1);
            hCoefficient = 10;
            qy = 2*pi/testCase.problem.originatingTransform.Ly;
            expected = testCase.problem.originatingTransform.f*audit.layout.k(iIn)*qy*hCoefficient;
            testCase.verifyEqual(audit.exactConvolution.bottomOperator(iOut,iIn),expected,"RelTol",1e-12)
        end

        function spectralEdgesUseProjectedRatherThanAliasedClosure(testCase)
            audit = testCase.problem.auditPeriodicCoupledPVOracle(polynomialDegree=12);
            sidebands = audit.exactConvolution.sidebands;
            testCase.verifyGreaterThan(sidebands.numberOfInteriorHorizontalModes,0)
            testCase.verifyGreaterThan(sidebands.numberOfEdgeHorizontalModes,0)
            testCase.verifyGreaterThan(sidebands.maximumDiscardedOperatorNorm,0)
            testCase.verifyEqual(sidebands.maximumMeanCoefficient,0,"AbsTol",1e-14)
            testCase.verifyLessThan(audit.fourier.projectedEdgeAPVDefect,1e-12)
            testCase.verifyLessThan(audit.dynamics.diagnostics.energyDefect,1e-12)
        end

        function flatTerrainProducesStationaryCoupledPVState(testCase)
            wvt = testCase.problem.originatingTransform;
            flat = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny));
            audit = flat.auditPeriodicCoupledPVOracle(polynomialDegree=8);
            testCase.verifyEqual(audit.status,"passed")
            testCase.verifyEqual(audit.exactConvolution.bottomOperator,zeros(size(audit.exactConvolution.bottomOperator)),"AbsTol",1e-14)
            testCase.verifyEqual(audit.dynamics.generator,zeros(size(audit.dynamics.generator)),"AbsTol",1e-14)
            testCase.verifyEqual(audit.dynamics.frequency,zeros(size(audit.dynamics.frequency)),"AbsTol",1e-14)
        end

        function balancedBasisBridgeConverges(testCase)
            resolution = [5 7 9];
            polynomialDegree = 2*resolution;
            defect = zeros(size(resolution));
            for iResolution = 1:numel(resolution)
                currentProblem = TestWVTerrainEnergyPeriodicCoupledPVOracle.sinusoidalProblem(resolution(iResolution),false,2,@(z)2e-5+0*z);
                audit = currentProblem.auditPeriodicCoupledPVOracle(polynomialDegree=polynomialDegree(iResolution));
                defect(iResolution) = audit.bridge.diagnostics.resolvedReconstructionDefect;
                testCase.verifyEqual(audit.bridge.diagnostics.rank,audit.bridge.diagnostics.expectedRank)
                testCase.verifyLessThan(audit.bridge.diagnostics.maximumBottomValueDefect,1e-12)
            end
            testCase.verifyTrue(all(diff(defect) < 0))
            testCase.verifyLessThan(defect(end),2e-2)
        end

        function arbitraryStratificationAndAntialiasingRetainClosure(testCase)
            currentProblem = TestWVTerrainEnergyPeriodicCoupledPVOracle.sinusoidalProblem(11,true,3,@(z)1.2e-5*exp(z/1800));
            audit = currentProblem.auditPeriodicCoupledPVOracle(polynomialDegree=16);
            testCase.verifyEqual(audit.status,"passed")
            testCase.verifyLessThan(audit.pseudospectral.diagnostics.exactConvolutionDefect,1e-11)
            testCase.verifyLessThan(audit.dynamics.diagnostics.energyDefect,1e-12)
            testCase.verifyLessThan(audit.dynamics.diagnostics.apvDefect,1e-12)
            testCase.verifyLessThan(audit.dynamics.diagnostics.bottomDefect,1e-12)
        end

        function auditPreservesLayoutAndRejectsUnsupportedInputs(testCase)
            horizontalBefore = testCase.problem.horizontalLayout;
            stateBefore = testCase.problem.stateLayout;
            testCase.problem.auditPeriodicCoupledPVOracle(polynomialDegree=12);
            testCase.verifyEqual(testCase.problem.horizontalLayout,horizontalBefore)
            testCase.verifyEqual(testCase.problem.stateLayout,stateBefore)
            testCase.verifyError(@()testCase.problem.auditPeriodicCoupledPVOracle(polynomialDegree=2), ...
                "WVTerrainEnergyGalerkin:InvalidPolynomialDegree")
            testCase.verifyError(@()testCase.problem.auditPeriodicCoupledPVOracle(polynomialDegree=8,quadratureOrder=8), ...
                "WVTerrainEnergyGalerkin:InsufficientQuadratureOrder")

            wvt = testCase.problem.originatingTransform;
            nyquistTerrain = repmat((-1).^(0:wvt.Nx-1)',1,wvt.Ny);
            nyquistProblem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=nyquistTerrain);
            testCase.verifyError(@()nyquistProblem.auditPeriodicCoupledPVOracle(polynomialDegree=8), ...
                "WVTerrainEnergyGalerkin:PeriodicPVOracleTerrainNyquist")
        end
    end

    methods (Static,Access=private)
        function problem = sinusoidalProblem(Nz,shouldAntialias,oversampling,N2)
            z = linspace(-1200,0,Nz)';
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[6 6 Nz], ...
                N2Function=N2,latitude=45,shouldAntialias=shouldAntialias,z=z);
            [~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
                (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
                topographicHeight=20*cos(2*pi*y/wvt.Ly), ...
                horizontalOversamplingFactor=oversampling);
        end
    end
end
