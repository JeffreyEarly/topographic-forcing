classdef TestWVTerrainEnergyHybridPrimitivePVOracle < matlab.unittest.TestCase
    % Verify the fixed-zonal hybrid primitive-PV oracle.

    properties
        problem
    end

    methods (TestClassSetup)
        function createProblem(testCase)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            testCase.problem = TestWVTerrainEnergyHybridPrimitivePVOracle.sinusoidalProblem(7,false,@(z)2e-5+0*z);
        end
    end

    methods (Test)
        function completeHybridRowsExposePrimitiveEquivalenceBlocker(testCase)
            audit = testCase.problem.auditHybridPrimitivePVOracle(polynomialDegree=4);
            block = audit.positiveBlock;
            testCase.verifyEqual(audit.status,"primitive-equivalence-blocker")
            testCase.verifyFalse(audit.isCompatible)
            testCase.verifyEqual(block.numberOfRows,block.numberOfCoordinates)
            testCase.verifyEqual(block.descriptorRank,block.numberOfCoordinates)
            testCase.verifyGreaterThan(block.descriptorReciprocalConditionNumber,1e-10)
            testCase.verifyLessThan(block.diagnostics.flatRecoveryDefect,1e-10)
            testCase.verifyLessThan(block.diagnostics.waveWeakDefect,1e-12)
            testCase.verifyLessThan(block.diagnostics.projectedAPVDefect,1e-12)
            testCase.verifyLessThan(block.diagnostics.bottomDefect,1e-12)
            testCase.verifyGreaterThan(block.diagnostics.fullWeakDefect,1e-4)
            testCase.verifyGreaterThan(block.diagnostics.energyDefect,1e-4)
            testCase.verifyGreaterThan(block.diagnostics.fullAPVDefect,0.1)
        end

        function independentTangentConjugacyAndQGClosurePass(testCase)
            audit = testCase.problem.auditHybridPrimitivePVOracle(polynomialDegree=4);
            testCase.verifyLessThan(audit.tangentAgreement.maximumFinalDefect,1e-9)
            testCase.verifyLessThan(audit.conjugacyDefect,1e-10)
            testCase.verifyTrue(audit.qgBridge.periodicOracleStatus == "passed" ...
                || audit.qgBridge.periodicOracleStatus == "basis-bridge-blocker")
            testCase.verifyLessThan(audit.qgBridge.periodicCoreMaximumDefect,1e-12)
            testCase.verifyLessThan(audit.qgBridge.projectedClosureDefect,1e-12)
            testCase.verifyEqual(audit.qgBridge.stationaryDimension, ...
                audit.qgBridge.pvCoordinateDimension)
            testCase.verifyEqual(audit.qgBridge.pvCoordinateRank, ...
                audit.qgBridge.pvCoordinateDimension)
        end

        function spectralEdgeDefectSurvivesVerticalRefinement(testCase)
            degree = [2 4 6];
            weak = zeros(size(degree));
            energy = zeros(size(degree));
            apv = zeros(size(degree));
            for iDegree = 1:numel(degree)
                audit = testCase.problem.auditHybridPrimitivePVOracle(polynomialDegree=degree(iDegree));
                testCase.verifyEqual(audit.status,"primitive-equivalence-blocker")
                weak(iDegree) = audit.positiveBlock.diagnostics.fullWeakDefect;
                energy(iDegree) = audit.positiveBlock.diagnostics.energyDefect;
                apv(iDegree) = audit.positiveBlock.diagnostics.fullAPVDefect;
                perMode = audit.positiveBlock.diagnostics.fullAPVDefectByHorizontalMode;
                testCase.verifyLessThan(max(perMode([1 2 5])),1e-10)
                testCase.verifyGreaterThan(min(perMode([3 4])),0.4)
            end
            testCase.verifyGreaterThan(min(weak),1e-5)
            testCase.verifyGreaterThan(min(energy),1e-5)
            testCase.verifyGreaterThan(min(apv),0.2)
        end

        function arbitraryStratificationRetainsBlocker(testCase)
            currentProblem = TestWVTerrainEnergyHybridPrimitivePVOracle.sinusoidalProblem(7,true,@(z)1.2e-5*exp(z/1800));
            audit = currentProblem.auditHybridPrimitivePVOracle(polynomialDegree=3);
            testCase.verifyEqual(audit.status,"primitive-equivalence-blocker")
            testCase.verifyLessThan(audit.tangentAgreement.maximumFinalDefect,1e-9)
            testCase.verifyLessThan(audit.positiveBlock.diagnostics.projectedAPVDefect,1e-12)
            testCase.verifyLessThan(audit.positiveBlock.diagnostics.bottomDefect,1e-12)
            testCase.verifyGreaterThan(audit.positiveBlock.diagnostics.fullAPVDefect,0.1)
        end

        function auditPreservesLayoutAndRejectsUnsupportedInputs(testCase)
            horizontalBefore = testCase.problem.horizontalLayout;
            stateBefore = testCase.problem.stateLayout;
            testCase.problem.auditHybridPrimitivePVOracle(polynomialDegree=2);
            testCase.verifyEqual(testCase.problem.horizontalLayout,horizontalBefore)
            testCase.verifyEqual(testCase.problem.stateLayout,stateBefore)
            testCase.verifyError(@()testCase.problem.auditHybridPrimitivePVOracle(zonalMode=0), ...
                "WVTerrainEnergyGalerkin:HybridOracleZeroZonalMode")
            testCase.verifyError(@()testCase.problem.auditHybridPrimitivePVOracle(zonalMode=20), ...
                "WVTerrainEnergyGalerkin:HybridOracleMissingZonalBlock")
            testCase.verifyError(@()testCase.problem.auditHybridPrimitivePVOracle(polynomialDegree=1), ...
                "WVTerrainEnergyGalerkin:InvalidPolynomialDegree")
            testCase.verifyError(@()testCase.problem.auditHybridPrimitivePVOracle( ...
                polynomialDegree=4,quadratureOrder=4), ...
                "WVTerrainEnergyGalerkin:InsufficientQuadratureOrder")

            wvt = testCase.problem.originatingTransform;
            [x,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
                (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
            twoDimensional = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
                topographicHeight=10*cos(2*pi*x/wvt.Lx)+10*cos(2*pi*y/wvt.Ly));
            testCase.verifyError(@()twoDimensional.auditHybridPrimitivePVOracle(), ...
                "WVTerrainEnergyGalerkin:HybridOracleRequiresZonallyInvariantTerrain")
        end
    end

    methods (Static,Access=private)
        function problem = sinusoidalProblem(Nz,shouldAntialias,N2)
            z = linspace(-1200,0,Nz)';
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[6 6 Nz], ...
                N2Function=N2,latitude=45,shouldAntialias=shouldAntialias,z=z);
            [~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
                (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
                topographicHeight=20*cos(2*pi*y/wvt.Ly));
        end
    end
end
