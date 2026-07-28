classdef TestWVTerrainEnergyBoundaryDescriptor < matlab.unittest.TestCase
    % Verify the pressure-retaining boundary-dynamical descriptor.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test)
        function flatDescriptorMatchesPublicOracle(testCase)
            problem = TestWVTerrainEnergyBoundaryDescriptor.problemWithHeight(0);
            audit = problem.auditBoundaryDynamicalDescriptor();
            testCase.verifyTrue(audit.isCompatible,string(join(audit.failedConditions,", ")))
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.descriptorResidual,1e-12)
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.continuityResidual,1e-12)
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.bottomResidual,1e-12)
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.surfaceResidual,1e-12)
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.pressureGaugeResidual,1e-12)
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.publicMapResidual,1e-12)
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.bottomCoordinateMapResidual,1e-12)
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.frequencyResidual,1e-11)
            for iK = 1:numel(audit.blocks)
                block = audit.blocks{iK};
                testCase.verifyEqual(block.diagnostics.numberOfFiniteModes,block.diagnostics.expectedFiniteDimension)
                testCase.verifyEqual(block.diagnostics.stationaryDimension,block.diagnostics.expectedStationaryDimension)
                testCase.verifyEqual(block.diagnostics.descriptorRank,block.diagnostics.expectedDescriptorRank)
                testCase.verifyEqual(block.diagnostics.constraintRank,block.diagnostics.expectedConstraintRank)
                testCase.verifyEqual(block.diagnostics.pressureNullity,0)
                testCase.verifyEqual(block.diagnostics.publicCoordinateRank,block.diagnostics.expectedFiniteDimension)
                testCase.verifyEqual(size(block.publicCoordinates,1),nnz(problem.stateLayout.horizontalIndex == iK))
            end
        end

        function uniformDepthMatchesAnalyticFrequencies(testCase)
            problem = TestWVTerrainEnergyBoundaryDescriptor.problemWithHeight(150);
            audit = problem.auditBoundaryDynamicalDescriptor();
            testCase.verifyTrue(audit.isCompatible,string(join(audit.failedConditions,", ")))
            testCase.verifyLessThanOrEqual(audit.diagnostics.maximum.frequencyResidual,1e-11)
        end

        function nonuniformTerrainIsDeferred(testCase)
            N2 = @(z) 2e-5+0*z;
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[4 4 5],N2Function=N2,latitude=45,shouldAntialias=false);
            [X,~] = ndgrid(wvt.x,wvt.y);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=20*cos(2*pi*X/wvt.Lx));
            testCase.verifyError(@()problem.auditBoundaryDynamicalDescriptor(), ...
                "WVTerrainEnergyGalerkin:BoundaryDescriptorRequiresUniformDepth")
        end
    end

    methods (Static)
        function problem = problemWithHeight(height)
            N2 = @(z) 2e-5+0*z;
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[4 4 5],N2Function=N2,latitude=45,shouldAntialias=false);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=height*ones(wvt.Nx,wvt.Ny));
        end
    end
end
