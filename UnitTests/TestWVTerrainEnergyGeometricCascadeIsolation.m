classdef TestWVTerrainEnergyGeometricCascadeIsolation < matlab.unittest.TestCase
    % Verify the Milestone-10.2.2 cascade and isolation audit.

    properties
        problem
        firstAudit
        cachedAudit
        cacheDirectory
        horizontalLayout
        stateLayout
    end

    methods (TestClassSetup)
        function createCascadeAudit(testCase)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(repositoryRoot);
            if ~isfolder(fullfile(repositoryRoot,"output"))
                mkdir(fullfile(repositoryRoot,"output"));
            end
            testCase.cacheDirectory = string(tempname( ...
                fullfile(repositoryRoot,"output")));
            testCase.problem = ...
                TestWVTerrainEnergyGeometricCascadeIsolation.createProblem;
            testCase.horizontalLayout = testCase.problem.horizontalLayout;
            testCase.stateLayout = testCase.problem.stateLayout;
            args = TestWVTerrainEnergyGeometricCascadeIsolation. ...
                reducedArguments(testCase.cacheDirectory);
            testCase.firstAudit = testCase.problem. ...
                auditGeometricCascadeIsolation(args{:});
            testCase.cachedAudit = testCase.problem. ...
                auditGeometricCascadeIsolation(args{:});
        end
    end

    methods (TestClassTeardown)
        function removeCache(testCase)
            if isfolder(testCase.cacheDirectory)
                rmdir(testCase.cacheDirectory,"s");
            end
        end
    end

    methods (Test)
        function outcomeUsesDeclaredClassification(testCase)
            allowed = ["cascade-resolved","cascade-slow", ...
                "resonant-block-required","cascade-nonconvergent"];
            testCase.verifyTrue(ismember( ...
                testCase.firstAudit.classification,allowed))
            testCase.verifyEqual(testCase.firstAudit.status, ...
                testCase.firstAudit.classification)
            testCase.verifyEqual(testCase.firstAudit.nextScope, ...
                "milestone-10.2.2-stop-for-analysis")
        end

        function sinusoidalShellsHaveExactFourierSelection(testCase)
            shells = testCase.firstAudit.zonal.shells;
            testCase.verifyEqual(shells.orders,(0:3).')
            testCase.verifyEqual(sortrows(shells.terrainModes),[0 -1;0 1])
            third = shells.shellModes{4};
            testCase.verifyTrue(any(third(:,2) == -3))
            testCase.verifyTrue(any(third(:,2) == 3))
            testCase.verifyLessThan( ...
                testCase.firstAudit.zonal.tail.maximumStructuralDefect,1e-11)
        end

        function exactFormsAndPhysicalDiagnosticsRemainAvailable(testCase)
            zonal = testCase.firstAudit.zonal;
            testCase.verifyLessThan( ...
                zonal.primary.energyHermitianDefect,1e-12)
            testCase.verifyLessThan( ...
                zonal.primary.exchangeSkewHermitianDefect,1e-12)
            testCase.verifyLessThan(zonal.physical.maximumRitzResidual,1e-10)
            testCase.verifySize(zonal.tail.amplitude,[4 4 5])
            testCase.verifySize(zonal.tail.horizontalShellEnergy,[4 5])
            testCase.verifyEqual(zonal.driftDecomposition.family, ...
                ["stationary";"declared-internal";"non-tangent-bottom"; ...
                "higher-flat-wave";"unresolved"])
            testCase.verifyLessThan( ...
                zonal.driftDecomposition.accountingDefect,1e-11)
        end

        function cacheIsContentAddressedAndResumable(testCase)
            first = testCase.firstAudit.cache;
            second = testCase.cachedAudit.cache;
            testCase.verifyTrue(first.enabled)
            testCase.verifyEqual(first.numberOfHits,0)
            testCase.verifyGreaterThan(first.numberOfMisses,0)
            testCase.verifyEqual(second.numberOfMisses,0)
            testCase.verifyEqual(second.numberOfHits,first.numberOfMisses)
            testCase.verifyGreaterThan(second.secondsSaved,0)
            testCase.verifyEqual(testCase.cachedAudit.classification, ...
                testCase.firstAudit.classification)
            paths = string({second.records.path});
            testCase.verifyTrue(all(startsWith(paths, ...
                testCase.cacheDirectory+filesep)))
            testCase.verifyTrue(all(isfile(paths)))
        end

        function cacheRejectsInvalidSignature(testCase)
            record = testCase.cachedAudit.cache.records(1);
            loaded = load(record.path,"entry");
            entry = loaded.entry;
            entry.schema = "invalid-test-schema";
            save(record.path,"entry","-v7.3");
            args = TestWVTerrainEnergyGeometricCascadeIsolation. ...
                reducedArguments(testCase.cacheDirectory);
            audit = testCase.problem.auditGeometricCascadeIsolation( ...
                args{:});
            testCase.verifyEqual(audit.cache.numberOfInvalidations,1)
            testCase.verifyGreaterThan(audit.cache.numberOfHits,0)
        end

        function publicLayoutsRemainUnchanged(testCase)
            testCase.verifyEqual(testCase.problem.horizontalLayout, ...
                testCase.horizontalLayout)
            testCase.verifyEqual(testCase.problem.stateLayout, ...
                testCase.stateLayout)
        end

        function validationRejectsIncompleteOrPersistentInputs(testCase)
            testCase.verifyError(@() testCase.problem. ...
                auditGeometricCascadeIsolation(scatteringOrders=[1;3;4]), ...
                "WVTerrainEnergyGalerkin:InvalidGeometricCascadeOrders")
            testCase.verifyError(@() testCase.problem. ...
                auditGeometricCascadeIsolation( ...
                primitivePolynomialDegrees=[3;4]), ...
                "WVTerrainEnergyGalerkin:InvalidGeometricCascadeDegrees")
            testCase.verifyError(@() testCase.problem. ...
                auditGeometricCascadeIsolation(cacheDirectory=tempdir), ...
                "WVTerrainEnergyGalerkin:InvalidGeometricCascadeCacheLocation")
        end
    end

    methods (Static)
        function args = reducedArguments(cacheDirectory)
            args = { ...
                "trustedModeBounds",[1 0], ...
                "waveModeIndices",1, ...
                "apvModeIndices",[0;1], ...
                "stationaryPolynomialDegree",2, ...
                "scatteringOrders",[1;2;3], ...
                "primitivePolynomialDegrees",[2;3;4], ...
                "comparisonPolynomialDegree",5, ...
                "paddingFactors",[2;3], ...
                "terrainScales",[0;1/8;1/4;1/2;1], ...
                "shouldRunTwoDimensionalControl",false, ...
                "cacheDirectory",cacheDirectory};
        end

        function problem = createProblem
            N0 = sqrt(2e-5);
            z = linspace(-1200,0,5)';
            wvt = WVTransformBoussinesq( ...
                [24e3 20e3 1200],[4 8 5], ...
                N2Function=@(z)N0^2+0*z,latitude=45, ...
                shouldAntialias=false,z=z);
            [~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
                (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
            problem = WVTerrainEnergyGalerkin.fromTopography( ...
                wvt,topographicHeight=2.5*cos(2*pi*y/wvt.Ly));
        end
    end
end
