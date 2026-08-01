classdef TestWVTerrainEnergyAnalyticGlobalTangentStationaryCompression ...
        < matlab.unittest.TestCase
    % Verify the Milestone-10.2.5.1 analytic attribution oracle.

    properties
        problem
        firstAudit
        cachedAudit
        uncachedAudit
        cacheDirectory
        providerPath
        horizontalLayout
        stateLayout
    end

    methods (TestClassSetup)
        function createAnalyticAttributionAudit(testCase)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(repositoryRoot);
            if ~isfolder(fullfile(repositoryRoot,"output"))
                mkdir(fullfile(repositoryRoot,"output"));
            end
            testCase.cacheDirectory = string(tempname( ...
                fullfile(repositoryRoot,"output")));
            testCase.providerPath = string(fullfile( ...
                repositoryRoot,"..","internal-modes-evp"));
            testCase.problem = createProblem(testCase.providerPath,true);
            testCase.horizontalLayout = testCase.problem.horizontalLayout;
            testCase.stateLayout = testCase.problem.stateLayout;
            withCache = reducedArguments(testCase.cacheDirectory);
            testCase.firstAudit = testCase.problem. ...
                auditAnalyticGlobalTangentStationaryCompression( ...
                withCache{:});
            testCase.cachedAudit = testCase.problem. ...
                auditAnalyticGlobalTangentStationaryCompression( ...
                withCache{:});
            withoutCache = reducedArguments("");
            testCase.uncachedAudit = testCase.problem. ...
                auditAnalyticGlobalTangentStationaryCompression( ...
                withoutCache{:});
        end
    end

    methods (TestClassTeardown)
        function removeCache(testCase)
            if isfolder(testCase.cacheDirectory)
                rmdir(testCase.cacheDirectory,"s");
            end
            if contains(path,testCase.providerPath)
                rmpath(testCase.providerPath)
            end
        end
    end

    methods (Test)
        function analyticFamilyPassesExactControls(testCase)
            control = testCase.firstAudit.training.analyticControl;
            testCase.verifyLessThanOrEqual( ...
                max(cell2mat(struct2cell(control))),1e-12)
            testCase.verifyTrue( ...
                testCase.firstAudit.training.analyticControlPasses)
        end

        function providerErrorDoesNotExplainTheModalDefect(testCase)
            training = testCase.firstAudit.training;
            testCase.verifyLessThan( ...
                max(training.verticalProviderProjectorDefect),1e-8)
            testCase.verifyLessThan( ...
                max(training.propagatedProviderProjectorDefect),1e-8)
            testCase.verifyGreaterThanOrEqual( ...
                training.stationaryProjectorImprovementFactor,4)
            testCase.verifyTrue(training.isMonotonic)
        end

        function approvedSeedClassificationIsRecorded(testCase)
            audit = testCase.firstAudit;
            testCase.verifyEqual(audit.classification, ...
                "analytic-modal-seed")
            testCase.verifyEqual(audit.status,audit.classification)
            testCase.verifyTrue(audit.validation.wasAttempted)
            testCase.verifyFalse(audit.validation.passes)
            testCase.verifyEqual(audit.nextScope, ...
                "milestone-10.2.5.1-stop-for-review")
        end

        function cacheIsOptInAndSignatureChecked(testCase)
            first = testCase.firstAudit.cache;
            second = testCase.cachedAudit.cache;
            disabled = testCase.uncachedAudit.cache;
            testCase.verifyTrue(first.enabled)
            testCase.verifyGreaterThan(first.numberOfMisses,0)
            testCase.verifyEqual(second.numberOfHits,1)
            testCase.verifyEqual(second.numberOfMisses,0)
            testCase.verifyGreaterThan(second.secondsSaved,0)
            testCase.verifyTrue(isfile(second.path))
            testCase.verifyFalse(disabled.enabled)
            testCase.verifyEqual(disabled.path,"")
        end

        function publicLayoutsRemainUnchanged(testCase)
            testCase.verifyEqual(testCase.problem.horizontalLayout, ...
                testCase.horizontalLayout)
            testCase.verifyEqual(testCase.problem.stateLayout, ...
                testCase.stateLayout)
        end

        function invalidInputsAreRejected(testCase)
            testCase.verifyError(@() testCase.problem. ...
                auditAnalyticGlobalTangentStationaryCompression( ...
                cacheDirectory=tempdir), ...
                "WVTerrainEnergyGalerkin:InvalidAnalyticGlobalTangentCache")
            testCase.verifyError(@() testCase.problem. ...
                auditAnalyticGlobalTangentStationaryCompression( ...
                supportModeBounds=[1 1;1 2]), ...
                "WVTerrainEnergyGalerkin:InvalidAnalyticGlobalTangentSupport")
            variableProblem = createProblem(testCase.providerPath,false);
            testCase.verifyError(@() variableProblem. ...
                auditAnalyticGlobalTangentStationaryCompression(), ...
                "WVTerrainEnergyGalerkin:AnalyticGlobalTangentRequiresConstantN")
        end
    end
end

function args = reducedArguments(cacheDirectory)
args = { ...
    "trustedModeBounds",[1 0], ...
    "supportModeBounds",[1 1;1 2;1 3], ...
    "targetWaveModeIndices",1, ...
    "interiorGeostrophicModeCounts",[1;2;3], ...
    "stationaryReferenceDegree",2, ...
    "primitiveReferenceDegrees",[2;3;4], ...
    "paddingFactors",[2;3], ...
    "terrainScales",[0;1/8;1/4;1/2;1], ...
    "comparisonEVPOrder",64, ...
    "cacheDirectory",cacheDirectory};
end

function problem = createProblem(providerPath,isConstant)
N0 = sqrt(2e-5);
z = linspace(-1200,0,5)';
if isConstant
    N2Function = @(zValue)N0^2+0*zValue;
else
    N2Function = @(zValue)N0^2*exp(zValue/600);
end
wvt = WVTransformBoussinesq( ...
    [24e3 20e3 1200],[4 8 5],N2Function=N2Function, ...
    latitude=45,shouldAntialias=false,z=z);
addpath(providerPath,"-end")
[~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
    (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
problem = WVTerrainEnergyGalerkin.fromTopography( ...
    wvt,topographicHeight=2.5*cos(2*pi*y/wvt.Ly));
end
