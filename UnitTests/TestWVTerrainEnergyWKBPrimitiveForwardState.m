classdef TestWVTerrainEnergyWKBPrimitiveForwardState < matlab.unittest.TestCase
    % Verify the Milestone-10.2.6 complete WKB primitive oracle.

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
        function createWKBPrimitiveAudit(testCase)
            repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(repositoryRoot);
            if ~isfolder(fullfile(repositoryRoot,"output"))
                mkdir(fullfile(repositoryRoot,"output"));
            end
            testCase.cacheDirectory = string(tempname( ...
                fullfile(repositoryRoot,"output")));
            testCase.providerPath = string(fullfile( ...
                repositoryRoot,"..","internal-modes-evp"));
            addpath(testCase.providerPath,"-end")
            testCase.problem = createProblem;
            testCase.horizontalLayout = testCase.problem.horizontalLayout;
            testCase.stateLayout = testCase.problem.stateLayout;
            withCache = reducedArguments(testCase.cacheDirectory);
            testCase.firstAudit = testCase.problem. ...
                auditWKBPrimitiveForwardState(withCache{:});
            testCase.cachedAudit = testCase.problem. ...
                auditWKBPrimitiveForwardState(withCache{:});
            withoutCache = reducedArguments("");
            testCase.uncachedAudit = testCase.problem. ...
                auditWKBPrimitiveForwardState(withoutCache{:});
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
        function wkbCoordinateAndWeakFormsPassExactControls(testCase)
            audit = testCase.firstAudit;
            wkb = audit.wkbDegreeConvergence(end);
            testCase.verifyTrue(wkb.isQualified)
            testCase.verifyLessThanOrEqual( ...
                wkb.verticalDiagnostics.forwardMapDefect,1e-12)
            testCase.verifyLessThanOrEqual( ...
                wkb.verticalDiagnostics.inverseMapDefect,1e-12)
            testCase.verifyLessThanOrEqual( ...
                wkb.verticalDiagnostics.jacobianDefect,1e-12)
            testCase.verifyLessThanOrEqual( ...
                wkb.verticalDiagnostics.endpointDefect,1e-12)
            testCase.verifyLessThanOrEqual( ...
                wkb.quadratureAdjointDefect,1e-12)
            testCase.verifyLessThanOrEqual(wkb.continuityDefect,1e-11)
            testCase.verifyLessThanOrEqual(wkb.conjugacyDefect,1e-11)
            testCase.verifyLessThanOrEqual( ...
                wkb.maximumStructuralDefect,1e-12)
            testCase.verifyLessThanOrEqual(wkb.maximumEnergyDefect,1e-10)
        end

        function constantNSpacesAreEquivalent(testCase)
            equivalence = testCase.firstAudit.constantNEquivalence;
            testCase.verifyTrue(equivalence.passes)
            testCase.verifyLessThanOrEqual( ...
                equivalence.maximumProjectorDefect,1e-10)
        end

        function auditReportsCompleteMatchedAccuracyEvidence(testCase)
            audit = testCase.firstAudit;
            testCase.verifyTrue(ismember(audit.classification,[ ...
                "wkb-spectral-acceleration","wkb-spectral-equivalent", ...
                "wkb-spectral-fallback","complete-spectral-blocker"]))
            testCase.verifyGreaterThan( ...
                audit.wkbDegreeConvergence(end).numberOfAdmissibleCoordinates,0)
            testCase.verifyGreaterThan( ...
                audit.legendreDegreeConvergence(end).numberOfAdmissibleCoordinates,0)
            testCase.verifyGreaterThan( ...
                audit.wkbDegreeConvergence(end).buildSeconds,0)
            testCase.verifyGreaterThan( ...
                audit.legendreDegreeConvergence(end).buildSeconds,0)
            testCase.verifyEqual(numel(audit.standardizedControls),4)
            testCase.verifyTrue(isfield(audit.twoDimensionalControls, ...
                "legendre"))
            testCase.verifyEqual(audit.nextScope, ...
                "milestone-10.3-only-after-review")
        end

        function cacheIsOptInContentAddressedAndResumable(testCase)
            first = testCase.firstAudit.cache;
            second = testCase.cachedAudit.cache;
            disabled = testCase.uncachedAudit.cache;
            testCase.verifyTrue(first.enabled)
            testCase.verifyGreaterThan(first.numberOfMisses,0)
            testCase.verifyEqual(second.numberOfMisses,0)
            testCase.verifyGreaterThan(second.numberOfHits,0)
            testCase.verifyGreaterThan(second.secondsSaved,0)
            testCase.verifyFalse(disabled.enabled)
            testCase.verifyTrue(all([disabled.records.status] == "disabled"))
            testCase.verifyTrue(all([disabled.records.path] == ""))
        end

        function publicLayoutsRemainUnchanged(testCase)
            testCase.verifyEqual(testCase.problem.horizontalLayout, ...
                testCase.horizontalLayout)
            testCase.verifyEqual(testCase.problem.stateLayout, ...
                testCase.stateLayout)
        end

        function invalidControlsAreRejected(testCase)
            testCase.verifyError(@() testCase.problem. ...
                auditWKBPrimitiveForwardState(cacheDirectory=tempdir), ...
                "WVTerrainEnergyGalerkin:InvalidWKBPrimitiveCache")
            testCase.verifyError(@() testCase.problem. ...
                auditWKBPrimitiveForwardState( ...
                supportModeBounds=[1 1]), ...
                "WVTerrainEnergyGalerkin:InvalidWKBPrimitiveSupport")
            testCase.verifyError(@() testCase.problem. ...
                auditWKBPrimitiveForwardState(verticalDegrees=[2;3;4], ...
                legendreReferenceDegrees=[5;6;7]), ...
                "WVTerrainEnergyGalerkin:MissingWKBPrimitiveCommonDegree")
        end
    end
end

function args = reducedArguments(cacheDirectory)
args = { ...
    "trustedModeBounds",[0 0], ...
    "supportModeBounds",[0 0;0 1;0 2], ...
    "verticalDegrees",[2;3;4], ...
    "legendreReferenceDegrees",[2;3;4], ...
    "paddingFactors",[2;3], ...
    "terrainScales",[0;1/2;1], ...
    "cacheDirectory",cacheDirectory};
end

function problem = createProblem
N0 = sqrt(2e-5);
z = linspace(-1200,0,5)';
wvt = WVTransformBoussinesq( ...
    [24e3 20e3 1200],[4 6 5],N2Function=@(zValue)N0^2+0*zValue, ...
    latitude=45,shouldAntialias=false,z=z);
[~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
    (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
problem = WVTerrainEnergyGalerkin.fromTopography( ...
    wvt,topographicHeight=0.2*cos(2*pi*y/wvt.Ly));
end
