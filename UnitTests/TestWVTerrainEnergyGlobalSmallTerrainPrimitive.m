classdef TestWVTerrainEnergyGlobalSmallTerrainPrimitive < matlab.unittest.TestCase
    % Verify the global small-terrain primitive compatibility oracle.

    properties
        problem
    end

    methods (TestClassSetup)
        function createProblem(testCase)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            N2 = @(z) 2e-5+0*z;
            z = linspace(-1200,0,5)';
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[6 6 5], ...
                N2Function=N2,latitude=45,shouldAntialias=false,z=z);
            [~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
                (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
            h = 20*cos(2*pi*y/wvt.Ly);
            testCase.problem = WVTerrainEnergyGalerkin.fromTopography( ...
                wvt,topographicHeight=h);
        end
    end

    methods (Test)
        function independentTangentsEstablishRepresentationBlocker(testCase)
            audit = testCase.problem.auditGlobalSmallTerrainPrimitive(polynomialDegree=2);
            testCase.verifyEqual(audit.status,"representation-blocker")
            testCase.verifyFalse(audit.isCompatible)
            testCase.verifyLessThan(audit.tangentAgreement.maximumRelativeDefect,1e-9)
            testCase.verifyLessThan(audit.compatibility.weakEvolutionDefect,1e-11)
            testCase.verifyLessThan(audit.compatibility.energyDefect,1e-11)
            testCase.verifyLessThan(audit.compatibility.bottomDefect,1e-11)
            testCase.verifyLessThan(audit.compatibility.conjugacyDefect,1e-11)
            testCase.verifyGreaterThan(audit.compatibility.apvDefect,1e-2)
            testCase.verifyGreaterThan(audit.compatibility.enstrophyDefect,1e-4)
        end

        function globalCouplingHasCorrectSelectionAndNontrivialCancellation(testCase)
            audit = testCase.problem.auditGlobalSmallTerrainPrimitive(polynomialDegree=2);
            testCase.verifyLessThan(audit.fourier.couplingLeakage,1e-12)
            testCase.verifyGreaterThan(audit.fourier.numberOfInteriorHorizontalModes,0)
            testCase.verifyLessThan(audit.fourier.interiorAPVDefect,audit.fourier.edgeAPVDefect)
            testCase.verifyGreaterThan(norm(audit.apvCancellation.tendencyCorrection,"fro"),0)
            testCase.verifyGreaterThan(norm(audit.apvCancellation.mapCorrection,"fro"),0)
            testCase.verifyGreaterThan(norm(audit.apvCancellation.total,"fro"),0)
            testCase.verifyEqual(sortrows(audit.fourier.terrainModes),[0 -1;0 1])
        end

        function verticalRefinementSeparatesInteriorAndTruncationDefects(testCase)
            degree = [2 4 6];
            interior = zeros(size(degree));
            edge = zeros(size(degree));
            dispersion = zeros(size(degree));
            for iDegree = 1:numel(degree)
                audit = testCase.problem.auditGlobalSmallTerrainPrimitive(polynomialDegree=degree(iDegree));
                testCase.verifyEqual(audit.status,"representation-blocker")
                interior(iDegree) = audit.fourier.interiorAPVDefect;
                edge(iDegree) = audit.fourier.edgeAPVDefect;
                dispersion(iDegree) = audit.flatDiagnostics.modeOneDispersionDefect;
            end
            testCase.verifyTrue(all(diff(interior) < 0))
            testCase.verifyGreaterThan(min(edge),0.1)
            testCase.verifyTrue(all(diff(dispersion) < 0))
            testCase.verifyLessThan(dispersion(end),1e-8)
        end

        function arbitraryStratificationRetainsClassification(testCase)
            N2 = @(z) 1.2e-5*exp(z/1800);
            z = linspace(-1200,0,5)';
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[4 4 5], ...
                N2Function=N2,latitude=45,shouldAntialias=false,z=z);
            [~,y] = ndgrid((0:wvt.Nx-1)'*wvt.Lx/wvt.Nx, ...
                (0:wvt.Ny-1)'*wvt.Ly/wvt.Ny);
            exponentialProblem = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
                topographicHeight=20*cos(2*pi*y/wvt.Ly));
            audit = exponentialProblem.auditGlobalSmallTerrainPrimitive(polynomialDegree=3);
            testCase.verifyEqual(audit.status,"representation-blocker")
            testCase.verifyLessThan(audit.tangentAgreement.maximumRelativeDefect,1e-9)
            testCase.verifyLessThan(audit.compatibility.energyDefect,1e-11)
            testCase.verifyLessThan(audit.compatibility.bottomDefect,1e-11)
            testCase.verifyGreaterThan(audit.compatibility.apvDefect,1e-2)
        end

        function auditPreservesLayoutAndRejectsInvalidInputs(testCase)
            horizontalBefore = testCase.problem.horizontalLayout;
            stateBefore = testCase.problem.stateLayout;
            conjugateBefore = testCase.problem.conjugateCoordinateIndex;
            testCase.problem.auditGlobalSmallTerrainPrimitive(polynomialDegree=2);
            testCase.verifyEqual(testCase.problem.horizontalLayout,horizontalBefore)
            testCase.verifyEqual(testCase.problem.stateLayout,stateBefore)
            testCase.verifyEqual(testCase.problem.conjugateCoordinateIndex,conjugateBefore)
            testCase.verifyError(@()testCase.problem.auditGlobalSmallTerrainPrimitive( ...
                polynomialDegree=1),"WVTerrainEnergyGalerkin:InvalidPolynomialDegree")
            testCase.verifyError(@()testCase.problem.auditGlobalSmallTerrainPrimitive( ...
                polynomialDegree=4,quadratureOrder=4), ...
                "WVTerrainEnergyGalerkin:InsufficientQuadratureOrder")
            testCase.verifyError(@()testCase.problem.auditGlobalSmallTerrainPrimitive( ...
                tangentStep=0.2),"WVTerrainEnergyGalerkin:InvalidTangentStep")

            wvt = testCase.problem.originatingTransform;
            flat = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
                topographicHeight=zeros(wvt.Nx,wvt.Ny));
            testCase.verifyError(@()flat.auditGlobalSmallTerrainPrimitive(), ...
                "WVTerrainEnergyGalerkin:GlobalPrimitiveAuditRequiresTerrainDirection")
        end
    end
end
