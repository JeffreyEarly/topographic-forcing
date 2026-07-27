classdef TestWVTerrainEnergyMixedBasis < matlab.unittest.TestCase
    % Verify mixed-coordinate reconstruction, conjugacy, and kinematics.

    properties
        wvt
        problem
    end

    methods (TestClassSetup)
        function createProblem(testCase)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            N2 = @(z) 2e-5+0*z;
            testCase.wvt = WVTransformBoussinesq([24e3 20e3 1200],[6 4 7],N2Function=N2,latitude=45,shouldAntialias=false);
            [X,Y] = ndgrid(testCase.wvt.x,testCase.wvt.y);
            h = 40*cos(2*pi*X/testCase.wvt.Lx)+20*sin(2*pi*Y/testCase.wvt.Ly);
            testCase.problem = WVTerrainEnergyGalerkin.fromTopography(testCase.wvt,topographicHeight=h);
        end
    end

    methods (Test)
        function layoutContainsExpectedCoordinates(testCase)
            problem = testCase.problem;
            Nj = numel(problem.verticalModeIndices);
            for iK = 1:height(problem.horizontalLayout)
                rows = problem.stateLayout.horizontalIndex == iK;
                if problem.horizontalLayout.k(iK)^2+problem.horizontalLayout.l(iK)^2 > 0
                    testCase.verifyEqual(nnz(rows),3*Nj-1)
                    testCase.verifyEqual(nnz(rows & problem.stateLayout.component == "etaB"),1)
                else
                    testCase.verifyEqual(nnz(rows),3*Nj)
                end
            end
            testCase.verifyEqual(problem.conjugateCoordinateIndex(problem.conjugateCoordinateIndex),(1:height(problem.stateLayout))')
        end

        function bottomFunctionAndVelocityConstraints(testCase)
            problem = testCase.problem;
            a = randn(height(problem.stateLayout),1)+1i*randn(height(problem.stateLayout),1);
            fields = problem.reconstructState(a);
            continuity = zeros(size(fields.uHat));
            for iK = 1:height(problem.horizontalLayout)
                continuity(:,iK) = 1i*problem.horizontalLayout.k(iK)*fields.uHat(:,iK)+1i*problem.horizontalLayout.l(iK)*fields.vHat(:,iK)+fields.wHatXi(:,iK);
                block = problem.basisBlocks{iK};
                rows = find(problem.stateLayout.horizontalIndex == iK & problem.stateLayout.component == "etaB");
                localRow = find(find(problem.stateLayout.horizontalIndex == iK) == rows);
                testCase.verifyEqual(block.etaHat(1,localRow),1,"AbsTol",2e-14)
                testCase.verifyEqual(block.etaHat(end,localRow),0,"AbsTol",2e-14)
            end
            scale = max(1,norm([fields.uHat(:);fields.vHat(:);fields.wHatXi(:)]));
            testCase.verifyLessThanOrEqual(norm(continuity(:))/scale,2e-12)
            testCase.verifyLessThanOrEqual(norm(fields.wHat([1 end],:),"fro")/max(1,norm(fields.wHat,"fro")),2e-12)
        end

        function packingRoundTripAndSpatialReality(testCase)
            [Ap,Am,A0,etaB] = TestWVTerrainEnergyMixedBasis.randomNativeState(testCase.problem);
            a = testCase.problem.packState(Ap=Ap,Am=Am,A0=A0,bottomDisplacement=etaB);
            [ApOut,AmOut,A0Out,etaBOut] = testCase.problem.unpackState(a);
            testCase.verifyEqual(ApOut,Ap,"AbsTol",1e-14)
            testCase.verifyEqual(AmOut,Am,"AbsTol",1e-14)
            testCase.verifyEqual(A0Out,A0,"AbsTol",1e-14)
            testCase.verifyEqual(etaBOut,etaB,"AbsTol",1e-14)
            fields = testCase.problem.reconstructState(a,outputDomain="spatial");
            names = ["uHat","vHat","wHat","etaHat","u","v","w"];
            for name = names
                testCase.verifyLessThanOrEqual(max(abs(imag(fields.(name))),[],"all"),1e-12*max(1,max(abs(fields.(name)),[],"all")))
            end

            bad = a;
            bad(1) = bad(1)+1i;
            testCase.verifyError(@()testCase.problem.unpackState(bad),"WVTerrainEnergyGalerkin:ConjugacyViolation")
        end

        function reconstructionAdjointIdentity(testCase)
            nState = height(testCase.problem.stateLayout);
            nK = height(testCase.problem.horizontalLayout);
            Nz = testCase.wvt.Nz;
            a = randn(nState,1)+1i*randn(nState,1);
            fields = testCase.problem.reconstructState(a);
            testFields = struct;
            names = ["uHat","vHat","wHat","etaHat"];
            for name = names
                testFields.(name) = randn(Nz,nK)+1i*randn(Nz,nK);
            end
            dual = testCase.problem.applyReconstructionAdjoint(testFields);
            lhs = 0;
            for name = names
                lhs = lhs+sum(conj(fields.(name)).*testFields.(name),"all");
            end
            rhs = a'*dual;
            testCase.verifyEqual(lhs,rhs,"RelTol",2e-12,"AbsTol",2e-12)
        end

        function interiorFlatEnergyMatchesDirectQuadrature(testCase)
            [Ap,Am,A0,etaB] = TestWVTerrainEnergyMixedBasis.randomNativeState(testCase.problem);
            a = testCase.problem.packState(Ap=Ap,Am=Am,A0=A0,bottomDisplacement=etaB);
            a(testCase.problem.stateLayout.component == "etaB") = 0;
            fields = testCase.problem.reconstructState(a,outputDomain="spatial");
            N2 = shiftdim(testCase.wvt.N2,-2);
            energyDensity = abs(fields.uHat).^2+abs(fields.vHat).^2+abs(fields.wHat).^2+N2.*abs(fields.etaHat).^2;
            directEnergy = testCase.wvt.rho0*sum(testCase.wvt.z_int.*squeeze(mean(mean(energyDensity,1),2)))/2;
            modalEnergy = 0;
            for iK = 1:height(testCase.problem.horizontalLayout)
                rows = find(testCase.problem.stateLayout.horizontalIndex == iK);
                modalEnergy = modalEnergy+real(a(rows)'*testCase.problem.flatModeBlocks{iK}.E*a(rows))/2;
            end
            testCase.verifyEqual(modalEnergy,directEnergy,"RelTol",2e-12)
        end

        function physicalBottomKinematicsAndDisplacementEquivalence(testCase)
            [Ap,Am,A0,etaB] = TestWVTerrainEnergyMixedBasis.randomNativeState(testCase.problem);
            a = testCase.problem.packState(Ap=Ap,Am=Am,A0=A0,bottomDisplacement=etaB);
            fields = testCase.problem.reconstructState(a,outputDomain="spatial");
            hx = testCase.wvt.diffX(testCase.problem.topographicHeight);
            hy = testCase.wvt.diffY(testCase.problem.topographicHeight);
            expectedBottomW = fields.u(:,:,1).*hx+fields.v(:,:,1).*hy;
            testCase.verifyLessThanOrEqual(norm(fields.w(:,:,1)-expectedBottomW,"fro")/max(1,norm(expectedBottomW,"fro")),2e-12)

            gamma = testCase.problem.gamma;
            gammaX = testCase.wvt.diffX(gamma);
            gammaY = testCase.wvt.diffY(gamma);
            xi = shiftdim(testCase.wvt.z,-2);
            etaITendency = (fields.wHat+xi.*fields.uHat.*(gammaX./gamma)+xi.*fields.vHat.*(gammaY./gamma));
            testCase.verifyLessThanOrEqual(norm(etaITendency-fields.w,"fro")/max(1,norm(fields.w,"fro")),2e-12)
            testCase.verifyLessThanOrEqual(norm(etaITendency(:,:,1)-expectedBottomW,"fro")/max(1,norm(expectedBottomW,"fro")),2e-12)
        end
    end

    methods (Static)
        function [Ap,Am,A0,etaB] = randomNativeState(problem)
            a = zeros(height(problem.stateLayout),1);
            primary = problem.stateLayout.isPrimary;
            a(primary) = randn(nnz(primary),1)+1i*randn(nnz(primary),1);
            a(~primary) = conj(a(problem.conjugateCoordinateIndex(~primary)));
            zeroK = find(problem.horizontalLayout.kMode == 0 & problem.horizontalLayout.lMode == 0,1);
            zeroRows = find(problem.stateLayout.horizontalIndex == zeroK);
            a(zeroRows) = 0;
            plusRows = zeroRows(problem.stateLayout.component(zeroRows) == "Ap");
            minusRows = zeroRows(problem.stateLayout.component(zeroRows) == "Am");
            a(plusRows) = randn(numel(plusRows),1)+1i*randn(numel(plusRows),1);
            a(minusRows) = conj(a(plusRows));
            realRows = setdiff(zeroRows,[plusRows;minusRows]);
            a(realRows) = randn(numel(realRows),1);
            [Ap,Am,A0,etaB] = problem.unpackState(a);
        end
    end
end
