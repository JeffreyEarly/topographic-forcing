classdef TestWVTerrainEnergyFourierGeometry < matlab.unittest.TestCase
    % Verify use of WaveVortexModel Fourier layouts and transforms.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test)
        function signedLayoutMatchesOriginatingGeometry(testCase)
            resolutions = {[6 4],[8 6]};
            for iResolution = 1:numel(resolutions)
                for shouldAntialias = [false true]
                    problem = TestWVTerrainEnergyFourierGeometry.createProblem(resolutions{iResolution},shouldAntialias);
                    wvt = problem.originatingTransform;
                    layout = problem.horizontalLayout;
                    primary = double(wvt.dftPrimaryIndices2D(:));
                    conjugate = double(wvt.dftConjugateIndices2D(:));
                    expectedDFTIndex = unique([primary;conjugate],"sorted");
                    testCase.verifyEqual(layout.dftIndex,expectedDFTIndex)

                    [expectedPrimary,expectedNativeIndex] = ismember(expectedDFTIndex,primary);
                    [expectedConjugate,conjugateNativeIndex] = ismember(expectedDFTIndex,conjugate);
                    expectedNativeIndex(~expectedPrimary) = conjugateNativeIndex(~expectedPrimary);
                    testCase.verifyTrue(all(expectedPrimary | expectedConjugate))
                    testCase.verifyEqual(layout.isPrimary,expectedPrimary)
                    testCase.verifyEqual(layout.isNative,expectedPrimary)
                    testCase.verifyEqual(layout.nativeIndex,expectedNativeIndex)

                    dftConjugate = double(WVGeometryDoublyPeriodic.indicesOfFourierConjugates(wvt.Nx,wvt.Ny));
                    [isRetained,expectedHorizontalConjugate] = ismember(dftConjugate(expectedDFTIndex),expectedDFTIndex);
                    testCase.verifyTrue(all(isRetained))
                    testCase.verifyEqual(expectedHorizontalConjugate(expectedHorizontalConjugate),(1:height(layout))')

                    actualHorizontalConjugate = zeros(height(layout),1);
                    for iK = 1:height(layout)
                        row = find(problem.stateLayout.horizontalIndex == iK & problem.stateLayout.component == "etaB",1);
                        conjugateRow = problem.conjugateCoordinateIndex(row);
                        actualHorizontalConjugate(iK) = problem.stateLayout.horizontalIndex(conjugateRow);
                    end
                    testCase.verifyEqual(actualHorizontalConjugate,expectedHorizontalConjugate)

                    nyquistMask = logical(WVGeometryDoublyPeriodic.maskForNyquistModes(wvt.Nx,wvt.Ny));
                    testCase.verifyFalse(any(nyquistMask(layout.dftIndex)))
                end
            end
        end

        function coordinateConjugacyIncludesWaveBranchExchange(testCase)
            problem = TestWVTerrainEnergyFourierGeometry.createProblem([6 4],false);
            layout = problem.stateLayout;
            conjugateRows = problem.conjugateCoordinateIndex;
            testCase.verifyEqual(conjugateRows(conjugateRows),(1:height(layout))')
            for row = 1:height(layout)
                conjugateRow = conjugateRows(row);
                expectedComponent = layout.component(row);
                if expectedComponent == "Ap"
                    expectedComponent = "Am";
                elseif expectedComponent == "Am"
                    expectedComponent = "Ap";
                end
                testCase.verifyEqual(layout.component(conjugateRow),expectedComponent)
                if isnan(layout.j(row))
                    testCase.verifyTrue(isnan(layout.j(conjugateRow)))
                else
                    testCase.verifyEqual(layout.j(conjugateRow),layout.j(row))
                end
            end
        end

        function realReconstructionUsesNativeWVOrdering(testCase)
            problem = TestWVTerrainEnergyFourierGeometry.createProblem([6 4],false);
            [Ap,Am,A0,etaB] = TestWVTerrainEnergyFourierGeometry.randomNativeState(problem);
            a = problem.packState(Ap=Ap,Am=Am,A0=A0,bottomDisplacement=etaB);
            spectral = problem.reconstructState(a);
            spatial = problem.reconstructState(a,outputDomain="spatial");
            primaryBlocks = find(problem.horizontalLayout.isPrimary);
            names = ["uHat","vHat","wHat","etaHat","uHatXi","vHatXi","wHatXi","etaHatXi"];
            for name = names
                native = zeros(problem.originatingTransform.Nz,problem.originatingTransform.Nkl);
                native(:,problem.horizontalLayout.nativeIndex(primaryBlocks)) = spectral.(name)(:,primaryBlocks);
                expected = problem.originatingTransform.transformToSpatialDomainWithFourier(native);
                testCase.verifyEqual(spatial.(name),expected,"RelTol",2e-13,"AbsTol",2e-13)
            end
        end

        function realDerivativesMatchOriginatingGeometry(testCase)
            problem = TestWVTerrainEnergyFourierGeometry.createProblem([6 4],false);
            [Ap,Am,A0,etaB] = TestWVTerrainEnergyFourierGeometry.randomNativeState(problem);
            a = problem.packState(Ap=Ap,Am=Am,A0=A0,bottomDisplacement=etaB);
            fields = problem.reconstructState(a,outputDomain="spatial");
            expected = problem.originatingTransform.diffX(fields.vHat)-problem.originatingTransform.diffY(fields.uHat)-problem.originatingTransform.f*fields.etaHatXi;
            actual = problem.evaluateAPV(a);
            testCase.verifyEqual(actual,expected,"RelTol",3e-12,"AbsTol",3e-12)
        end

        function complexTrialModeRetainsPhaseAndDerivative(testCase)
            problem = TestWVTerrainEnergyFourierGeometry.createProblem([6 4],false);
            iK = find(problem.horizontalLayout.kMode == 1 & problem.horizontalLayout.lMode == 0,1);
            rows = find(problem.stateLayout.horizontalIndex == iK);
            row = rows(find(problem.stateLayout.component(rows) == "Ap" & problem.stateLayout.j(rows) == 1,1));
            localRow = find(rows == row);
            a = zeros(height(problem.stateLayout),1);
            a(row) = 1;

            spectral = problem.reconstructState(a);
            spatial = problem.reconstructState(a,outputDomain="spatial");
            [X,Y] = ndgrid(problem.originatingTransform.x,problem.originatingTransform.y);
            phase = exp(1i*(problem.horizontalLayout.k(iK)*X+problem.horizontalLayout.l(iK)*Y));
            expectedU = phase.*shiftdim(problem.basisBlocks{iK}.uHat(:,localRow),-2);
            testCase.verifyEqual(spatial.uHat,expectedU,"RelTol",2e-13,"AbsTol",2e-13)
            testCase.verifyGreaterThan(max(abs(imag(spatial.uHat)),[],"all"),1e-6*max(abs(spatial.uHat),[],"all"))

            qProfile = 1i*problem.horizontalLayout.k(iK)*spectral.vHat(:,iK)-1i*problem.horizontalLayout.l(iK)*spectral.uHat(:,iK)-problem.originatingTransform.f*spectral.etaHatXi(:,iK);
            expectedQ = phase.*shiftdim(qProfile,-2);
            actualQ = problem.evaluateAPV(a);
            testCase.verifyEqual(actualQ,expectedQ,"RelTol",3e-12,"AbsTol",3e-12)
        end
    end

    methods (Static)
        function problem = createProblem(Nxy,shouldAntialias)
            N2 = @(z) 2e-5+0*z;
            wvt = WVTransformBoussinesq([24e3 20e3 1200],[Nxy 7],N2Function=N2,latitude=45,shouldAntialias=shouldAntialias);
            problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(Nxy));
        end

        function [Ap,Am,A0,etaB] = randomNativeState(problem)
            a = zeros(height(problem.stateLayout),1);
            primary = problem.stateLayout.isPrimary;
            a(primary) = randn(nnz(primary),1)+1i*randn(nnz(primary),1);
            a(~primary) = conj(a(problem.conjugateCoordinateIndex(~primary)));
            zeroK = find(problem.horizontalLayout.kMode == 0 & problem.horizontalLayout.lMode == 0,1);
            zeroRows = find(problem.stateLayout.horizontalIndex == zeroK);
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
