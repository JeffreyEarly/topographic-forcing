classdef TestWVTerrainEnergyGalerkinContract < matlab.unittest.TestCase
    % Verify the public terrain-energy Galerkin contract.

    methods (TestClassSetup)
        function addRepositoryToPath(~)
            addpath(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (Test)
        function validTerrainAndAntialiasConventions(testCase)
            for shouldAntialias = [false true]
                wvt = TestWVTerrainEnergyGalerkinContract.createTransform(shouldAntialias);
                [X,Y] = ndgrid(wvt.x,wvt.y);
                terrains = {
                    zeros(wvt.Nx,wvt.Ny)
                    50*ones(wvt.Nx,wvt.Ny)
                    25*cos(2*pi*X/wvt.Lx)
                    15*cos(2*pi*X/wvt.Lx)+10*sin(2*pi*Y/wvt.Ly)
                    };
                for iTerrain = 1:numel(terrains)
                    Ap = wvt.Ap;
                    Am = wvt.Am;
                    A0 = wvt.A0;
                    t = wvt.t;
                    forcing = wvt.forcing;
                    problem = WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=terrains{iTerrain});
                    testCase.verifyClass(problem,"WVTerrainEnergyGalerkin")
                    testCase.verifyEqual(problem.topographicHeight,terrains{iTerrain})
                    testCase.verifyEqual(problem.gamma,1-terrains{iTerrain}/wvt.Lz)
                    testCase.verifyEqual(problem.hydrostaticTransform.shouldAntialias,shouldAntialias)
                    testCase.verifyEqual(problem.originatingTransform,wvt)
                    testCase.verifyEqual(wvt.Ap,Ap)
                    testCase.verifyEqual(wvt.Am,Am)
                    testCase.verifyEqual(wvt.A0,A0)
                    testCase.verifyEqual(wvt.t,t)
                    testCase.verifyEqual(wvt.forcing,forcing)
                    testCase.verifyFalse(isa(problem,"WVForcing"))
                end
            end
        end

        function validationErrorsAreStructured(testCase)
            wvt = TestWVTerrainEnergyGalerkinContract.createTransform(false);
            h = zeros(wvt.Nx,wvt.Ny);
            hydro = WVTransformHydrostatic([wvt.Lx wvt.Ly wvt.Lz],[wvt.Nx wvt.Ny wvt.Nz],N2Function=wvt.N2Function,latitude=wvt.latitude,shouldAntialias=false);
            testCase.verifyError(@()WVTerrainEnergyGalerkin.fromTopography(hydro,topographicHeight=h),"WVTerrainEnergyGalerkin:UnsupportedTransform")
            testCase.verifyError(@()WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=zeros(wvt.Nx-1,wvt.Ny)),"WVTerrainEnergyGalerkin:InvalidTopographySize")
            testCase.verifyError(@()WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=1i*h+1i),"WVTerrainEnergyGalerkin:ComplexTopography")
            hNonfinite = h;
            hNonfinite(1) = NaN;
            testCase.verifyError(@()WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=hNonfinite),"WVTerrainEnergyGalerkin:NonfiniteTopography")
            testCase.verifyError(@()WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=wvt.Lz*ones(size(h))),"WVTerrainEnergyGalerkin:NonpositiveGamma")
            testCase.verifyError(@()WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h,verticalModeIndices=[0; 2; 1]),"WVTerrainEnergyGalerkin:InvalidVerticalModeIndices")
            testCase.verifyError(@()WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h,verticalModeIndices=[1; 2]),"WVTerrainEnergyGalerkin:IncompleteVerticalModeSet")
            testCase.verifyError(@()WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h,verticalModeIndices=0),"WVTerrainEnergyGalerkin:IncompleteVerticalModeSet")
            testCase.verifyError(@()WVTerrainEnergyGalerkin.fromTopography(wvt,topographicHeight=h,horizontalOversamplingFactor=1.5),"WVTerrainEnergyGalerkin:InvalidOversamplingFactor")
        end
    end

    methods (Static)
        function wvt = createTransform(shouldAntialias)
            N2 = @(z) 2e-5+0*z;
            wvt = WVTransformBoussinesq([20e3 20e3 1000],[4 4 5],N2Function=N2,latitude=45,shouldAntialias=shouldAntialias);
        end
    end
end
