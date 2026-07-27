classdef WVTerrainEnergyGalerkin < handle
    % Construct a pressure-free terrain-energy Galerkin problem.
    %
    % `WVTerrainEnergyGalerkin` uses hydrostatic wave-vortex modes as
    % coordinates for the complete nonhydrostatic weak equations. One
    % additional displacement coordinate per horizontal wavenumber carries
    % the bottom value. The class is a standalone linear scientific system;
    % it is not a `WVForcing` and does not modify its originating transform.
    % The finite-terrain energy, exchange, and APV forms are evaluated by
    % common oversampled quadrature.
    %
    % ```matlab
    % problem = WVTerrainEnergyGalerkin.fromTopography(wvt, ...
    %     topographicHeight=h);
    % ```
    %
    % - Topic: Create a Galerkin problem
    % - Topic: Inspect the mixed basis
    % - Topic: Transform Galerkin states
    % - Topic: Inspect flat modes
    % - Topic: Inspect terrain forms
    % - Topic: Audit constrained evolution
    % - Declaration: classdef WVTerrainEnergyGalerkin < handle

    properties (SetAccess=private)
        % Originating nonhydrostatic wave-vortex transform.
        %
        % - Topic: Inspect the mixed basis
        originatingTransform

        % Matched hydrostatic transform supplying the reference coordinates.
        %
        % - Topic: Inspect the mixed basis
        hydrostaticTransform

        % Upward-positive stationary topographic height in meters.
        %
        % - Topic: Inspect the mixed basis
        topographicHeight

        % Terrain coordinate scale factor, $$\gamma=1-h/D$$.
        %
        % - Topic: Inspect the mixed basis
        gamma

        % Retained hydrostatic vertical-mode labels.
        %
        % - Topic: Inspect the mixed basis
        verticalModeIndices

        % Requested horizontal oversampling factor for terrain products.
        %
        % The finite-terrain forms use this factor in both horizontal
        % directions for their common quadrature grid.
        %
        % - Topic: Inspect the mixed basis
        horizontalOversamplingFactor

        % Metadata for the retained full-complex horizontal Fourier grid.
        %
        % - Topic: Inspect the mixed basis
        horizontalLayout

        % Metadata for every mixed Galerkin coordinate.
        %
        % - Topic: Inspect the mixed basis
        stateLayout

        % Coordinate index of the real-field conjugate of each coordinate.
        %
        % - Topic: Inspect the mixed basis
        conjugateCoordinateIndex

        % Per-wavenumber hydrostatic coordinate blocks.
        %
        % Each entry contains vertical profiles for the hatted fields and
        % their required vertical derivatives.
        %
        % - Topic: Inspect the mixed basis
        basisBlocks

        % Per-wavenumber flat nonhydrostatic oracle results.
        %
        % - Topic: Inspect flat modes
        flatModeBlocks

        % Dense finite-terrain energy, exchange, and APV forms.
        %
        % `finiteTerrainForms.energyMatrix` is the positive Hermitian
        % representation of the finite-terrain energy inner product,
        % `exchangeMatrix` is the skew-Hermitian Coriolis--buoyancy form,
        % and `apvMatrix` maps mixed coefficients to quadrature-weighted
        % finite-terrain APV samples.
        %
        % - Topic: Inspect terrain forms
        finiteTerrainForms

        % Construction and scientific-gate diagnostics.
        %
        % - Topic: Inspect flat modes
        constructionDiagnostics
    end

    properties (Access=private)
        % Horizontal index of the retained Fourier conjugate.
        horizontalConjugateIndex
    end

    methods (Static)
        function self = fromTopography(wvt,options)
            % Construct a terrain-energy Galerkin problem from topography.
            %
            % The retained vertical modes must include mode zero and at
            % least one internal mode. The source transform's antialiasing
            % convention is inherited without alteration.
            %
            % - Topic: Create a Galerkin problem
            % - Declaration: problem = WVTerrainEnergyGalerkin.fromTopography(wvt,options)
            % - Parameter wvt: originating `WVTransformBoussinesq`
            % - Parameter topographicHeight: real `Nx`-by-`Ny` terrain in meters
            % - Parameter verticalModeIndices: retained labels from `wvt.j`
            % - Parameter horizontalOversamplingFactor: positive integer, default 2
            % - Returns problem: constructed `WVTerrainEnergyGalerkin`
            arguments
                wvt (1,1) WVTransform
                options.topographicHeight (:,:) double
                options.verticalModeIndices (:,1) double = wvt.j
                options.horizontalOversamplingFactor (1,1) double = 2
            end

            if ~isa(wvt,"WVTransformBoussinesq")
                error("WVTerrainEnergyGalerkin:UnsupportedTransform", ...
                    "The originating transform must be a WVTransformBoussinesq.")
            end

            h = options.topographicHeight;
            if ~isequal(size(h),[wvt.Nx wvt.Ny])
                error("WVTerrainEnergyGalerkin:InvalidTopographySize", ...
                    "topographicHeight must have size [%d %d].",wvt.Nx,wvt.Ny)
            end
            if ~isreal(h)
                error("WVTerrainEnergyGalerkin:ComplexTopography", ...
                    "topographicHeight must be real.")
            end
            if any(~isfinite(h),"all")
                error("WVTerrainEnergyGalerkin:NonfiniteTopography", ...
                    "topographicHeight must contain only finite values.")
            end

            gamma = 1-h/wvt.Lz;
            if any(gamma <= 0,"all")
                error("WVTerrainEnergyGalerkin:NonpositiveGamma", ...
                    "The terrain must satisfy gamma=1-h/D>0 at every grid point.")
            end

            j = options.verticalModeIndices;
            if any(~isfinite(j)) || any(j ~= round(j)) || numel(unique(j)) ~= numel(j) || ~issorted(j) || any(~ismember(j,wvt.j))
                error("WVTerrainEnergyGalerkin:InvalidVerticalModeIndices", ...
                    "verticalModeIndices must be sorted, unique integer labels present in wvt.j.")
            end
            if ~ismember(0,j) || ~any(j > 0)
                error("WVTerrainEnergyGalerkin:IncompleteVerticalModeSet", ...
                    "verticalModeIndices must include mode zero and at least one internal mode.")
            end

            oversampling = options.horizontalOversamplingFactor;
            if ~isfinite(oversampling) || oversampling ~= round(oversampling) || oversampling < 1
                error("WVTerrainEnergyGalerkin:InvalidOversamplingFactor", ...
                    "horizontalOversamplingFactor must be a positive integer.")
            end

            self = WVTerrainEnergyGalerkin(wvt,h,gamma,j,oversampling);
        end
    end

    methods
        function a = packState(self,options)
            % Pack WaveVortexModel coefficients into the full-complex state.
            %
            % `bottomDisplacement` uses the originating transform's
            % nonredundant horizontal layout and has units of meters.
            %
            % - Topic: Transform Galerkin states
            % - Declaration: a = packState(Ap=Ap,Am=Am,A0=A0,bottomDisplacement=etaB)
            % - Parameter Ap: positive-wave coefficients in native layout
            % - Parameter Am: negative-wave coefficients in native layout
            % - Parameter A0: balanced coefficients in native layout
            % - Parameter bottomDisplacement: bottom values in native horizontal layout
            % - Returns a: full-complex mixed Galerkin vector
            arguments
                self (1,1) WVTerrainEnergyGalerkin
                options.Ap (:,:) double = zeros(self.originatingTransform.spectralMatrixSize)
                options.Am (:,:) double = zeros(self.originatingTransform.spectralMatrixSize)
                options.A0 (:,:) double = zeros(self.originatingTransform.spectralMatrixSize)
                options.bottomDisplacement (:,:) double = zeros(1,self.originatingTransform.Nkl)
            end

            wvt = self.originatingTransform;
            if ~isequal(size(options.Ap),wvt.spectralMatrixSize) || ~isequal(size(options.Am),wvt.spectralMatrixSize) || ~isequal(size(options.A0),wvt.spectralMatrixSize)
                error("WVTerrainEnergyGalerkin:InvalidNativeCoefficientSize", ...
                    "Ap, Am, and A0 must have size [%d %d].",wvt.spectralMatrixSize(1),wvt.spectralMatrixSize(2))
            end
            if ~isequal(size(options.bottomDisplacement),[1 wvt.Nkl])
                error("WVTerrainEnergyGalerkin:InvalidBottomCoefficientSize", ...
                    "bottomDisplacement must have size [1 %d].",wvt.Nkl)
            end

            a = zeros(height(self.stateLayout),1);
            primaryBlocks = find(self.horizontalLayout.isPrimary);
            for iBlock = reshape(primaryBlocks,1,[])
                iNative = self.horizontalLayout.nativeIndex(iBlock);
                rows = find(self.stateLayout.horizontalIndex == iBlock);
                labels = self.stateLayout.component(rows);
                for iRow = 1:numel(rows)
                    row = rows(iRow);
                    jLabel = self.stateLayout.j(row);
                    if labels(iRow) == "Ap"
                        a(row) = options.Ap(wvt.j == jLabel,iNative);
                    elseif labels(iRow) == "Am"
                        a(row) = options.Am(wvt.j == jLabel,iNative);
                    elseif labels(iRow) == "A0"
                        a(row) = options.A0(wvt.j == jLabel,iNative);
                    else
                        a(row) = options.bottomDisplacement(iNative);
                    end
                end
            end
            unfilled = ~self.stateLayout.isPrimary;
            a(unfilled) = conj(a(self.conjugateCoordinateIndex(unfilled)));

            self.validateConjugacy(a);
        end

        function [Ap,Am,A0,etaB] = unpackState(self,a)
            % Unpack a conjugate-symmetric state into WaveVortexModel arrays.
            %
            % - Topic: Transform Galerkin states
            % - Declaration: [Ap,Am,A0,etaB] = unpackState(a)
            % - Parameter a: full-complex mixed Galerkin vector
            % - Returns Ap: positive-wave coefficients
            % - Returns Am: negative-wave coefficients
            % - Returns A0: balanced coefficients
            % - Returns etaB: bottom displacement in native horizontal layout
            arguments
                self (1,1) WVTerrainEnergyGalerkin
                a (:,1) double
            end
            self.validateStateSize(a);
            self.validateConjugacy(a);

            wvt = self.originatingTransform;
            Ap = zeros(wvt.spectralMatrixSize);
            Am = zeros(wvt.spectralMatrixSize);
            A0 = zeros(wvt.spectralMatrixSize);
            etaB = zeros(1,wvt.Nkl);
            nativeBlocks = find(self.horizontalLayout.isNative);
            for iBlock = reshape(nativeBlocks,1,[])
                iNative = self.horizontalLayout.nativeIndex(iBlock);
                rows = find(self.stateLayout.horizontalIndex == iBlock);
                for row = reshape(rows,1,[])
                    jLabel = self.stateLayout.j(row);
                    component = self.stateLayout.component(row);
                    if component == "Ap"
                        Ap(wvt.j == jLabel,iNative) = a(row);
                    elseif component == "Am"
                        Am(wvt.j == jLabel,iNative) = a(row);
                    elseif component == "A0"
                        A0(wvt.j == jLabel,iNative) = a(row);
                    else
                        etaB(iNative) = a(row);
                    end
                end
            end
        end

        function fields = reconstructState(self,a,options)
            % Reconstruct hatted or physical fields from mixed coefficients.
            %
            % Spectral output has size `Nz`-by-`nK`. Spatial output has size
            % `Nx`-by-`Ny`-by-`Nz` and may be complex unless `a` satisfies
            % the real-field conjugacy relation.
            %
            % - Topic: Transform Galerkin states
            % - Declaration: fields = reconstructState(a,outputDomain=domain)
            % - Parameter a: full-complex mixed Galerkin vector
            % - Parameter outputDomain: `"spectral"` or `"spatial"`
            % - Returns fields: structure containing hatted and physical fields
            arguments
                self (1,1) WVTerrainEnergyGalerkin
                a (:,1) double
                options.outputDomain (1,1) string {mustBeMember(options.outputDomain,["spectral","spatial"])} = "spectral"
            end
            self.validateStateSize(a);
            nK = height(self.horizontalLayout);
            Nz = self.originatingTransform.Nz;
            names = ["uHat","vHat","wHat","etaHat","uHatXi","vHatXi","wHatXi","etaHatXi"];
            fields = struct;
            for name = names
                fields.(name) = zeros(Nz,nK);
            end
            for iK = 1:nK
                rows = self.stateLayout.horizontalIndex == iK;
                block = self.basisBlocks{iK};
                for name = names
                    fields.(name)(:,iK) = block.(name)*a(rows);
                end
            end

            if options.outputDomain == "spatial"
                isHermitianState = self.isConjugateSymmetricState(a);
                for name = names
                    fields.(name) = self.spectralToSpatial(fields.(name),isHermitianState);
                end
                gamma3 = self.gamma;
                gradLnGammaX = self.originatingTransform.diffX(self.gamma)./self.gamma;
                gradLnGammaY = self.originatingTransform.diffY(self.gamma)./self.gamma;
                xi3 = shiftdim(self.originatingTransform.z,-2);
                fields.u = fields.uHat./gamma3;
                fields.v = fields.vHat./gamma3;
                fields.w = fields.wHat + xi3.*(fields.uHat.*gradLnGammaX + fields.vHat.*gradLnGammaY);
            end
        end

        function dual = applyReconstructionAdjoint(self,fields)
            % Apply the Euclidean adjoint of spectral reconstruction.
            %
            % This is the adjoint map, not the terrain-energy inverse
            % projection introduced in later milestones.
            %
            % - Topic: Transform Galerkin states
            % - Declaration: dual = applyReconstructionAdjoint(fields)
            % - Parameter fields: structure with spectral hatted fields
            % - Returns dual: dual mixed-coordinate vector
            arguments
                self (1,1) WVTerrainEnergyGalerkin
                fields (1,1) struct
            end
            required = ["uHat","vHat","wHat","etaHat"];
            for name = required
                if ~isfield(fields,name) || ~isequal(size(fields.(name)),[self.originatingTransform.Nz height(self.horizontalLayout)])
                    error("WVTerrainEnergyGalerkin:InvalidSpectralField", ...
                        "%s must have size [%d %d].",name,self.originatingTransform.Nz,height(self.horizontalLayout))
                end
            end
            dual = zeros(height(self.stateLayout),1);
            for iK = 1:height(self.horizontalLayout)
                rows = self.stateLayout.horizontalIndex == iK;
                block = self.basisBlocks{iK};
                dual(rows) = block.uHat'*fields.uHat(:,iK) + block.vHat'*fields.vHat(:,iK) + block.wHat'*fields.wHat(:,iK) + block.etaHat'*fields.etaHat(:,iK);
            end
        end

        function q = evaluateAPV(self,a)
            % Evaluate finite-terrain linear APV on the native spatial grid.
            %
            % - Topic: Transform Galerkin states
            % - Declaration: q = evaluateAPV(a)
            % - Parameter a: full-complex mixed Galerkin vector
            % - Returns q: finite-terrain APV on the native grid
            arguments
                self (1,1) WVTerrainEnergyGalerkin
                a (:,1) double
            end
            spectral = self.reconstructState(a);
            isHermitianState = self.isConjugateSymmetricState(a);
            uHat = self.spectralToSpatial(spectral.uHat,isHermitianState);
            vHat = self.spectralToSpatial(spectral.vHat,isHermitianState);
            uHatXi = self.spectralToSpatial(spectral.uHatXi,isHermitianState);
            vHatXi = self.spectralToSpatial(spectral.vHatXi,isHermitianState);
            etaXi = self.spectralToSpatial(spectral.etaHatXi,isHermitianState);
            gamma3 = self.gamma;
            u = uHat./gamma3;
            v = vHat./gamma3;
            uXi = uHatXi./gamma3;
            vXi = vHatXi./gamma3;
            gradLnGammaX = self.originatingTransform.diffX(self.gamma)./self.gamma;
            gradLnGammaY = self.originatingTransform.diffY(self.gamma)./self.gamma;
            xi3 = shiftdim(self.originatingTransform.z,-2);
            vx = self.horizontalDerivative(v,1);
            uy = self.horizontalDerivative(u,2);
            q = vx-xi3.*gradLnGammaX.*vXi-uy+xi3.*gradLnGammaY.*uXi-(self.originatingTransform.f./gamma3).*etaXi;
        end

        function block = flatModesForHorizontalMode(self,kMode,lMode)
            % Return the flat dense-oracle block for a horizontal mode.
            %
            % - Topic: Inspect flat modes
            % - Declaration: block = flatModesForHorizontalMode(kMode,lMode)
            % - Parameter kMode: integer zonal mode number
            % - Parameter lMode: integer meridional mode number
            % - Returns block: flat generalized-eigenproblem diagnostics
            arguments
                self (1,1) WVTerrainEnergyGalerkin
                kMode (1,1) double {mustBeInteger}
                lMode (1,1) double {mustBeInteger}
            end
            iK = find(self.horizontalLayout.kMode == kMode & self.horizontalLayout.lMode == lMode,1);
            if isempty(iK)
                error("WVTerrainEnergyGalerkin:UnknownHorizontalMode", ...
                    "The horizontal mode (%d,%d) is not retained.",kMode,lMode)
            end
            block = self.flatModeBlocks{iK};
        end

        function audit = auditConstrainedClosure(self)
            % Audit an energy-, APV-, and bottom-compatible dense closure.
            %
            % The audit preserves the raw finite-terrain Galerkin forms and
            % asks whether a skew-Hermitian generator in energy coordinates
            % can also satisfy zero APV tendency and the resolved strong
            % bottom-displacement equation. Incompatibility is scientific
            % output rather than an exception. `constrainedGenerator` and
            % `constrainedExchangeMatrix` are empty when no algebraically
            % feasible closure exists.
            %
            % ```matlab
            % audit = problem.auditConstrainedClosure();
            % ```
            %
            % - Topic: Audit constrained evolution
            % - Declaration: audit = auditConstrainedClosure()
            % - Returns audit: dense closure matrices and compatibility diagnostics
            audit = buildConstrainedClosureAudit(self);
        end
    end

    methods (Access=private)
        function self = WVTerrainEnergyGalerkin(wvt,h,gamma,j,oversampling)
            self.originatingTransform = wvt;
            self.topographicHeight = h;
            self.gamma = gamma;
            self.verticalModeIndices = j;
            self.horizontalOversamplingFactor = oversampling;
            self.hydrostaticTransform = self.constructMatchedHydrostaticTransform;
            [self.horizontalLayout,self.stateLayout,self.horizontalConjugateIndex] = self.constructLayouts;
            self.basisBlocks = self.constructBasisBlocks;
            self.conjugateCoordinateIndex = self.constructConjugateCoordinateMap;
            [self.flatModeBlocks,self.constructionDiagnostics] = self.constructFlatOracle;
            [self.finiteTerrainForms,finiteTerrainDiagnostics] = buildFiniteTerrainForms(self);
            self.constructionDiagnostics.finiteTerrain = finiteTerrainDiagnostics;
        end

        function hydro = constructMatchedHydrostaticTransform(self)
            wvt = self.originatingTransform;
            hydro = WVTransformHydrostatic([wvt.Lx wvt.Ly wvt.Lz],[wvt.Nx wvt.Ny wvt.Nz], ...
                shouldAntialias=wvt.shouldAntialias,z=wvt.z,j=wvt.j,Nj=wvt.Nj,N2Function=wvt.N2Function, ...
                rho0=wvt.rho0,planetaryRadius=wvt.planetaryRadius,rotationRate=wvt.rotationRate,latitude=wvt.latitude,g=wvt.g, ...
                dLnN2=wvt.dLnN2,PF0inv=wvt.PF0inv,QG0inv=wvt.QG0inv,PF0=wvt.PF0,QG0=wvt.QG0,h_0=wvt.h_0,P0=wvt.P0,Q0=wvt.Q0,z_int=wvt.z_int);
            hydro.removeAllForcing;
        end

        function [horizontal,state,horizontalConjugateIndex] = constructLayouts(self)
            wvt = self.originatingTransform;
            primary = double(wvt.dftPrimaryIndices2D(:));
            conjugate = double(wvt.dftConjugateIndices2D(:));
            dftIndex = unique([primary;conjugate],"sorted");
            [isPrimary,nativeIndex] = ismember(dftIndex,primary);
            [isConjugate,conjugateNativeIndex] = ismember(dftIndex,conjugate);
            if ~all(isPrimary | isConjugate)
                error("WVTerrainEnergyGalerkin:InvalidHorizontalLayout", ...
                    "The retained signed layout contains an index absent from the originating geometry.")
            end
            nativeIndex(~isPrimary) = conjugateNativeIndex(~isPrimary);

            dftConjugateIndex = double(WVGeometryDoublyPeriodic.indicesOfFourierConjugates(wvt.Nx,wvt.Ny));
            [isConjugateRetained,horizontalConjugateIndex] = ismember(dftConjugateIndex(dftIndex),dftIndex);
            if ~all(isConjugateRetained) || any(horizontalConjugateIndex(horizontalConjugateIndex) ~= (1:numel(dftIndex))')
                error("WVTerrainEnergyGalerkin:InvalidHorizontalConjugacy", ...
                    "The originating geometry does not provide a complete retained conjugate pairing.")
            end

            nyquistMask = logical(WVGeometryDoublyPeriodic.maskForNyquistModes(wvt.Nx,wvt.Ny));
            if any(nyquistMask(dftIndex))
                error("WVTerrainEnergyGalerkin:RetainedNyquistMode", ...
                    "The originating geometry must exclude horizontal Nyquist modes.")
            end

            [kModeDFT,lModeDFT] = ndgrid(wvt.kMode_dft,wvt.lMode_dft);
            [kDFT,lDFT] = ndgrid(wvt.k_dft,wvt.l_dft);
            kMode = kModeDFT(dftIndex);
            lMode = lModeDFT(dftIndex);
            isNative = isPrimary;
            horizontal = table(dftIndex,kMode,lMode,kDFT(dftIndex),lDFT(dftIndex),isPrimary,isNative,nativeIndex, ...
                'VariableNames',["dftIndex","kMode","lMode","k","l","isPrimary","isNative","nativeIndex"]);

            horizontalIndex = zeros(0,1);
            component = strings(0,1);
            jLabel = zeros(0,1);
            coordinateIsPrimary = false(0,1);
            internalJ = self.verticalModeIndices(self.verticalModeIndices > 0);
            for iK = 1:height(horizontal)
                if horizontal.k(iK)^2+horizontal.l(iK)^2 > 0
                    components = [repmat("Ap",numel(internalJ),1);repmat("Am",numel(internalJ),1);repmat("A0",numel(self.verticalModeIndices),1);"etaB"];
                    js = [internalJ;internalJ;self.verticalModeIndices;NaN];
                else
                    components = [repmat("Ap",numel(self.verticalModeIndices),1);repmat("Am",numel(self.verticalModeIndices),1);repmat("A0",numel(internalJ),1);"etaB"];
                    js = [self.verticalModeIndices;self.verticalModeIndices;internalJ;NaN];
                end
                n = numel(components);
                horizontalIndex(end+(1:n),1) = iK;
                component(end+(1:n),1) = components;
                jLabel(end+(1:n),1) = js;
                coordinateIsPrimary(end+(1:n),1) = horizontal.isPrimary(iK);
            end
            state = table((1:numel(component))',horizontalIndex,component,jLabel,coordinateIsPrimary, ...
                'VariableNames',["index","horizontalIndex","component","j","isPrimary"]);
        end

        function blocks = constructBasisBlocks(self)
            hydro = self.hydrostaticTransform;
            wvt = self.originatingTransform;
            Nz = wvt.Nz;
            nK = height(self.horizontalLayout);
            blocks = cell(nK,1);
            F = hydro.FinvMatrix;
            G = hydro.GinvMatrix;
            DzF = -(hydro.N2/hydro.g).*hydro.QG0inv*(squeeze(hydro.Q0./hydro.P0).*hydro.PF0);
            DzG = hydro.PF0inv*(squeeze(hydro.P0./(hydro.Q0.*hydro.h_0)).*hydro.QG0);
            FXi = DzF*F;
            GXi = DzG*G;

            for iK = 1:nK
                rows = find(self.stateLayout.horizontalIndex == iK);
                n = numel(rows);
                block = struct;
                names = ["uHat","vHat","wHat","etaHat","uHatXi","vHatXi","wHatXi","etaHatXi"];
                for name = names
                    block.(name) = zeros(Nz,n);
                end
                k = self.horizontalLayout.k(iK);
                l = self.horizontalLayout.l(iK);
                kappa = hypot(k,l);
                for i = 1:n
                    row = rows(i);
                    component = self.stateLayout.component(row);
                    jLabel = self.stateLayout.j(row);
                    if component == "etaB"
                        [chi,chiXi] = self.bottomFunction(kappa);
                        block.etaHat(:,i) = chi;
                        block.etaHatXi(:,i) = chiXi;
                        continue
                    end
                    ij = find(wvt.j == jLabel,1);
                    Fj = F(:,ij);
                    Gj = G(:,ij);
                    if kappa > 0 && (component == "Ap" || component == "Am")
                        omega = sqrt(wvt.f^2+wvt.g*hydro.h_0(ij)*kappa^2);
                        alpha = atan2(l,k);
                        if component == "Ap"
                            uFactor = cos(alpha)-1i*(wvt.f/omega)*sin(alpha);
                            vFactor = sin(alpha)+1i*(wvt.f/omega)*cos(alpha);
                            etaFactor = -kappa*hydro.h_0(ij)/omega;
                        else
                            uFactor = cos(alpha)+1i*(wvt.f/omega)*sin(alpha);
                            vFactor = sin(alpha)-1i*(wvt.f/omega)*cos(alpha);
                            etaFactor = kappa*hydro.h_0(ij)/omega;
                        end
                        block.uHat(:,i) = uFactor*Fj;
                        block.vHat(:,i) = vFactor*Fj;
                        block.wHat(:,i) = -1i*kappa*hydro.h_0(ij)*Gj;
                        block.wHatXi(:,i) = -1i*kappa*hydro.h_0(ij)*GXi(:,ij);
                        block.etaHat(:,i) = etaFactor*Gj;
                        block.uHatXi(:,i) = uFactor*FXi(:,ij);
                        block.vHatXi(:,i) = vFactor*FXi(:,ij);
                        block.etaHatXi(:,i) = etaFactor*GXi(:,ij);
                    elseif kappa > 0
                        Lr2inv = 0;
                        if jLabel > 0
                            Lr2inv = wvt.f^2/(wvt.g*hydro.h_0(ij));
                        end
                        denominator = kappa^2+Lr2inv;
                        uFactor = 1i*l/denominator;
                        vFactor = -1i*k/denominator;
                        etaFactor = 0;
                        if jLabel > 0
                            etaFactor = -(wvt.f/wvt.g)/denominator;
                        end
                        block.uHat(:,i) = uFactor*Fj;
                        block.vHat(:,i) = vFactor*Fj;
                        block.etaHat(:,i) = etaFactor*Gj;
                        block.uHatXi(:,i) = uFactor*FXi(:,ij);
                        block.vHatXi(:,i) = vFactor*FXi(:,ij);
                        block.etaHatXi(:,i) = etaFactor*GXi(:,ij);
                    elseif component == "Ap" || component == "Am"
                        signV = 1;
                        if component == "Am"
                            signV = -1;
                        end
                        block.uHat(:,i) = Fj;
                        block.vHat(:,i) = signV*1i*Fj;
                        block.uHatXi(:,i) = FXi(:,ij);
                        block.vHatXi(:,i) = signV*1i*FXi(:,ij);
                    else
                        block.etaHat(:,i) = Gj;
                        block.etaHatXi(:,i) = GXi(:,ij);
                    end
                end
                blocks{iK} = block;
            end

            for iK = find(~self.horizontalLayout.isPrimary).'
                iPartner = self.horizontalConjugateIndex(iK);
                partner = blocks{iPartner};
                rows = find(self.stateLayout.horizontalIndex == iK);
                partnerRows = find(self.stateLayout.horizontalIndex == iPartner);
                for i = 1:numel(rows)
                    component = self.stateLayout.component(rows(i));
                    jLabel = self.stateLayout.j(rows(i));
                    partnerComponent = component;
                    if component == "Ap"
                        partnerComponent = "Am";
                    elseif component == "Am"
                        partnerComponent = "Ap";
                    end
                    ip = find(self.stateLayout.component(partnerRows) == partnerComponent & (self.stateLayout.j(partnerRows) == jLabel | (isnan(self.stateLayout.j(partnerRows)) & isnan(jLabel))),1);
                    names = ["uHat","vHat","wHat","etaHat","uHatXi","vHatXi","wHatXi","etaHatXi"];
                    for name = names
                        blocks{iK}.(name)(:,i) = conj(partner.(name)(:,ip));
                    end
                end
            end
        end

        function map = constructConjugateCoordinateMap(self)
            n = height(self.stateLayout);
            map = zeros(n,1);
            for row = 1:n
                iK = self.stateLayout.horizontalIndex(row);
                partnerK = self.horizontalConjugateIndex(iK);
                component = self.stateLayout.component(row);
                if component == "Ap"
                    partnerComponent = "Am";
                elseif component == "Am"
                    partnerComponent = "Ap";
                else
                    partnerComponent = component;
                end
                jLabel = self.stateLayout.j(row);
                candidates = find(self.stateLayout.horizontalIndex == partnerK & self.stateLayout.component == partnerComponent);
                if isnan(jLabel)
                    partner = candidates(isnan(self.stateLayout.j(candidates)));
                else
                    partner = candidates(self.stateLayout.j(candidates) == jLabel);
                end
                map(row) = partner;
            end
        end

        function [blocks,diagnostics] = constructFlatOracle(self)
            wvt = self.originatingTransform;
            zWeight = wvt.z_int(:);
            N2Weight = zWeight.*wvt.N2(:);
            nK = height(self.horizontalLayout);
            blocks = cell(nK,1);
            hermitianDefect = zeros(nK,1);
            skewDefect = zeros(nK,1);
            scaledRcond = zeros(nK,1);
            maximumResidual = zeros(nK,1);
            for iK = 1:nK
                basis = self.basisBlocks{iK};
                U = basis.uHat;
                V = basis.vHat;
                W = basis.wHat;
                Eta = basis.etaHat;
                E = wvt.rho0*(U'*(zWeight.*U)+V'*(zWeight.*V)+W'*(zWeight.*W)+Eta'*(N2Weight.*Eta));
                exchange = wvt.f*U'*(zWeight.*V)+Eta'*(N2Weight.*W);
                J = wvt.rho0*(exchange-exchange');
                k = self.horizontalLayout.k(iK);
                l = self.horizontalLayout.l(iK);
                Q = 1i*k*V-1i*l*U-wvt.f*basis.etaHatXi;
                hermitianDefect(iK) = norm(E-E',"fro")/max(norm(E,"fro"),realmin);
                skewDefect(iK) = norm(J+J',"fro")/max(norm(J,"fro"),realmin);

                scale = 1./sqrt(real(diag(E)));
                D = diag(scale);
                Es = D*E*D;
                Js = D*J*D;
                scaledRcond(iK) = rcond(Es);
                if hermitianDefect(iK) > 1e-13 || skewDefect(iK) > 1e-13
                    error("WVTerrainEnergyGalerkin:FlatStructureFailure", ...
                        "Raw flat forms failed their structural gate at horizontal mode (%d,%d): E %.3g, J %.3g.",self.horizontalLayout.kMode(iK),self.horizontalLayout.lMode(iK),hermitianDefect(iK),skewDefect(iK))
                end
                if scaledRcond(iK) <= 1e-12
                    error("WVTerrainEnergyGalerkin:IllConditionedFlatEnergy", ...
                        "The scaled flat energy form is ill-conditioned at horizontal mode (%d,%d).",self.horizontalLayout.kMode(iK),self.horizontalLayout.lMode(iK))
                end
                Es = (Es+Es')/2;
                Js = (Js-Js')/2;
                chol(Es,"lower");
                [C,frequency] = eig(1i*J,E,"vector");
                [frequency,order] = sort(real(frequency));
                C = C(:,order);
                for iMode = 1:size(C,2)
                    C(:,iMode) = C(:,iMode)/sqrt(real(C(:,iMode)'*E*C(:,iMode)));
                end
                frequencyGroupTolerance = 100*eps*size(E,1)*max([abs(wvt.f);abs(frequency)]);
                iFirst = 1;
                while iFirst <= numel(frequency)
                    iLast = find(abs(frequency-frequency(iFirst)) <= frequencyGroupTolerance,1,"last");
                    group = iFirst:iLast;
                    gram = (C(:,group)'*E*C(:,group));
                    gram = (gram+gram')/2;
                    R = chol(gram);
                    C(:,group) = C(:,group)/R;
                    iFirst = iLast+1;
                end
                residual = zeros(numel(frequency),1);
                qgpvNorm = zeros(numel(frequency),1);
                for iMode = 1:numel(frequency)
                    residual(iMode) = norm(1i*J*C(:,iMode)-frequency(iMode)*E*C(:,iMode))/(max(norm(E*C(:,iMode))*max(abs(frequency(iMode)),abs(wvt.f)),realmin));
                    qgpvNorm(iMode) = sqrt(real((Q*C(:,iMode))'*(zWeight.*(Q*C(:,iMode)))));
                end
                maximumResidual(iK) = max(residual);
                blocks{iK} = struct("E",E,"J",J,"Q",Q,"coordinateScale",scale,"scaledE",Es,"scaledJ",Js, ...
                    "frequency",frequency,"eigenvectors",C,"residual",residual,"qgpvNorm",qgpvNorm, ...
                    "horizontalIndex",iK,"kMode",self.horizontalLayout.kMode(iK),"lMode",self.horizontalLayout.lMode(iK));
            end
            diagnostics = struct("maximumHermitianDefect",max(hermitianDefect),"maximumSkewHermitianDefect",max(skewDefect), ...
                "minimumScaledEnergyRcond",min(scaledRcond),"maximumEigenResidual",max(maximumResidual), ...
                "hermitianDefect",hermitianDefect,"skewHermitianDefect",skewDefect,"scaledEnergyRcond",scaledRcond);
            internal = find(self.verticalModeIndices > 0);
            verticalColumns = arrayfun(@(j)find(wvt.j == j,1),self.verticalModeIndices(internal));
            G = self.hydrostaticTransform.GinvMatrix(:,verticalColumns);
            DzG = self.hydrostaticTransform.PF0inv*(squeeze(self.hydrostaticTransform.P0./(self.hydrostaticTransform.Q0.*self.hydrostaticTransform.h_0)).*self.hydrostaticTransform.QG0);
            GXi = DzG*G;
            effectiveVerticalWavenumber = sqrt(sum(abs(GXi).^2.*zWeight,1)./sum(abs(G).^2.*zWeight,1)).';
            horizontalWavenumber = hypot(self.horizontalLayout.k,self.horizontalLayout.l).';
            diagnostics.effectiveVerticalWavenumber = effectiveVerticalWavenumber;
            diagnostics.nonhydrostaticIndicator = horizontalWavenumber./effectiveVerticalWavenumber;
        end

        function [chi,chiXi] = bottomFunction(self,kappa)
            xi = self.originatingTransform.z(:);
            D = self.originatingTransform.Lz;
            if kappa == 0
                chi = -xi/D;
                chiXi = -ones(size(xi))/D;
                chi(1) = 1;
                chi(end) = 0;
                return
            end
            x = -xi;
            denominator = -expm1(-2*kappa*D);
            chi = exp(kappa*(x-D)).*(-expm1(-2*kappa*x))/denominator;
            chiXi = -kappa*exp(kappa*(x-D)).*(1+exp(-2*kappa*x))/denominator;
            chi(1) = 1;
            chi(end) = 0;
        end

        function spatial = spectralToSpatial(self,spectral,isHermitian)
            wvt = self.originatingTransform;
            if isHermitian
                native = zeros(size(spectral,1),wvt.Nkl);
                primaryBlocks = find(self.horizontalLayout.isPrimary);
                native(:,self.horizontalLayout.nativeIndex(primaryBlocks)) = spectral(:,primaryBlocks);
                spatial = wvt.transformToSpatialDomainWithFourier(native);
                return
            end

            dft = zeros(wvt.Nx*wvt.Ny,size(spectral,1));
            dft(self.horizontalLayout.dftIndex,:) = spectral.';
            dft = reshape(dft,[wvt.Nx wvt.Ny size(spectral,1)]);
            spatial = self.transformComplexDFTGridToSpatial(dft);
        end

        function derivative = horizontalDerivative(self,field,dimension)
            wvt = self.originatingTransform;
            if isreal(field)
                if dimension == 1
                    derivative = wvt.diffX(field);
                else
                    derivative = wvt.diffY(field);
                end
                return
            end

            if ismatrix(field)
                field = reshape(field,[wvt.Nx wvt.Ny 1]);
            end
            spectral = wvt.transformFromSpatialDomainToDFTGrid(field);
            nyquistMask = logical(WVGeometryDoublyPeriodic.maskForNyquistModes(wvt.Nx,wvt.Ny));
            if dimension == 1
                factor = 1i*wvt.k_dft;
                factor(nyquistMask(:,1)) = 0;
            else
                factor = 1i*reshape(wvt.l_dft,1,wvt.Ny,1);
                factor(nyquistMask(1,:,1)) = 0;
            end
            derivative = self.transformComplexDFTGridToSpatial(factor.*spectral);
        end

        function spatial = transformComplexDFTGridToSpatial(self,dft)
            wvt = self.originatingTransform;
            spatial = ifft(ifft(dft,wvt.Nx,1),wvt.Ny,2)*(wvt.Nx*wvt.Ny);
        end

        function bool = isConjugateSymmetricState(self,a)
            defect = norm(a-conj(a(self.conjugateCoordinateIndex)))/max(1,norm(a));
            bool = defect <= 1e-12;
        end

        function validateStateSize(self,a)
            if numel(a) ~= height(self.stateLayout)
                error("WVTerrainEnergyGalerkin:InvalidStateSize", ...
                    "The state must contain %d coordinates.",height(self.stateLayout))
            end
        end

        function validateConjugacy(self,a)
            self.validateStateSize(a);
            defect = norm(a-conj(a(self.conjugateCoordinateIndex)))/max(1,norm(a));
            if defect > 1e-12
                error("WVTerrainEnergyGalerkin:ConjugacyViolation", ...
                    "The state does not satisfy real-field conjugacy; relative defect is %.3g.",defect)
            end
        end
    end
end
