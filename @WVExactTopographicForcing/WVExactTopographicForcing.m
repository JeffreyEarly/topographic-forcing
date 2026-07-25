classdef WVExactTopographicForcing < WVForcing
    % Apply the exact flow-linear terrain terms in mapped coordinates.
    %
    % `WVExactTopographicForcing` evaluates the terrain terms on the native
    % spatial grid of a `WVTransformBoussinesq`. The transform coordinate
    % `z` is interpreted as the mapped coordinate $$\xi$$, while the
    % transform fields are interpreted as the hatted projection-ready
    % variables. The terrain is exact in height and the forcing is linear
    % in flow amplitude.
    %
    % The initial implementation supports discretely constant
    % stratification and performs no wave--vortex projection itself.
    % `WVTransformBoussinesq.nonlinearFlux` performs the single projection
    % after accumulating all spatial forcings.
    %
    % ```matlab
    % force = WVExactTopographicForcing(wvt,topographicHeight=h);
    % wvt.removeAllForcing();
    % wvt.addForcing(force);
    % ```
    %
    % - Topic: Create the forcing
    % - Topic: Inspect mapped geometry
    % - Topic: Evaluate the forcing
    % - Topic: CAAnnotatedClass requirement
    % - Declaration: classdef WVExactTopographicForcing < WVForcing

    properties (SetAccess = private)
        % Upward-positive topographic height $$h(x,y)$$ in meters.
        %
        % The field is interpreted periodically on the transform's
        % horizontal Fourier grid.
        %
        % - Topic: Inspect mapped geometry
        topographicHeight (:,:) double

        % Mapped-depth factor $$\gamma=1-h/D$$.
        %
        % - Topic: Inspect mapped geometry
        gamma (:,:) double

        % Physical vertical coordinate $$z=\gamma\xi$$ in meters.
        %
        % The array has size $$N_x\times N_y\times N_z$$.
        %
        % - Topic: Inspect mapped geometry
        physicalZ (:,:,:) double
    end

    properties (Access = private)
        xi
        dGammaDx
        dGammaDy
        dLogGammaDx
        dLogGammaDy
    end

    methods
        function self = WVExactTopographicForcing(wvt,options)
            % Create an exact linear topographic forcing.
            %
            % - Topic: Create the forcing
            % - Declaration: force = WVExactTopographicForcing(wvt,options)
            % - Parameter wvt: `WVTransformBoussinesq` receiving the forcing
            % - Parameter options.topographicHeight: real upward-positive terrain of size $$N_x\times N_y$$ in meters
            % - Parameter options.name: forcing name registered with the transform
            % - Returns force: configured `WVExactTopographicForcing`
            arguments (Input)
                wvt WVTransform {mustBeNonempty}
                options.topographicHeight double
                options.name (1,1) string = "exact topographic forcing"
            end

            if ~isa(wvt,"WVTransformBoussinesq")
                error("WVExactTopographicForcing:UnsupportedTransform", "WVExactTopographicForcing currently supports only WVTransformBoussinesq.")
            end
            if wvt.shouldAntialias
                error("WVExactTopographicForcing:AntialiasingUnsupported", "Construct the initial proof of concept with shouldAntialias=false. Antialiasing support is deferred to Milestone 8.")
            end
            if ~isequal(size(options.topographicHeight),[wvt.Nx wvt.Ny])
                error("WVExactTopographicForcing:InvalidTopographicHeightSize", "topographicHeight must have size [%d %d], matching the transform horizontal grid.", wvt.Nx, wvt.Ny)
            end
            if ~isreal(options.topographicHeight) || any(~isfinite(options.topographicHeight),"all")
                error("WVExactTopographicForcing:InvalidTopographicHeight", "topographicHeight must be real and finite.")
            end
            if any(wvt.N2 ~= wvt.N2(1))
                error("WVExactTopographicForcing:NonconstantStratificationUnsupported", "The initial proof of concept requires discretely constant N2. General stationary stratification is deferred to Milestone 9.")
            end

            gamma = 1-options.topographicHeight/wvt.Lz;
            if any(gamma <= 0,"all")
                error("WVExactTopographicForcing:NonpositiveMappedDepth", "topographicHeight must satisfy gamma=1-h/D>0 everywhere.")
            end

            self@WVForcing(wvt,options.name,WVForcingType("NonhydrostaticSpatial"));
            self.topographicHeight = options.topographicHeight;
            self.gamma = gamma;
            self.xi = reshape(wvt.z,1,1,[]);
            self.physicalZ = gamma.*self.xi;
            self.dGammaDx = wvt.diffX(gamma);
            self.dGammaDy = wvt.diffY(gamma);
            self.dLogGammaDx = self.dGammaDx./gamma;
            self.dLogGammaDy = self.dGammaDy./gamma;
        end

        function [u,v,w] = physicalVelocity(self)
            % Reconstruct physical velocity from the hatted transform fields.
            %
            % - Topic: Inspect mapped geometry
            % - Declaration: [u,v,w] = physicalVelocity()
            % - Returns u: physical zonal velocity in meters per second
            % - Returns v: physical meridional velocity in meters per second
            % - Returns w: physical vertical velocity in meters per second
            [hatU,hatV,hatW] = self.wvt.variableWithName('u','v','w');
            u = hatU./self.gamma;
            v = hatV./self.gamma;
            w = hatW+self.xi.*(hatU.*self.dLogGammaDx+hatV.*self.dLogGammaDy);
        end

        function residual = bottomKinematicResidual(self)
            % Evaluate the physical sloping-bottom kinematic residual.
            %
            % The returned $$N_x\times N_y$$ field is
            % $$w_b-u_b\partial_xh-v_b\partial_yh$$.
            %
            % - Topic: Inspect mapped geometry
            % - Declaration: residual = bottomKinematicResidual()
            % - Returns residual: bottom-normal velocity residual in meters per second
            [u,v,w] = self.physicalVelocity();
            [~,iBottom] = min(self.wvt.z);
            dHdx = -self.wvt.Lz*self.dGammaDx;
            dHdy = -self.wvt.Lz*self.dGammaDy;
            residual = w(:,:,iBottom)-u(:,:,iBottom).*dHdx-v(:,:,iBottom).*dHdy;
        end

        function value = mappedVolumeIntegral(self,field)
            % Integrate a mapped field with the physical Jacobian.
            %
            % The horizontal integral is area-normalized, matching
            % `WVTransformBoussinesq.volumeIntegral`.
            %
            % - Topic: Inspect mapped geometry
            % - Declaration: value = mappedVolumeIntegral(field)
            % - Parameter field: array of size $$N_x\times N_y\times N_z$$
            % - Returns value: horizontal-mean physical volume integral
            arguments (Input)
                self WVExactTopographicForcing
                field double
            end
            if ~isequal(size(field),self.wvt.spatialMatrixSize)
                expectedSize = self.wvt.spatialMatrixSize;
                error("WVExactTopographicForcing:InvalidSpatialFieldSize", "field must have size [%d %d %d].", expectedSize(1), expectedSize(2), expectedSize(3))
            end
            value = sum(shiftdim(self.wvt.z_int,-2).*mean(mean(self.gamma.*field,1),2),"all");
        end

        function [Fu,Fv,Fw,Feta] = linearTerrainTendency(self)
            % Evaluate the projection-ready exact linear terrain tendency.
            %
            % The returned fields are
            % $$-(\mathcal T_u,\mathcal T_v,\mathcal T_w,\mathcal T_\eta)$$
            % for discretely constant $$N^2$$. Pressure is reconstructed
            % directly from the transform's `p` field.
            %
            % - Topic: Evaluate the forcing
            % - Declaration: [Fu,Fv,Fw,Feta] = linearTerrainTendency()
            % - Returns Fu: zonal momentum tendency
            % - Returns Fv: meridional momentum tendency
            % - Returns Fw: vertical momentum tendency
            % - Returns Feta: displacement tendency
            [hatU,hatV,hatP] = self.wvt.variableWithName('u','v','p');
            dPdx = self.wvt.diffX(hatP);
            dPdy = self.wvt.diffY(hatP);
            dPdxi = self.wvt.diffZF(hatP);
            inverseRho0 = 1/self.wvt.rho0;

            Tu = inverseRho0*((self.gamma-1).*dPdx-self.xi.*self.dGammaDx.*dPdxi);
            Tv = inverseRho0*((self.gamma-1).*dPdy-self.xi.*self.dGammaDy.*dPdxi);
            HuL = self.wvt.f*hatV-inverseRho0*dPdx-Tu;
            HvL = -self.wvt.f*hatU-inverseRho0*dPdy-Tv;
            Tw = self.xi.*(HuL.*self.dLogGammaDx+HvL.*self.dLogGammaDy) ...
                +(1./self.gamma-1).*inverseRho0.*dPdxi;
            Teta = -self.xi.*(hatU.*self.dLogGammaDx+hatV.*self.dLogGammaDy);

            Fu = -Tu;
            Fv = -Tv;
            Fw = -Tw;
            Feta = -Teta;
        end

        function [Fu,Fv,Fw,Feta] = addNonhydrostaticSpatialForcing(self,wvt,Fu,Fv,Fw,Feta)
            % Add the exact linear terrain tendency to spatial forcing.
            %
            % - Topic: Evaluate the forcing
            % - Declaration: [Fu,Fv,Fw,Feta] = addNonhydrostaticSpatialForcing(wvt,Fu,Fv,Fw,Feta)
            % - Parameter wvt: transform at the current model time
            % - Parameter Fu: accumulated zonal momentum tendency
            % - Parameter Fv: accumulated meridional momentum tendency
            % - Parameter Fw: accumulated vertical momentum tendency
            % - Parameter Feta: accumulated displacement tendency
            % - Returns Fu: zonal tendency including this forcing
            % - Returns Fv: meridional tendency including this forcing
            % - Returns Fw: vertical tendency including this forcing
            % - Returns Feta: displacement tendency including this forcing
            self.requireOriginatingTransform(wvt);
            [terrainFu,terrainFv,terrainFw,terrainFeta] = self.linearTerrainTendency();
            Fu = Fu+terrainFu;
            Fv = Fv+terrainFv;
            Fw = Fw+terrainFw;
            Feta = Feta+terrainFeta;
        end

        function force = forcingWithResolutionOfTransform(~,~)
            % Reject transform resolution changes until Milestone 8.
            %
            % - Topic: Create the forcing
            % - Declaration: force = forcingWithResolutionOfTransform(wvtX2)
            % - Returns force: no value; this method always throws
            force = WVForcing.empty(0,0); %#ok<NASGU>
            error("WVExactTopographicForcing:ResolutionChangeUnsupported", "Resolution conversion is deferred to Milestone 8, where terrain resampling and geometry rebuilding will be implemented.")
        end
    end

    methods (Access = private)
        function requireOriginatingTransform(self,wvt)
            if wvt ~= self.wvt
                error("WVExactTopographicForcing:TransformMismatch", "The forcing can only be evaluated with the WVTransformBoussinesq instance used during construction.")
            end
        end
    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            % Return required persisted property names.
            %
            % Persistence is intentionally deferred for this research
            % proof of concept.
            %
            % - Topic: CAAnnotatedClass requirement
            % - Declaration: vars = classRequiredPropertyNames()
            % - Returns vars: empty property-name cell array
            vars = {};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            % Return property annotations defined by this class.
            %
            % Persistence is intentionally deferred for this research
            % proof of concept.
            %
            % - Topic: CAAnnotatedClass requirement
            % - Declaration: propertyAnnotations = classDefinedPropertyAnnotations()
            % - Returns propertyAnnotations: empty annotation array
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
        end
    end
end
