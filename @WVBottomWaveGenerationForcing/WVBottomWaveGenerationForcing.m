classdef WVBottomWaveGenerationForcing < WVForcing
    % Generate internal waves from prescribed barotropic flow over topography.
    %
    % `WVBottomWaveGenerationForcing` projects the first-order bottom
    % velocity
    %
    % $$
    % g_b=\boldsymbol U_{\mathrm{bt}}(t)\boldsymbol{\cdot}\nabla_Hh
    % $$
    %
    % onto the rigid-lid wave modes using their bottom pressure. The
    % projection is precomputed, so ordinary forcing calls add spectral
    % wave tendencies without a pressure solve or spatial transform. The
    % incoming balanced tendency is left unchanged.
    %
    % ```matlab
    % forcing = WVBottomWaveGenerationForcing(wvt,topographicHeight=h,barotropicVelocityAmplitude=[0.05; 0]);
    % wvt.removeAllForcing();
    % wvt.addForcing(forcing);
    % ```
    %
    % - Topic: Create the forcing
    % - Topic: Inspect the forcing
    % - Topic: Evaluate the forcing
    % - Topic: CAAnnotatedClass requirement
    % - Declaration: classdef WVBottomWaveGenerationForcing < WVForcing

    properties (SetAccess = private)
        % Upward-positive topographic height $$h(x,y)$$ in meters.
        %
        % The field is stationary and periodic on the transform's
        % horizontal grid.
        %
        % - Topic: Inspect the forcing
        topographicHeight (:,:) double

        % Complex barotropic velocity amplitude in meters per second.
        %
        % The two entries are the zonal and meridional amplitudes in
        % $$\boldsymbol U_{\mathrm{bt}}=R(t)\operatorname{Re}
        % \{\widehat{\boldsymbol U}_{\mathrm{bt}}e^{-i\omega(t-t_0)}\}$$.
        %
        % - Topic: Inspect the forcing
        barotropicVelocityAmplitude (2,1) double

        % Barotropic angular frequency $$\omega$$ in radians per second.
        %
        % - Topic: Inspect the forcing
        frequency (1,1) double

        % Duration of the half-cosine startup ramp in seconds.
        %
        % - Topic: Inspect the forcing
        rampDuration (1,1) double

        % Time at which the prescribed barotropic forcing begins, in seconds.
        %
        % - Topic: Inspect the forcing
        startTime (1,1) double
    end

    properties (Access = private)
        dHdx
        dHdy
        responsePlusX
        responsePlusY
        responseMinusX
        responseMinusY
    end

    methods
        function self = WVBottomWaveGenerationForcing(wvt,options)
            % Create a prescribed bottom wave-generation forcing.
            %
            % The default frequency is the M2 tidal frequency. A zero ramp
            % duration activates the harmonic current immediately at
            % `startTime`.
            %
            % - Topic: Create the forcing
            % - Declaration: forcing = WVBottomWaveGenerationForcing(wvt,options)
            % - Parameter wvt: `WVTransformBoussinesq` receiving the forcing
            % - Parameter options.topographicHeight: real stationary terrain of size $$N_x\times N_y$$ in meters
            % - Parameter options.barotropicVelocityAmplitude: finite complex two-component velocity amplitude in meters per second
            % - Parameter options.frequency: positive angular frequency in radians per second
            % - Parameter options.rampDuration: nonnegative startup-ramp duration in seconds
            % - Parameter options.startTime: finite forcing start time in seconds
            % - Parameter options.name: forcing name registered with the transform
            % - Returns forcing: configured `WVBottomWaveGenerationForcing`
            arguments (Input)
                wvt WVTransform {mustBeNonempty}
                options.topographicHeight double
                options.barotropicVelocityAmplitude double
                options.frequency double = 2*pi/(12.4206012*3600)
                options.rampDuration double = 0
                options.startTime double = wvt.t
                options.name (1,1) string = "bottom wave generation"
            end

            if ~isa(wvt,"WVTransformBoussinesq")
                error("WVBottomWaveGenerationForcing:UnsupportedTransform", "WVBottomWaveGenerationForcing currently supports only WVTransformBoussinesq.")
            end
            if any(wvt.N2 ~= wvt.N2(1))
                error("WVBottomWaveGenerationForcing:NonconstantStratificationUnsupported", "The initial implementation requires discretely constant N2.")
            end
            if ~isequal(size(options.topographicHeight),[wvt.Nx wvt.Ny])
                error("WVBottomWaveGenerationForcing:InvalidTopographicHeightSize", "topographicHeight must have size [%d %d], matching the transform horizontal grid.", wvt.Nx, wvt.Ny)
            end
            if ~isreal(options.topographicHeight) || any(~isfinite(options.topographicHeight),"all")
                error("WVBottomWaveGenerationForcing:InvalidTopographicHeight", "topographicHeight must be real and finite.")
            end
            if ~isequal(size(options.barotropicVelocityAmplitude),[2 1]) || any(~isfinite(options.barotropicVelocityAmplitude),"all")
                error("WVBottomWaveGenerationForcing:InvalidBarotropicVelocityAmplitude", "barotropicVelocityAmplitude must be a finite complex 2-by-1 vector.")
            end
            if ~isscalar(options.frequency) || ~isreal(options.frequency) || ~isfinite(options.frequency) || options.frequency <= 0
                error("WVBottomWaveGenerationForcing:InvalidFrequency", "frequency must be a finite positive real scalar.")
            end
            if ~isscalar(options.rampDuration) || ~isreal(options.rampDuration) || ~isfinite(options.rampDuration) || options.rampDuration < 0
                error("WVBottomWaveGenerationForcing:InvalidRampDuration", "rampDuration must be a finite nonnegative real scalar.")
            end
            if ~isscalar(options.startTime) || ~isreal(options.startTime) || ~isfinite(options.startTime)
                error("WVBottomWaveGenerationForcing:InvalidStartTime", "startTime must be a finite real scalar.")
            end
            if strlength(options.name) == 0
                error("WVBottomWaveGenerationForcing:InvalidName", "name must be a nonempty string.")
            end

            self@WVForcing(wvt,options.name,WVForcingType("Spectral"));
            self.topographicHeight = options.topographicHeight;
            self.barotropicVelocityAmplitude = options.barotropicVelocityAmplitude;
            self.frequency = options.frequency;
            self.rampDuration = options.rampDuration;
            self.startTime = options.startTime;
            self.dHdx = wvt.diffX(self.topographicHeight);
            self.dHdy = wvt.diffY(self.topographicHeight);

            terrainFourier = wvt.transformFromSpatialDomainWithFourier(repmat(self.topographicHeight,1,1,wvt.Nz));
            [self.responsePlusX,self.responsePlusY,self.responseMinusX,self.responseMinusY] = self.buildResponses(wvt,terrainFourier(1,:));
        end

        function velocity = barotropicVelocityAtTime(self,t)
            % Evaluate the prescribed horizontally uniform current.
            %
            % - Topic: Evaluate the forcing
            % - Declaration: velocity = barotropicVelocityAtTime(t)
            % - Parameter t: finite scalar time in seconds
            % - Returns velocity: real two-component velocity in meters per second
            arguments (Input)
                self WVBottomWaveGenerationForcing
                t (1,1) double {mustBeFinite}
            end
            arguments (Output)
                velocity (2,1) double
            end

            elapsed = t-self.startTime;
            if elapsed < 0
                velocity = zeros(2,1);
                return
            end
            if self.rampDuration == 0 || elapsed >= self.rampDuration
                ramp = 1;
            else
                ramp = 0.5*(1-cos(pi*elapsed/self.rampDuration));
            end
            velocity = ramp*real(self.barotropicVelocityAmplitude*exp(-1i*self.frequency*elapsed));
        end

        function gBottom = bottomVelocityAtTime(self,t)
            % Evaluate $$g_b=\boldsymbol U_{\mathrm{bt}}\boldsymbol{\cdot}\nabla_Hh$$.
            %
            % - Topic: Evaluate the forcing
            % - Declaration: gBottom = bottomVelocityAtTime(t)
            % - Parameter t: finite scalar time in seconds
            % - Returns gBottom: real bottom-normal velocity on the horizontal grid
            velocity = self.barotropicVelocityAtTime(t);
            gBottom = velocity(1)*self.dHdx+velocity(2)*self.dHdy;
        end

        function [Fp,Fm,F0] = addSpectralForcing(self,wvt,Fp,Fm,F0)
            % Add the precomputed wave-generation tendency.
            %
            % The physical modal tendencies are converted componentwise to
            % WaveVortexModel's stored interaction representation. `F0` is
            % returned without modification.
            %
            % - Topic: Evaluate the forcing
            % - Declaration: [Fp,Fm,F0] = addSpectralForcing(wvt,Fp,Fm,F0)
            % - Parameter wvt: transform at the current model time
            % - Parameter Fp: accumulated positive-wave tendency
            % - Parameter Fm: accumulated negative-wave tendency
            % - Parameter F0: accumulated balanced tendency
            % - Returns Fp: positive-wave tendency including this forcing
            % - Returns Fm: negative-wave tendency including this forcing
            % - Returns F0: unchanged incoming balanced tendency
            self.requireOriginatingTransform(wvt);
            velocity = self.barotropicVelocityAtTime(wvt.t);
            Fpt = velocity(1)*self.responsePlusX+velocity(2)*self.responsePlusY;
            Fmt = velocity(1)*self.responseMinusX+velocity(2)*self.responseMinusY;
            Fp = Fp+Fpt.*wvt.conjPhase;
            Fm = Fm+Fmt.*wvt.phase;
        end

        function forcing = forcingWithResolutionOfTransform(~,~)
            % Reject explicit transform-resolution changes until Milestone 8.
            %
            % Native transforms with either value of `shouldAntialias` are
            % supported. This method concerns rebuilding for a different
            % transform resolution.
            %
            % - Topic: Create the forcing
            % - Declaration: forcing = forcingWithResolutionOfTransform(wvtX2)
            % - Returns forcing: no value; this method always throws
            forcing = WVForcing.empty(0,0); %#ok<NASGU>
            error("WVBottomWaveGenerationForcing:ResolutionChangeUnsupported", "Explicit transform-resolution conversion is deferred to Milestone 8.")
        end
    end

    methods (Access = private)
        function requireOriginatingTransform(self,wvt)
            if wvt ~= self.wvt
                error("WVBottomWaveGenerationForcing:TransformMismatch", "The forcing can only be evaluated with the WVTransformBoussinesq instance used during construction.")
            end
        end
    end

    methods (Static, Access = private)
        function [responsePlusX,responsePlusY,responseMinusX,responseMinusY] = buildResponses(wvt,terrainFourier)
            [~,iBottom] = min(wvt.z);
            piPlus = complex(zeros(size(wvt.Apm_TE_factor)));
            piMinus = complex(zeros(size(wvt.Apm_TE_factor)));
            for iK = 1:numel(wvt.K2unique)
                indices = wvt.K2uniqueK2Map{iK};
                bottomF = reshape(wvt.PFpmInv(iBottom,:,iK),[],1).*wvt.Ppm(:,iK);
                piPlus(:,indices) = wvt.g*bottomF.*wvt.NAp(:,indices);
                piMinus(:,indices) = wvt.g*bottomF.*wvt.NAm(:,indices);
            end

            dHdxFourier = 1i*wvt.K.*terrainFourier;
            dHdyFourier = 1i*wvt.L.*terrainFourier;
            responsePlusX = complex(zeros(size(piPlus)));
            responsePlusY = complex(zeros(size(piPlus)));
            responseMinusX = complex(zeros(size(piMinus)));
            responseMinusY = complex(zeros(size(piMinus)));
            maskPlus = logical(wvt.waveComponent.maskAp);
            maskMinus = logical(wvt.waveComponent.maskAm);
            responsePlusX(maskPlus) = conj(piPlus(maskPlus)).*dHdxFourier(maskPlus)./wvt.Apm_TE_factor(maskPlus);
            responsePlusY(maskPlus) = conj(piPlus(maskPlus)).*dHdyFourier(maskPlus)./wvt.Apm_TE_factor(maskPlus);
            responseMinusX(maskMinus) = conj(piMinus(maskMinus)).*dHdxFourier(maskMinus)./wvt.Apm_TE_factor(maskMinus);
            responseMinusY(maskMinus) = conj(piMinus(maskMinus)).*dHdyFourier(maskMinus)./wvt.Apm_TE_factor(maskMinus);
        end
    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            % Return required persisted property names.
            %
            % Persistence is deferred until the scientific forcing has
            % passed its prescribed-generation validation.
            %
            % - Topic: CAAnnotatedClass requirement
            % - Declaration: vars = classRequiredPropertyNames()
            % - Returns vars: empty property-name cell array
            vars = {};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            % Return property annotations defined by this class.
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
