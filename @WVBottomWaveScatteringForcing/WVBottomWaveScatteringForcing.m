classdef WVBottomWaveScatteringForcing < WVForcing
    % Scatter an evolving internal-wave field from stationary topography.
    %
    % `WVBottomWaveScatteringForcing` evaluates the first-order mean-depth
    % bottom velocity
    %
    % $$
    % g_b=\boldsymbol u_{H,d}\boldsymbol{\cdot}\nabla_Hh
    % -h\,\partial_z w_d
    % $$
    %
    % from the instantaneous wave coefficients and projects it onto the
    % rigid-lid wave modes with their bottom pressure. The forcing consumes
    % and produces only the wave branches; the incoming balanced tendency
    % is returned unchanged.
    %
    % ```matlab
    % forcing = WVBottomWaveScatteringForcing(wvt,topographicHeight=h);
    % wvt.addForcing(forcing);
    % ```
    %
    % - Topic: Create the forcing
    % - Topic: Inspect the forcing
    % - Topic: Evaluate the forcing
    % - Topic: Diagnose bottom displacement
    % - Topic: Restart persistence
    % - Declaration: classdef WVBottomWaveScatteringForcing < WVForcing

    properties (SetAccess = private)
        % Upward-positive topographic height $$h(x,y)$$ in meters.
        %
        % The field is stationary and periodic on the transform's
        % horizontal grid.
        %
        % - Topic: Inspect the forcing
        topographicHeight (:,:) double
    end

    properties (Access = private)
        dHdx
        dHdy
        bottomF
        bottomGz
        pressureProjectionPlus
        pressureProjectionMinus
        primaryDFTIndex
        conjugateDFTIndex
        conjugateWVIndex
    end

    methods
        function self = WVBottomWaveScatteringForcing(wvt,options)
            % Create an autonomous first-order wave-scattering forcing.
            %
            % The transform must contain a wave component and implement
            % `waveModeVerticalStructureAtIndex`.
            %
            % - Topic: Create the forcing
            % - Declaration: forcing = WVBottomWaveScatteringForcing(wvt,options)
            % - Parameter wvt: supported wave-bearing `WVTransform` receiving the forcing
            % - Parameter options.topographicHeight: real stationary terrain of size $$N_x\times N_y$$ in meters
            % - Parameter options.name: forcing name registered with the transform
            % - Returns forcing: configured `WVBottomWaveScatteringForcing`
            arguments (Input)
                wvt WVTransform {mustBeNonempty}
                options.topographicHeight double
                options.name (1,1) string = "bottom wave scattering"
            end

            if ~WVBottomWaveScatteringForcing.isSupportedTransform(wvt)
                error("WVBottomWaveScatteringForcing:UnsupportedTransform", "WVBottomWaveScatteringForcing requires a wave-bearing transform that implements waveModeVerticalStructureAtIndex.")
            end
            if ~isequal(size(options.topographicHeight),[wvt.Nx wvt.Ny])
                error("WVBottomWaveScatteringForcing:InvalidTopographicHeightSize", "topographicHeight must have size [%d %d], matching the transform horizontal grid.", wvt.Nx, wvt.Ny)
            end
            if ~isreal(options.topographicHeight) || any(~isfinite(options.topographicHeight),"all")
                error("WVBottomWaveScatteringForcing:InvalidTopographicHeight", "topographicHeight must be real and finite.")
            end
            if strlength(options.name) == 0
                error("WVBottomWaveScatteringForcing:InvalidName", "name must be a nonempty string.")
            end

            self@WVForcing(wvt,options.name,WVForcingType("Spectral"));
            self.topographicHeight = options.topographicHeight;
            self.dHdx = wvt.diffX(self.topographicHeight);
            self.dHdy = wvt.diffY(self.topographicHeight);
            [self.bottomF,self.bottomGz,self.pressureProjectionPlus,self.pressureProjectionMinus] = self.buildEndpointFactors(wvt);
            [self.primaryDFTIndex,self.conjugateDFTIndex,self.conjugateWVIndex] = self.horizontalFourierIndices(wvt);
        end

        function [gBottom,bottomFields] = bottomVelocityFromWaveState(self,wvt)
            % Reconstruct the instantaneous first-order bottom velocity.
            %
            % `bottomFields` contains the wave-only fields `u`, `v`, and
            % `dWdz` evaluated at the mean bottom. Balanced coefficients
            % are deliberately ignored.
            %
            % - Topic: Evaluate the forcing
            % - Declaration: [gBottom,bottomFields] = bottomVelocityFromWaveState(wvt)
            % - Parameter wvt: originating transform at the requested model time
            % - Returns gBottom: bottom velocity $$g_b$$ of size $$N_x\times N_y$$
            % - Returns bottomFields: structure containing `u`, `v`, and `dWdz`
            arguments (Input)
                self WVBottomWaveScatteringForcing {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            arguments (Output)
                gBottom (:,:) double
                bottomFields (1,1) struct
            end

            self.requireOriginatingTransform(wvt);
            [gBottom,bottomFields] = self.bottomVelocityFromCoefficients(wvt.Ap,wvt.Am,wvt.t);
        end

        function [Fp,Fm,F0] = addSpectralForcing(self,wvt,Fp,Fm,F0)
            % Add the autonomous first-order wave-scattering tendency.
            %
            % - Topic: Evaluate the forcing
            % - Declaration: [Fp,Fm,F0] = addSpectralForcing(wvt,Fp,Fm,F0)
            % - Parameter wvt: transform at the current model time
            % - Parameter Fp: accumulated positive-wave tendency
            % - Parameter Fm: accumulated negative-wave tendency
            % - Parameter F0: accumulated balanced tendency
            % - Returns Fp: positive-wave tendency including scattering
            % - Returns Fm: negative-wave tendency including scattering
            % - Returns F0: unchanged incoming balanced tendency
            self.requireOriginatingTransform(wvt);
            [~,~,gBottomFourier] = self.bottomVelocityFromCoefficients(wvt.Ap,wvt.Am,wvt.t);
            Fpt = self.pressureProjectionPlus.*gBottomFourier;
            Fmt = self.pressureProjectionMinus.*gBottomFourier;
            Fp = Fp+Fpt.*wvt.conjPhase;
            Fm = Fm+Fmt.*wvt.phase;
        end

        function forcing = forcingWithResolutionOfTransform(self,wvtX2)
            % Rebuild the scattering forcing at another resolution.
            %
            % - Topic: Create the forcing
            % - Declaration: forcing = forcingWithResolutionOfTransform(wvtX2)
            % - Parameter wvtX2: compatible supported transform at the target resolution
            % - Returns forcing: equivalent forcing rebuilt for `wvtX2`
            arguments (Input)
                self WVBottomWaveScatteringForcing {mustBeNonempty}
                wvtX2 WVTransform {mustBeNonempty}
            end
            arguments (Output)
                forcing WVBottomWaveScatteringForcing
            end

            if ~WVBottomWaveScatteringForcing.isSupportedTransform(wvtX2)
                error("WVBottomWaveScatteringForcing:UnsupportedTransform", "WVBottomWaveScatteringForcing requires a wave-bearing transform that implements waveModeVerticalStructureAtIndex.")
            end
            if ~isequal([self.wvt.Lx self.wvt.Ly self.wvt.Lz],[wvtX2.Lx wvtX2.Ly wvtX2.Lz])
                error("WVBottomWaveScatteringForcing:IncompatibleDomain", "Resolution conversion requires transforms with identical domain dimensions.")
            end

            terrainFourier = self.wvt.transformFromSpatialDomainWithFourier(repmat(self.topographicHeight,1,1,self.wvt.Nz));
            terrainFourierX2 = self.wvt.spectralVariableWithResolution(wvtX2,terrainFourier);
            terrainX2 = wvtX2.transformToSpatialDomainWithFourier(repmat(terrainFourierX2(1,:),wvtX2.Nz,1));
            forcing = WVBottomWaveScatteringForcing(wvtX2,topographicHeight=real(terrainX2(:,:,1)),name=string(self.name));
        end

        function writeToGroup(self,group,propertyAnnotations,attributes)
            % Write the scattering forcing to a transform-owned group.
            %
            % - Topic: Restart persistence
            % - Declaration: writeToGroup(group,propertyAnnotations,attributes)
            % - Parameter group: transform-owned NetCDF group for this forcing
            % - Parameter propertyAnnotations: forcing properties to persist
            % - Parameter attributes: additional NetCDF attributes
            arguments (Input)
                self WVBottomWaveScatteringForcing {mustBeNonempty}
                group NetCDFGroup {mustBeNonempty}
                propertyAnnotations CAPropertyAnnotation = CAPropertyAnnotation.empty(0,0)
                attributes = configureDictionary("string","string")
            end

            isTopographicHeight = string({propertyAnnotations.name}) == "topographicHeight";
            writeToGroup@CAAnnotatedClass(self,group,propertyAnnotations(~isTopographicHeight),attributes);
            if any(isTopographicHeight)
                annotation = propertyAnnotations(find(isTopographicHeight,1));
                variableAttributes = annotation.attributes;
                variableAttributes('units') = annotation.units;
                variableAttributes('long_name') = annotation.description;
                group.addVariable(annotation.name,annotation.dimensions,self.topographicHeight,isComplex=annotation.isComplex,attributes=variableAttributes);
            end
        end
    end

    methods (Access = private)
        function [gBottom,bottomFields,gBottomFourier] = bottomVelocityFromCoefficients(self,Ap,Am,t)
            phase = exp(self.wvt.iOmega*(t-self.wvt.t0));
            Apt = self.wvt.waveComponent.maskAp.*Ap.*phase;
            Amt = self.wvt.waveComponent.maskAm.*Am.*conj(phase);
            uBottomFourier = sum(self.bottomF.*(self.wvt.UAp.*Apt+self.wvt.UAm.*Amt),1);
            vBottomFourier = sum(self.bottomF.*(self.wvt.VAp.*Apt+self.wvt.VAm.*Amt),1);
            dWdzBottomFourier = sum(self.bottomGz.*(self.wvt.WAp.*Apt+self.wvt.WAm.*Amt),1);
            uBottom = self.spatialFieldFromHorizontalFourier(uBottomFourier);
            vBottom = self.spatialFieldFromHorizontalFourier(vBottomFourier);
            dWdzBottom = self.spatialFieldFromHorizontalFourier(dWdzBottomFourier);
            gBottom = uBottom.*self.dHdx+vBottom.*self.dHdy-self.topographicHeight.*dWdzBottom;
            gBottomFourier = self.horizontalFourierFromSpatialField(gBottom);
            bottomFields = struct(u=uBottom,v=vBottom,dWdz=dWdzBottom);
        end

        function field = spatialFieldFromHorizontalFourier(self,coefficients)
            dft = complex(zeros(self.wvt.Nx,self.wvt.Ny));
            dft(self.primaryDFTIndex) = coefficients;
            dft(self.conjugateDFTIndex) = conj(coefficients(self.conjugateWVIndex));
            field = ifft2(dft,"symmetric")*(self.wvt.Nx*self.wvt.Ny);
        end

        function coefficients = horizontalFourierFromSpatialField(self,field)
            dft = fft2(field)/(self.wvt.Nx*self.wvt.Ny);
            coefficients = reshape(dft(self.primaryDFTIndex),1,[]);
        end

        function requireOriginatingTransform(self,wvt)
            if wvt ~= self.wvt
                error("WVBottomWaveScatteringForcing:TransformMismatch", "The forcing can only be evaluated with the WVTransform instance used during construction.")
            end
        end
    end

    methods (Static, Access = private)
        function [bottomF,bottomGz,pressureProjectionPlus,pressureProjectionMinus] = buildEndpointFactors(wvt)
            [~,iBottom] = min(wvt.z);
            [bottomF,bottomGz] = wvt.waveModeVerticalStructureAtIndex(iBottom);

            pressurePlus = wvt.g*bottomF.*wvt.NAp;
            pressureMinus = wvt.g*bottomF.*wvt.NAm;
            pressureProjectionPlus = complex(zeros(size(wvt.Ap)));
            pressureProjectionMinus = complex(zeros(size(wvt.Am)));
            maskPlus = logical(wvt.waveComponent.maskAp);
            maskMinus = logical(wvt.waveComponent.maskAm);
            pressureProjectionPlus(maskPlus) = conj(pressurePlus(maskPlus))./wvt.Apm_TE_factor(maskPlus);
            pressureProjectionMinus(maskMinus) = conj(pressureMinus(maskMinus))./wvt.Apm_TE_factor(maskMinus);
        end

        function tf = isSupportedTransform(wvt)
            tf = wvt.hasWaveComponent && ismethod(wvt,"waveModeVerticalStructureAtIndex");
        end

        function [primaryDFTIndex,conjugateDFTIndex,conjugateWVIndex] = horizontalFourierIndices(wvt)
            numberOfHorizontalPoints = wvt.Nx*wvt.Ny;
            primaryIndex = reshape(wvt.dftPrimaryIndex,wvt.Nz,[]);
            conjugateIndex = reshape(wvt.dftConjugateIndex,wvt.Nz,[]);
            wvConjugateIndex = reshape(wvt.wvConjugateIndex,wvt.Nz,[]);
            primaryDFTIndex = mod(primaryIndex(1,:)-1,numberOfHorizontalPoints)+1;
            conjugateDFTIndex = mod(conjugateIndex(1,:)-1,numberOfHorizontalPoints)+1;
            conjugateWVIndex = ceil(wvConjugateIndex(1,:)/wvt.Nz);
        end
    end

    methods (Static)
        function diagnostics = bottomDisplacementFromFile(path,options)
            % Reconstruct bottom displacement from standard model output.
            %
            % The method reads accepted `Ap` and `Am` samples from the
            % `wave-vortex` output group, recomputes $$g_b$$, and applies
            % cumulative trapezoidal quadrature to
            % $$\partial_t\eta_d=g_b$$.
            %
            % - Topic: Diagnose bottom displacement
            % - Declaration: diagnostics = bottomDisplacementFromFile(path,options)
            % - Parameter path: standard WaveVortexModel NetCDF output file
            % - Parameter options.forcingName: persisted scattering forcing to use; empty requires exactly one
            % - Parameter options.iTime: unique increasing output indices; empty selects all
            % - Parameter options.initialBottomDisplacement: real bottom displacement at the first selected time
            % - Returns diagnostics: time, bottom velocity, bottom displacement, and quadrature metadata
            arguments (Input)
                path {mustBeFile}
                options.forcingName (1,1) string = ""
                options.iTime double {mustBeInteger,mustBePositive} = double.empty(0,1)
                options.initialBottomDisplacement double = []
            end
            arguments (Output)
                diagnostics (1,1) struct
            end

            [wvt,ncfile] = WVTransform.waveVortexTransformFromFile(char(path),iTime=1,shouldReadOnly=true);
            fileCleanup = onCleanup(@()ncfile.close());
            if ~ncfile.hasGroupWithName("wave-vortex")
                error("WVBottomWaveScatteringForcing:MissingOutputGroup", "The file does not contain the standard 'wave-vortex' output group.")
            end
            group = ncfile.groupWithName("wave-vortex");
            requiredVariables = ["t" "Ap" "Am"];
            for variableName = requiredVariables
                if ~group.hasVariableWithName(variableName)
                    error("WVBottomWaveScatteringForcing:MissingOutputVariable", "The 'wave-vortex' group does not contain the required variable '%s'.", variableName)
                end
            end

            forcingNames = string(wvt.forcingNames());
            isScattering = false(size(forcingNames));
            for iForcing = 1:numel(forcingNames)
                isScattering(iForcing) = isa(wvt.forcingWithName(char(forcingNames(iForcing))),"WVBottomWaveScatteringForcing");
            end
            if strlength(options.forcingName) == 0
                matchingNames = forcingNames(isScattering);
                if numel(matchingNames) ~= 1
                    error("WVBottomWaveScatteringForcing:AmbiguousForcing", "The output must contain exactly one WVBottomWaveScatteringForcing when forcingName is omitted; found %d.", numel(matchingNames))
                end
                forcingName = matchingNames;
            else
                forcingName = options.forcingName;
                if ~any(forcingNames == forcingName)
                    error("WVBottomWaveScatteringForcing:UnknownForcing", "The output does not contain a forcing named '%s'.", forcingName)
                end
                if ~isa(wvt.forcingWithName(char(forcingName)),"WVBottomWaveScatteringForcing")
                    error("WVBottomWaveScatteringForcing:IncorrectForcingType", "The forcing named '%s' is not a WVBottomWaveScatteringForcing.", forcingName)
                end
            end
            forcing = wvt.forcingWithName(char(forcingName));

            allTime = reshape(group.readVariables("t"),[],1);
            if isempty(options.iTime)
                timeIndices = reshape(1:numel(allTime),[],1);
            else
                if ~isvector(options.iTime)
                    error("WVBottomWaveScatteringForcing:InvalidTimeIndices", "iTime must be a vector of output indices.")
                end
                timeIndices = reshape(options.iTime,[],1);
            end
            if any(timeIndices > numel(allTime))
                error("WVBottomWaveScatteringForcing:InvalidTimeIndices", "iTime contains an index larger than the %d saved output times.", numel(allTime))
            end
            if any(diff(timeIndices) <= 0)
                error("WVBottomWaveScatteringForcing:InvalidTimeIndices", "iTime must contain unique strictly increasing indices.")
            end
            time = allTime(timeIndices);
            if any(diff(time) <= 0)
                error("WVBottomWaveScatteringForcing:InvalidOutputTime", "Selected output times must be strictly increasing.")
            end

            if isempty(options.initialBottomDisplacement)
                initialBottomDisplacement = zeros(wvt.Nx,wvt.Ny);
            else
                initialBottomDisplacement = options.initialBottomDisplacement;
                if ~isequal(size(initialBottomDisplacement),[wvt.Nx wvt.Ny]) || ~isreal(initialBottomDisplacement) || any(~isfinite(initialBottomDisplacement),"all")
                    error("WVBottomWaveScatteringForcing:InvalidInitialBottomDisplacement", "initialBottomDisplacement must be a real finite [%d %d] field.", wvt.Nx, wvt.Ny)
                end
            end

            numberOfTimes = numel(time);
            bottomVelocity = zeros(wvt.Nx,wvt.Ny,numberOfTimes);
            bottomDisplacement = zeros(wvt.Nx,wvt.Ny,numberOfTimes);
            bottomVelocityFourier = complex(zeros(wvt.Nkl,numberOfTimes));
            bottomDisplacementFourier = complex(zeros(wvt.Nkl,numberOfTimes));
            bottomDisplacementFourier(:,1) = forcing.horizontalFourierFromSpatialField(initialBottomDisplacement);
            bottomDisplacement(:,:,1) = forcing.spatialFieldFromHorizontalFourier(reshape(bottomDisplacementFourier(:,1),1,[]));
            for iOutput = 1:numberOfTimes
                [Ap,Am] = group.readVariablesAtIndexAlongDimension("t",timeIndices(iOutput),"Ap","Am");
                [bottomVelocity(:,:,iOutput),~,gBottomFourier] = forcing.bottomVelocityFromCoefficients(Ap,Am,time(iOutput));
                bottomVelocityFourier(:,iOutput) = reshape(gBottomFourier,[],1);
                if iOutput > 1
                    deltaT = time(iOutput)-time(iOutput-1);
                    bottomDisplacementFourier(:,iOutput) = bottomDisplacementFourier(:,iOutput-1)+0.5*deltaT*(bottomVelocityFourier(:,iOutput-1)+bottomVelocityFourier(:,iOutput));
                    bottomDisplacement(:,:,iOutput) = forcing.spatialFieldFromHorizontalFourier(reshape(bottomDisplacementFourier(:,iOutput),1,[]));
                end
            end

            diagnostics = struct(time=time,timeIndices=timeIndices,bottomVelocity=bottomVelocity,bottomVelocityFourier=bottomVelocityFourier,bottomDisplacement=bottomDisplacement,bottomDisplacementFourier=bottomDisplacementFourier,forcingName=forcingName,outputInterval=diff(time),quadratureMethod="trapezoidal");
            clear fileCleanup
        end

        function requiredPropertyNames = classRequiredPropertyNames()
            % Return the forcing properties required for restart.
            %
            % - Topic: Restart persistence
            % - Declaration: requiredPropertyNames = classRequiredPropertyNames()
            % - Returns requiredPropertyNames: properties required to reconstruct the forcing
            arguments (Output)
                requiredPropertyNames cell
            end
            requiredPropertyNames = {'topographicHeight','name'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            % Return metadata used to persist the forcing configuration.
            %
            % - Topic: Restart persistence
            % - Declaration: propertyAnnotations = classDefinedPropertyAnnotations()
            % - Returns propertyAnnotations: annotated forcing properties
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CANumericProperty('topographicHeight',{'x','y'},'m','upward-positive topographic height');
            propertyAnnotations(end+1) = CAPropertyAnnotation('name','name of the forcing');
        end

        function forcing = forcingFromGroup(group,wvt)
            % Reconstruct a scattering forcing from its annotated group.
            %
            % - Topic: Restart persistence
            % - Declaration: forcing = forcingFromGroup(group,wvt)
            % - Parameter group: NetCDF group containing the forcing state
            % - Parameter wvt: restored transform receiving the forcing
            % - Returns forcing: reconstructed `WVBottomWaveScatteringForcing`
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            arguments (Output)
                forcing WVBottomWaveScatteringForcing
            end

            requiredPropertyNames = WVBottomWaveScatteringForcing.classRequiredPropertyNames();
            missingPropertyNames = string.empty(1,0);
            for iProperty = 1:numel(requiredPropertyNames)
                propertyName = requiredPropertyNames{iProperty};
                isPresent = group.hasVariableWithName(propertyName) || group.hasGroupWithName(propertyName) || isKey(group.attributes,propertyName);
                if ~isPresent
                    missingPropertyNames(end+1) = string(propertyName); %#ok<AGROW>
                end
            end
            if ~isempty(missingPropertyNames)
                error("WVBottomWaveScatteringForcing:IncompleteRestart", "The restart group is missing required wave-scattering properties: %s.", join(missingPropertyNames,", "))
            end

            options = CAAnnotatedClass.propertyValuesFromGroup(group,requiredPropertyNames);
            options.name = string(options.name);
            optionArguments = namedargs2cell(options);
            forcing = WVBottomWaveScatteringForcing(wvt,optionArguments{:});
        end
    end
end
