function [wvt,configuration] = gaussianRidgeScatteringConfiguration(resolution,domainWavelengths,shouldAntialias)
% Create the constant-stratification Gaussian-ridge scattering transform.

depth = 5e3;
M2Frequency = 2*pi/(12.4206012*3600);
coriolisFrequency = M2Frequency/2;
buoyancyFrequency = 10*coriolisFrequency;
rotationRate = 7.2921e-5;
latitude = asind(coriolisFrequency/(2*rotationRate));
verticalWavenumber = pi/depth;
mu = sqrt((M2Frequency^2-coriolisFrequency^2)/(buoyancyFrequency^2-M2Frequency^2));
firstModeWavelength = 2*depth/mu;
domainSize = [domainWavelengths*firstModeWavelength 20e3 depth];
carrierMode = domainWavelengths;

wvt = WVTransformBoussinesq(domainSize,resolution,N2=@(z)buoyancyFrequency^2*ones(size(z)),latitude=latitude,rotationRate=rotationRate,shouldAntialias=shouldAntialias);
wvt.t0 = 0;
wvt.t = 0;
maximumRetainedMode = max(round(abs(wvt.K(:))*wvt.Lx/(2*pi)));
if carrierMode > maximumRetainedMode
    error("GaussianRidgeWaveScatteringExample:UnresolvedCarrier", "The transform retains horizontal modes only through k=%d, below the requested carrier k=%d.",maximumRetainedMode,carrierMode)
end
carrierIndex = wvt.indexFromModeNumber(carrierMode,0,1);
horizontalWavenumber = abs(wvt.K(carrierIndex));
actualFrequency = wvt.Omega(carrierIndex);
groupVelocity = horizontalWavenumber*verticalWavenumber^2*(buoyancyFrequency^2-coriolisFrequency^2)/(actualFrequency*(horizontalWavenumber^2+verticalWavenumber^2)^2);

configuration = struct(domainSize=domainSize,resolution=resolution,depth=depth,latitude=latitude,rotationRate=rotationRate,coriolisFrequency=coriolisFrequency,buoyancyFrequency=buoyancyFrequency,targetFrequency=M2Frequency,actualFrequency=actualFrequency,verticalWavenumber=verticalWavenumber,horizontalWavenumber=horizontalWavenumber,mu=mu,firstModeWavelength=firstModeWavelength,domainWavelengths=domainWavelengths,carrierMode=carrierMode,maximumRetainedMode=maximumRetainedMode,groupVelocity=groupVelocity,shouldAntialias=shouldAntialias);
end
