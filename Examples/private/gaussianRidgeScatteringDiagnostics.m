function diagnostics = gaussianRidgeScatteringDiagnostics(wvt,forcing,Ap,Am,A0,time)
% Diagnose modal partitions and first-order energy from sampled coefficients.

numberOfTimes = numel(time);
plusEnergy = zeros([size(wvt.Ap) numberOfTimes]);
minusEnergy = zeros([size(wvt.Am) numberOfTimes]);
balancedEnergy = zeros(numberOfTimes,1);
qgpvNorm = zeros(numberOfTimes,1);
firstOrderEnergy = zeros(numberOfTimes,1);
bottomCorrection = zeros(numberOfTimes,1);
for iTime = 1:numberOfTimes
    wvt.Ap = Ap(:,:,iTime);
    wvt.Am = Am(:,:,iTime);
    wvt.A0 = A0(:,:,iTime);
    wvt.t = time(iTime);
    plusEnergy(:,:,iTime) = wvt.Apm_TE_factor.*abs(wvt.Ap).^2;
    minusEnergy(:,:,iTime) = wvt.Apm_TE_factor.*abs(wvt.Am).^2;
    balancedEnergy(iTime) = sum(wvt.A0_TE_factor(:).*abs(wvt.A0(:)).^2);
    qgpvNorm(iTime) = norm(wvt.A0_QGPV_factor(:).*wvt.A0(:));
    [~,bottomFields] = forcing.bottomVelocityFromWaveState(wvt);
    bottomCorrection(iTime) = 0.5*mean(forcing.topographicHeight.*(bottomFields.u.^2+bottomFields.v.^2),"all");
    firstOrderEnergy(iTime) = sum(plusEnergy(:,:,iTime),"all")+sum(minusEnergy(:,:,iTime),"all")-bottomCorrection(iTime);
end

positiveK = wvt.K > 0;
negativeK = wvt.K < 0;
firstMode = wvt.J == 1;
higherMode = wvt.J >= 2;
rightwardFirstModeEnergy = squeeze(sum(plusEnergy.*(negativeK & firstMode),[1 2])+sum(minusEnergy.*(positiveK & firstMode),[1 2]));
leftwardFirstModeEnergy = squeeze(sum(plusEnergy.*(positiveK & firstMode),[1 2])+sum(minusEnergy.*(negativeK & firstMode),[1 2]));
rightwardHigherModeEnergy = squeeze(sum(plusEnergy.*(negativeK & higherMode),[1 2])+sum(minusEnergy.*(positiveK & higherMode),[1 2]));
leftwardHigherModeEnergy = squeeze(sum(plusEnergy.*(positiveK & higherMode),[1 2])+sum(minusEnergy.*(negativeK & higherMode),[1 2]));
waveEnergy = reshape(sum(plusEnergy+minusEnergy,[1 2]),[],1);
scatteredCoefficientNorm = zeros(numberOfTimes,1);
for iTime = 1:numberOfTimes
    scatteredCoefficientNorm(iTime) = sqrt(sum(wvt.Apm_TE_factor.*(abs(Ap(:,:,iTime)-Ap(:,:,1)).^2+abs(Am(:,:,iTime)-Am(:,:,1)).^2),"all"));
end

initialEnergy = waveEnergy(1);
diagnostics = struct(time=reshape(time,[],1),waveEnergy=waveEnergy,firstOrderEnergy=firstOrderEnergy,bottomEnergyCorrection=bottomCorrection,relativeFlatEnergyChange=(waveEnergy-initialEnergy)/initialEnergy,relativeFirstOrderEnergyChange=(firstOrderEnergy-firstOrderEnergy(1))/initialEnergy,rightwardFirstModeEnergy=rightwardFirstModeEnergy,leftwardFirstModeEnergy=leftwardFirstModeEnergy,rightwardHigherModeEnergy=rightwardHigherModeEnergy,leftwardHigherModeEnergy=leftwardHigherModeEnergy,higherModeEnergy=rightwardHigherModeEnergy+leftwardHigherModeEnergy,balancedEnergy=balancedEnergy,qgpvNorm=qgpvNorm,scatteredCoefficientNorm=scatteredCoefficientNorm,initialEnergy=initialEnergy);
end
