function audit = buildBoundaryGreenIdentityAudit(problem,bottomSlope)
% Classify the invariant branch of the local constant-slope descriptor.

descriptor = buildBoundaryDynamicalDescriptorAudit(problem,bottomSlope);
flatDescriptor = buildBoundaryDynamicalDescriptorAudit(problem,[0 0]);
slopeMagnitude = norm(bottomSlope);
tangentStep = 1e-4;
tangentConvergence = 0;
if slopeMagnitude > 0
    direction = bottomSlope/slopeMagnitude;
    plusDescriptor = buildBoundaryDynamicalDescriptorAudit(problem,tangentStep*direction);
    minusDescriptor = buildBoundaryDynamicalDescriptorAudit(problem,-tangentStep*direction);
    plusHalfDescriptor = buildBoundaryDynamicalDescriptorAudit(problem,(tangentStep/2)*direction);
    minusHalfDescriptor = buildBoundaryDynamicalDescriptorAudit(problem,-(tangentStep/2)*direction);
    for iK = 1:numel(descriptor.blocks)
        derivative = (plusDescriptor.blocks{iK}.publicGenerator-minusDescriptor.blocks{iK}.publicGenerator)/(2*tangentStep);
        derivativeHalf = (plusHalfDescriptor.blocks{iK}.publicGenerator-minusHalfDescriptor.blocks{iK}.publicGenerator)/tangentStep;
        tangentConvergence = max(tangentConvergence, ...
            norm(derivative-derivativeHalf,"fro")/max(norm(derivativeHalf,"fro"),realmin));
        L = flatDescriptor.blocks{iK}.localForms.generator+slopeMagnitude*derivativeHalf;
        descriptor.blocks{iK}.localForms.rawGenerator = descriptor.blocks{iK}.localForms.generator;
        descriptor.blocks{iK}.localForms.firstOrderGenerator = L;
        descriptor.blocks{iK}.localForms.generator = L;
        descriptor.blocks{iK}.localForms.diagnostics = updateGeneratorDiagnostics( ...
            descriptor.blocks{iK}.localForms,descriptor.blocks{iK}.localForms.diagnostics);
    end
else
    for iK = 1:numel(descriptor.blocks)
        descriptor.blocks{iK}.localForms.rawGenerator = descriptor.blocks{iK}.publicGenerator;
        descriptor.blocks{iK}.localForms.firstOrderGenerator = descriptor.blocks{iK}.localForms.generator;
    end
end

names = ["energyHermitianDefect","exchangeSkewDefect","weakEvolutionDefect", ...
    "energySkewDefect","apvTendencyDefect","potentialEnstrophyDefect", ...
    "bottomEvolutionDefect","pressureGreenDefect","maximumRealGrowthRate"];
maximum = struct;
for name = names
    maximum.(name) = max(cellfun(@(block)block.localForms.diagnostics.(name),descriptor.blocks));
end
minimumEnergyEigenvalue = min(cellfun(@(block)block.localForms.diagnostics.minimumEnergyEigenvalue,descriptor.blocks));
maximumEnergyConditionNumber = max(cellfun(@(block)block.localForms.diagnostics.energyConditionNumber,descriptor.blocks));

activeBoundaryModeCount = 0;
maximumBoundaryParticipation = 0;
flatBoundaryInvariantDefect = 0;
slopeBoundaryInvariantDefect = 0;
for iK = 1:numel(descriptor.blocks)
    block = descriptor.blocks{iK};
    forms = block.localForms;
    [C,lambda] = eig(forms.generator,"vector");
    frequency = real(1i*lambda);
    B = forms.bottomValue;
    frequencyScale = max([abs(problem.originatingTransform.f);abs(frequency)]);
    oscillatory = abs(frequency) > 1000*eps*max(frequencyScale,1);
    participation = abs(B*C(:,oscillatory))./max(vecnorm(C(:,oscillatory)),realmin);
    activeBoundaryModeCount = activeBoundaryModeCount+nnz(participation > 1e-10);
    if ~isempty(participation)
        maximumBoundaryParticipation = max(maximumBoundaryParticipation,max(participation));
    end

    Gbottom = B'*B;
    slopeBoundaryInvariantDefect = max(slopeBoundaryInvariantDefect, ...
        invariantDefect(forms.generator,Gbottom));
    flatForms = flatDescriptor.blocks{iK}.localForms;
    flatBoundaryInvariantDefect = max(flatBoundaryInvariantDefect, ...
        invariantDefect(flatForms.generator,Gbottom));
end

metricNames = names(1:end-1);
failedConditions = metricNames(arrayfun(@(name)maximum.(name) > 1e-11,metricNames)).';
if ~descriptor.isCompatible
    failedConditions = [failedConditions;descriptor.failedConditions];
end
frequencyScale = max([abs(problem.originatingTransform.f); ...
    cellfun(@(block)max(abs(block.frequency)),descriptor.blocks)]);
if maximum.maximumRealGrowthRate > 1e-11*max(frequencyScale,1)
    failedConditions(end+1,1) = "complexFrequency";
end
if minimumEnergyEigenvalue <= 0
    failedConditions(end+1,1) = "nonpositivePhysicalEnergy";
end
if any(bottomSlope ~= 0) && activeBoundaryModeCount == 0
    failedConditions(end+1,1) = "inactiveBoundaryCoordinate";
end
failedConditions = unique(failedConditions,"stable");

branch = "P";
description = "physical-energy-only";
if ~isempty(failedConditions)
    branch = "incompatible";
    description = "incompatible";
end

audit = struct;
audit.status = description;
audit.branch = branch;
audit.continuumBranch = "P";
audit.continuumBranchDescription = "physical-energy-only";
audit.isCompatible = isempty(failedConditions);
audit.failedConditions = failedConditions;
audit.bottomSlope = bottomSlope;
audit.descriptor = descriptor;
audit.diagnostics = struct("maximum",maximum, ...
    "minimumEnergyEigenvalue",minimumEnergyEigenvalue, ...
    "maximumEnergyConditionNumber",maximumEnergyConditionNumber, ...
    "flatBottomValueQuadraticDefect",flatBoundaryInvariantDefect, ...
    "slopeBottomValueQuadraticDefect",slopeBoundaryInvariantDefect, ...
    "activeBoundaryModeCount",activeBoundaryModeCount, ...
    "maximumBoundaryParticipation",maximumBoundaryParticipation, ...
    "tangentStep",tangentStep, ...
    "tangentConvergence",tangentConvergence, ...
    "generalizedBoundaryMetricRank",0, ...
    "generalizedBoundaryMetricDefiniteness","absent");
end

function diagnostics = updateGeneratorDiagnostics(forms,diagnostics)
L = forms.generator;
E = forms.energyMatrix;
J = forms.exchangeMatrix;
Q = forms.apvMatrix;
Z = forms.potentialEnstrophyMatrix;
B = forms.bottomValue;
R = forms.bottomTendency;
diagnostics.weakEvolutionDefect = norm(E*L-J,"fro")/max(norm(E*L,"fro")+norm(J,"fro"),realmin);
diagnostics.energySkewDefect = norm(L'*E+E*L,"fro")/max(norm(E,"fro")*norm(L,"fro"),realmin);
diagnostics.apvTendencyDefect = norm(Q*L,"fro")/max(norm(Q,"fro")*norm(L,"fro"),realmin);
diagnostics.potentialEnstrophyDefect = norm(L'*Z+Z*L,"fro")/max(norm(Z,"fro")*norm(L,"fro"),realmin);
diagnostics.bottomEvolutionDefect = norm(B*L-R,"fro")/max(norm(B,"fro")*norm(L,"fro")+norm(R,"fro"),realmin);
diagnostics.maximumRealGrowthRate = max(abs(real(eig(L))));
end

function value = invariantDefect(L,G)
value = norm(L'*G+G*L,"fro")/max(norm(G,"fro")*norm(L,"fro"),realmin);
end
