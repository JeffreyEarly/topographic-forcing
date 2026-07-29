function audit = buildDealiasedProjectedPrimitiveTangentAudit(problem,trustedBounds,supportBounds,degrees,paddingFactors,quadratureOrder,tangentStep)
% Build the dealiased projected primitive tangent refinement oracle.

layouts = retainedLayouts(problem,supportBounds);
terrain = terrainSupport(problem);
validateTrustedGuard(layouts{1},trustedBounds,terrain);

nLevel = size(supportBounds,1);
diagonal = repmat(emptySummary,0,1);
padding = repmat(emptyPaddingComparison,0,1);
for iLevel = 1:nLevel
    [auditTwo,summaryTwo] = runAudit(problem,layouts{iLevel},trustedBounds, ...
        degrees(iLevel),paddingFactors(1),quadratureOrder,tangentStep);
    diagonal(end+1,1) = summaryTwo; %#ok<AGROW>

    comparison = emptyPaddingComparison;
    comparison.supportModeBounds = supportBounds(iLevel,:);
    comparison.polynomialDegree = degrees(iLevel);
    for iPadding = 2:numel(paddingFactors)
        [auditOther,summaryOther] = runAudit(problem,layouts{iLevel},trustedBounds, ...
            degrees(iLevel),paddingFactors(iPadding),quadratureOrder,tangentStep);
        comparison.paddingFactors(end+1,1) = paddingFactors(iPadding);
        action = paddingActionDefects(auditTwo,auditOther);
        comparison.actionDefects(end+1,1) = action.maximumDefect;
        comparison.actionDefectsByField(end+1,1) = action;
        comparison.projectedResidualDefects(end+1,1) = paddingResidualDefect( ...
            auditTwo,auditOther);
        comparison.maximumDefect = max([comparison.maximumDefect; ...
            comparison.actionDefects(end);comparison.projectedResidualDefects(end)]);
        comparison.summaries(end+1,1) = summaryOther;
    end
    padding(end+1,1) = comparison; %#ok<AGROW>
end

horizontal = repmat(emptySummary,nLevel,1);
vertical = repmat(emptySummary,nLevel,1);
horizontal(end) = diagonal(end);
vertical(end) = diagonal(end);
for iLevel = 1:nLevel-1
    [~,horizontal(iLevel)] = runAudit(problem,layouts{iLevel},trustedBounds, ...
        degrees(end),paddingFactors(1),quadratureOrder,tangentStep);
    [~,vertical(iLevel)] = runAudit(problem,layouts{end},trustedBounds, ...
        degrees(iLevel),paddingFactors(1),quadratureOrder,tangentStep);
end

structural = structuralDiagnostics(diagonal,horizontal,vertical,padding);
convergence = convergenceDiagnostics(diagonal,horizontal,vertical);
tolerance = struct("projection",1e-12,"tangent",1e-9,"weak",1e-11, ...
    "energy",1e-11,"bottom",1e-11,"conjugacy",1e-11,"padding",1e-10, ...
    "exactProjected",1e-10,"convergentProjected",1e-8,"refinementFactor",4);
structurePasses = structural.maximumProjectionDefect <= tolerance.projection ...
    && structural.maximumTangentDefect <= tolerance.tangent ...
    && structural.maximumWeakDefect <= tolerance.weak ...
    && structural.maximumEnergyDefect <= tolerance.energy ...
    && structural.maximumBottomDefect <= tolerance.bottom ...
    && structural.maximumConjugacyDefect <= tolerance.conjugacy ...
    && structural.maximumPaddingDefect <= tolerance.padding;
finest = diagonal(end);
exactPasses = finest.projectedAPVDefect <= tolerance.exactProjected ...
    && finest.projectedEnstrophyDefect <= tolerance.exactProjected;
convergencePasses = all(convergence.diagonalAPVReduction >= tolerance.refinementFactor) ...
    && all(convergence.diagonalEnstrophyReduction >= tolerance.refinementFactor) ...
    && finest.trustedAPVDefect <= tolerance.convergentProjected ...
    && finest.trustedEnstrophyDefect <= tolerance.convergentProjected;

if ~structurePasses
    status = "implementation-unresolved";
    diagnosis = "The common projected primitive construction failed an adjointness, convolution, tangent, weak, energy, bottom, conjugacy, or padding-independence gate.";
elseif exactPasses
    status = "exact-projected";
    diagnosis = "The unmodified primitive tangent conserves the declared projected APV and potential enstrophy to the exact-projection tolerance.";
elseif convergencePasses
    status = "convergent-projected";
    diagnosis = "The unmodified primitive tangent is not exactly APV closed at finite support, but its trusted-band APV and potential-enstrophy defects meet the required refinement rate and final tolerance.";
else
    status = "nonconvergent-projected";
    diagnosis = "The common dealiased primitive projection preserves weak evolution, physical energy, bottom evolution, and conjugacy, but trusted-band APV or potential enstrophy does not meet the required convergence gate.";
end

audit = struct;
audit.scope = "dealiased-projected-primitive-tangent-only";
audit.status = status;
audit.isCompatible = status == "exact-projected" || status == "convergent-projected";
audit.diagnosis = diagnosis;
audit.trustedModeBounds = trustedBounds;
audit.supportModeBounds = supportBounds;
audit.polynomialDegrees = degrees;
audit.paddingFactors = paddingFactors;
audit.quadratureOrder = quadratureOrder;
audit.tangentStep = tangentStep;
audit.terrain = terrain;
audit.diagonalRefinement = diagonal;
audit.horizontalRefinement = horizontal;
audit.verticalRefinement = vertical;
audit.paddingComparisons = padding;
audit.structural = structural;
audit.convergence = convergence;
audit.requiredTolerance = tolerance;
audit.externalSidebands = externalSidebandDiagnostics(diagonal);
audit.nextScope = "milestone-7-only-if-exact-or-convergent-projected";
end

function [result,summary] = runAudit(problem,layout,trustedBounds,degree,paddingFactor,quadratureOrder,tangentStep)
if isempty(quadratureOrder)
    order = max(2*degree+7,18);
else
    order = quadratureOrder;
end
result = buildGlobalSmallTerrainPrimitiveAudit(problem,degree,order,tangentStep, ...
    horizontalLayout=layout,paddingFactor=paddingFactor, ...
    trustedModeBounds=trustedBounds,rejectTerrainNyquist=true);
summary = summarizeAudit(result);
end

function summary = summarizeAudit(audit)
summary = emptySummary;
summary.supportModeBounds = [max(abs(audit.layout.horizontalLayout.kMode)) ...
    max(abs(audit.layout.horizontalLayout.lMode))];
summary.numberOfHorizontalModes = audit.layout.numberOfHorizontalModes;
summary.numberOfTrustedHorizontalModes = nnz(audit.layout.trustedHorizontalModes);
summary.polynomialDegree = audit.polynomialDegree;
summary.quadratureOrder = audit.quadratureOrder;
summary.paddingFactor = audit.layout.paddingFactor;
summary.projectionDefect = max([audit.projection.adjointDefect ...
    audit.projection.restrictionIdentityDefect ...
    audit.projection.maximumExactConvolutionDefect]);
summary.tangentDefect = audit.tangentAgreement.maximumRelativeDefect;
summary.weakDefect = audit.compatibility.weakEvolutionDefect;
summary.energyDefect = audit.compatibility.energyDefect;
summary.bottomDefect = audit.compatibility.bottomDefect;
summary.conjugacyDefect = audit.compatibility.conjugacyDefect;
summary.projectedAPVDefect = audit.projectedCompatibility.apvDefect;
summary.trustedAPVDefect = audit.projectedCompatibility.trustedAPVDefect;
summary.projectedEnstrophyDefect = audit.projectedCompatibility.enstrophyDefect;
summary.trustedEnstrophyDefect = audit.projectedCompatibility.trustedEnstrophyDefect;
summary.externalAPVDefect = audit.projectedCompatibility.externalAPVDefect;
summary.numberOfEdgeHorizontalModes = audit.projection.numberOfEdgeHorizontalModes;
summary.maximumDiscardedOperatorNorm = audit.projection.maximumDiscardedOperatorNorm;
summary.exactConvolutionDefects = audit.projection.exactConvolutionDefects;
end

function summary = emptySummary
summary = struct("supportModeBounds",[0 0],"numberOfHorizontalModes",0, ...
    "numberOfTrustedHorizontalModes",0,"polynomialDegree",0,"quadratureOrder",0, ...
    "paddingFactor",0,"projectionDefect",NaN,"tangentDefect",NaN, ...
    "weakDefect",NaN,"energyDefect",NaN,"bottomDefect",NaN, ...
    "conjugacyDefect",NaN,"projectedAPVDefect",NaN,"trustedAPVDefect",NaN, ...
    "projectedEnstrophyDefect",NaN,"trustedEnstrophyDefect",NaN, ...
    "externalAPVDefect",NaN,"numberOfEdgeHorizontalModes",0, ...
    "maximumDiscardedOperatorNorm",0,"exactConvolutionDefects",nan(1,3));
end

function comparison = emptyPaddingComparison
comparison = struct("supportModeBounds",[0 0],"polynomialDegree",0, ...
    "paddingFactors",2,"actionDefects",zeros(0,1), ...
    "actionDefectsByField",repmat(emptyActionDefects,0,1), ...
    "projectedResidualDefects",zeros(0,1),"maximumDefect",0, ...
    "summaries",repmat(emptySummary,0,1));
end

function diagnostics = paddingActionDefects(first,second)
trusted = first.layout.trustedColumns;
names = ["flatL","L1","H1","J1","Q1","R1"];
pairs = {
    first.flatReference.L(:,trusted),second.flatReference.L(:,trusted)
    first.analyticTangent.L(:,trusted),second.analyticTangent.L(:,trusted)
    first.analyticTangent.E(:,trusted),second.analyticTangent.E(:,trusted)
    first.analyticTangent.J(:,trusted),second.analyticTangent.J(:,trusted)
    first.analyticTangent.QProjected(:,trusted),second.analyticTangent.QProjected(:,trusted)
    first.analyticTangent.R(:,trusted),second.analyticTangent.R(:,trusted)
    };
values = zeros(numel(names),1);
for iPair = 1:size(pairs,1)
    values(iPair) = matrixPairDefect(pairs{iPair,1},pairs{iPair,2});
end
diagnostics = cell2struct(num2cell(values),cellstr(names),1);
diagnostics.maximumDefect = max(values);
end

function diagnostics = emptyActionDefects
diagnostics = struct("flatL",NaN,"L1",NaN,"H1",NaN,"J1",NaN, ...
    "Q1",NaN,"R1",NaN,"maximumDefect",NaN);
end

function value = matrixPairDefect(first,second)
value = norm(first-second,"fro")/max(norm(first,"fro")+norm(second,"fro"),realmin);
end

function value = paddingResidualDefect(first,second)
trusted = first.layout.trustedColumns;
firstResidual = first.projectedCompatibility.apvResidual(:,trusted);
secondResidual = second.projectedCompatibility.apvResidual(:,trusted);
firstTerms = first.projectedCompatibility.apvTerms;
secondTerms = second.projectedCompatibility.apvTerms;
scale = 0;
for iTerm = 1:numel(firstTerms)
    scale = scale+norm(firstTerms{iTerm}(:,trusted),"fro") ...
        +norm(secondTerms{iTerm}(:,trusted),"fro");
end
value = norm(firstResidual-secondResidual,"fro")/max(scale,realmin);
end

function diagnostics = structuralDiagnostics(diagonal,horizontal,vertical,padding)
allSummaries = [diagonal;horizontal;vertical];
diagnostics = struct( ...
    "maximumProjectionDefect",max([allSummaries.projectionDefect]), ...
    "maximumTangentDefect",max([allSummaries.tangentDefect]), ...
    "maximumWeakDefect",max([allSummaries.weakDefect]), ...
    "maximumEnergyDefect",max([allSummaries.energyDefect]), ...
    "maximumBottomDefect",max([allSummaries.bottomDefect]), ...
    "maximumConjugacyDefect",max([allSummaries.conjugacyDefect]), ...
    "maximumPaddingDefect",max([padding.maximumDefect]));
end

function diagnostics = convergenceDiagnostics(diagonal,horizontal,vertical)
diagnostics = struct( ...
    "diagonalAPVReduction",reduction([diagonal.trustedAPVDefect]), ...
    "diagonalEnstrophyReduction",reduction([diagonal.trustedEnstrophyDefect]), ...
    "horizontalAPVReduction",reduction([horizontal.trustedAPVDefect]), ...
    "horizontalEnstrophyReduction",reduction([horizontal.trustedEnstrophyDefect]), ...
    "verticalAPVReduction",reduction([vertical.trustedAPVDefect]), ...
    "verticalEnstrophyReduction",reduction([vertical.trustedEnstrophyDefect]));
end

function ratio = reduction(values)
ratio = values(1:end-1)./max(values(2:end),realmin);
end

function diagnostics = externalSidebandDiagnostics(diagonal)
diagnostics = struct( ...
    "isPresent",any([diagonal.numberOfEdgeHorizontalModes] > 0) ...
        && any([diagonal.maximumDiscardedOperatorNorm] > 0), ...
    "numberOfEdgeHorizontalModes",[diagonal.numberOfEdgeHorizontalModes].', ...
    "maximumDiscardedOperatorNorm",[diagonal.maximumDiscardedOperatorNorm].', ...
    "projectedComplementAPVDefect",[diagonal.externalAPVDefect].');
end

function layouts = retainedLayouts(problem,bounds)
source = problem.horizontalLayout;
layouts = cell(size(bounds,1),1);
for iLevel = 1:size(bounds,1)
    mask = abs(source.kMode) <= bounds(iLevel,1) ...
        & abs(source.lMode) <= bounds(iLevel,2);
    layout = source(mask,:);
    if isempty(layout)
        error("WVTerrainEnergyGalerkin:EmptyProjectedPrimitiveSupport", ...
            "A requested supportModeBounds row retains no horizontal modes.")
    end
    expected = (2*bounds(iLevel,1)+1)*(2*bounds(iLevel,2)+1);
    if height(layout) ~= expected
        error("WVTerrainEnergyGalerkin:UnavailableProjectedPrimitiveSupport", ...
            "The originating transform does not retain the complete signed rectangle [%d %d].", ...
            bounds(iLevel,1),bounds(iLevel,2))
    end
    layouts{iLevel} = layout;
end
end

function terrain = terrainSupport(problem)
wvt = problem.originatingTransform;
spectrum = fft2(problem.topographicHeight)/(wvt.Nx*wvt.Ny);
nyquist = logical(WVGeometryDoublyPeriodic.maskForNyquistModes(wvt.Nx,wvt.Ny));
scale = max(abs(spectrum),[],"all");
tolerance = 100*eps*max(scale,1);
if max(abs(spectrum(nyquist)),[],"all") > tolerance
    error("WVTerrainEnergyGalerkin:ProjectedPrimitiveTerrainNyquist", ...
        "The dealiased projected primitive audit requires terrain with no material Nyquist coefficient.")
end
[kMode,lMode] = ndgrid(wvt.kMode_dft,wvt.lMode_dft);
active = abs(spectrum) > tolerance & ~nyquist;
terrain = struct("kMode",kMode(active),"lMode",lMode(active), ...
    "coefficient",spectrum(active),"spectralTolerance",tolerance);
end

function validateTrustedGuard(smallestLayout,trustedBounds,terrain)
trusted = abs(smallestLayout.kMode) <= trustedBounds(1) ...
    & abs(smallestLayout.lMode) <= trustedBounds(2);
if ~any(trusted)
    error("WVTerrainEnergyGalerkin:EmptyProjectedPrimitiveTrustedBand", ...
        "trustedModeBounds retains no horizontal modes.")
end
for iInput = find(trusted).'
    destinations = [smallestLayout.kMode(iInput)+terrain.kMode, ...
        smallestLayout.lMode(iInput)+terrain.lMode];
    for iDestination = 1:size(destinations,1)
        if ~any(smallestLayout.kMode == destinations(iDestination,1) ...
                & smallestLayout.lMode == destinations(iDestination,2))
            error("WVTerrainEnergyGalerkin:InsufficientProjectedPrimitiveGuard", ...
                "The smallest support does not contain every first terrain sideband of the trusted band.")
        end
    end
end
end
