function [W,diagnostics] = buildRealConjugacyCoordinates(conjugateIndex)
% Build a unitary map from independent real to full complex coordinates.

conjugateIndex = conjugateIndex(:);
n = numel(conjugateIndex);
if any(conjugateIndex < 1 | conjugateIndex > n | conjugateIndex ~= round(conjugateIndex)) || any(conjugateIndex(conjugateIndex) ~= (1:n)')
    error("WVTerrainEnergyGalerkin:InvalidConjugacyMap", ...
        "The conjugate-coordinate map must be a complete involution.")
end

W = zeros(n);
visited = false(n,1);
column = 0;
for i = 1:n
    if visited(i)
        continue
    end
    partner = conjugateIndex(i);
    if partner == i
        column = column+1;
        W(i,column) = 1;
        visited(i) = true;
    else
        column = column+1;
        W([i partner],column) = 1/sqrt(2);
        column = column+1;
        W(i,column) = 1i/sqrt(2);
        W(partner,column) = -1i/sqrt(2);
        visited([i partner]) = true;
    end
end
if column ~= n
    error("WVTerrainEnergyGalerkin:IncompleteRealCoordinateMap", ...
        "The real-coordinate map contains %d columns instead of %d.",column,n)
end

C = sparse((1:n)',conjugateIndex,1,n,n);
diagnostics = struct;
diagnostics.unitaryDefect = norm(W'*W-eye(n),"fro")/sqrt(n);
diagnostics.conjugacyDefect = norm(C*conj(W)-W,"fro")/sqrt(n);
diagnostics.numberOfSelfConjugateCoordinates = nnz(conjugateIndex == (1:n)');
diagnostics.numberOfConjugatePairs = (n-diagnostics.numberOfSelfConjugateCoordinates)/2;
end
