function [P,Pr] = legendreValuesForPrimitiveAudit(r,degree)
% Evaluate Legendre polynomials and reference-coordinate derivatives.
P = zeros(numel(r),degree+1);
Pr = zeros(numel(r),degree+1);
P(:,1) = 1;
if degree == 0
    return
end
P(:,2) = r;
Pr(:,2) = 1;
for n = 2:degree
    P(:,n+1) = ((2*n-1)*r.*P(:,n)-(n-1)*P(:,n-1))/n;
    Pr(:,n+1) = ((2*n-1)*(P(:,n)+r.*Pr(:,n)) ...
        -(n-1)*Pr(:,n-1))/n;
end
end
