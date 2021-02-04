%%%%% Weighted L2 norm factorization %%%%%
function [OutV] = EfficientMCL2_V(Matrix, W, InU,InV, MaxIt, tol)
%%%%% Matrix - data matrix for factorization (d*n) %%%%%
%%%%% W - weight matrix (d*n) %%%%%
%%%%% InU,InV - initialization of U (d*k dimensional) and V (n*k dimensional)     %%%%%
%%%%% MaxIt - Maximal iteration number      %%%%%
%%%%% tol - Tolerance for RobustL1 algorithm       %%%%%
if nargin<5
    MaxIt=20;
end
if nargin<6
    tol=0.000001;
end



[k] = size(InU,2);
OutU = InU;
OutV = InV;
for i = 1:MaxIt
     ind = randperm(k);
%    ind = 1:k;
    for j = ind
        TX = Matrix - OutU*OutV' + OutU(:,j)*OutV(:,j)';
        u = InU(:,j);
        OutV(:,j) = optimMCL2(TX,W,u);    
%         OutU(:,j) = optimMCL2(TX',W',OutV(:,j));     
    end
    if norm(InV - OutV) < tol
        break;
    else
        InU = OutU;
    end
end
 
% Nu = sqrt(sum(OutU.^2))';
% Nv = sqrt(sum(OutV.^2))';
% No = diag(Nu.*Nv);
% OutU = OutU*diag(1./Nu)*sqrt(No);
% OutV = OutV*diag(1./Nv)*sqrt(No);

