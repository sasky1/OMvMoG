function [model] = randn_warmstart(X,r,k)
% Use reweighted L2 matrix factorization to calculate U
% Use mog to calculatethe model
 disp('Warm-start...');

 U=[];
 if r>1
[L,~] = inexact_alm_rpca(X,0.5/sqrt(max(size(X,1),size(X,2))));
[U,~] = PCA(L,r-1);
 end
 U=[U,median(X,2)];
m=size(X,1);
A=repmat(eye(r),[1,1,m])*0.1;
B=U'/0.1;
Sigma=zeros(1,k);Sigma(1)=0.025;
for i=2:k
 Sigma(i)=Sigma(i-1)/15;
end

model.Sigma=Sigma;
model.weight=ones(1,k)/k;
model.mu=zeros(1,k);
model.A=A;
model.B=B;
model.U=U;
end