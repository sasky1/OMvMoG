function [model] = t_warmstart(X,k)
% Use reweighted L2 matrix factorization to calculate U
% Use mog to calculatethe model
 disp('Warm-start...');
 if nargin<4
 ref=1;
 end
 disp('Align image...');
 tau = zeros(6,size(X,3));
 if size(X,3)>1
 [X] = preAlign(X,ref,tau); 
X=reshape(X,[size(X,1)*size(X,2),size(X,3)]);
 else
     X=X(:);
 end
m=size(X,1);
A=repmat(eye(size(X,2)),[1,1,m])/20;
B=X'*20;
Sigma=zeros(1,k);Sigma(1)=0.025;
for i=2:k
 Sigma(i)=Sigma(i-1)/15;
end

model.Sigma=Sigma;
model.weight=ones(1,k)/k;
model.mu=zeros(1,k);
model.A=A;
model.B=B;
model.U=X;
end