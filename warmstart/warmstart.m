function [model] = warmstart(X,r,k)
% Use reweighted L2 matrix factorization to calculate U
% Use mog to calculatethe model
 disp('Warm-start...');
 disp('Calculating U...');
W=ones(size(X));
%RPCA for warmstart
U=[];
% if r>1
[L,~] = inexact_alm_rpca(X,0.5/sqrt(max(size(X,1),size(X,2))));
 [U,V] = PCA(L,r);%U=[U,median(X,2)];V=X'*U/(U'*U+0.0001*eye(r));
%PCA for warmstart
% [U,~] = PCA(X,r-1);
% end
% iter=10;
%  disp('Warm-start...');
%      disp('Calculating U...');
% for i=1:iter
% E=abs(X-U*V');
% [U,V] = EfficientMCL2(X, W, U,V, 100, 0.00000001);
% W=1./(abs(X-U*V')+0.001);
% if sum(sum(abs(E-abs(X-U*V'))))/sum(E(:))<0.000001
% break;
% end
% end
 disp('Calculating the model...');
XX=U*V';
noi=X-XX;
ind=randperm(size(noi,1)*size(noi,2));
noi_sample=noi(ind(1:fix(1*size(noi,1)*size(noi,2))))';
[nmodel,label,R] =mogcluster(noi_sample,k); 
 Sigma=nmodel.Sigma;
 W=ones(size(X));
for j = 1:k
    W(:) = W(:) + R(:,j)/(2*Sigma(j));
end
   disp('Warm start is over. ');
m=size(X,1);
A=zeros(r,r,m);B=zeros(r,m);
for i=1:m
    A(:,:,i)=(V'*diag(W(i,:))*V+0.001*eye(r))^-1*0.2;%%%%%%%%%%%%%%%%%%update sudu
    B(:,i)=V'*diag(W(i,:))*(X(i,:))'/0.2;%%%%%%%%%%%%%%%%%%%%
end
% A=repmat(eye(r),[1,1,m]);
% B=U';
model.Sigma=nmodel.Sigma;
model.weight=nmodel.weight;
model.mu=nmodel.mu;
model.A=A;
model.B=B;
model.U=U;
end