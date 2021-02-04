function [model] = warmstart_MV(X_start,r,k)
% Use reweighted L2 matrix factorization to calculate U
% Use mog to calculatethe model
 disp('Warm-start...');
 disp('Calculating U...');
 X=[];
 for i=1:size(X_start,2)
 X=[X;X_start{i}];
 end
 
W=ones(size(X));
%RPCA for warmstart
U=[];
[L,~] = inexact_alm_rpca(X,0.5/sqrt(max(size(X,1),size(X,2))));
 [U,V] = PCA(L,r);%U=[U,median(X,2)];V=X'*U/(U'*U+0.0001*eye(r));
%PCA for warmstart
% [U,~] = PCA(X,r-1);
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
% ind=randperm(size(noi,1)*size(noi,2));
% noi_sample=noi(ind(1:fix(1*size(noi,1)*size(noi,2))))';
% [nmodel,label,R] =mogcluster(noi,k); 
%  Sigma=nmodel.Sigma;
 W=ones(size(X));
% for j = 1:k
%     W(:) = W(:) + R(:,j)/(2*Sigma(j));
% end

s=1;
for j=1:size(X_start,2)
 m=size(X_start{j},1);  
A{j}=zeros(r,r,m);B{j}=zeros(r,m);
for i=1:m
    A{j}(:,:,i)=(V'*diag(W(i+s-1,:))*V+0.001*eye(r))^-1*0.2;
    B{j}(:,i)=V'*diag(W(i+s-1,:))*(X(i+s-1,:))'/0.2;
end
UU{j}=U(s:(s+m-1),:);noise=noi(s:(s+m-1),:);
[nmodel,~,R] =mogcluster(noise(:),k); 
model{j}.Sigma=nmodel.Sigma;
model{j}.weight=nmodel.weight;
model{j}.mu=nmodel.mu;
model{j}.A=A{j};
model{j}.B=B{j};
model{j}.U=UU{j};
s=s+m;
end
   disp('Warm start is over. ');
% A=repmat(eye(r),[1,1,m]);
% B=U';

% for j=1:size(X_start,2)
% model{j}.Sigma=nmodel.Sigma;
% model{j}.weight=nmodel.weight;
% model{j}.mu=nmodel.mu;
% model{j}.A=A{j};
% model{j}.B=B{j};
% model{j}.U=UU{j};
% end

end