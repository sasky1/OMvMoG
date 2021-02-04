function [label,model,W,U,V,A,B] = onlinestart(InW,X,r,k,mod)

if nargin<5
    mod=1;
end

[U V] = RPMF(X, r, 1, 1, 1e-2);
V=V';
XX=U*V';
noi=X-XX;

if mod==1%soft cluster
    [model label,R] =mogcluster(noi(:),k); 
end

if mod==2%hard cluster
   [label,mu,sigma] = maxguass(noi(:)',k);
for i=1:k
nk(i)=length(find(label==i))/length(label);
end
model.Sigma=sigma;model.weight=nk;model.mu=mu;
end

[m,n]=size(X);W=ones(m,n);%set the weight to a ones matrix;

A=zeros(r,r,m);B=zeros(r,m);
for i=1:m
    A(:,:,i)=(V'*diag(W(i,:))*V)^-1;
    B(:,i)=V'*diag(W(i,:))*(X(i,:))';
end
