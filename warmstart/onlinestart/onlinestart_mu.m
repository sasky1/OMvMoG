function [label,model,W,U,V,A,B] = onlinestart_mu(InW,X,r,k,mod)

if nargin<5
    mod=1;
end

[U V] = RPMF(X, r, 1, 1, 1e-2);
V=V';
XX=U*V';
noi=X-XX;

if mod==1%soft cluster
    [model label,R] =mogcluster_mu(noi(:),k); 
end
for i=1:k
RR(:,:,k)=reshape(R(:,k),size(X,1),size(X,2));
end
if mod==2%hard cluster
   [label,mu,sigma] = maxguass(noi(:)',k);
for i=1:k
nk(i)=length(find(label==i))/length(label);
end
model.Sigma=sigma;model.weight=nk;model.mu=mu;
end
sigma=model.Sigma;nk=model.weight;mu=model.mu;
[m,n]=size(X);W=ones(m,n);%set the weight to a ones matrix;

A=zeros(r,r,m);B=zeros(r,m);
for i=1:m
    for j=1:k
%         size(V')
%         size(repmat(RR(:,i,k)',r,1))
      A(:,:,i)=A(:,:,i)+V'.*repmat(RR(i,:,j),r,1)*V/(2*sigma(j)) ;
         B(:,i)=B(:,i)+V'.*repmat(RR(i,:,j),r,1)*(bsxfun(@minus,X(i,:),mu(j)))'/(2*sigma(j));
    end
    A(:,:,i)=A(:,:,i)^-1;
end
