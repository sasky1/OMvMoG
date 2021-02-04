function [label,model,W,U,V,A,B] = onlinestart1(InW,X,r,k)


[~ ,U ,V]=RobustApproximation_M_UV_TraceNormReg(X,InW,r,1e-3,1.2,1,0);
V=V';
iter=3;eps=0.2;
for i=1:iter
    sum(sum(abs(X-U*V')));
W=1./(abs(X-U*V').^3+eps);
% [U,V] = EfficientMCL2(X, W, U,V, 100, 0.00000001);
[~,U,V] = RobustApproximation_M_UV_TraceNormReg(X,W,r,1e-3,1.2,1,0);
V=V';
end

XX=U*V';

noi=X-XX;

[label,mu,sigma] = maxguass(noi(:)',k);

for i=1:k
nk(i)=length(find(label==i))/length(label);
end
model.Sigma=sigma;model.weight=nk;model.mu=mu;
[m,n]=size(X);W=ones(m,n);
lab=reshape(label,m,n);
% for i=1:m
%     for j=1:n
% W(i,j)=1/(2*sigma(lab(i,j)));
%     end
% end
A=zeros(r,r,m);B=zeros(r,m);
for i=1:m
    A(:,:,i)=(V'*diag(W(i,:))*V)^-1;
    B(:,i)=V'*diag(W(i,:))*(X(i,:))';
end






