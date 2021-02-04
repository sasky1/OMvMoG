function [L,E,F,model,H,VV,w_x]= MV_OMoGMF_subsample(model,X,lamda,ro,H,VV)
%online mog matrix factorization
%param.N controls the speed of updating model
%param.lamda controls the speed of updating U
% if nargin<4
%   w_X=ones(size(X));
% end
% siz=model.siz;
% U=model.U;
% if nargin<4
%     
% end

V=size(X,2);
for i=1:V
%  siz{i}=model{i}.siz;
[m{i},n{i}]=size(X{i});
E{i}=zeros(m{i},n{i});F{i}=zeros(m{i},n{i});
W_X{i} = [];
end
tic
 for i=1:size(X{1},2)   
if mod(i,20)==0||i==1
      disp(['Calculating the model of the ',num2str(i),'th frame']);
end
for j=1:V
 Ind=randperm(size(X{j},1));
 ind{j}=Ind(1:fix(ro*size(X{j},1)));
x{j}=X{j}(ind{j},i);
h{j}=H{j}(ind{j},i);
w_x{j}=ones(fix(ro*size(X{j},1)),1);
n_model{j}=model{j};
n_model{j}.U=model{j}.U(ind{j},:);
n_model{j}.A=model{j}.A(:,:,ind{j});
n_model{j}.B=model{j}.B(:,ind{j});
end
v=VV(:,i);
[n_model,v,~,~,vv,h] =MV_onlinemogmf(n_model,x,lamda,w_x,v,h);  

for j=1:V
model{j}.U(ind{j},:)=n_model{j}.U;
model{j}.A(:,:,ind{j})=n_model{j}.A;
model{j}.B(:,ind{j})=n_model{j}.B;       
model{j}.weight=n_model{j}.weight;
model{j}.mu=n_model{j}.mu;
model{j}.Sigma=n_model{j}.Sigma;   
 L{j}(:,i)=model{j}.U*v{j};
 E{j}(:,i)=X{j}(:,i)-L{j}(:,i);
%  [~,a]=sort(model{i}.Sigma);
%  F0=(X(:,i)-U*v).*(label0==a(end));
%  FF=reshape(E(:,i),siz);
% F(:,i)=F0(:);
H{j}(ind{j},i)=h{j};
end
VV(:,i)=vv;
 end
 end