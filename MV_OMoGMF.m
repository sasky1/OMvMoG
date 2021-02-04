function[L,E,F,model,label,H,VV,Onoise,Ow]= MV_OMoGMF(model,X,lamda,w_X,H,VV)
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
E{i}=zeros(m{i},n{i});label{i}=zeros(m{i},n{i});F{i}=zeros(m{i},n{i});
end
tic
 for i=1:size(X{1},2)   
if mod(i,20)==0||i==1
      disp(['Calculating the model of the ',num2str(i),'th frame']);
end
for j=1:V
x{j}=X{j}(:,i);
w_x{j}=w_X{j}(:,i);
h{j}=H{j}(:,i);
end
v=VV(:,i);
[model,v,label0,~,vv,h] =MV_onlinemogmf(model,x,lamda,w_x,v,h);  

for j=1:V
 label{j}(:,i)=label0{j};
 L{j}(:,i)=model{j}.U*v{j};
 E{j}(:,i)=X{j}(:,i)-L{j}(:,i);
 
 Onoise{j,i}=model{j}.Sigma;
 Ow{j,i}=model{j}.weight;
%  [~,a]=sort(model{i}.Sigma);
%  F0=(X(:,i)-U*v).*(label0==a(end));
%  FF=reshape(E(:,i),siz);
% F(:,i)=F0(:);
H{j}(:,i)=h{j};
end
VV(:,i)=vv;
 end


 end