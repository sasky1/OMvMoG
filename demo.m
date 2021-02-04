% % Contact Email: lijinxing158@gmail.com or cshyong@comp.polyu.edu.hk
% %

clear all;clc;
load ('data.mat')% a Cell Maxtrix X: X{i} is D*N, D is the dimension and N is the number
lambda=5;
Num_train = size(X{1},2);
Num_rank  = 30;
Num_view  = numel(X);
Num_iter  = 2;

for i=1:Num_view
    Y{i}=X{i};
    W_X{i}=ones(size(X{i}));
end

ind=randperm(Num_train);
for i=1:Num_view
    [model{i}]  = warmstart(X{i}(:,ind(1:100)),Num_rank,4);
    model{i}.N=50;model{i}.lamda=0.97;
    H{i}=ones(size(X{i}));
end
VV=zeros(Num_rank,Num_train);
for i=1:Num_iter
    [L,E,F,model,label,H,VV]= MV_OMoGMF(model,X,lambda,W_X,H,VV); %L is the reconstructed data
%     ro = 0.01;% The preserved rate
%     [L,E,F,model,H,VV,w_x]= MV_OMoGMF_subsample(model,X,lambda,ro,H,VV);
end