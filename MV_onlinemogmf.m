function [nmodel,v,label,R,vv,h] =MV_onlinemogmf(model,x,lamda,w_x,vv,h)
V=size(x,2);
if nargin <4
for i=1:V
w_x{i}=ones(size(x{i}));
end
end
[nmodel,label,R,h,v,vv]=MV_update_MoGparameters(x,model,lamda,w_x,vv,h);
nmodel=MV_update_subspace(x,h,v,nmodel);
